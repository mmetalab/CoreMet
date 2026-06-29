"""
REST API routes for CoreMet Web Application.

Provides programmatic access to interaction prediction, database search,
species listing, and job result retrieval.
"""

import io
import json
import logging
import tempfile
import zipfile
from pathlib import Path
from threading import Thread

import pandas as pd
from flask import Blueprint, Response, after_this_request, jsonify, request, send_file

from flask_limiter import Limiter
from flask_limiter.util import get_remote_address

from app.config import Config
from app.services.data_service import DataService
from app.services.job_service import create_job, get_job, update_job

# ---------------------------------------------------------------------------
# Blueprint & rate-limiter
# ---------------------------------------------------------------------------

api_bp = Blueprint("api", __name__, url_prefix="/api/v1")

limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["100 per hour"],
    storage_uri="memory://",
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Lazy singletons — initialised on first request so the import stays light
# ---------------------------------------------------------------------------

_config: Config | None = None
_data_service: DataService | None = None


def _get_config() -> Config:
    global _config
    if _config is None:
        _config = Config()
    return _config


def _get_data_service() -> DataService:
    global _data_service
    if _data_service is None:
        _data_service = DataService()
    return _data_service


def _data_dir() -> Path:
    return _get_config().DATA_DIR


def _release_dir() -> Path:
    return _data_dir() / "databases" / "release"


_DOWNLOAD_FILES = {
    "mpi": "coremetdb_mpi.csv",
    "mei": "coremetdb_mei.csv",
    "mdi": "coremetdb_mdi.csv",
    "mmi": "coremetdb_mmi.csv",
    "mdri": "coremetdb_mdri.csv",
    "mgi": "coremetdb_mgi.csv",
    "mgwas": "coremetdb_mgwas.csv",
}


def _send_existing_file(path: Path, download_name: str, mimetype: str):
    if not path.exists():
        return jsonify({"error": f"{download_name} is not available on this deployment"}), 404
    return send_file(path, as_attachment=True, download_name=download_name, mimetype=mimetype)


# ---------------------------------------------------------------------------
# GET /api/v1/health and /api/v1/stats
# ---------------------------------------------------------------------------

@api_bp.route("/stats", methods=["GET"])
@limiter.limit("100 per hour")
def stats():
    """Return canonical CoreMet database statistics."""
    stats_path = _data_dir() / "coremetdb_stats.json"
    if not stats_path.exists():
        return jsonify({"error": "Database statistics are not available"}), 404
    try:
        return jsonify(json.loads(stats_path.read_text())), 200
    except json.JSONDecodeError:
        return jsonify({"error": "Database statistics file is invalid"}), 500


# ---------------------------------------------------------------------------
# POST /api/v1/predict
# ---------------------------------------------------------------------------

@api_bp.route("/predict", methods=["POST"])
@limiter.limit("100 per hour")
def predict():
    """Accept metabolite + protein lists and launch an async prediction job.

    Request JSON
    -------------
    {
        "metabolites": [
            {"name": "...", "hmdb_id": "...", "smiles": "..."}
        ],
        "proteins": [
            {"uniprot_id": "...", "name": "...", "gene": "...",
             "organism": "...", "sequence": "..."}
        ],
        "organism": "Homo sapiens"          # optional, default "All"
    }

    Response JSON
    --------------
    {"job_id": "<id>", "status": "running"}
    """
    data = request.get_json(silent=True)
    if data is None:
        return jsonify({"error": "Request body must be valid JSON"}), 400

    # --- validate required fields ----------------------------------------
    metabolites = data.get("metabolites")
    proteins = data.get("proteins")

    if not metabolites or not isinstance(metabolites, list):
        return jsonify({"error": "'metabolites' must be a non-empty array"}), 400
    if not proteins or not isinstance(proteins, list):
        return jsonify({"error": "'proteins' must be a non-empty array"}), 400

    # Validate each metabolite entry
    for idx, met in enumerate(metabolites):
        for field in ("name", "hmdb_id", "smiles"):
            if not met.get(field):
                return jsonify({
                    "error": f"metabolites[{idx}] is missing required field '{field}'"
                }), 400

    # Validate each protein entry
    for idx, prot in enumerate(proteins):
        for field in ("uniprot_id", "name", "gene", "organism", "sequence"):
            if not prot.get(field):
                return jsonify({
                    "error": f"proteins[{idx}] is missing required field '{field}'"
                }), 400

    organism = data.get("organism", "All")

    # --- build DataFrames matching internal column conventions -------------
    met_df = pd.DataFrame([
        {
            "Metabolite Name": m["name"],
            "HMDB ID": m["hmdb_id"],
            "SMILES": m["smiles"],
        }
        for m in metabolites
    ])

    prot_df = pd.DataFrame([
        {
            "UniprotID": p["uniprot_id"],
            "Protein Name": p["name"],
            "Gene Name": p["gene"],
            "Organism": p["organism"],
            "Sequence": p["sequence"],
        }
        for p in proteins
    ])

    # --- persist the job and kick off background prediction ----------------
    job_id = create_job(
        metabolites_json=met_df.to_json(orient="records"),
        proteins_json=prot_df.to_json(orient="records"),
        organism=organism,
    )

    thread = Thread(
        target=_run_prediction,
        args=(job_id, met_df, prot_df, organism),
        daemon=True,
    )
    thread.start()

    return jsonify({"job_id": job_id, "status": "running"}), 202


def _run_prediction(
    job_id: str,
    met_df: pd.DataFrame,
    prot_df: pd.DataFrame,
    organism: str,
):
    """Execute prediction in a background thread and persist the results."""
    try:
        # Import here to avoid heavy DGL/torch load at module level
        from app.services.prediction_service import PredictionService

        service = PredictionService()
        results_df = service.predict_interactions(met_df, prot_df, organism)
        update_job(job_id, "completed", results_df.to_json(orient="records"))
    except Exception as exc:
        logger.exception("Prediction failed for job %s", job_id)
        update_job(job_id, "failed", json.dumps({"error": str(exc)}))


# ---------------------------------------------------------------------------
# GET /api/v1/database/search
# ---------------------------------------------------------------------------

@api_bp.route("/database/search", methods=["GET"])
@limiter.limit("100 per hour")
def database_search():
    """Search the CoreMet database for metabolite-protein interactions.

    Query parameters
    -----------------
    metabolite : str   — partial match on Metabolite Name or HMDB ID
    protein    : str   — partial match on Protein Name, Uniprot ID, or Gene Name
    organism   : str   — exact match on Species column
    pathway    : str   — partial match on Pathway_ID (if column exists)
    limit      : int   — max rows to return (default 100, max 1000)

    Returns JSON array of matching interaction records.
    """
    ds = _get_data_service()
    mpi_db = ds.mpi_db

    if mpi_db is None or mpi_db.empty:
        return jsonify({"error": "Interaction database is not loaded"}), 503

    # --- query length guard ------------------------------------------------
    for _param in ("metabolite", "protein", "organism", "pathway"):
        if len(request.args.get(_param, "")) > _MAX_QUERY_LEN:
            return jsonify({"error": f"'{_param}' exceeds {_MAX_QUERY_LEN} characters"}), 400

    mask = pd.Series(True, index=mpi_db.index)

    # --- metabolite filter ------------------------------------------------
    metabolite_q = request.args.get("metabolite", "").strip()
    if metabolite_q:
        met_mask = pd.Series(False, index=mpi_db.index)
        for col in ("Metabolite Name", "HMDB ID"):
            if col in mpi_db.columns:
                met_mask |= mpi_db[col].astype(str).str.contains(
                    metabolite_q, case=False, na=False, regex=False
                )
        mask &= met_mask

    # --- protein filter ---------------------------------------------------
    protein_q = request.args.get("protein", "").strip()
    if protein_q:
        prot_mask = pd.Series(False, index=mpi_db.index)
        for col in ("Protein Name", "Uniprot ID", "Gene Name"):
            if col in mpi_db.columns:
                prot_mask |= mpi_db[col].astype(str).str.contains(
                    protein_q, case=False, na=False, regex=False
                )
        mask &= prot_mask

    # --- organism filter --------------------------------------------------
    organism_q = request.args.get("organism", "").strip()
    if organism_q and "Species" in mpi_db.columns:
        mask &= mpi_db["Species"].str.lower() == organism_q.lower()

    # --- pathway filter ---------------------------------------------------
    pathway_q = request.args.get("pathway", "").strip()
    if pathway_q and "Pathway_ID" in mpi_db.columns:
        mask &= mpi_db["Pathway_ID"].astype(str).str.contains(
            pathway_q, case=False, na=False, regex=False
        )

    # --- apply limit ------------------------------------------------------
    try:
        limit = min(int(request.args.get("limit", 100)), 1000)
    except (TypeError, ValueError):
        limit = 100

    results = mpi_db.loc[mask].head(limit)

    return jsonify(results.to_dict(orient="records")), 200


# ---------------------------------------------------------------------------
# GET /api/v1/species
# ---------------------------------------------------------------------------

@api_bp.route("/species", methods=["GET"])
@limiter.limit("100 per hour")
def list_species():
    """Return supported species together with their interaction counts.

    Response JSON
    --------------
    {
        "species": [
            {"name": "Homo sapiens", "mpi_count": 12345},
            ...
        ]
    }
    """
    ds = _get_data_service()
    mpi_db = ds.mpi_db

    if mpi_db is None or mpi_db.empty:
        return jsonify({"error": "Interaction database is not loaded"}), 503

    if "Species" not in mpi_db.columns:
        return jsonify({"error": "Species column not found in database"}), 500

    counts = mpi_db["Species"].value_counts()
    species_list = [
        {"name": species, "mpi_count": int(count)}
        for species, count in counts.items()
    ]

    return jsonify({"species": species_list}), 200


# ---------------------------------------------------------------------------
# GET /api/v1/results/<job_id>
# ---------------------------------------------------------------------------

@api_bp.route("/results/<job_id>", methods=["GET"])
@limiter.limit("100 per hour")
def get_results(job_id: str):
    """Retrieve prediction results for a previously submitted job.

    Response JSON (job completed)
    ------------------------------
    {
        "job_id": "...",
        "status": "completed",
        "organism": "...",
        "created_at": "...",
        "expires_at": "...",
        "results": [ ... ]
    }
    """
    job = get_job(job_id)
    if job is None:
        return jsonify({"error": f"Job '{job_id}' not found"}), 404

    response = {
        "job_id": job["job_id"],
        "status": job["status"],
        "organism": job["organism"],
        "created_at": job["created_at"],
        "expires_at": job["expires_at"],
    }

    if job["status"] == "completed" and job["results"]:
        try:
            response["results"] = json.loads(job["results"])
        except (json.JSONDecodeError, TypeError):
            response["results"] = []
    elif job["status"] == "failed" and job["results"]:
        try:
            response["error_detail"] = json.loads(job["results"])
        except (json.JSONDecodeError, TypeError):
            response["error_detail"] = job["results"]

    return jsonify(response), 200


# ---------------------------------------------------------------------------
# GET /api/v1/health
# ---------------------------------------------------------------------------

@api_bp.route("/health", methods=["GET"])
@limiter.limit("100 per hour")
def health_check():
    """Lightweight health check endpoint for deployment monitoring."""
    return jsonify({"status": "ok"}), 200


# ---------------------------------------------------------------------------
# GET /api/v1/mmi/stats
# ---------------------------------------------------------------------------

@api_bp.route("/mmi/stats", methods=["GET"])
@limiter.limit("100 per hour")
def mmi_stats():
    """Return summary statistics for the MMI (Metabolite-Microbe) database.

    Response JSON
    --------------
    {
        "total": 83149,
        "metabolites": ...,
        "microbes": ...,
        "relationship_types": ...,
        "organisms": {...},
        "available": true
    }
    """
    try:
        from app.services.mmi_service import get_mmi_stats
        stats = get_mmi_stats()
        return jsonify(stats), 200
    except Exception:
        logger.exception("Failed to retrieve MMI stats")
        return jsonify({"error": "Service temporarily unavailable", "status": 500}), 500


# ---------------------------------------------------------------------------
# GET /api/v1/mmi/search
# ---------------------------------------------------------------------------

@api_bp.route("/mmi/search", methods=["GET"])
@limiter.limit("100 per hour")
def mmi_search():
    """Search the MMI database for metabolite-microbe interactions.

    Query parameters
    -----------------
    q          : str   — text search across metabolite/microbe names and IDs
    metabolite : str   — partial match on Metabolite_Name or HMDB_ID
    microbe    : str   — partial match on Microbe_Name or Taxonomy_ID
    organism   : str   — exact match: 'human' or 'mouse'
    relationship : str — exact match: 'causal' or 'correlative'
    limit      : int   — max rows (default 100, max 1000)

    Returns JSON array of matching interaction records.
    """
    try:
        from app.services.mmi_service import get_mmi_db
        mmi_db = get_mmi_db()
    except Exception:
        return jsonify({"error": "MMI database is not loaded"}), 503

    if mmi_db.empty:
        return jsonify([]), 200

    for _param in ("q", "metabolite", "microbe"):
        if len(request.args.get(_param, "")) > _MAX_QUERY_LEN:
            return jsonify({"error": f"'{_param}' exceeds {_MAX_QUERY_LEN} characters"}), 400

    mask = pd.Series(True, index=mmi_db.index)

    # Free-text search
    q = request.args.get("q", "").strip()
    if q:
        q_mask = pd.Series(False, index=mmi_db.index)
        for col in ("Metabolite_Name", "HMDB_ID", "KEGG_ID", "PubChem_CID",
                     "Microbe_Name", "Taxonomy_ID"):
            if col in mmi_db.columns:
                q_mask |= mmi_db[col].astype(str).str.contains(q, case=False, na=False, regex=False)
        mask &= q_mask

    # Metabolite filter
    met_q = request.args.get("metabolite", "").strip()
    if met_q:
        met_mask = pd.Series(False, index=mmi_db.index)
        for col in ("Metabolite_Name", "HMDB_ID"):
            if col in mmi_db.columns:
                met_mask |= mmi_db[col].astype(str).str.contains(met_q, case=False, na=False, regex=False)
        mask &= met_mask

    # Microbe filter
    mic_q = request.args.get("microbe", "").strip()
    if mic_q:
        mic_mask = pd.Series(False, index=mmi_db.index)
        for col in ("Microbe_Name", "Taxonomy_ID"):
            if col in mmi_db.columns:
                mic_mask |= mmi_db[col].astype(str).str.contains(mic_q, case=False, na=False, regex=False)
        mask &= mic_mask

    # Organism filter
    org_q = request.args.get("organism", "").strip()
    if org_q and "Organism" in mmi_db.columns:
        mask &= mmi_db["Organism"].str.lower() == org_q.lower()

    # Relationship type filter
    rel_q = request.args.get("relationship", "").strip()
    if rel_q and "Relationship_Type" in mmi_db.columns:
        mask &= mmi_db["Relationship_Type"].str.lower() == rel_q.lower()

    try:
        limit = min(int(request.args.get("limit", 100)), 1000)
    except (TypeError, ValueError):
        limit = 100

    results = mmi_db.loc[mask].head(limit)
    return jsonify(results.to_dict(orient="records")), 200


# ---------------------------------------------------------------------------
# GET /api/v1/mdri/stats
# ---------------------------------------------------------------------------

@api_bp.route("/mdri/stats", methods=["GET"])
@limiter.limit("100 per hour")
def mdri_stats():
    """Return summary statistics for the MDrI (Metabolite-Drug) database.

    Response JSON
    --------------
    {
        "total": 3500,
        "metabolites": 1200,
        "drugs": 350,
        "interaction_types": {"PK": 1800, "PD": 1700},
        "sources": {"CoreMet_curated": 97, "DrugBank_cross_ref": 3283, "enzyme_bridge": 120},
        "available": true
    }
    """
    try:
        from app.services.mdri_service import get_mdri_stats
        stats = get_mdri_stats()
        return jsonify(stats), 200
    except Exception:
        logger.exception("Failed to retrieve MDrI stats")
        return jsonify({"error": "Service temporarily unavailable", "status": 500}), 500


# ---------------------------------------------------------------------------
# GET /api/v1/mdri/search
# ---------------------------------------------------------------------------

@api_bp.route("/mdri/search", methods=["GET"])
@limiter.limit("100 per hour")
def mdri_search():
    """Search the MDrI database for metabolite-drug interactions.

    Query parameters
    -----------------
    q          : str   — text search across metabolite/drug names and IDs
    metabolite : str   — partial match on Metabolite_Name or HMDB_ID
    drug       : str   — partial match on Drug_Name or DrugBank_ID
    interaction_type : str — exact match: 'PK' or 'PD'
    evidence   : str   — exact match on Evidence_Level
    limit      : int   — max rows (default 100, max 1000)

    Returns JSON array of matching interaction records.
    """
    try:
        from app.services.mdri_service import get_mdri_db
        mdri_db = get_mdri_db()
    except Exception:
        return jsonify({"error": "MDrI database is not loaded"}), 503

    if mdri_db.empty:
        return jsonify([]), 200

    for _param in ("q", "metabolite", "drug"):
        if len(request.args.get(_param, "")) > _MAX_QUERY_LEN:
            return jsonify({"error": f"'{_param}' exceeds {_MAX_QUERY_LEN} characters"}), 400

    mask = pd.Series(True, index=mdri_db.index)

    # Free-text search
    q = request.args.get("q", "").strip()
    if q:
        q_mask = pd.Series(False, index=mdri_db.index)
        for col in ("Metabolite_Name", "HMDB_ID", "Drug_Name", "DrugBank_ID"):
            if col in mdri_db.columns:
                q_mask |= mdri_db[col].astype(str).str.contains(
                    q, case=False, na=False, regex=False
                )
        mask &= q_mask

    # Metabolite filter
    met_q = request.args.get("metabolite", "").strip()
    if met_q:
        met_mask = pd.Series(False, index=mdri_db.index)
        for col in ("Metabolite_Name", "HMDB_ID"):
            if col in mdri_db.columns:
                met_mask |= mdri_db[col].astype(str).str.contains(
                    met_q, case=False, na=False, regex=False
                )
        mask &= met_mask

    # Drug filter
    drug_q = request.args.get("drug", "").strip()
    if drug_q:
        drug_mask = pd.Series(False, index=mdri_db.index)
        for col in ("Drug_Name", "DrugBank_ID"):
            if col in mdri_db.columns:
                drug_mask |= mdri_db[col].astype(str).str.contains(
                    drug_q, case=False, na=False, regex=False
                )
        mask &= drug_mask

    # Interaction type filter
    itype_q = request.args.get("interaction_type", "").strip()
    if itype_q and "Interaction_Type" in mdri_db.columns:
        mask &= mdri_db["Interaction_Type"].str.upper() == itype_q.upper()

    # Evidence level filter
    ev_q = request.args.get("evidence", "").strip()
    if ev_q and "Evidence_Level" in mdri_db.columns:
        mask &= mdri_db["Evidence_Level"].str.lower() == ev_q.lower()

    try:
        limit = min(int(request.args.get("limit", 100)), 1000)
    except (TypeError, ValueError):
        limit = 100

    results = mdri_db.loc[mask].head(limit)
    return jsonify(results.to_dict(orient="records")), 200


# ---------------------------------------------------------------------------
# GET /api/v1/mgwas/stats
# ---------------------------------------------------------------------------

@api_bp.route("/mgwas/stats", methods=["GET"])
@limiter.limit("100 per hour")
def mgwas_stats():
    """Return summary statistics for the mGWAS (Metabolite-SNP) database.

    Response JSON
    -------------
    total, metabolites, snps, genes, chromosomes, sources, available
    """
    try:
        from app.services.mgwas_service import get_mgwas_stats
        stats = get_mgwas_stats()
        return jsonify(stats), 200
    except Exception:
        logger.exception("Failed to retrieve mGWAS stats")
        return jsonify({"error": "Service temporarily unavailable", "status": 500}), 500


# ---------------------------------------------------------------------------
# GET /api/v1/mgwas/search
# ---------------------------------------------------------------------------

@api_bp.route("/mgwas/search", methods=["GET"])
@limiter.limit("100 per hour")
def mgwas_search():
    """Search the mGWAS database for metabolite-SNP associations.

    Query parameters
    ----------------
    q          – free-text search across all columns
    metabolite – filter by metabolite name or HMDB ID
    snp        – filter by rsID
    gene       – filter by mapped gene
    chromosome – filter by chromosome (1-22, X)
    limit      – max results (default 100, max 1000)
    """
    try:
        from app.services.mgwas_service import get_mgwas_db
        mgwas_db = get_mgwas_db()
    except Exception:
        return jsonify({"error": "mGWAS database is not loaded"}), 503

    if mgwas_db.empty:
        return jsonify([]), 200

    for _param in ("q", "metabolite", "snp", "gene"):
        if len(request.args.get(_param, "")) > _MAX_QUERY_LEN:
            return jsonify({"error": f"'{_param}' exceeds {_MAX_QUERY_LEN} characters"}), 400

    mask = pd.Series(True, index=mgwas_db.index)

    # Free-text search
    q = request.args.get("q", "").strip()
    if q:
        q_mask = pd.Series(False, index=mgwas_db.index)
        for col in ("Metabolite_Name", "HMDB_ID", "rsID", "Mapped_Gene", "Trait"):
            if col in mgwas_db.columns:
                q_mask |= mgwas_db[col].astype(str).str.contains(q, case=False, na=False, regex=False)
        mask &= q_mask

    # Metabolite filter
    met_q = request.args.get("metabolite", "").strip()
    if met_q:
        met_mask = pd.Series(False, index=mgwas_db.index)
        for col in ("Metabolite_Name", "HMDB_ID"):
            if col in mgwas_db.columns:
                met_mask |= mgwas_db[col].astype(str).str.contains(met_q, case=False, na=False, regex=False)
        mask &= met_mask

    # SNP filter
    snp_q = request.args.get("snp", "").strip()
    if snp_q:
        mask &= mgwas_db["rsID"].astype(str).str.contains(snp_q, case=False, na=False, regex=False)

    # Gene filter
    gene_q = request.args.get("gene", "").strip()
    if gene_q:
        mask &= mgwas_db["Mapped_Gene"].astype(str).str.contains(gene_q, case=False, na=False, regex=False)

    # Chromosome filter
    chr_q = request.args.get("chromosome", "").strip()
    if chr_q and "Chromosome" in mgwas_db.columns:
        mask &= mgwas_db["Chromosome"] == chr_q

    try:
        limit = min(int(request.args.get("limit", 100)), 1000)
    except (TypeError, ValueError):
        limit = 100

    results = mgwas_db.loc[mask].head(limit)
    return jsonify(results.to_dict(orient="records")), 200


# ---------------------------------------------------------------------------
# GET /api/v1/mgi/stats
# ---------------------------------------------------------------------------

@api_bp.route("/mgi/stats", methods=["GET"])
@limiter.limit("100 per hour")
def mgi_stats():
    """Return summary statistics for the MGI (Metabolite-Gene) database.

    Response JSON
    -------------
    total, metabolites, genes, organisms, interaction_types, sources, available
    """
    try:
        from app.services.mgi_service import get_mgi_stats
        stats = get_mgi_stats()
        return jsonify(stats), 200
    except Exception:
        logger.exception("Failed to retrieve MGI stats")
        return jsonify({"error": "Service temporarily unavailable", "status": 500}), 500


# ---------------------------------------------------------------------------
# GET /api/v1/mgi/search
# ---------------------------------------------------------------------------

@api_bp.route("/mgi/search", methods=["GET"])
@limiter.limit("100 per hour")
def mgi_search():
    """Search the MGI database for metabolite-gene interactions.

    Query parameters
    ----------------
    q          – free-text search across all columns
    metabolite – filter by metabolite name or HMDB ID
    gene       – filter by gene symbol or gene ID
    organism   – filter by organism name
    limit      – max results (default 100, max 1000)
    """
    try:
        from app.services.mgi_service import get_mgi_db
        mgi_db = get_mgi_db()
    except Exception:
        return jsonify({"error": "MGI database is not loaded"}), 503

    if mgi_db.empty:
        return jsonify([]), 200

    for _param in ("q", "metabolite", "gene"):
        if len(request.args.get(_param, "")) > _MAX_QUERY_LEN:
            return jsonify({"error": f"'{_param}' exceeds {_MAX_QUERY_LEN} characters"}), 400

    mask = pd.Series(True, index=mgi_db.index)

    # Free-text search
    q = request.args.get("q", "").strip()
    if q:
        q_mask = pd.Series(False, index=mgi_db.index)
        for col in ("Metabolite_Name", "HMDB_ID", "Gene_Symbol", "Gene_ID", "Organism", "Interaction_Type"):
            if col in mgi_db.columns:
                q_mask |= mgi_db[col].astype(str).str.contains(q, case=False, na=False, regex=False)
        mask &= q_mask

    # Metabolite filter
    met_q = request.args.get("metabolite", "").strip()
    if met_q:
        met_mask = pd.Series(False, index=mgi_db.index)
        for col in ("Metabolite_Name", "HMDB_ID"):
            if col in mgi_db.columns:
                met_mask |= mgi_db[col].astype(str).str.contains(met_q, case=False, na=False, regex=False)
        mask &= met_mask

    # Gene filter
    gene_q = request.args.get("gene", "").strip()
    if gene_q:
        gene_mask = pd.Series(False, index=mgi_db.index)
        for col in ("Gene_Symbol", "Gene_ID"):
            if col in mgi_db.columns:
                gene_mask |= mgi_db[col].astype(str).str.contains(gene_q, case=False, na=False, regex=False)
        mask &= gene_mask

    # Organism filter
    org_q = request.args.get("organism", "").strip()
    if org_q and "Organism" in mgi_db.columns:
        mask &= mgi_db["Organism"].astype(str).str.contains(org_q, case=False, na=False, regex=False)

    try:
        limit = min(int(request.args.get("limit", 100)), 1000)
    except (TypeError, ValueError):
        limit = 100

    results = mgi_db.loc[mask].head(limit)
    return jsonify(results.to_dict(orient="records")), 200


# ---------------------------------------------------------------------------
# GET /api/v1/autocomplete
# ---------------------------------------------------------------------------

_autocomplete_index = None  # lazy-loaded list of (label, hmdb_id)


def _get_autocomplete_index():
    """Build or return cached autocomplete index from MPI database."""
    global _autocomplete_index
    if _autocomplete_index is not None:
        return _autocomplete_index
    try:
        cfg = _get_config()
        mpi = pd.read_csv(cfg.MPI_DB_PATH, usecols=["HMDB ID", "Metabolite Name"], dtype=str)
        mpi = mpi.dropna(subset=["Metabolite Name"]).drop_duplicates(subset=["HMDB ID"])
        _autocomplete_index = [
            {"label": f"{row['Metabolite Name']} ({row['HMDB ID']})",
             "name": row["Metabolite Name"],
             "hmdb_id": row["HMDB ID"]}
            for _, row in mpi.iterrows()
            if pd.notna(row["HMDB ID"])
        ]
    except Exception:
        _autocomplete_index = []
    return _autocomplete_index


@api_bp.route("/autocomplete", methods=["GET"])
@limiter.limit("300 per hour")
def autocomplete():
    """Return top-10 metabolite matches for a query string.

    Query parameters
    -----------------
    q     : str — partial name or HMDB ID prefix (min 2 chars)
    limit : int — max results (default 10, max 30)
    """
    q = request.args.get("q", "").strip()
    if len(q) < 2:
        return jsonify([]), 200
    if len(q) > _MAX_QUERY_LEN:
        return jsonify({"error": f"Query exceeds {_MAX_QUERY_LEN} characters"}), 400

    try:
        limit = min(int(request.args.get("limit", 10)), 30)
    except (TypeError, ValueError):
        limit = 10

    q_lower = q.lower()
    index = _get_autocomplete_index()

    # Score: exact start > word-start > contains
    scored = []
    for entry in index:
        name_l = entry["name"].lower()
        hmdb_l = entry["hmdb_id"].lower()
        if name_l.startswith(q_lower) or hmdb_l.startswith(q_lower):
            scored.append((0, entry))
        elif q_lower in name_l or q_lower in hmdb_l:
            scored.append((1, entry))
        if len(scored) >= limit * 3:  # gather enough, then sort
            break

    scored.sort(key=lambda x: x[0])
    results = [s[1] for s in scored[:limit]]
    return jsonify(results), 200


# ---------------------------------------------------------------------------
# Dataset downloads
# ---------------------------------------------------------------------------

@api_bp.route("/download/<dataset>", methods=["GET"])
@limiter.limit("30 per hour")
def download_dataset(dataset):
    """Download release CSVs and lightweight metadata files."""
    dataset = dataset.strip().lower()

    if dataset in _DOWNLOAD_FILES:
        filename = _DOWNLOAD_FILES[dataset]
        return _send_existing_file(_release_dir() / filename, filename, "text/csv")

    if dataset == "node-metadata":
        return _send_existing_file(
            _data_dir() / "coremetdb_entity_registry.json",
            "coremetdb_entity_registry.json",
            "application/json",
        )

    if dataset == "schema":
        schema_path = _get_config().BASE_DIR / "DATA_README.md"
        return _send_existing_file(schema_path, "CoreMet_DATA_README.md", "text/markdown")

    if dataset in {"full-edges", "full", "all"}:
        missing = [
            filename for filename in _DOWNLOAD_FILES.values()
            if not (_release_dir() / filename).exists()
        ]
        if missing:
            return jsonify({
                "error": "Full database release is not available on this deployment",
                "missing": missing,
            }), 404

        tmp = tempfile.NamedTemporaryFile(prefix="coremetdb_release_", suffix=".zip", delete=False)
        tmp_path = Path(tmp.name)
        tmp.close()
        with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for filename in _DOWNLOAD_FILES.values():
                archive.write(_release_dir() / filename, arcname=f"CoreMet-DB_v1/{filename}")
            stats_path = _data_dir() / "coremetdb_stats.json"
            if stats_path.exists():
                archive.write(stats_path, arcname="CoreMet-DB_v1/coremetdb_stats.json")
            readme_path = _get_config().BASE_DIR / "DATA_README.md"
            if readme_path.exists():
                archive.write(readme_path, arcname="CoreMet-DB_v1/DATA_README.md")

        @after_this_request
        def cleanup(response):
            tmp_path.unlink(missing_ok=True)
            return response

        return send_file(
            tmp_path,
            as_attachment=True,
            download_name="CoreMet-DB_v1_release.zip",
            mimetype="application/zip",
        )

    if dataset in {"embeddings", "model"}:
        return jsonify({
            "error": f"{dataset} is an optional machine-learning resource and is not bundled with the database deployment"
        }), 404

    return jsonify({"error": f"Unknown dataset '{dataset}'"}), 404


# ---------------------------------------------------------------------------
# Metabolite CSV export
# ---------------------------------------------------------------------------

_MAX_QUERY_LEN = 200
_VALID_DBS = {"mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi", "all"}


@api_bp.route("/export/metabolite", methods=["GET"])
@limiter.limit("100 per hour")
def export_metabolite_csv():
    """Download interaction data for a metabolite as CSV.

    Query params:
        id   – HMDB ID  (e.g. HMDB0000191)
        name – metabolite name (fallback)
        db   – one of mpi, mei, mdi, mmi, mdri, all  (default: all)
    """
    hmdb_id = request.args.get("id", "").strip()
    name = request.args.get("name", "").strip()
    db = request.args.get("db", "all").strip().lower()

    if not hmdb_id and not name:
        return jsonify({"error": "Provide 'id' or 'name' parameter"}), 400
    if db not in _VALID_DBS:
        return jsonify({"error": f"Invalid db. Choose from: {sorted(_VALID_DBS)}"}), 400

    from pages.metabolite_detail import _lookup_metabolite
    data = _lookup_metabolite(hmdb_id=hmdb_id, name=name)

    db_keys = ["mpi", "mei", "mdi", "mmi", "mdri"] if db == "all" else [db]
    frames = []
    for key in db_keys:
        df = data.get(key, pd.DataFrame())
        if not df.empty:
            df_copy = df.copy()
            if db == "all":
                df_copy.insert(0, "Database", key.upper())
            frames.append(df_copy)

    if not frames:
        return jsonify({"error": "No data found"}), 404

    combined = pd.concat(frames, ignore_index=True)
    buf = io.StringIO()
    combined.to_csv(buf, index=False)

    met_label = (data.get("name") or hmdb_id or name).replace(" ", "_")
    filename = f"{met_label}_{db}.csv"

    return Response(
        buf.getvalue(),
        mimetype="text/csv",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


# ---------------------------------------------------------------------------
# CoreMet-AI endpoint
# ---------------------------------------------------------------------------

@api_bp.route("/ai/query", methods=["POST"])
@limiter.limit("30 per hour")
def ai_query():
    """
    POST /api/v1/ai/query

    Execute a CoreMet-AI natural language query.

    Request body (JSON):
        {"query": "How does butyrate influence colorectal cancer?"}

    Returns the full AI pipeline result including query plan,
    ranked paths, evidence, grounded summary, and subgraph.
    """
    body = request.get_json(silent=True) or {}
    query = body.get("query", "").strip()
    if not query:
        return jsonify({"error": "Missing 'query' field"}), 400
    if len(query) > 500:
        return jsonify({"error": "Query too long (max 500 characters)"}), 400

    try:
        from app.services.ai_orchestrator import execute_ai_query
        result = execute_ai_query(query)
        return jsonify(result)
    except Exception as e:
        logger.error(f"AI query API error: {e}", exc_info=True)
        return jsonify({"error": "Internal server error"}), 500
