"""
CoreMet-AI — Entity Resolver

Deterministic resolution of free-text entity names to internal database IDs.
Uses synonym dictionaries built from all 7 databases and applies fuzzy matching
when exact match fails. Returns confidence scores.
"""

import logging
import threading
from difflib import SequenceMatcher
from typing import Dict, List, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

_synonym_index: Optional[Dict] = None
_build_lock = threading.Lock()


# ── Entity type detection heuristics ──────────────────────────────────
_TYPE_HINTS = {
    "metabolite": ["hmdb", "kegg", "chebi", "pubchem", "metabolite"],
    "disease": ["cancer", "diabetes", "syndrome", "disorder", "disease",
                 "carcinoma", "adenoma", "obesity", "asthma", "infection"],
    "gene": ["gene", "cyp", "slc", "abc", "hla-"],
    "protein": ["protein", "enzyme", "kinase", "transferase", "synthase",
                "receptor", "transporter", "albumin"],
    "drug": ["drug", "metformin", "aspirin", "statin", "inhibitor",
             "atorvastatin", "ibuprofen"],
    "microbe": ["bacterium", "lactobacillus", "bifidobacterium", "bacteroides",
                "clostridium", "streptococcus", "escherichia", "microbe"],
    "snp": ["rs"],
}

# ── Common synonym mappings (query term → canonical DB name) ──────────
_COMMON_SYNONYMS = {
    # Metabolites — common name → DB name
    "butyrate": "butyric acid",
    "butanoate": "butyric acid",
    "butanoic acid": "butyric acid",
    "propionate": "propionic acid",
    "acetate": "acetic acid",
    "lactate": "lactic acid",
    "pyruvate": "pyruvic acid",
    "succinate": "succinic acid",
    "fumarate": "fumaric acid",
    "malate": "malic acid",
    "citrate": "citric acid",
    "oxalate": "oxalic acid",
    "glutamate": "glutamic acid",
    "aspartate": "aspartic acid",
    "formate": "formic acid",
    "valerate": "valeric acid",
    "caproate": "caproic acid",
    "palmitate": "palmitic acid",
    "stearate": "stearic acid",
    "oleate": "oleic acid",
    "linoleate": "linoleic acid",
    "arachidonate": "arachidonic acid",
    "tryptophan": "l-tryptophan",
    "glucose": "d-glucose",
    "scfa": "butyric acid",
    "short-chain fatty acid": "butyric acid",
    "short-chain fatty acids": "butyric acid",
    # Disease synonyms
    "colorectal cancer": "colorectal neoplasms",
    "colon cancer": "colorectal neoplasms",
    "breast cancer": "breast neoplasms",
    "lung cancer": "lung neoplasms",
    "metabolic syndrome": "metabolic syndrome x",
    "alzheimer": "alzheimer disease",
    "alzheimers": "alzheimer disease",
    "parkinson": "parkinson disease",
    "parkinsons": "parkinson disease",
    "depression": "depressive disorder",
    "ibd": "inflammatory bowel diseases",
    "crohn": "crohn disease",
    "crohns": "crohn disease",
    "t2d": "diabetes mellitus, type 2",
    "type 2 diabetes": "diabetes mellitus, type 2",
    "t1d": "diabetes mellitus, type 1",
    "nafld": "non-alcoholic fatty liver disease",
    "nash": "non-alcoholic fatty liver disease",
    "schizophrenia": "schizophrenia",
    # Protein/gene synonyms
    "albumin": "ALB",
    "serum albumin": "ALB",
    "hsa": "ALB",
    "hemoglobin": "hemoglobin subunit alpha",
    "insulin": "insulin",
    "transferrin": "serotransferrin",
}


def _build_synonym_index() -> Dict:
    """Build a unified synonym → (canonical_name, entity_type, internal_id) index."""
    global _synonym_index
    if _synonym_index is not None:
        return _synonym_index

    with _build_lock:
        if _synonym_index is not None:
            return _synonym_index

        index: Dict[str, List[dict]] = {}

        def _add(name: str, etype: str, eid: str, canonical: str):
            key = name.strip().lower()
            if not key:
                return
            entry = {"name": canonical, "type": etype, "id": eid}
            index.setdefault(key, []).append(entry)

        # ── Metabolites (from ALL databases for maximum coverage) ─────
        _met_seen = set()
        for _db_name, _db_loader, _met_name_col, _hmdb_col in [
            ("MEI", "app.services.mei_service.get_mei_db", "Metabolite_Name", "HMDB_ID"),
            ("MDI", "app.services.mdi_service.get_mdi_db", "Metabolite_Name", "HMDB_ID"),
            ("MMI", "app.services.mmi_service.get_mmi_db", "Metabolite_Name", "HMDB_ID"),
            ("MDrI", "app.services.mdri_service.get_mdri_db", "Metabolite_Name", "HMDB_ID"),
            ("MGI", "app.services.mgi_service.get_mgi_db", "Metabolite_Name", "HMDB_ID"),
            ("mGWAS", "app.services.mgwas_service.get_mgwas_db", "Metabolite_Name", "HMDB_ID"),
        ]:
            try:
                import importlib
                _mod_path, _func = _db_loader.rsplit(".", 1)
                _mod = importlib.import_module(_mod_path)
                _df = getattr(_mod, _func)()
                if not _df.empty and _met_name_col in _df.columns and _hmdb_col in _df.columns:
                    _mets = _df[[_met_name_col, _hmdb_col]].drop_duplicates()
                    for _, r in _mets.iterrows():
                        name = str(r.get(_met_name_col, ""))
                        hmdb = str(r.get(_hmdb_col, ""))
                        _mkey = (name.lower(), hmdb)
                        if name and _mkey not in _met_seen:
                            _met_seen.add(_mkey)
                            _add(name, "metabolite", hmdb, name)
            except Exception as e:
                logger.warning(f"Entity resolver: {_db_name} metabolite index failed: {e}")

        # ── Diseases ─────────────────────────────────────────────────
        try:
            from app.services.mdi_service import get_mdi_db
            mdi = get_mdi_db()
            if not mdi.empty:
                dis = mdi[["Disease_Name"]].drop_duplicates()
                for _, r in dis.iterrows():
                    name = str(r.get("Disease_Name", ""))
                    if name:
                        _add(name, "disease", name, name)
        except Exception as e:
            logger.warning(f"Entity resolver: MDI load failed: {e}")

        # ── Genes ────────────────────────────────────────────────────
        try:
            from app.services.mgi_service import get_mgi_db
            mgi = get_mgi_db()
            if not mgi.empty:
                genes = mgi[["Gene_Symbol", "Gene_ID"]].drop_duplicates()
                for _, r in genes.iterrows():
                    sym = str(r.get("Gene_Symbol", ""))
                    gid = str(r.get("Gene_ID", ""))
                    if sym:
                        _add(sym, "gene", gid, sym)
        except Exception as e:
            logger.warning(f"Entity resolver: MGI load failed: {e}")

        # ── Drugs ────────────────────────────────────────────────────
        try:
            from app.services.mdri_service import get_mdri_db
            mdri = get_mdri_db()
            if not mdri.empty:
                drugs = mdri[["Drug_Name", "DrugBank_ID"]].drop_duplicates()
                for _, r in drugs.iterrows():
                    name = str(r.get("Drug_Name", ""))
                    dbid = str(r.get("DrugBank_ID", ""))
                    if name:
                        _add(name, "drug", dbid, name)
        except Exception as e:
            logger.warning(f"Entity resolver: MDrI load failed: {e}")

        # ── Microbes ────────────────────────────────────────────────
        try:
            from app.services.mmi_service import get_mmi_db
            mmi = get_mmi_db()
            if not mmi.empty:
                mics = mmi[["Microbe_Name", "Taxonomy_ID"]].drop_duplicates()
                for _, r in mics.iterrows():
                    name = str(r.get("Microbe_Name", ""))
                    tid = str(r.get("Taxonomy_ID", ""))
                    if name:
                        _add(name, "microbe", tid, name)
        except Exception as e:
            logger.warning(f"Entity resolver: MMI load failed: {e}")

        # ── SNPs ─────────────────────────────────────────────────────
        try:
            from app.services.mgwas_service import get_mgwas_db
            mgwas = get_mgwas_db()
            if not mgwas.empty:
                snps = mgwas[["rsID"]].drop_duplicates()
                for _, r in snps.iterrows():
                    rsid = str(r.get("rsID", ""))
                    if rsid:
                        _add(rsid, "snp", rsid, rsid)
        except Exception as e:
            logger.warning(f"Entity resolver: mGWAS load failed: {e}")

        # ── Proteins (from MEI) ──────────────────────────────────────
        try:
            from app.services.mei_service import get_mei_db
            mei = get_mei_db()
            if not mei.empty:
                prots = mei[["Enzyme_Name", "Uniprot_ID"]].drop_duplicates()
                for _, r in prots.iterrows():
                    name = str(r.get("Enzyme_Name", ""))
                    uid = str(r.get("Uniprot_ID", ""))
                    if name:
                        _add(name, "protein", uid, name)
                    if uid:
                        _add(uid, "protein", uid, name or uid)
        except Exception as e:
            logger.warning(f"Entity resolver: MEI proteins load failed: {e}")

        # ── Proteins from MPI (Protein Name, Uniprot ID) ─────────────
        try:
            from app.config import Config
            cfg_obj = Config()
            mpi_path = cfg_obj.MPI_DB_PATH
            if mpi_path and mpi_path.exists():
                mpi_df = pd.read_csv(mpi_path, low_memory=False)
                prot_col = None
                uid_col = None
                for c in mpi_df.columns:
                    cl = c.lower()
                    if 'protein' in cl and 'name' in cl:
                        prot_col = c
                    elif 'uniprot' in cl:
                        uid_col = c
                if prot_col and uid_col:
                    mpi_prots = mpi_df[[prot_col, uid_col]].drop_duplicates()
                    for _, r in mpi_prots.iterrows():
                        name = str(r.get(prot_col, ""))
                        uid = str(r.get(uid_col, ""))
                        if name and name != "nan" and len(name) > 1:
                            _add(name, "protein", uid, name)
        except Exception as e:
            logger.warning(f"Entity resolver: MPI protein index failed: {e}")

        _synonym_index = index
        logger.info(f"Entity resolver: built index with {len(index):,} unique keys")
        return _synonym_index


def _guess_type(query: str) -> Optional[str]:
    """Heuristic entity type guess from the query string."""
    q = query.lower().strip()
    # SNP is easy — starts with rs followed by digits
    if q.startswith("rs") and len(q) > 2 and q[2:].isdigit():
        return "snp"
    # HMDB ID
    if q.startswith("hmdb"):
        return "metabolite"
    for etype, hints in _TYPE_HINTS.items():
        for h in hints:
            if h in q:
                return etype
    return None


def resolve_entity(
    query: str,
    expected_type: Optional[str] = None,
    top_k: int = 5,
    min_confidence: float = 0.3,
) -> List[dict]:
    """
    Resolve a free-text entity name to database records.

    Returns a ranked list of candidates:
        [{"name": ..., "type": ..., "id": ..., "confidence": 0.0-1.0}, ...]
    """
    index = _build_synonym_index()
    q = query.strip().lower()
    if not q:
        return []

    candidates = []

    # ── Stage 0: common synonym lookup ────────────────────────────────
    canonical = _COMMON_SYNONYMS.get(q)
    if canonical:
        canon_key = canonical.lower()
        if canon_key in index:
            for entry in index[canon_key]:
                if expected_type and entry["type"] != expected_type:
                    continue
                candidates.append({**entry, "confidence": 0.98})

    # ── Stage 1: exact match ──────────────────────────────────────────
    if q in index:
        for entry in index[q]:
            if expected_type and entry["type"] != expected_type:
                continue
            candidates.append({**entry, "confidence": 1.0})

    # ── Stage 2: prefix match ─────────────────────────────────────────
    if len(candidates) < top_k:
        for key, entries in index.items():
            if key.startswith(q) and key != q:
                for entry in entries:
                    if expected_type and entry["type"] != expected_type:
                        continue
                    score = len(q) / len(key)
                    candidates.append({**entry, "confidence": round(score, 3)})

    # ── Stage 3: substring match ──────────────────────────────────────
    if len(candidates) < top_k:
        for key, entries in index.items():
            if q in key and not key.startswith(q):
                for entry in entries:
                    if expected_type and entry["type"] != expected_type:
                        continue
                    score = len(q) / len(key) * 0.8
                    if score >= min_confidence:
                        candidates.append({**entry, "confidence": round(score, 3)})

    # ── Stage 4: fuzzy match (expensive, only if needed) ──────────────
    if len(candidates) < 2:
        best_fuzzy = []
        for key, entries in index.items():
            ratio = SequenceMatcher(None, q, key).ratio()
            if ratio >= 0.6:
                for entry in entries:
                    if expected_type and entry["type"] != expected_type:
                        continue
                    best_fuzzy.append({**entry, "confidence": round(ratio * 0.9, 3)})
        best_fuzzy.sort(key=lambda x: x["confidence"], reverse=True)
        candidates.extend(best_fuzzy[:top_k])

    # ── Deduplicate, prioritize metabolite type, and rank ───────────
    # In CoreMet, metabolites are central. When an entity appears as
    # both metabolite and drug (e.g., Butyric acid), prefer metabolite.
    _TYPE_PRIORITY = {
        "metabolite": 0, "disease": 1, "gene": 2, "protein": 3,
        "microbe": 4, "drug": 5, "snp": 6,
    }

    def _sort_key(c):
        type_prio = _TYPE_PRIORITY.get(c["type"], 9)
        # Prefer entries with non-empty IDs
        id_penalty = 0 if (c.get("id") and c["id"] != "nan") else 1
        return (-c["confidence"], id_penalty, type_prio)

    # Try to fill missing IDs from other candidates of the same type
    id_map = {}  # (name_lower, type) -> id
    for c in candidates:
        if c.get("id") and c["id"] != "nan" and c["id"] != "":
            nkey = (c["type"],)
            if nkey not in id_map:
                id_map[nkey] = c["id"]
    for c in candidates:
        if not c.get("id") or c["id"] == "nan" or c["id"] == "":
            nkey = (c["type"],)
            if nkey in id_map:
                c["id"] = id_map[nkey]

    seen = set()
    unique = []
    for c in sorted(candidates, key=_sort_key):
        # Deduplicate by canonical name (across types, keep best type)
        name_key = c["name"].lower()
        type_key = (c["name"], c["type"])
        if type_key not in seen:
            seen.add(type_key)
            unique.append(c)
    # Re-sort: metabolite first when confidence is close
    unique.sort(key=_sort_key)
    return unique[:top_k]


def resolve_pair(
    source_text: str,
    target_text: str,
    source_type: Optional[str] = None,
    target_type: Optional[str] = None,
) -> Tuple[List[dict], List[dict]]:
    """Resolve a source-target entity pair for AI queries."""
    src = resolve_entity(source_text, expected_type=source_type)
    tgt = resolve_entity(target_text, expected_type=target_type)
    return src, tgt
