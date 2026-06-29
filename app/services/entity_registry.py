"""
CoreMet Entity ID Registry.

Assigns and manages unique CoreMet IDs for all entities across the database.

ID Format:  CM{type_code}{6-digit zero-padded number}
  CMM000001  — Metabolite
  CMG000001  — Gene
  CMP000001  — Protein
  CMD000001  — Disease
  CMB000001  — Microbe (Bacteria)
  CMR000001  — Drug (Rx/pharmaceutical)
  CMS000001  — SNP

The registry is built at startup from all loaded databases and cached in memory.
A JSON export is also saved to data/coremetdb_entity_registry.json for persistence.
"""

import json
import os
import logging
from pathlib import Path
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# ── ID Prefixes ──────────────────────────────────────────────────────────
TYPE_PREFIX = {
    "metabolite": "CMM",
    "gene":       "CMG",
    "protein":    "CMP",
    "disease":    "CMD",
    "microbe":    "CMB",
    "drug":       "CMR",
    "snp":        "CMS",
}

PREFIX_TYPE = {v: k for k, v in TYPE_PREFIX.items()}

_REGISTRY_FILE = Path(__file__).parent.parent.parent / "data" / "coremetdb_entity_registry.json"

# In-memory registry:  {entity_type: {canonical_name_lower: {id, name, external_ids}}}
_registry: dict = {}
# Reverse index:  {coremetdb_id: {type, name, external_ids}}
_reverse: dict = {}
# Name alias index: {alias_lower: coremetdb_id}
_aliases: dict = {}
_built = False


def _make_id(prefix: str, num: int) -> str:
    return f"{prefix}{num:06d}"


def build_registry():
    """Scan all databases and assign CoreMet IDs to every unique entity."""
    global _registry, _reverse, _aliases, _built

    if _built:
        return

    # Try loading from disk first
    if _REGISTRY_FILE.exists():
        try:
            data = json.loads(_REGISTRY_FILE.read_text())
            _registry = data.get("registry", {})
            _reverse = data.get("reverse", {})
            _aliases = data.get("aliases", {})
            _built = True
            logger.info("CoreMet ID registry loaded from %s (%d entities)",
                        _REGISTRY_FILE, len(_reverse))
            return
        except Exception as e:
            logger.warning("Failed to load registry from disk: %s", e)

    logger.info("Building CoreMet ID registry from databases...")

    _registry = {t: {} for t in TYPE_PREFIX}
    _reverse = {}
    _aliases = {}
    counters = {t: 0 for t in TYPE_PREFIX}

    def _register(entity_type, name, external_ids=None):
        """Register an entity and return its CoreMet ID."""
        if not name or not str(name).strip():
            return None
        name = str(name).strip()
        key = name.lower()
        if key in _registry[entity_type]:
            # Update external_ids if new ones provided
            if external_ids:
                existing = _registry[entity_type][key]
                for k, v in external_ids.items():
                    if v and (k not in existing["external_ids"] or not existing["external_ids"][k]):
                        existing["external_ids"][k] = v
            return _registry[entity_type][key]["id"]

        counters[entity_type] += 1
        cid = _make_id(TYPE_PREFIX[entity_type], counters[entity_type])
        entry = {
            "id": cid,
            "name": name,
            "external_ids": external_ids or {},
        }
        _registry[entity_type][key] = entry
        _reverse[cid] = {"type": entity_type, "name": name, "external_ids": external_ids or {}}

        # Register aliases
        _aliases[key] = cid
        if external_ids:
            for v in external_ids.values():
                if v and str(v).strip():
                    _aliases[str(v).strip().lower()] = cid

        return cid

    # ── Load from databases (vectorized for speed) ────────────────────
    def _bulk_register(entity_type, names, ext_id_col=None, ext_ids=None):
        """Register multiple entities at once using vectorized approach."""
        for i, name in enumerate(names):
            name = str(name).strip()
            if not name:
                continue
            ext = {}
            if ext_id_col and ext_ids is not None and i < len(ext_ids):
                val = str(ext_ids[i]).strip()
                if val:
                    ext[ext_id_col] = val
            _register(entity_type, name, ext if ext else None)

    try:
        from app.config import Config
        cfg = Config()

        # MPI: metabolites, proteins, genes
        mpi = pd.read_csv(cfg.MPI_DB_PATH, dtype=str).fillna("")

        # Metabolites from MPI
        met_df = mpi[["HMDB ID", "Metabolite Name"]].drop_duplicates(subset=["HMDB ID"])
        met_df = met_df[met_df["HMDB ID"].str.startswith("HMDB", na=False)]
        _bulk_register("metabolite", met_df["Metabolite Name"].values,
                       "hmdb_id", met_df["HMDB ID"].values)

        # Proteins from MPI
        prot_df = mpi[["Uniprot ID", "Protein Name"]].drop_duplicates(subset=["Uniprot ID"])
        prot_df = prot_df[prot_df["Uniprot ID"].str.len() > 0]
        _bulk_register("protein", prot_df["Protein Name"].values,
                       "uniprot_id", prot_df["Uniprot ID"].values)

        # Genes from MPI
        gene_names = mpi["Gene Name"].dropna().unique()
        gene_names = [g for g in gene_names if str(g).strip()]
        _bulk_register("gene", gene_names)

    except Exception as e:
        logger.warning("Failed to index MPI: %s", e)

    # MDI: metabolites, diseases
    try:
        from app.services.mdi_service import get_mdi_db
        mdi = get_mdi_db()
        if mdi is not None:
            mdi_met = mdi[["HMDB_ID", "Metabolite_Name"]].drop_duplicates(subset=["HMDB_ID"])
            mdi_met = mdi_met[mdi_met["HMDB_ID"].str.startswith("HMDB", na=False)]
            _bulk_register("metabolite", mdi_met["Metabolite_Name"].values,
                           "hmdb_id", mdi_met["HMDB_ID"].values)

            dis_df = mdi[["Disease_Name"]].copy()
            if "Disease_ID" in mdi.columns:
                dis_df["Disease_ID"] = mdi["Disease_ID"]
            else:
                dis_df["Disease_ID"] = ""
            dis_df = dis_df.drop_duplicates(subset=["Disease_Name"])
            dis_df = dis_df[dis_df["Disease_Name"].str.len() > 0]
            _bulk_register("disease", dis_df["Disease_Name"].values,
                           "disease_id", dis_df["Disease_ID"].values)
    except Exception as e:
        logger.warning("Failed to index MDI: %s", e)

    # MMI: metabolites, microbes
    try:
        from app.services.mmi_service import get_mmi_db
        mmi = get_mmi_db()
        if mmi is not None:
            mmi_met = mmi["Metabolite_Name"].dropna().unique()
            _bulk_register("metabolite", [str(n).strip() for n in mmi_met if str(n).strip()])

            mic_df = mmi[["Microbe_Name"]].copy()
            if "Taxonomy_ID" in mmi.columns:
                mic_df["Taxonomy_ID"] = mmi["Taxonomy_ID"]
            else:
                mic_df["Taxonomy_ID"] = ""
            mic_df = mic_df.drop_duplicates(subset=["Microbe_Name"])
            mic_df = mic_df[mic_df["Microbe_Name"].str.len() > 0]
            _bulk_register("microbe", mic_df["Microbe_Name"].values,
                           "taxonomy_id", mic_df["Taxonomy_ID"].values)
    except Exception as e:
        logger.warning("Failed to index MMI: %s", e)

    # MDrI: drugs, metabolites
    try:
        from app.services.mdri_service import get_mdri_db
        mdri = get_mdri_db()
        if mdri is not None:
            drug_df = mdri[["Drug_Name"]].copy()
            if "DrugBank_ID" in mdri.columns:
                drug_df["DrugBank_ID"] = mdri["DrugBank_ID"]
            else:
                drug_df["DrugBank_ID"] = ""
            drug_df = drug_df.drop_duplicates(subset=["Drug_Name"])
            drug_df = drug_df[drug_df["Drug_Name"].str.len() > 0]
            _bulk_register("drug", drug_df["Drug_Name"].values,
                           "drugbank_id", drug_df["DrugBank_ID"].values)

            mdri_met = mdri[["HMDB_ID", "Metabolite_Name"]].drop_duplicates(subset=["HMDB_ID"])
            mdri_met = mdri_met[mdri_met["HMDB_ID"].str.startswith("HMDB", na=False)]
            _bulk_register("metabolite", mdri_met["Metabolite_Name"].values,
                           "hmdb_id", mdri_met["HMDB_ID"].values)
    except Exception as e:
        logger.warning("Failed to index MDrI: %s", e)

    # MGI: genes (vectorized — 1.6M rows)
    try:
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if mgi is not None:
            gene_syms = mgi["Gene_Symbol"].dropna().unique()
            gene_syms = [str(g).strip() for g in gene_syms if str(g).strip()]
            _bulk_register("gene", gene_syms)

            if "HMDB_ID" in mgi.columns:
                mgi_met = mgi[["HMDB_ID", "Metabolite_Name"]].drop_duplicates(subset=["HMDB_ID"])
                mgi_met = mgi_met[mgi_met["HMDB_ID"].str.startswith("HMDB", na=False)]
                _bulk_register("metabolite", mgi_met["Metabolite_Name"].values,
                               "hmdb_id", mgi_met["HMDB_ID"].values)
    except Exception as e:
        logger.warning("Failed to index MGI: %s", e)

    # mGWAS: SNPs, genes
    try:
        from app.services.mgwas_service import get_mgwas_db
        mgwas = get_mgwas_db()
        if mgwas is not None:
            snp_ids = mgwas["rsID"].dropna().unique()
            snp_ids = [str(s).strip() for s in snp_ids if str(s).strip()]
            _bulk_register("snp", snp_ids, "dbsnp_id", snp_ids)

            if "Mapped_Gene" in mgwas.columns:
                mgenes = mgwas["Mapped_Gene"].dropna().unique()
                mgenes = [str(g).strip() for g in mgenes if str(g).strip()]
                _bulk_register("gene", mgenes)
    except Exception as e:
        logger.warning("Failed to index mGWAS: %s", e)

    _built = True
    logger.info("CoreMet ID registry built: %s",
                ", ".join(f"{t}={len(d)}" for t, d in _registry.items()))

    # Save to disk
    try:
        _REGISTRY_FILE.parent.mkdir(parents=True, exist_ok=True)
        data = {"registry": _registry, "reverse": _reverse, "aliases": _aliases}
        _REGISTRY_FILE.write_text(json.dumps(data, indent=1))
        logger.info("Registry saved to %s", _REGISTRY_FILE)
    except Exception as e:
        logger.warning("Failed to save registry: %s", e)


def lookup_id(name: str, entity_type: str = None) -> Optional[str]:
    """Look up the CoreMet ID for an entity by name or external ID."""
    if not _built:
        build_registry()

    key = str(name).strip().lower()

    # Direct alias match
    if key in _aliases:
        cid = _aliases[key]
        if entity_type:
            if _reverse.get(cid, {}).get("type") == entity_type:
                return cid
        else:
            return cid

    # Type-specific registry lookup
    if entity_type and entity_type in _registry:
        entry = _registry[entity_type].get(key)
        if entry:
            return entry["id"]

    # Search all types
    for t in TYPE_PREFIX:
        entry = _registry.get(t, {}).get(key)
        if entry:
            return entry["id"]

    return None


def get_entity(coremetdb_id: str) -> Optional[dict]:
    """Get entity info by CoreMet ID."""
    if not _built:
        build_registry()
    return _reverse.get(coremetdb_id)


def get_entity_by_name(name: str, entity_type: str = None) -> Optional[dict]:
    """Get full entity info by name."""
    cid = lookup_id(name, entity_type)
    if cid:
        info = _reverse.get(cid, {}).copy()
        info["coremetdb_id"] = cid
        return info
    return None


def get_type_stats() -> dict:
    """Return count of entities per type."""
    if not _built:
        build_registry()
    return {t: len(d) for t, d in _registry.items()}


def format_id(coremetdb_id: str) -> str:
    """Format a CoreMet ID for display."""
    return coremetdb_id  # Already formatted as CMM000001
