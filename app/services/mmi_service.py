"""
CoreMet — Metabolite–Microbe Interaction (MMI) Service.

Provides MMI database access, search, and annotation functions.
Data source: gutMGene v2.0 (experimentally validated microbe-metabolite interactions).

Schema: Metabolite_Name, HMDB_ID, KEGG_ID, PubChem_CID, ChEBI_ID, SMILES,
        Microbe_Name, Taxonomy_ID, Rank, Substrate, Substrate_PubChem_CID,
        Relationship_Type, Tissue, Organism, Evidence_Level,
        Experimental_Method, PMID, Source
"""

import logging
import threading
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
_rel = _DATA_DIR / "databases" / "release" / "coremetdb_mmi.csv"   # deduplicated, paper-matching
_v3 = _DATA_DIR / "databases" / "mmi_database_v3.csv"
MMI_DB_PATH = _rel if _rel.exists() else (_v3 if _v3.exists() else _DATA_DIR / "databases" / "mmi_database.csv")

_mmi_cache = None
_mmi_lock = threading.Lock()


def _invalidate_cache():
    """Clear the cached MMI database (useful after rebuilds)."""
    global _mmi_cache
    _mmi_cache = None


def get_mmi_db() -> pd.DataFrame:
    """Load and cache the MMI database (thread-safe). Returns empty DataFrame if not available."""
    global _mmi_cache
    if _mmi_cache is not None:
        return _mmi_cache
    with _mmi_lock:
        if _mmi_cache is not None:
            return _mmi_cache
        try:
            df = pd.read_csv(MMI_DB_PATH, dtype=str).fillna("")
            logger.info(f"MMI database loaded: {len(df):,} records")
            _mmi_cache = df
        except FileNotFoundError:
            logger.warning(f"MMI database not found at {MMI_DB_PATH}")
            _mmi_cache = pd.DataFrame()
    return _mmi_cache


def get_mmi_stats() -> dict:
    """Return summary statistics for the MMI database."""
    df = get_mmi_db()
    if df.empty:
        return {"total": 0, "metabolites": 0, "microbes": 0,
                "relationship_types": 0, "available": False}
    return {
        "total": len(df),
        "metabolites": df["Metabolite_Name"].nunique() if "Metabolite_Name" in df.columns else 0,
        "microbes": df["Microbe_Name"].nunique() if "Microbe_Name" in df.columns else 0,
        "relationship_types": df["Relationship_Type"].nunique() if "Relationship_Type" in df.columns else 0,
        "organisms": df["Organism"].value_counts().to_dict() if "Organism" in df.columns else {},
        "tissues": df["Tissue"].nunique() if "Tissue" in df.columns else 0,
        "available": True,
    }


def search_mmi(query: str = "", filters: dict | None = None) -> pd.DataFrame:
    """Search the MMI database by text query and/or column filters."""
    df = get_mmi_db()
    if df.empty:
        return df

    if query:
        q = query.lower()
        mask = pd.Series(False, index=df.index)
        for col in ["Metabolite_Name", "HMDB_ID", "KEGG_ID", "PubChem_CID",
                     "ChEBI_ID", "Microbe_Name", "Taxonomy_ID"]:
            if col in df.columns:
                mask |= df[col].astype(str).str.lower().str.contains(q, na=False)
        df = df[mask]

    if filters:
        for col, vals in filters.items():
            if col in df.columns and vals:
                if isinstance(vals, list):
                    df = df[df[col].isin(vals)]
                else:
                    df = df[df[col] == vals]

    return df


def get_microbes_for_metabolite(metabolite: str) -> pd.DataFrame:
    """Given a metabolite name or HMDB ID, return associated microbes."""
    df = get_mmi_db()
    if df.empty:
        return df
    q = metabolite.lower()
    mask = pd.Series(False, index=df.index)
    for col in ["Metabolite_Name", "HMDB_ID", "KEGG_ID", "PubChem_CID"]:
        if col in df.columns:
            mask |= df[col].astype(str).str.lower().str.contains(q, na=False)
    return df[mask]


def get_metabolites_for_microbe(microbe: str) -> pd.DataFrame:
    """Given a microbe name or taxonomy ID, return associated metabolites."""
    df = get_mmi_db()
    if df.empty:
        return df
    q = microbe.lower()
    mask = pd.Series(False, index=df.index)
    for col in ["Microbe_Name", "Taxonomy_ID"]:
        if col in df.columns:
            mask |= df[col].astype(str).str.lower().str.contains(q, na=False)
    return df[mask]


def annotate_metabolites_with_microbes(metabolite_names: list[str]) -> pd.DataFrame:
    """Annotate a list of metabolites with their known microbe associations."""
    df = get_mmi_db()
    if df.empty or not metabolite_names:
        return pd.DataFrame()

    names_lower = {n.lower() for n in metabolite_names if n}
    mask = df["Metabolite_Name"].astype(str).str.lower().isin(names_lower)
    return df[mask].copy()


def get_unique_microbes() -> list[str]:
    """Return sorted list of unique microbe names."""
    df = get_mmi_db()
    if df.empty or "Microbe_Name" not in df.columns:
        return []
    return sorted(df["Microbe_Name"].dropna().unique().tolist())


def get_unique_metabolites() -> list[str]:
    """Return sorted list of unique metabolite names."""
    df = get_mmi_db()
    if df.empty or "Metabolite_Name" not in df.columns:
        return []
    return sorted(df["Metabolite_Name"].dropna().unique().tolist())
