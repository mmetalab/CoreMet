"""
CoreMet — Metabolite–Drug Interaction (MDrI) Service.

Provides MDrI database access, search, and annotation functions.
Data sources: HMDB drug associations, CTD, CoreMet curated literature.

Schema: Metabolite_Name, HMDB_ID, SMILES, Drug_Name, DrugBank_ID,
        Interaction_Type, Tissue, Cell_Location, Biospecimen,
        Evidence_Level, Description, PMID, Source
"""

import logging
import threading
from pathlib import Path
import pandas as pd

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).resolve().parent.parent.parent / "data"
_rel = _DATA_DIR / "databases" / "release" / "coremetdb_mdri.csv"   # deduplicated, paper-matching
_v3 = _DATA_DIR / "databases" / "mdri_database_v3.csv"
MDRI_DB_PATH = _rel if _rel.exists() else (_v3 if _v3.exists() else _DATA_DIR / "databases" / "mdri_database.csv")

_mdri_cache = None
_mdri_lock = threading.Lock()


def _invalidate_cache():
    """Clear the cached MDrI database (useful after rebuilds)."""
    global _mdri_cache
    _mdri_cache = None


def get_mdri_db() -> pd.DataFrame:
    """Load and cache the MDrI database (thread-safe). Returns empty DataFrame if not available."""
    global _mdri_cache
    if _mdri_cache is not None:
        return _mdri_cache
    with _mdri_lock:
        if _mdri_cache is not None:
            return _mdri_cache
        try:
            df = pd.read_csv(MDRI_DB_PATH, dtype=str).fillna("")
            logger.info(f"MDrI database loaded: {len(df):,} records")
            _mdri_cache = df
        except FileNotFoundError:
            logger.warning(f"MDrI database not found at {MDRI_DB_PATH}")
            _mdri_cache = pd.DataFrame()
    return _mdri_cache


def get_mdri_stats() -> dict:
    """Return summary statistics for the MDrI database."""
    df = get_mdri_db()
    if df.empty:
        return {"total": 0, "metabolites": 0, "drugs": 0,
                "interaction_types": 0, "available": False}
    return {
        "total": len(df),
        "metabolites": df["Metabolite_Name"].nunique() if "Metabolite_Name" in df.columns else 0,
        "drugs": df["Drug_Name"].nunique() if "Drug_Name" in df.columns else 0,
        "interaction_types": df["Interaction_Type"].nunique() if "Interaction_Type" in df.columns else 0,
        "tissues": df["Tissue"].nunique() if "Tissue" in df.columns else 0,
        "evidence_levels": df["Evidence_Level"].value_counts().to_dict() if "Evidence_Level" in df.columns else {},
        "sources": df["Source"].value_counts().to_dict() if "Source" in df.columns else {},
        "available": True,
    }


def search_mdri(query: str = "", filters: dict | None = None) -> pd.DataFrame:
    """Search the MDrI database by text query and/or column filters."""
    df = get_mdri_db()
    if df.empty:
        return df

    if query:
        q = query.lower()
        mask = pd.Series(False, index=df.index)
        for col in ["Metabolite_Name", "HMDB_ID", "Drug_Name", "DrugBank_ID",
                     "Tissue", "Interaction_Type"]:
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

    return df.reset_index(drop=True)
