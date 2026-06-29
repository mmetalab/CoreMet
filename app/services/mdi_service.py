"""
Metabolite–Disease Interaction (MDI) service.

Loads the curated MDI database and provides query/search functionality.
"""

import logging
import threading
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)

_DB_DIR = Path(__file__).parent.parent.parent / "data" / "databases"
# Prefer v3 (expanded with HMDB diseases, backward-compatible columns)
_rel = _DB_DIR / "release" / "coremetdb_mdi.csv"   # deduplicated, paper-matching
_v3 = _DB_DIR / "mdi_database_v3.csv"
MDI_DB_PATH = _rel if _rel.exists() else (_v3 if _v3.exists() else _DB_DIR / "mdi_database.csv")

_mdi_db: Optional[pd.DataFrame] = None
_mdi_lock = threading.Lock()


def get_mdi_db() -> pd.DataFrame:
    """Load and return the MDI database (cached on first call, thread-safe)."""
    global _mdi_db
    if _mdi_db is not None:
        return _mdi_db
    with _mdi_lock:
        if _mdi_db is not None:
            return _mdi_db
        try:
            _mdi_db = pd.read_csv(MDI_DB_PATH, dtype=str).fillna("")
            logger.info(f"MDI database loaded: {len(_mdi_db):,} records")
        except FileNotFoundError:
            logger.warning(f"MDI database not found at {MDI_DB_PATH}")
            _mdi_db = pd.DataFrame()
    return _mdi_db


def get_mdi_stats() -> dict:
    """Return summary statistics for the MDI database."""
    df = get_mdi_db()
    if df.empty:
        return {"total": 0, "metabolites": 0, "diseases": 0, "categories": 0}
    return {
        "total": len(df),
        "metabolites": df["HMDB_ID"].nunique(),
        "diseases": df["Disease_Name"].nunique(),
        "categories": df["Category"].nunique(),
    }


def search_mdi(
    metabolite: str = "",
    disease: str = "",
    category: str = "",
    association_type: str = "",
    evidence_level: str = "",
    limit: int = 500,
) -> pd.DataFrame:
    """Search the MDI database with filters."""
    df = get_mdi_db()
    if df.empty:
        return df

    mask = pd.Series(True, index=df.index)

    if metabolite:
        q = metabolite.lower()
        met_mask = pd.Series(False, index=df.index)
        for col in ["Metabolite_Name", "HMDB_ID"]:
            if col in df.columns:
                met_mask |= df[col].str.lower().str.contains(q, na=False)
        mask &= met_mask

    if disease:
        q = disease.lower()
        dis_mask = pd.Series(False, index=df.index)
        for col in ["Disease_Name", "Disease_ID", "MeSH_ID"]:
            if col in df.columns:
                dis_mask |= df[col].str.lower().str.contains(q, na=False)
        mask &= dis_mask

    if category:
        mask &= df["Category"].str.lower() == category.lower()

    if association_type:
        mask &= df["Association_Type"].str.lower() == association_type.lower()

    if evidence_level:
        mask &= df["Evidence_Level"].str.lower() == evidence_level.lower()

    return df.loc[mask].head(limit)


def get_metabolites_for_disease(disease_name: str) -> pd.DataFrame:
    """Get all metabolites associated with a specific disease."""
    df = get_mdi_db()
    if df.empty:
        return df
    mask = df["Disease_Name"].str.lower() == disease_name.lower()
    return df.loc[mask]


def get_diseases_for_metabolite(hmdb_id: str = "", metabolite_name: str = "") -> pd.DataFrame:
    """Get all diseases associated with a specific metabolite."""
    df = get_mdi_db()
    if df.empty:
        return df
    mask = pd.Series(False, index=df.index)
    if hmdb_id:
        mask |= df["HMDB_ID"] == hmdb_id
    if metabolite_name:
        mask |= df["Metabolite_Name"].str.lower() == metabolite_name.lower()
    return df.loc[mask]


def annotate_metabolites_with_diseases(hmdb_ids: List[str]) -> Dict[str, dict]:
    """
    Batch lookup: given a list of HMDB IDs, return disease annotations.

    Returns dict mapping HMDB_ID → {diseases: [...], categories: [...], count: N}
    """
    df = get_mdi_db()
    if df.empty:
        return {}

    result = {}
    subset = df[df["HMDB_ID"].isin(hmdb_ids)]
    for hmdb_id, group in subset.groupby("HMDB_ID"):
        diseases = sorted(group["Disease_Name"].unique().tolist())
        categories = sorted(group["Category"].unique().tolist())
        result[hmdb_id] = {
            "diseases": diseases,
            "categories": categories,
            "count": len(diseases),
            "disease_list": "; ".join(diseases),
        }
    return result
