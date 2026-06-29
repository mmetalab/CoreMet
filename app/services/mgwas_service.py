"""
Metabolite–SNP / mGWAS Interaction service.

Loads the curated mGWAS database and provides query/search functionality.
Data source: GWAS Catalog (EBI) + curated mGWAS publications.
"""

import logging
import threading
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)

_DB_DIR = Path(__file__).parent.parent.parent / "data" / "databases"
_rel = _DB_DIR / "release" / "coremetdb_mgwas.csv"   # deduplicated, paper-matching
_v3 = _DB_DIR / "mgwas_database_v3.csv"
MGWAS_DB_PATH = _rel if _rel.exists() else (_v3 if _v3.exists() else _DB_DIR / "mgwas_database.csv")

_mgwas_db: Optional[pd.DataFrame] = None
_mgwas_lock = threading.Lock()


def get_mgwas_db() -> pd.DataFrame:
    """Load and return the mGWAS database (cached on first call, thread-safe)."""
    global _mgwas_db
    if _mgwas_db is not None:
        return _mgwas_db
    with _mgwas_lock:
        if _mgwas_db is not None:
            return _mgwas_db
        try:
            from app.services.csv_loader import load_optimized
            _mgwas_db = load_optimized(MGWAS_DB_PATH)
            logger.info(f"mGWAS database loaded: {len(_mgwas_db):,} records")
        except FileNotFoundError:
            logger.warning(f"mGWAS database not found at {MGWAS_DB_PATH}")
            _mgwas_db = pd.DataFrame()
    return _mgwas_db


def get_mgwas_stats() -> dict:
    """Return summary statistics for the mGWAS database."""
    df = get_mgwas_db()
    if df.empty:
        return {"total": 0, "metabolites": 0, "snps": 0, "genes": 0,
                "chromosomes": 0, "available": False}
    return {
        "total": len(df),
        "metabolites": df["HMDB_ID"].nunique(),
        "snps": df["rsID"].nunique(),
        "genes": df["Mapped_Gene"].nunique(),
        "chromosomes": df["Chromosome"].nunique(),
        "sources": df["Source"].value_counts().to_dict() if "Source" in df.columns else {},
        "available": True,
    }


def search_mgwas(
    metabolite: str = "",
    snp: str = "",
    gene: str = "",
    chromosome: str = "",
    source: str = "",
    limit: int = 500,
) -> pd.DataFrame:
    """Search the mGWAS database with filters."""
    df = get_mgwas_db()
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

    if snp:
        q = snp.lower()
        mask &= df["rsID"].str.lower().str.contains(q, na=False)

    if gene:
        q = gene.lower()
        mask &= df["Mapped_Gene"].str.lower().str.contains(q, na=False)

    if chromosome:
        mask &= df["Chromosome"] == chromosome

    if source:
        mask &= df["Source"].str.lower() == source.lower()

    return df.loc[mask].head(limit)


def get_snps_for_metabolite(hmdb_id: str = "", metabolite_name: str = "") -> pd.DataFrame:
    """Get all SNPs associated with a specific metabolite."""
    df = get_mgwas_db()
    if df.empty:
        return df
    mask = pd.Series(False, index=df.index)
    if hmdb_id:
        mask |= df["HMDB_ID"] == hmdb_id
    if metabolite_name:
        mask |= df["Metabolite_Name"].str.lower() == metabolite_name.lower()
    return df.loc[mask]


def get_metabolites_for_snp(rsid: str) -> pd.DataFrame:
    """Get all metabolites associated with a specific SNP."""
    df = get_mgwas_db()
    if df.empty:
        return df
    mask = df["rsID"].str.lower() == rsid.lower()
    return df.loc[mask]


def annotate_metabolites_with_snps(hmdb_ids: List[str]) -> Dict[str, dict]:
    """Batch lookup: given HMDB IDs, return mGWAS annotations."""
    df = get_mgwas_db()
    if df.empty:
        return {}

    result = {}
    subset = df[df["HMDB_ID"].isin(hmdb_ids)]
    for hmdb_id, group in subset.groupby("HMDB_ID"):
        snps = sorted(group["rsID"].unique().tolist())
        genes = sorted(group["Mapped_Gene"].unique().tolist())
        result[hmdb_id] = {
            "snps": snps,
            "genes": genes,
            "count": len(snps),
            "snp_list": "; ".join(snps[:10]),  # Cap at 10 for display
        }
    return result
