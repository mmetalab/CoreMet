"""
Metabolite–Gene Interaction (MGI) service.

Loads the curated MGI database and provides query/search functionality.
Data source: CTD (Comparative Toxicogenomics Database) chemical–gene interactions.
"""

import logging
import threading
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)

_DB_DIR = Path(__file__).parent.parent.parent / "data" / "databases"
_rel = _DB_DIR / "release" / "coremetdb_mgi.csv"   # deduplicated, paper-matching
_v3 = _DB_DIR / "mgi_database_v3.csv"
MGI_DB_PATH = _rel if _rel.exists() else (_v3 if _v3.exists() else _DB_DIR / "mgi_database.csv")

_mgi_db: Optional[pd.DataFrame] = None
_mgi_lock = threading.Lock()


def get_mgi_db() -> pd.DataFrame:
    """Load and return the MGI database (cached on first call, thread-safe)."""
    global _mgi_db
    if _mgi_db is not None:
        return _mgi_db
    with _mgi_lock:
        if _mgi_db is not None:
            return _mgi_db
        try:
            from app.services.csv_loader import load_optimized
            _mgi_db = load_optimized(MGI_DB_PATH)
            logger.info(f"MGI database loaded: {len(_mgi_db):,} records")
        except FileNotFoundError:
            logger.warning(f"MGI database not found at {MGI_DB_PATH}")
            _mgi_db = pd.DataFrame()
    return _mgi_db


def get_mgi_stats() -> dict:
    """Return summary statistics for the MGI database."""
    df = get_mgi_db()
    if df.empty:
        return {"total": 0, "metabolites": 0, "genes": 0, "organisms": 0,
                "interaction_types": 0, "available": False}
    return {
        "total": len(df),
        "metabolites": df["HMDB_ID"].nunique(),
        "genes": df["Gene_Symbol"].nunique(),
        "organisms": df["Organism"].nunique(),
        "interaction_types": df["Interaction_Type"].nunique(),
        "sources": df["Source"].value_counts().to_dict() if "Source" in df.columns else {},
        "available": True,
    }


def search_mgi(
    metabolite: str = "",
    gene: str = "",
    organism: str = "",
    interaction_type: str = "",
    source: str = "",
    limit: int = 500,
) -> pd.DataFrame:
    """Search the MGI database with filters."""
    df = get_mgi_db()
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

    if gene:
        q = gene.lower()
        gene_mask = pd.Series(False, index=df.index)
        for col in ["Gene_Symbol", "Gene_ID"]:
            if col in df.columns:
                gene_mask |= df[col].str.lower().str.contains(q, na=False)
        mask &= gene_mask

    if organism:
        mask &= df["Organism"].str.lower() == organism.lower()

    if interaction_type:
        mask &= df["Interaction_Type"].str.lower().str.contains(interaction_type.lower(), na=False)

    if source:
        mask &= df["Source"].str.lower() == source.lower()

    return df.loc[mask].head(limit)


def get_genes_for_metabolite(hmdb_id: str = "", metabolite_name: str = "") -> pd.DataFrame:
    """Get all genes associated with a specific metabolite."""
    df = get_mgi_db()
    if df.empty:
        return df
    mask = pd.Series(False, index=df.index)
    if hmdb_id:
        mask |= df["HMDB_ID"] == hmdb_id
    if metabolite_name:
        mask |= df["Metabolite_Name"].str.lower() == metabolite_name.lower()
    return df.loc[mask]


def get_metabolites_for_gene(gene_symbol: str) -> pd.DataFrame:
    """Get all metabolites associated with a specific gene."""
    df = get_mgi_db()
    if df.empty:
        return df
    mask = df["Gene_Symbol"].str.lower() == gene_symbol.lower()
    return df.loc[mask]


def annotate_metabolites_with_genes(hmdb_ids: List[str]) -> Dict[str, dict]:
    """
    Batch lookup: given a list of HMDB IDs, return gene annotations.

    Returns dict mapping HMDB_ID → {genes: [...], organisms: [...], count: N}
    """
    df = get_mgi_db()
    if df.empty:
        return {}

    result = {}
    subset = df[df["HMDB_ID"].isin(hmdb_ids)]
    for hmdb_id, group in subset.groupby("HMDB_ID"):
        genes = sorted(group["Gene_Symbol"].unique().tolist())
        organisms = sorted(group["Organism"].unique().tolist())
        result[hmdb_id] = {
            "genes": genes[:50],  # Cap for memory
            "organisms": organisms,
            "count": len(genes),
            "gene_list": "; ".join(genes[:10]),  # Cap at 10 for display
        }
    return result
