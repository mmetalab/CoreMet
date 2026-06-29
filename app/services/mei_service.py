"""
Metabolite–Enzyme Interaction (MEI) service.

Loads the curated MEI database and provides query/search functionality.
"""

import logging
import threading
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List

logger = logging.getLogger(__name__)

_MEI_RELEASE = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mei.csv"
MEI_DB_PATH = _MEI_RELEASE if _MEI_RELEASE.exists() else (
    Path(__file__).parent.parent.parent / "data" / "databases" / "mei_database.csv")

_mei_db: Optional[pd.DataFrame] = None
_mei_lock = threading.Lock()
_uniprot_to_ec: Optional[Dict[str, dict]] = None


def get_mei_db() -> pd.DataFrame:
    """Load and return the MEI database (cached on first call, thread-safe)."""
    global _mei_db
    if _mei_db is not None:
        return _mei_db
    with _mei_lock:
        if _mei_db is not None:
            return _mei_db
        try:
            from app.services.csv_loader import load_optimized
            _mei_db = load_optimized(MEI_DB_PATH)
            logger.info(f"MEI database loaded: {len(_mei_db):,} records")
        except FileNotFoundError:
            logger.warning(f"MEI database not found at {MEI_DB_PATH}")
            _mei_db = pd.DataFrame()
    return _mei_db


def get_mei_stats() -> dict:
    """Return summary statistics for the MEI database."""
    df = get_mei_db()
    if df.empty:
        return {"total": 0, "metabolites": 0, "enzymes": 0, "organisms": 0, "ecs": 0}
    return {
        "total": len(df),
        "metabolites": df["Metabolite_Name"].nunique(),
        "enzymes": df["Uniprot_ID"].nunique(),
        "organisms": df["Species"].nunique(),
        "ecs": df["EC_Number"].nunique(),
    }


def search_mei(
    metabolite: str = "",
    enzyme: str = "",
    ec_number: str = "",
    organism: str = "",
    pathway: str = "",
    limit: int = 500,
) -> pd.DataFrame:
    """Search the MEI database with filters."""
    df = get_mei_db()
    if df.empty:
        return df

    mask = pd.Series(True, index=df.index)

    if metabolite:
        q = metabolite.lower()
        met_mask = pd.Series(False, index=df.index)
        for col in ["Metabolite_Name", "HMDB_ID", "KEGG_Compound"]:
            if col in df.columns:
                met_mask |= df[col].str.lower().str.contains(q, na=False)
        mask &= met_mask

    if enzyme:
        q = enzyme.lower()
        enz_mask = pd.Series(False, index=df.index)
        for col in ["Enzyme_Name", "Uniprot_ID", "Gene_Name"]:
            if col in df.columns:
                enz_mask |= df[col].str.lower().str.contains(q, na=False)
        mask &= enz_mask

    if ec_number:
        mask &= df["EC_Number"].str.contains(ec_number, na=False)

    if organism:
        mask &= df["Species"].str.lower() == organism.lower()

    if pathway:
        q = pathway.lower()
        pw_mask = pd.Series(False, index=df.index)
        for col in ["Pathway_ID", "Pathway_Name"]:
            if col in df.columns:
                pw_mask |= df[col].str.lower().str.contains(q, na=False)
        mask &= pw_mask

    return df.loc[mask].head(limit)


def _build_uniprot_ec_lookup() -> Dict[str, dict]:
    """Build a UniProt → EC/Enzyme lookup from the MEI database (cached)."""
    global _uniprot_to_ec
    if _uniprot_to_ec is not None:
        return _uniprot_to_ec
    
    df = get_mei_db()
    _uniprot_to_ec = {}
    if df.empty:
        return _uniprot_to_ec
    
    # Group by UniProt ID, collect unique EC numbers and enzyme names
    for uid, group in df.groupby("Uniprot_ID"):
        if not uid:
            continue
        ecs = sorted(group["EC_Number"].unique().tolist())
        names = group["Enzyme_Name"].unique().tolist()
        # Pick the most common enzyme name
        name = names[0] if names else ""
        _uniprot_to_ec[uid] = {
            "EC_Number": "; ".join([e for e in ecs if e]),
            "Enzyme_Name": name,
            "Is_Enzyme": True,
        }
    
    logger.info(f"MEI UniProt lookup built: {len(_uniprot_to_ec):,} enzymes")
    return _uniprot_to_ec


def annotate_predictions_with_mei(predictions_df: pd.DataFrame, uniprot_col: str = "Protein") -> pd.DataFrame:
    """
    Annotate prediction results with MEI enzyme info.
    
    Adds columns: Is_Enzyme, EC_Number, Enzyme_Name, Interaction_Type (MPI/MEI).
    Works on the predictions DataFrame returned from PredictionService.
    
    Args:
        predictions_df: DataFrame with prediction results (must have a UniProt ID column)
        uniprot_col: Name of the column containing UniProt IDs
    
    Returns:
        DataFrame with additional enzyme annotation columns
    """
    df = predictions_df.copy()
    lookup = _build_uniprot_ec_lookup()
    
    if not lookup or uniprot_col not in df.columns:
        df["Is_Enzyme"] = False
        df["EC_Number"] = ""
        df["Enzyme_Name"] = ""
        df["Interaction_Type"] = "MPI"
        return df
    
    # Vectorized lookup
    ec_numbers = []
    enzyme_names = []
    is_enzyme = []
    
    for uid in df[uniprot_col]:
        info = lookup.get(str(uid), None)
        if info:
            ec_numbers.append(info["EC_Number"])
            enzyme_names.append(info["Enzyme_Name"])
            is_enzyme.append(True)
        else:
            ec_numbers.append("")
            enzyme_names.append("")
            is_enzyme.append(False)
    
    df["Is_Enzyme"] = is_enzyme
    df["EC_Number"] = ec_numbers
    df["Enzyme_Name"] = enzyme_names
    df["Interaction_Type"] = df["Is_Enzyme"].map({True: "MEI", False: "MPI"})
    
    n_mei = sum(is_enzyme)
    logger.info(f"MEI annotation: {n_mei}/{len(df)} predictions are enzyme interactions")
    
    return df
