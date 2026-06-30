"""
Metabolite–Protein Interaction (MPI) service.

Loads the core MPI database for graph traversal.
"""

import logging
import threading
import pandas as pd
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

_mpi_db: Optional[pd.DataFrame] = None
_mpi_lock = threading.Lock()


def get_mpi_db() -> pd.DataFrame:
    """Load and return the MPI database (cached, thread-safe)."""
    global _mpi_db
    if _mpi_db is not None:
        return _mpi_db
    with _mpi_lock:
        if _mpi_db is not None:
            return _mpi_db
        try:
            from app.config import Config
            cfg = Config()
            path = cfg.MPI_DB_PATH
            if path and path.exists():
                from app.services.csv_loader import load_optimized
                df = load_optimized(path)
                # Standardize column names for graph traversal compatibility
                col_map = {
                    "HMDB ID": "HMDB_ID",
                    "Uniprot ID": "Uniprot_ID",
                    "Metabolite Name": "Metabolite_Name",
                    "Protein Name": "Protein_Name",
                    "Gene Name": "Gene_Name",
                }
                df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})
                # Keep string operations stable even when columns are categorical.
                df["Metabolite_Name"] = df["Metabolite_Name"].astype("string").fillna("")
                df["Protein_Name"] = df["Protein_Name"].astype("string").fillna("")
                df["Gene_Name"] = df["Gene_Name"].astype("string").fillna("")
                _mpi_db = df
                logger.info(f"MPI service: loaded {len(df):,} interactions")
                return _mpi_db
            else:
                logger.warning("MPI service: database file not found")
                _mpi_db = pd.DataFrame()
                return _mpi_db
        except Exception as e:
            logger.error(f"MPI service: failed to load: {e}")
            _mpi_db = pd.DataFrame()
            return _mpi_db
