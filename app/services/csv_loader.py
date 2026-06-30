"""
Memory-efficient CSV loader for CoreMet database tables.

The interaction tables (especially MGI, 1.66M rows) are dominated by repeated
strings (HMDB IDs, SMILES, organism, interaction type, source). Loading them as
plain object columns costs ~1.5 GB for MGI alone; converting low-cardinality
columns to pandas `category` dtype cuts that by ~95% with no change to read,
filter, search, or display behavior. This keeps the whole app well under the
2 GB Render instance limit.
"""
from __future__ import annotations

import pandas as pd


def load_optimized(path, **read_csv_kwargs) -> pd.DataFrame:
    """Read a CoreMet interaction CSV with minimal memory.

    Every column is read directly as pandas `category` dtype, which keeps both the
    transient load spike and the resident footprint low (MGI: ~1.5 GB as object ->
    ~0.07 GB as category). Empty cells are read as "" (na_filter off), matching the
    prior `.fillna("")` behavior, so downstream string filters/searches are unchanged.
    """
    read_csv_kwargs.setdefault("dtype", "category")
    read_csv_kwargs.setdefault("keep_default_na", False)  # blanks -> "" (a valid category), no NaN
    read_csv_kwargs.setdefault("na_filter", False)
    read_csv_kwargs.pop("low_memory", None)               # irrelevant with an explicit dtype
    return pd.read_csv(path, **read_csv_kwargs)
