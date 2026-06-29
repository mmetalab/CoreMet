#!/usr/bin/env python3
"""Compute metabolite cross-reference coverage; append to coremetdb_stats.json."""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
STATS = DATA / "coremetdb_stats.json"
CHUNK = 200_000

# (file, hmdb_col, {xref_name: col})
SOURCES = [
    (DATA / "mpidatabase" / "MPIDB_v3.csv",            "HMDB ID", {"SMILES": "SMILES"}),
    (DATA / "databases" / "mei_database_v2_enriched.csv", "HMDB_ID", {"SMILES": "SMILES", "KEGG": "KEGG_Compound"}),
    (DATA / "databases" / "mdi_database_v3.csv",       "HMDB_ID", {"SMILES": "SMILES"}),
    (DATA / "databases" / "mmi_database_v3.csv",       "HMDB_ID", {"SMILES": "SMILES", "KEGG": "KEGG_ID",
                                                                    "PubChem": "PubChem_CID", "ChEBI": "ChEBI_ID"}),
    (DATA / "databases" / "mdri_database_v3.csv",      "HMDB_ID", {"SMILES": "SMILES"}),
    (DATA / "databases" / "mgwas_database_v3.csv",     "HMDB_ID", {"SMILES": "SMILES"}),
    (DATA / "databases" / "mgi_database_v3.csv",       "HMDB_ID", {}),
]
NULL = {"", "nan", "None", "-", "NA", "n/a"}


def main() -> None:
    all_mets: set[str] = set()
    has: dict[str, set[str]] = {k: set() for k in ["SMILES", "KEGG", "PubChem", "ChEBI"]}

    for path, hcol, xrefs in SOURCES:
        cols = [hcol] + list(xrefs.values())
        for ch in pd.read_csv(path, dtype=str, usecols=lambda c: c in cols,
                              chunksize=CHUNK, on_bad_lines="skip", low_memory=False):
            h = ch[hcol].astype(str).str.strip()
            valid = ~h.isin(NULL) & h.notna()
            all_mets.update(h[valid])
            for name, col in xrefs.items():
                if col in ch.columns:
                    v = ch[col].astype(str).str.strip()
                    ok = valid & v.notna() & ~v.isin(NULL)
                    has[name].update(h[ok])

    n = len(all_mets)
    cov = {"unique_metabolites": n,
           "HMDB": {"count": n, "pct": 100.0}}
    for name, s in has.items():
        c = len(s & all_mets)
        cov[name] = {"count": c, "pct": round(100 * c / n, 1) if n else 0.0}

    stats = json.loads(STATS.read_text())
    stats["cross_references"] = cov
    STATS.write_text(json.dumps(stats, indent=2))
    print(json.dumps(cov, indent=2))
    print("Updated", STATS)


if __name__ == "__main__":
    main()
