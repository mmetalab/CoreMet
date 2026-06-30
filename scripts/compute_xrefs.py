#!/usr/bin/env python3
"""Compute release metabolite cross-reference coverage; append to coremetdb_stats.json."""
from __future__ import annotations
import json
import re
from pathlib import Path
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
STATS = DATA / "coremetdb_stats.json"
CHUNK = 200_000

# (layer, file, hmdb_col, {xref_name: col})
SOURCES = [
    ("MPI", DATA / "databases" / "release" / "coremetdb_mpi.csv", "HMDB ID", {"SMILES": "SMILES"}),
    ("MEI", DATA / "databases" / "release" / "coremetdb_mei.csv", "HMDB_ID", {"SMILES": "SMILES", "KEGG": "KEGG_Compound"}),
    ("MDI", DATA / "databases" / "release" / "coremetdb_mdi.csv", "HMDB_ID", {"SMILES": "SMILES"}),
    ("MMI", DATA / "databases" / "release" / "coremetdb_mmi.csv", "HMDB_ID", {"SMILES": "SMILES", "KEGG": "KEGG_ID",
                                                                                  "PubChem": "PubChem_CID", "ChEBI": "ChEBI_ID"}),
    ("MDrI", DATA / "databases" / "release" / "coremetdb_mdri.csv", "HMDB_ID", {"SMILES": "SMILES"}),
    ("mGWAS", DATA / "databases" / "release" / "coremetdb_mgwas.csv", "HMDB_ID", {"SMILES": "SMILES"}),
    ("MGI", DATA / "databases" / "release" / "coremetdb_mgi.csv", "HMDB_ID", {}),
]
NULL = {"", "nan", "None", "-", "NA", "n/a"}
HMDB_RE = re.compile(r"^HMDB\d+$")


def main() -> None:
    all_mets: set[str] = set()
    valid_hmdb: set[str] = set()
    has: dict[str, set[str]] = {k: set() for k in ["SMILES", "KEGG", "PubChem", "ChEBI"]}
    row_coverage: dict[str, dict[str, float | int]] = {}

    for layer, path, hcol, xrefs in SOURCES:
        cols = [hcol] + list(xrefs.values())
        total_rows = 0
        valid_rows = 0
        for ch in pd.read_csv(path, dtype=str, usecols=lambda c: c in cols,
                              chunksize=CHUNK, on_bad_lines="skip", low_memory=False):
            h = ch[hcol].astype(str).str.strip()
            present = ~h.isin(NULL) & h.notna()
            valid = present & h.map(lambda x: bool(HMDB_RE.match(x)))
            total_rows += len(ch)
            valid_rows += int(valid.sum())
            all_mets.update(h[present])
            valid_hmdb.update(h[valid])
            for name, col in xrefs.items():
                if col in ch.columns:
                    v = ch[col].astype(str).str.strip()
                    ok = valid & v.notna() & ~v.isin(NULL)
                    has[name].update(h[ok])
        row_coverage[layer] = {
            "valid_rows": valid_rows,
            "total_rows": total_rows,
            "pct": round(100 * valid_rows / total_rows, 1) if total_rows else 0.0,
        }

    n = len(all_mets)
    h_n = len(valid_hmdb)
    cov = {
        "note": (
            "HMDB identifiers are provided where available. Unique-metabolite HMDB coverage is near-complete, "
            "but row-level HMDB coverage varies by layer because KEGG/Rhea MPI and MEI records can retain "
            "source identifiers when no unambiguous HMDB mapping is available."
        ),
        "unique_metabolites": n,
        "unique_valid_hmdb_ids": h_n,
        "HMDB": {"count": h_n, "pct": round(100 * h_n / n, 2) if n else 0.0},
    }
    for name, s in has.items():
        c = len(s & valid_hmdb)
        cov[name] = {"count": c, "pct": round(100 * c / h_n, 1) if h_n else 0.0}
    cov["per_layer_hmdb_row_coverage"] = row_coverage

    stats = json.loads(STATS.read_text())
    stats["cross_references"] = cov
    STATS.write_text(json.dumps(stats, indent=2))
    print(json.dumps(cov, indent=2))
    print("Updated", STATS)


if __name__ == "__main__":
    main()
