#!/usr/bin/env python3
"""
build_runtime_bundle.py — Assemble the CoreMet runtime data bundle for Zenodo.

Produces ../CoreMetDB/coremet_runtime_core.zip containing exactly the data the live
web app loads (no raw sources, no v3 duplicates). Archive paths are relative to data/
so scripts/fetch_data.py can extract straight into CoreMet/data/.

Core (database app: browse, search, profile, network, download, API):
  coremetdb_stats.json, coremetdb_use_cases.json, coremetdb_entity_registry.json,
  databases/release/coremetdb_*.csv, mpidatabase/MPIDB_v3.csv,
  module_summaries/*.json

    conda run -n mpi-vgae python scripts/build_runtime_bundle.py
"""
from __future__ import annotations
import zipfile
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT_DIR = DATA.parent.parent / "CoreMetDB"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "coremet_runtime_core.zip"

MEMBERS = [
    "coremetdb_stats.json",
    "coremetdb_use_cases.json",
    "coremetdb_entity_registry.json",
    "mpidatabase/MPIDB_v3.csv",
    "module_summaries/mgi_summary.json",
    "module_summaries/mmi_summary.json",
] + [f"databases/release/coremetdb_{k}.csv"
     for k in ("mpi", "mei", "mdi", "mmi", "mdri", "mgi", "mgwas")]


def main() -> None:
    missing = [m for m in MEMBERS if not (DATA / m).exists()]
    if missing:
        print("WARNING missing (skipped):", missing)
    total = 0
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for m in MEMBERS:
            p = DATA / m
            if p.exists():
                z.write(p, arcname=m)         # arcname relative to data/
                total += p.stat().st_size
    print(f"raw {total/1e6:.0f} MB -> {OUT} ({OUT.stat().st_size/1e6:.0f} MB)")
    print("Upload to Zenodo; set DATA_BUNDLE_URL to its direct download link.")


if __name__ == "__main__":
    main()
