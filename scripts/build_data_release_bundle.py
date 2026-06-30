#!/usr/bin/env python3
"""
build_data_release_bundle.py — Package the CoreMet deposit for Zenodo/Figshare.

Bundles the 7 deduplicated release CSVs + stats JSON + DATA_README into a single
zip whose contents match the manuscript exactly. Upload this zip to Zenodo/Figshare
to obtain a citable DOI for the data-availability statement.

    conda run -n mpi-vgae python scripts/build_data_release_bundle.py
"""
from __future__ import annotations
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REL = ROOT / "data" / "databases" / "release"
STATS = ROOT / "data" / "coremetdb_stats.json"
README = ROOT / "DATA_README.md"
OUT_DIR = ROOT.parent / "CoreMetDB"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT = OUT_DIR / "coremet_dataset_v1.zip"

FILES = [REL / f"coremetdb_{k}.csv" for k in
         ("mpi", "mei", "mdi", "mmi", "mdri", "mgi", "mgwas")] + [STATS, README]


def main() -> None:
    missing = [f for f in FILES if not f.exists()]
    if missing:
        raise FileNotFoundError(f"missing: {missing}")
    total_rows = 0
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in FILES:
            z.write(f, arcname=f"CoreMet-DB_v1/{f.name}")
            if f.suffix == ".csv":
                with open(f) as fh:
                    total_rows += sum(1 for _ in fh) - 1
    size_mb = OUT.stat().st_size / 1e6
    print(f"Bundled {len(FILES)} files, {total_rows:,} interaction rows -> {OUT} ({size_mb:.0f} MB)")
    print("Upload this zip to Zenodo or Figshare to mint a DOI, then fill [URL]/DOI in the manuscript.")


if __name__ == "__main__":
    main()
