#!/usr/bin/env python3
"""
fetch_data.py — Download the CoreMet runtime data at deploy time.

Large data files are NOT committed to git (see .gitignore); they are hosted on
Zenodo/Figshare and fetched here during the Render build. Idempotent: skips
download if the data already exists (e.g., on a persistent disk).

Set DATA_BUNDLE_URL to a direct-download archive (.zip or .tar.gz) containing the
runtime `data/` subset. For Zenodo: use the file's direct URL
(https://zenodo.org/records/<id>/files/<file>?download=1).

    DATA_BUNDLE_URL=https://zenodo.org/records/XXXX/files/coremetdb_runtime.zip?download=1 \
        python scripts/fetch_data.py
"""
from __future__ import annotations
import os
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
SENTINEL = DATA / "coremetdb_stats.json"   # present => data already in place
URL = os.environ.get("DATA_BUNDLE_URL", "").strip()


def main() -> None:
    # If the runtime data is already present (committed sentinel or persistent disk), skip.
    if (DATA / "databases" / "release" / "coremetdb_mgi.csv").exists():
        print("fetch_data: runtime data already present — skipping download.")
        return
    if not URL:
        print("fetch_data: DATA_BUNDLE_URL not set and data absent. "
              "Set it to the Zenodo/Figshare archive URL, or attach a persistent disk "
              "with the data/ directory. See DEPLOY_RENDER.md.", file=sys.stderr)
        # Do not hard-fail the build; the app can still start its lighter pages.
        return

    DATA.mkdir(parents=True, exist_ok=True)
    archive = DATA / "_bundle"
    print(f"fetch_data: downloading {URL} ...")
    urllib.request.urlretrieve(URL, archive)

    print("fetch_data: extracting ...")
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as z:
            z.extractall(DATA)
    elif tarfile.is_tarfile(archive):
        with tarfile.open(archive) as t:
            t.extractall(DATA)
    else:
        print("fetch_data: unrecognized archive format", file=sys.stderr)
        sys.exit(1)
    archive.unlink(missing_ok=True)
    print("fetch_data: done.")


if __name__ == "__main__":
    main()
