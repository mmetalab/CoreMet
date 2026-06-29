#!/usr/bin/env python3
"""
build_release_csvs.py — Clean, deduplicated download files for CoreMet.

Produces data/databases/release/*.csv whose row counts EXACTLY match the numbers
reported in the manuscript / coremetdb_stats.json, so a user who downloads a layer
counts the same number of interactions stated in the paper.

  - MPI : curated-only (drop unsourced pathway_participant) + dedup
  - MEI / MMI : dedup by edge key (collapses pathway/role multiplicity)
  - MDI / MDrI / MGI / mGWAS : already unique -> copied verbatim

    conda run -n mpi-vgae python scripts/build_release_csvs.py
"""
from __future__ import annotations
import json
import shutil
from pathlib import Path
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
REL = DATA / "databases" / "release"
REL.mkdir(parents=True, exist_ok=True)
STATS = json.loads((DATA / "coremetdb_stats.json").read_text())

# name -> (source path, dedup-key cols, needs_processing)
SRC = {
    "MPI":  (DATA / "mpidatabase" / "MPIDB_v3.csv",            ["HMDB ID", "Uniprot ID", "Species"], True),
    "MEI":  (DATA / "databases" / "mei_database_v2_enriched.csv", ["HMDB_ID", "EC_Number", "Species"], True),
    "MMI":  (DATA / "databases" / "mmi_database_v3.csv",       ["HMDB_ID", "Microbe_Name", "Relationship_Type"], True),
    "MDI":  (DATA / "databases" / "mdi_database_v3.csv",       None, False),
    "MDrI": (DATA / "databases" / "mdri_database_v3.csv",      None, False),
    "MGI":  (DATA / "databases" / "mgi_database_v3.csv",       None, False),
    "mGWAS":(DATA / "databases" / "mgwas_database_v3.csv",     None, False),
}
OUTNAME = {"MPI": "mpi", "MEI": "mei", "MMI": "mmi", "MDI": "mdi",
           "MDrI": "mdri", "MGI": "mgi", "mGWAS": "mgwas"}

# Empty entity-name fields to backfill so every browsed record shows a name.
# (name_col, id_col, kind): kind "met" resolves the HMDB id to a metabolite name
# via the entity registry (real name), falling back to the id; "id" uses the id column.
NAME_FILL = {
    "MPI":  [("Metabolite Name", "HMDB ID", "met"), ("Protein Name", "Uniprot ID", "id")],
    "MEI":  [("Metabolite_Name", "HMDB_ID", "met"), ("Enzyme_Name", "Uniprot_ID", "id")],
}


def _hmdb2name() -> dict:
    """HMDB id -> metabolite name, from the CoreMet entity registry."""
    try:
        reg = json.loads((DATA / "coremetdb_entity_registry.json").read_text())
    except Exception:
        return {}
    out = {}
    for _id, info in reg.get("reverse", {}).items():
        if info.get("type") == "metabolite":
            hmdb = (info.get("external_ids") or {}).get("hmdb_id")
            if hmdb and info.get("name"):
                out[hmdb] = info["name"]
    return out


def _backfill_names(name: str, df, h2n: dict):
    for name_col, id_col, kind in NAME_FILL.get(name, []):
        if name_col not in df.columns or id_col not in df.columns:
            continue
        empty = df[name_col].isna() | (df[name_col].astype(str).str.strip() == "")
        if not empty.any():
            continue
        ids = df.loc[empty, id_col].astype(str)
        if kind == "met":
            df.loc[empty, name_col] = ids.map(lambda h: h2n.get(h, h))
        else:
            df.loc[empty, name_col] = ids
    return df


def main() -> None:
    h2n = _hmdb2name()
    for name, (path, keys, process) in SRC.items():
        out = REL / f"coremetdb_{OUTNAME[name]}.csv"
        target = STATS["databases"][name]["interactions"]
        if not process:
            shutil.copyfile(path, out)
            n = sum(1 for _ in open(out)) - 1
            print(f"{name:6s} copied            -> {n:>9,} rows  (target {target:,})")
            continue

        df = pd.read_csv(path, dtype=str, low_memory=False)
        if name == "MPI":
            df = df[df["Evidence_Source"].notna() | (df["interaction_subtype"] == "reaction_participant")]
            # Rhea reaction_participant rows carry a RHEA reaction_id but no Evidence_Source;
            # label their provenance so every released record has a source.
            _no_src = df["Evidence_Source"].isna() | (df["Evidence_Source"].astype(str).str.strip() == "")
            df.loc[_no_src, "Evidence_Source"] = "Rhea"
        df = df.drop_duplicates(subset=keys)
        df = _backfill_names(name, df, h2n)
        df.to_csv(out, index=False)
        flag = "OK" if len(df) == target else "MISMATCH!"
        print(f"{name:6s} curated/deduped   -> {len(df):>9,} rows  (target {target:,})  {flag}")

    print("Release files in", REL)


if __name__ == "__main__":
    main()
