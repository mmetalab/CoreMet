#!/usr/bin/env python3
"""
compute_use_cases.py — Cross-layer traversal counts for the four CoreMet
manuscript use cases, computed from the canonical v3 CSVs (MPI curated-only).

Writes data/coremetdb_use_cases.json. Feeds Manuscript "Use Cases" + Figure 3.

    conda activate mpi-vgae
    python scripts/compute_use_cases.py
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = DATA / "coremetdb_use_cases.json"
CHUNK = 200_000

FILES = {
    "MPI":  (DATA / "mpidatabase" / "MPIDB_v3.csv",            "HMDB ID",  "Uniprot ID"),
    "MEI":  (DATA / "databases" / "mei_database_v2_enriched.csv", "HMDB_ID", "EC_Number"),
    "MDI":  (DATA / "databases" / "mdi_database_v3.csv",       "HMDB_ID",  "Disease_Name"),
    "MMI":  (DATA / "databases" / "mmi_database_v3.csv",       "HMDB_ID",  "Microbe_Name"),
    "MDrI": (DATA / "databases" / "mdri_database_v3.csv",      "HMDB_ID",  "Drug_Name"),
    "MGI":  (DATA / "databases" / "mgi_database_v3.csv",       "HMDB_ID",  "Gene_Symbol"),
    "mGWAS":(DATA / "databases" / "mgwas_database_v3.csv",     "HMDB_ID",  "rsID"),
}


def _norm(s: pd.Series) -> pd.Series:
    s = s.dropna().astype(str).str.strip()
    return s[~s.isin({"", "nan", "None", "-", "NA"})]


def load(db: str, usecols: list[str]) -> pd.DataFrame:
    path, *_ = FILES[db]
    frames = []
    for ch in pd.read_csv(path, dtype=str, usecols=lambda c: c in usecols,
                          chunksize=CHUNK, on_bad_lines="skip", low_memory=False):
        if db == "MPI" and "Evidence_Source" in ch and "interaction_subtype" in ch:
            ch = ch[ch["Evidence_Source"].notna() | (ch["interaction_subtype"] == "reaction_participant")]
        frames.append(ch)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=usecols)


def targets_for_metabolites(db: str, hmdb_ids: set[str]) -> int:
    """Unique partner entities connected to the given metabolites in `db`."""
    _, mcol, tcol = FILES[db]
    cols = [mcol, tcol]
    if db == "MPI":
        cols += ["Evidence_Source", "interaction_subtype"]
    found: set[str] = set()
    path = FILES[db][0]
    for ch in pd.read_csv(path, dtype=str, usecols=lambda c: c in cols,
                          chunksize=CHUNK, on_bad_lines="skip", low_memory=False):
        if db == "MPI":
            ch = ch[ch["Evidence_Source"].notna() | (ch["interaction_subtype"] == "reaction_participant")]
        hit = ch[ch[mcol].isin(hmdb_ids)]
        found.update(_norm(hit[tcol]))
    return len(found)


def metabolites_from(db: str, key_col: str, match) -> set[str]:
    """HMDB ids whose `key_col` matches `match` (callable or contains-string)."""
    _, mcol, _ = FILES[db]
    path = FILES[db][0]
    out: set[str] = set()
    for ch in pd.read_csv(path, dtype=str, usecols=lambda c: c in [mcol, key_col],
                          chunksize=CHUNK, on_bad_lines="skip", low_memory=False):
        col = ch[key_col].fillna("")
        mask = col.map(match) if callable(match) else col.str.contains(match, case=False, regex=False)
        out.update(_norm(ch.loc[mask, mcol]))
    return out


def trace(anchor_label: str, hmdb_ids: set[str], layers: list[str]) -> dict:
    res = {"anchor": anchor_label, "metabolites": len(hmdb_ids)}
    for db in layers:
        res[db] = targets_for_metabolites(db, hmdb_ids)
    return res


def main() -> None:
    out = {}

    # Use Case 3: Butyrate (single metabolite) — direct profile across layers.
    but = {"HMDB0000039"}
    out["butyrate"] = trace("Butyrate (HMDB0000039)", but,
                            ["MPI", "MEI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"])

    # Use Case 1: HCC disease tracing.
    hcc_mets = metabolites_from("MDI", "Disease_Name", "hepatocellular carcinoma")
    out["hcc"] = trace("Hepatocellular carcinoma", hcc_mets,
                       ["MPI", "MGI", "MMI", "MDrI", "mGWAS"])

    # Use Case 2: Metformin drug bridge.
    met_mets = metabolites_from("MDrI", "Drug_Name", "metformin")
    out["metformin"] = trace("Metformin", met_mets, ["MDI", "MMI", "MGI", "MPI"])

    # Use Case 4: GCKR locus (rs1260326) GWAS mechanistic linking.
    gckr_mets = metabolites_from("mGWAS", "rsID", lambda v: v.strip() == "rs1260326")
    out["gckr"] = trace("GCKR locus (rs1260326)", gckr_mets,
                       ["MPI", "MDI", "MGI", "MMI"])

    OUT.write_text(json.dumps(out, indent=2))
    for k, v in out.items():
        print(k, "->", v)
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
