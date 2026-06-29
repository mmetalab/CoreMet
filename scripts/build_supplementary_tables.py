#!/usr/bin/env python3
"""Build CoreMet supplementary tables (CSV) from the canonical stats JSON."""
from __future__ import annotations
import csv
import json
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"
OUT = DATA.parent.parent / "Submission_NAR" / "Supplementary"
OUT.mkdir(parents=True, exist_ok=True)
S = json.loads((DATA / "coremetdb_stats.json").read_text())
U = json.loads((DATA / "coremetdb_use_cases.json").read_text())


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print("wrote", OUT / name)


# Table S1 — database content & provenance
rows = []
for k, d in S["databases"].items():
    pm = "" if d["pmid_coverage"] is None else f"{d['pmid_coverage']*100:.1f}%"
    rows.append([k, d["label"], d["interactions"], d["metabolites"],
                 d["targets"], d["target_label"], d["source"], pm])
rows.append(["Total", "7 types", S["totals"]["interactions"],
             S["totals"]["unique_metabolites_union"], "", "", "", ""])
write("Table_S1_database_content.csv",
      ["Database", "Type", "Interactions", "Metabolites", "Targets",
       "Target_type", "Source", "PubMed_coverage"], rows)

# Table S2 — metabolite cross-reference coverage
xr = S["cross_references"]
write("Table_S2_cross_reference_coverage.csv",
      ["Identifier", "Metabolites", "Coverage_pct"],
      [[k, xr[k]["count"], xr[k]["pct"]] for k in ["HMDB", "SMILES", "KEGG", "PubChem", "ChEBI"]])

# Table S3 — cross-type metabolite overlap
ov = S["cross_type_overlap"]["metabolites_by_db_count"]
write("Table_S3_cross_type_overlap.csv",
      ["N_interaction_layers", "N_metabolites"],
      [[k, ov[k]] for k in sorted(ov, key=int)])

# Table S4 — use-case cross-layer traversal
layers = ["MPI", "MEI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"]
rows = []
for key in ["hcc", "metformin", "butyrate", "gckr"]:
    uc = U[key]
    rows.append([uc["anchor"], uc.get("metabolites", "")] + [uc.get(L, "") for L in layers])
write("Table_S4_use_case_traversal.csv",
      ["Use_case_anchor", "Seed_metabolites"] + layers, rows)

print("Done.")
