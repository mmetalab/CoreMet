#!/usr/bin/env python3
"""
Expand MDI (Metabolite–Disease Interactions) using CTD Chemical–Disease data.

Strategy:
  1. Parse HMDB XML to build CAS → (HMDB_ID, name, SMILES) lookup.
  2. Stream CTD Chemical–Disease TSV (852 MB), filtering for:
     - DirectEvidence = "therapeutic" or "marker/mechanism" (curated, not inferred).
     - Chemical has a CAS number matching an HMDB metabolite.
  3. Merge with existing MDI database, deduplicate.

Output: data/databases/mdi_database.csv (updated in place, backup created).
"""
import csv
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path

from lxml import etree
import pandas as pd

BASE        = Path(__file__).resolve().parent.parent
CTD_FILE    = BASE / "data" / "raw" / "ctd" / "CTD_chemicals_diseases.tsv"
HMDB_XML    = BASE / "data" / "raw" / "hmdb" / "hmdb_metabolites.xml"
MDI_DB      = BASE / "data" / "databases" / "mdi_database.csv"
MDI_BACKUP  = BASE / "data" / "databases" / "mdi_database.csv.bak"

HMDB_NS = "{http://www.hmdb.ca}"

# MeSH disease-ID → Category mapping (top-level MeSH tree)
MESH_CATEGORY_MAP = {
    "C01": "Infections",        "C04": "Neoplasms",
    "C05": "Musculoskeletal",   "C06": "Digestive",
    "C07": "Stomatognathic",    "C08": "Respiratory",
    "C09": "Otorhinolaryngologic", "C10": "Nervous System",
    "C11": "Eye",               "C12": "Urogenital",
    "C13": "Reproductive",      "C14": "Cardiovascular",
    "C15": "Hematologic",       "C16": "Congenital",
    "C17": "Skin",              "C18": "Nutritional/Metabolic",
    "C19": "Endocrine",         "C20": "Immune",
    "C23": "Pathological Conditions", "C25": "Substance-Related",
    "C26": "Wounds",            "F01": "Behavioral",
    "F03": "Mental",
}

def _guess_category(disease_id: str) -> str:
    """Rough category from MeSH disease ID prefix."""
    if disease_id.startswith("MESH:"):
        return "CTD_curated"
    return "CTD_curated"


def parse_hmdb_cas():
    """Build CAS → (HMDB_ID, name, SMILES) from HMDB XML."""
    print("[1/3] Parsing HMDB XML for CAS → metabolite mapping …")
    t0 = time.time()
    cas_to_metab = {}  # CAS → (hmdb_id, name, smiles)
    n = 0
    context = etree.iterparse(str(HMDB_XML), events=("end",), tag=f"{HMDB_NS}metabolite")
    for _, mel in context:
        parent = mel.getparent()
        if parent is None or parent.tag != f"{HMDB_NS}hmdb":
            mel.clear()
            continue
        n += 1
        hmdb_id = mel.findtext(f"{HMDB_NS}accession", "")
        name    = mel.findtext(f"{HMDB_NS}name", "")
        cas     = mel.findtext(f"{HMDB_NS}cas_registry_number", "") or ""
        smiles  = mel.findtext(f"{HMDB_NS}smiles", "") or ""
        if cas and cas != "0":
            cas_to_metab[cas] = (hmdb_id, name, smiles)
        mel.clear()
        if n % 50000 == 0:
            print(f"  … {n} metabolites")
    print(f"  Done: {n} metabolites, {len(cas_to_metab)} with CAS  ({time.time()-t0:.1f}s)")
    return cas_to_metab


def stream_ctd(cas_to_metab):
    """Stream CTD Chemical-Disease, keep only direct-evidence rows matching HMDB."""
    print("[2/3] Streaming CTD Chemical–Disease (direct evidence only) …")
    t0 = time.time()
    col_names = [
        "ChemicalName", "ChemicalID", "CasRN",
        "DiseaseName", "DiseaseID", "DirectEvidence",
        "InferenceGeneSymbol", "InferenceScore", "OmimIDs", "PubMedIDs",
    ]
    new_rows = []
    seen = set()
    line_count = 0

    with open(CTD_FILE, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            line_count += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            row = dict(zip(col_names, parts))
            direct = row.get("DirectEvidence", "").strip()
            if not direct:
                continue  # Skip inferred rows (99%+ of the file)
            cas = row.get("CasRN", "").strip()
            if not cas or cas not in cas_to_metab:
                continue

            hmdb_id, met_name, smiles = cas_to_metab[cas]
            disease_name = row.get("DiseaseName", "")
            disease_id = row.get("DiseaseID", "")
            pmids = row.get("PubMedIDs", "")

            # Map DirectEvidence to Association_Type
            if "therapeutic" in direct.lower():
                assoc = "Therapeutic"
            elif "marker" in direct.lower():
                assoc = "Biomarker"
            else:
                assoc = direct.capitalize()

            key = (hmdb_id, disease_id)
            if key in seen:
                continue
            seen.add(key)

            mesh_id = disease_id.replace("MESH:", "").replace("OMIM:", "")

            new_rows.append({
                "Metabolite_Name": met_name,
                "HMDB_ID": hmdb_id,
                "SMILES": smiles,
                "Disease_Name": disease_name,
                "Disease_ID": disease_id,
                "MeSH_ID": mesh_id,
                "Category": "CTD_curated",
                "Association_Type": assoc,
                "Evidence_Level": "CTD-curated",
                "Avg_Network_Score": "",
                "Source": "CTD",
                "PMID": pmids,
            })

            if line_count % 2000000 == 0:
                print(f"  … scanned {line_count:,} lines, {len(new_rows)} hits so far")

    print(f"  Done: scanned {line_count:,} lines, {len(new_rows)} new MDI from CTD  ({time.time()-t0:.1f}s)")
    return new_rows


def merge_and_save(new_rows):
    """Merge with existing MDI DB."""
    print("[3/3] Merging with existing MDI database …")
    # Read existing
    existing = []
    fields = [
        "Metabolite_Name", "HMDB_ID", "SMILES",
        "Disease_Name", "Disease_ID", "MeSH_ID", "Category",
        "Association_Type", "Evidence_Level", "Avg_Network_Score", "Source",
    ]
    if MDI_DB.exists():
        with open(MDI_DB) as f:
            reader = csv.DictReader(f)
            existing = list(reader)
            # Use the actual header from file
            if reader.fieldnames:
                fields = reader.fieldnames
        print(f"  Existing entries: {len(existing)}")
        shutil.copy2(MDI_DB, MDI_BACKUP)
        print(f"  Backup: {MDI_BACKUP}")

    # Dedup key
    existing_keys = set()
    for r in existing:
        key = (r.get("HMDB_ID", ""), r.get("Disease_ID", ""))
        existing_keys.add(key)

    added = 0
    for r in new_rows:
        key = (r["HMDB_ID"], r["Disease_ID"])
        if key not in existing_keys:
            existing.append(r)
            existing_keys.add(key)
            added += 1

    # Sort
    existing.sort(key=lambda r: (r.get("Metabolite_Name", "").lower(),
                                  r.get("Disease_Name", "").lower()))

    with open(MDI_DB, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(existing)

    total = len(existing)
    print(f"\n{'='*60}")
    print(f"MDI database expanded: {total - added} existing + {added} new = {total} total")
    print(f"Saved to {MDI_DB}")

    from collections import Counter
    src = Counter(r.get("Source", "") for r in existing)
    for s, c in src.most_common():
        print(f"  {s}: {c}")

    return total


if __name__ == "__main__":
    t_start = time.time()
    for p, label in [(CTD_FILE, "CTD"), (HMDB_XML, "HMDB")]:
        if not p.exists():
            print(f"ERROR: {label} file not found at {p}")
            sys.exit(1)

    cas_map = parse_hmdb_cas()
    new_rows = stream_ctd(cas_map)
    merge_and_save(new_rows)
    print(f"\nTotal time: {time.time() - t_start:.1f}s")
