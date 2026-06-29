#!/usr/bin/env python3
"""
CoreMet — Build Metabolite–Microbe Interaction (MMI) Database.

Data sources:
  1. gutMGene v2.0 — gut microbiota–metabolite interactions (experimentally validated)
  2. AGORA2 — genome-scale metabolic models (future expansion)

Output: data/databases/mmi_database.csv

Schema:
  Metabolite_Name, HMDB_ID, KEGG_ID, PubChem_CID, ChEBI_ID, SMILES,
  Microbe_Name, Taxonomy_ID, Rank,
  Substrate, Substrate_PubChem_CID,
  Relationship_Type (causal/correlative),
  Tissue, Organism, Evidence_Level (experimental/predicted),
  Experimental_Method, PMID, Source (gutMGene/AGORA2)
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data"
DB_DIR = DATA_DIR / "databases"

# gutMGene source (placed by user outside workspace)
# Script: .../mpi-vgae-web/scripts/  -> parent.parent = mpi-vgae-web
# gutMGene: .../MPI_web_server/gutMGene/  -> parent.parent.parent = MPI_web_server
GUTMGENE_DIR = Path(__file__).resolve().parent.parent.parent / "gutMGene"

# MEI / MPI databases for SMILES cross-reference
MEI_DB_PATH = DB_DIR / "mei_database.csv"
MPI_RAW_PATH = DATA_DIR / "raw" / "mpidatabase" / "MPIDB_May2024.csv"

OUTPUT_PATH = DB_DIR / "mmi_database.csv"

# Output columns (extended schema)
OUTPUT_COLS = [
    "Metabolite_Name", "HMDB_ID", "KEGG_ID", "PubChem_CID", "ChEBI_ID", "SMILES",
    "Microbe_Name", "Taxonomy_ID", "Rank",
    "Substrate", "Substrate_PubChem_CID",
    "Relationship_Type", "Tissue", "Organism",
    "Evidence_Level", "Experimental_Method", "PMID", "Source",
]


def _build_smiles_lookup() -> dict:
    """Build HMDB_ID -> SMILES lookup from MEI and MPI databases."""
    lookup = {}

    # From MEI database
    if MEI_DB_PATH.exists():
        mei = pd.read_csv(MEI_DB_PATH, usecols=["HMDB_ID", "SMILES"])
        for _, row in mei.dropna(subset=["HMDB_ID", "SMILES"]).iterrows():
            hmdb = str(row["HMDB_ID"]).strip()
            smiles = str(row["SMILES"]).strip()
            if hmdb and smiles and smiles != "nan":
                lookup[hmdb] = smiles
        print(f"  MEI SMILES lookup: {len(lookup)} unique HMDB->SMILES mappings")

    # From MPI raw database (uses spaces in column names)
    if MPI_RAW_PATH.exists():
        mpi = pd.read_csv(MPI_RAW_PATH, usecols=["HMDB ID", "SMILES"])
        added = 0
        for _, row in mpi.dropna(subset=["HMDB ID", "SMILES"]).iterrows():
            hmdb = str(row["HMDB ID"]).strip()
            smiles = str(row["SMILES"]).strip()
            if hmdb and smiles and smiles != "nan" and hmdb not in lookup:
                lookup[hmdb] = smiles
                added += 1
        print(f"  MPI SMILES lookup: +{added} new mappings (total {len(lookup)})")

    return lookup


def _clean_str(val) -> str:
    """Convert a value to clean string, handling NaN/float."""
    if pd.isna(val):
        return ""
    s = str(val).strip()
    return "" if s == "nan" else s


def _clean_int_str(val) -> str:
    """Convert a numeric value to clean integer string (no .0 suffix)."""
    s = _clean_str(val)
    if s and "." in s:
        try:
            s = str(int(float(s)))
        except (ValueError, OverflowError):
            pass
    return s


def _classify_relationship(associative_mode) -> str:
    """Map gutMGene 'Associative mode' to a cleaner relationship type."""
    mode = _clean_str(associative_mode).lower()
    if mode == "causally":
        return "causal"
    elif mode == "correlatively":
        return "correlative"
    return mode if mode else "unknown"


def _classify_tissue(sample) -> str:
    """Map gutMGene 'Sample' field to tissue category."""
    s = _clean_str(sample).lower()
    if not s:
        return "gut"
    if any(x in s for x in ["fecal", "feces", "stool"]):
        return "fecal"
    elif any(x in s for x in ["cecal", "cecum", "caecal"]):
        return "cecal"
    elif "colon" in s:
        return "colon"
    elif "plasma" in s:
        return "plasma"
    elif "serum" in s:
        return "serum"
    elif "urine" in s:
        return "urine"
    elif "brain" in s:
        return "brain"
    elif "liver" in s:
        return "liver"
    elif "intestin" in s:
        return "intestinal"
    return s


def parse_gutmgene() -> pd.DataFrame:
    """Parse gutMGene Microbe-Metabolite interactions."""
    fpath = GUTMGENE_DIR / "Gut Microbe-Microbial metabolite.csv"
    if not fpath.exists():
        print(f"  [SKIP] gutMGene not found: {fpath}")
        return pd.DataFrame()

    print(f"  Reading: {fpath}")
    gm = pd.read_csv(fpath)
    print(f"  Raw rows: {len(gm)}")

    # Build SMILES lookup
    smiles_lookup = _build_smiles_lookup()

    records = []
    for _, row in gm.iterrows():
        metabolite = _clean_str(row.get("Metabolite"))
        if not metabolite:
            continue
        microbe = _clean_str(row.get("Gut Microbiota"))
        if not microbe:
            continue

        hmdb = _clean_str(row.get("Metabolite HMDB"))
        kegg = _clean_str(row.get("Metabolite KEGG"))
        pubchem = _clean_int_str(row.get("Metabolite PubChem CID"))
        chebi = _clean_str(row.get("Metabolite ChEBI"))
        tax_id = _clean_int_str(row.get("Gut Microbiota NCBI ID"))
        rank = _clean_str(row.get("Rank"))
        substrate = _clean_str(row.get("Substrate"))
        sub_pubchem = _clean_int_str(row.get("Substrate PubChem CID"))
        organism = _clean_str(row.get("human/mouse"))
        exp_method = _clean_str(row.get("Experimental method"))
        pmid = _clean_int_str(row.get("PMID"))

        rel_type = _classify_relationship(row.get("Associative mode"))
        tissue = _classify_tissue(row.get("Sample"))
        smiles = smiles_lookup.get(hmdb, "") if hmdb else ""

        records.append({
            "Metabolite_Name": metabolite,
            "HMDB_ID": hmdb,
            "KEGG_ID": kegg,
            "PubChem_CID": pubchem,
            "ChEBI_ID": chebi,
            "SMILES": smiles,
            "Microbe_Name": microbe,
            "Taxonomy_ID": tax_id,
            "Rank": rank,
            "Substrate": substrate,
            "Substrate_PubChem_CID": sub_pubchem,
            "Relationship_Type": rel_type,
            "Tissue": tissue,
            "Organism": organism,
            "Evidence_Level": "experimental",
            "Experimental_Method": exp_method,
            "PMID": pmid,
            "Source": "gutMGene",
        })

    df = pd.DataFrame(records, columns=OUTPUT_COLS)
    print(f"  Parsed: {len(df)} records")
    return df


def parse_gutmgene_metabolite_gene() -> pd.DataFrame:
    """Parse gutMGene Metabolite-Host Gene data (supplementary info)."""
    fpath = GUTMGENE_DIR / "Microbial metabolite-Host Gene.csv"
    if not fpath.exists():
        return pd.DataFrame()
    mg = pd.read_csv(fpath)
    print(f"  Supplementary metabolite-gene data: {len(mg)} records "
          f"({mg['Metabolite'].nunique()} metabolites)")
    return mg


def build_mmi_database():
    """Build the full MMI database."""
    print("=" * 60)
    print("CoreMet — Build MMI Database")
    print("=" * 60)

    gm_avail = (GUTMGENE_DIR / "Gut Microbe-Microbial metabolite.csv").exists()
    print(f"\nData sources:")
    print(f"  gutMGene Microbe-Metabolite: {'OK' if gm_avail else 'NOT FOUND'}")
    print(f"  gutMGene dir: {GUTMGENE_DIR}")
    print(f"  AGORA2: DEFERRED (5.8 GB SBML models — future expansion)")

    if not gm_avail:
        print("\nNo data sources available. Creating empty schema CSV.")
        empty = pd.DataFrame(columns=OUTPUT_COLS)
        DB_DIR.mkdir(parents=True, exist_ok=True)
        empty.to_csv(OUTPUT_PATH, index=False)
        return

    # Step 1: Parse gutMGene
    print("\n--- Step 1: Parse gutMGene Microbe-Metabolite ---")
    gm_df = parse_gutmgene()

    # Step 2: supplementary info
    print("\n--- Step 2: Supplementary metabolite-gene data ---")
    parse_gutmgene_metabolite_gene()

    if gm_df.empty:
        print("\nNo records extracted.")
        return

    # Step 3: Deduplicate
    print("\n--- Step 3: Deduplicate ---")
    before = len(gm_df)
    gm_df = gm_df.drop_duplicates(
        subset=["Metabolite_Name", "Microbe_Name", "Relationship_Type", "Tissue"]
    )
    after = len(gm_df)
    print(f"  {before} -> {after} records (removed {before - after} duplicates)")

    # Step 4: Summary
    print("\n--- Step 4: Summary ---")
    print(f"  Total MMI records: {len(gm_df)}")
    print(f"  Unique metabolites: {gm_df['Metabolite_Name'].nunique()}")
    print(f"  Unique microbes: {gm_df['Microbe_Name'].nunique()}")
    print(f"  Relationship types: {gm_df['Relationship_Type'].value_counts().to_dict()}")
    print(f"  Organisms: {gm_df['Organism'].value_counts().to_dict()}")
    for col in ["HMDB_ID", "KEGG_ID", "PubChem_CID", "ChEBI_ID", "SMILES", "Taxonomy_ID"]:
        filled = (gm_df[col] != "").sum()
        print(f"  {col}: {filled}/{len(gm_df)} ({100*filled/len(gm_df):.1f}%)")

    # Write
    DB_DIR.mkdir(parents=True, exist_ok=True)
    gm_df.to_csv(OUTPUT_PATH, index=False)
    print(f"\n  Written: {OUTPUT_PATH}")
    print("  MMI database built successfully!")


if __name__ == "__main__":
    build_mmi_database()
