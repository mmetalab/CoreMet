#!/usr/bin/env python3
"""
Build Metabolite–Disease Interaction (MDI) Database.

Phase 2a: Seed MDI database from existing curated disease case studies.
- Extracts metabolite-disease associations from disease_mpi/ folders
- Maps HMDB IDs to metabolite names and SMILES from MPIDB_v2
- Assigns Disease IDs (mapped to DOID where possible)
- Classifies association types (biomarker, therapeutic target, etc.)

Phase 2b (future): Extend with HMDB XML, CTD, DisGeNET downloads.
"""

import json
import os
import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent.parent
DISEASE_MPI_DIR = BASE_DIR / "data" / "mpidatabase" / "disease_mpi"
MPIDB_PATH = BASE_DIR / "data" / "mpidatabase" / "MPIDB_v2.csv"
MEI_DB_PATH = BASE_DIR / "data" / "databases" / "mei_database.csv"
OUTPUT_PATH = BASE_DIR / "data" / "databases" / "mdi_database.csv"
CASE_STUDY_DIR = BASE_DIR / "data" / "case_studies"

# ---------------------------------------------------------------------------
# Disease metadata — curated mapping to DOID and MeSH IDs
# ---------------------------------------------------------------------------
DISEASE_METADATA = {
    # ── Existing (6) ──
    "alzheimers": {
        "Disease_Name": "Alzheimer's disease",
        "Disease_ID": "DOID:10652",
        "MeSH_ID": "D000544",
        "Category": "Neurodegenerative",
    },
    "breast_cancer": {
        "Disease_Name": "Breast cancer",
        "Disease_ID": "DOID:1612",
        "MeSH_ID": "D001943",
        "Category": "Cancer",
    },
    "hcc": {
        "Disease_Name": "Hepatocellular carcinoma",
        "Disease_ID": "DOID:684",
        "MeSH_ID": "D006528",
        "Category": "Cancer",
    },
    "schizophrenia": {
        "Disease_Name": "Schizophrenia",
        "Disease_ID": "DOID:5419",
        "MeSH_ID": "D012559",
        "Category": "Neuropsychiatric",
    },
    "t2_diabetes": {
        "Disease_Name": "Type 2 diabetes mellitus",
        "Disease_ID": "DOID:9352",
        "MeSH_ID": "D003924",
        "Category": "Metabolic",
    },
    "thyroid_cancer": {
        "Disease_Name": "Thyroid cancer",
        "Disease_ID": "DOID:1781",
        "MeSH_ID": "D013964",
        "Category": "Cancer",
    },
    # ── Cancers (new) ──
    "lung_cancer": {
        "Disease_Name": "Lung cancer",
        "Disease_ID": "DOID:1324",
        "MeSH_ID": "D008175",
        "Category": "Cancer",
    },
    "colorectal_cancer": {
        "Disease_Name": "Colorectal cancer",
        "Disease_ID": "DOID:9256",
        "MeSH_ID": "D015179",
        "Category": "Cancer",
    },
    "pancreatic_cancer": {
        "Disease_Name": "Pancreatic cancer",
        "Disease_ID": "DOID:1793",
        "MeSH_ID": "D010190",
        "Category": "Cancer",
    },
    "prostate_cancer": {
        "Disease_Name": "Prostate cancer",
        "Disease_ID": "DOID:10283",
        "MeSH_ID": "D011471",
        "Category": "Cancer",
    },
    "ovarian_cancer": {
        "Disease_Name": "Ovarian cancer",
        "Disease_ID": "DOID:2394",
        "MeSH_ID": "D010051",
        "Category": "Cancer",
    },
    "gastric_cancer": {
        "Disease_Name": "Gastric cancer",
        "Disease_ID": "DOID:10534",
        "MeSH_ID": "D013274",
        "Category": "Cancer",
    },
    "bladder_cancer": {
        "Disease_Name": "Bladder cancer",
        "Disease_ID": "DOID:11054",
        "MeSH_ID": "D001749",
        "Category": "Cancer",
    },
    "renal_cell_carcinoma": {
        "Disease_Name": "Renal cell carcinoma",
        "Disease_ID": "DOID:4450",
        "MeSH_ID": "D002292",
        "Category": "Cancer",
    },
    "melanoma": {
        "Disease_Name": "Melanoma",
        "Disease_ID": "DOID:1909",
        "MeSH_ID": "D008545",
        "Category": "Cancer",
    },
    "glioblastoma": {
        "Disease_Name": "Glioblastoma",
        "Disease_ID": "DOID:3068",
        "MeSH_ID": "D005909",
        "Category": "Cancer",
    },
    "esophageal_cancer": {
        "Disease_Name": "Esophageal cancer",
        "Disease_ID": "DOID:5041",
        "MeSH_ID": "D004938",
        "Category": "Cancer",
    },
    "endometrial_cancer": {
        "Disease_Name": "Endometrial cancer",
        "Disease_ID": "DOID:1380",
        "MeSH_ID": "D016889",
        "Category": "Cancer",
    },
    "aml_leukemia": {
        "Disease_Name": "Acute myeloid leukemia",
        "Disease_ID": "DOID:9119",
        "MeSH_ID": "D015470",
        "Category": "Cancer",
    },
    "lymphoma": {
        "Disease_Name": "Lymphoma",
        "Disease_ID": "DOID:0060058",
        "MeSH_ID": "D008223",
        "Category": "Cancer",
    },
    "cervical_cancer": {
        "Disease_Name": "Cervical cancer",
        "Disease_ID": "DOID:4362",
        "MeSH_ID": "D002583",
        "Category": "Cancer",
    },
    "oral_cancer": {
        "Disease_Name": "Oral squamous cell carcinoma",
        "Disease_ID": "DOID:0050866",
        "MeSH_ID": "D009062",
        "Category": "Cancer",
    },
    "cholangiocarcinoma": {
        "Disease_Name": "Cholangiocarcinoma",
        "Disease_ID": "DOID:4947",
        "MeSH_ID": "D018281",
        "Category": "Cancer",
    },
    # ── Neurodegenerative / Neuropsychiatric (new) ──
    "parkinsons": {
        "Disease_Name": "Parkinson's disease",
        "Disease_ID": "DOID:14330",
        "MeSH_ID": "D010300",
        "Category": "Neurodegenerative",
    },
    "als": {
        "Disease_Name": "Amyotrophic lateral sclerosis",
        "Disease_ID": "DOID:332",
        "MeSH_ID": "D000690",
        "Category": "Neurodegenerative",
    },
    "multiple_sclerosis": {
        "Disease_Name": "Multiple sclerosis",
        "Disease_ID": "DOID:2377",
        "MeSH_ID": "D009103",
        "Category": "Autoimmune",
    },
    "depression": {
        "Disease_Name": "Major depressive disorder",
        "Disease_ID": "DOID:1596",
        "MeSH_ID": "D003865",
        "Category": "Neuropsychiatric",
    },
    # ── Metabolic (new) ──
    "obesity": {
        "Disease_Name": "Obesity",
        "Disease_ID": "DOID:9970",
        "MeSH_ID": "D009765",
        "Category": "Metabolic",
    },
    "nafld": {
        "Disease_Name": "Non-alcoholic fatty liver disease",
        "Disease_ID": "DOID:0080208",
        "MeSH_ID": "D065626",
        "Category": "Metabolic",
    },
    "metabolic_syndrome": {
        "Disease_Name": "Metabolic syndrome",
        "Disease_ID": "DOID:14688",
        "MeSH_ID": "D024821",
        "Category": "Metabolic",
    },
    # ── Cardiovascular (new) ──
    "heart_failure": {
        "Disease_Name": "Heart failure",
        "Disease_ID": "DOID:6000",
        "MeSH_ID": "D006333",
        "Category": "Cardiovascular",
    },
    "atherosclerosis": {
        "Disease_Name": "Atherosclerosis",
        "Disease_ID": "DOID:1936",
        "MeSH_ID": "D050197",
        "Category": "Cardiovascular",
    },
    "hypertension": {
        "Disease_Name": "Hypertension",
        "Disease_ID": "DOID:10763",
        "MeSH_ID": "D006973",
        "Category": "Cardiovascular",
    },
    # ── Autoimmune / Inflammatory (new) ──
    "ibd": {
        "Disease_Name": "Inflammatory bowel disease",
        "Disease_ID": "DOID:0050589",
        "MeSH_ID": "D015212",
        "Category": "Autoimmune",
    },
    "rheumatoid_arthritis": {
        "Disease_Name": "Rheumatoid arthritis",
        "Disease_ID": "DOID:7148",
        "MeSH_ID": "D001172",
        "Category": "Autoimmune",
    },
    "lupus": {
        "Disease_Name": "Systemic lupus erythematosus",
        "Disease_ID": "DOID:9074",
        "MeSH_ID": "D008180",
        "Category": "Autoimmune",
    },
    # ── Respiratory (new) ──
    "copd": {
        "Disease_Name": "Chronic obstructive pulmonary disease",
        "Disease_ID": "DOID:3083",
        "MeSH_ID": "D029424",
        "Category": "Respiratory",
    },
    "asthma": {
        "Disease_Name": "Asthma",
        "Disease_ID": "DOID:2841",
        "MeSH_ID": "D001249",
        "Category": "Respiratory",
    },
    # ── Other (new) ──
    "chronic_kidney_disease": {
        "Disease_Name": "Chronic kidney disease",
        "Disease_ID": "DOID:784",
        "MeSH_ID": "D051436",
        "Category": "Renal",
    },
    "covid19": {
        "Disease_Name": "COVID-19",
        "Disease_ID": "DOID:0080600",
        "MeSH_ID": "D000086382",
        "Category": "Infectious",
    },
}


def load_mpidb_metabolite_info():
    """Build HMDB → (name, SMILES) lookup from MPIDB_v2."""
    if not MPIDB_PATH.exists():
        return {}
    df = pd.read_csv(MPIDB_PATH, dtype=str).fillna("")
    lookup = {}
    for _, row in df.iterrows():
        hmdb = row.get("HMDB ID", "")
        if hmdb:
            name = row.get("Metabolite Name", "")
            smiles = row.get("SMILES", "")
            if hmdb not in lookup or (not lookup[hmdb]["SMILES"] and smiles):
                lookup[hmdb] = {"Metabolite_Name": name, "SMILES": smiles}
    return lookup


def load_mei_metabolite_info():
    """Supplement metabolite info from MEI database."""
    if not MEI_DB_PATH.exists():
        return {}
    df = pd.read_csv(MEI_DB_PATH, dtype=str).fillna("")
    lookup = {}
    for _, row in df.iterrows():
        hmdb = row.get("HMDB_ID", "")
        if hmdb:
            name = row.get("Metabolite_Name", "")
            smiles = row.get("SMILES", "")
            if hmdb not in lookup or (not lookup[hmdb]["SMILES"] and smiles):
                lookup[hmdb] = {"Metabolite_Name": name, "SMILES": smiles}
    return lookup


def extract_disease_metabolites():
    """Extract metabolite-disease pairs from all disease_mpi subdirectories."""
    records = []
    
    for disease_dir in sorted(DISEASE_MPI_DIR.iterdir()):
        if not disease_dir.is_dir():
            continue
        
        disease_key = disease_dir.name
        if disease_key not in DISEASE_METADATA:
            print(f"  Warning: Unknown disease '{disease_key}', skipping")
            continue
        
        meta = DISEASE_METADATA[disease_key]
        met_path = disease_dir / "metabolites.csv"
        pred_path = disease_dir / "predictions.csv"
        stats_path = disease_dir / "network_stats.json"
        
        if not met_path.exists():
            print(f"  Warning: No metabolites.csv for {disease_key}, skipping")
            continue
        
        # Load metabolite seed list (curated from literature)
        df_met = pd.read_csv(met_path, dtype=str).fillna("")
        
        # Load prediction data for additional context (confidence scores)
        pred_scores = {}
        if pred_path.exists():
            df_pred = pd.read_csv(pred_path, dtype=str).fillna("")
            if "HMDB_ID" in df_pred.columns and "Prediction Score" in df_pred.columns:
                # Average prediction score per metabolite as a measure of network involvement
                df_pred["Prediction Score"] = pd.to_numeric(df_pred["Prediction Score"], errors="coerce")
                for hmdb, group in df_pred.groupby("HMDB_ID"):
                    if hmdb:
                        pred_scores[hmdb] = round(group["Prediction Score"].mean(), 4)
        
        # Load network stats for hub information
        hub_metabolites = set()
        if stats_path.exists():
            with open(stats_path) as f:
                stats = json.load(f)
            for hub in stats.get("top_hubs", []):
                # Try to identify if hub is a metabolite (starts with HMDB)
                hub_name = hub.get("name", "") if isinstance(hub, dict) else str(hub)
                if hub_name.startswith("HMDB"):
                    hub_metabolites.add(hub_name)
        
        print(f"  {disease_key}: {len(df_met)} metabolites, "
              f"{len(pred_scores)} with prediction scores, "
              f"{len(hub_metabolites)} hub metabolites")
        
        for _, row in df_met.iterrows():
            hmdb_id = row.get("hmdb_id", "")
            name = row.get("name", "")
            
            # Determine association type
            if hmdb_id in hub_metabolites:
                assoc_type = "Hub biomarker"
            else:
                assoc_type = "Disease-associated biomarker"
            
            # Evidence level based on curation source
            evidence = "Literature-curated"
            
            # Average prediction score (network involvement)
            avg_score = pred_scores.get(hmdb_id, None)
            
            records.append({
                "Metabolite_Name": name,
                "HMDB_ID": hmdb_id,
                "SMILES": "",  # Will be filled later
                "Disease_Name": meta["Disease_Name"],
                "Disease_ID": meta["Disease_ID"],
                "MeSH_ID": meta["MeSH_ID"],
                "Category": meta["Category"],
                "Association_Type": assoc_type,
                "Evidence_Level": evidence,
                "Avg_Network_Score": avg_score if avg_score else "",
                "Source": f"CoreMet_disease_mpi/{disease_key}",
            })
    
    return pd.DataFrame(records)


def extract_case_study_metabolites():
    """Extract additional metabolite-disease pairs from case_studies/ directory."""
    records = []
    
    if not CASE_STUDY_DIR.exists():
        return pd.DataFrame(records)
    
    # Check for schizophrenia case study
    for subdir in sorted(CASE_STUDY_DIR.iterdir()):
        if not subdir.is_dir():
            continue
        
        disease_key = subdir.name
        pred_path = subdir / "predictions.csv"
        
        if not pred_path.exists():
            continue
        
        # Map case study folder names to disease metadata
        case_to_disease = {
            "schizophrenia": "schizophrenia",
            "thyroid_cancer": "thyroid_cancer",
            "hcc": "hcc",
        }
        
        mapped_key = case_to_disease.get(disease_key)
        if not mapped_key or mapped_key not in DISEASE_METADATA:
            continue
        
        meta = DISEASE_METADATA[mapped_key]
        
        df_pred = pd.read_csv(pred_path, dtype=str).fillna("")
        
        # Extract unique metabolites from the predictions
        # These have high-confidence predictions in the case study
        met_col = "Metabolite" if "Metabolite" in df_pred.columns else None
        hmdb_col = "HMDB_ID" if "HMDB_ID" in df_pred.columns else None
        
        if met_col is None:
            continue
        
        seen_mets = set()
        for _, row in df_pred.iterrows():
            met = row.get(met_col, "")
            hmdb = row.get(hmdb_col, "") if hmdb_col else ""
            
            if met and met not in seen_mets:
                seen_mets.add(met)
                records.append({
                    "Metabolite_Name": met,
                    "HMDB_ID": hmdb,
                    "SMILES": "",
                    "Disease_Name": meta["Disease_Name"],
                    "Disease_ID": meta["Disease_ID"],
                    "MeSH_ID": meta["MeSH_ID"],
                    "Category": meta["Category"],
                    "Association_Type": "Case-study predicted",
                    "Evidence_Level": "Predicted",
                    "Avg_Network_Score": "",
                    "Source": f"CoreMet_case_study/{disease_key}",
                })
        
        print(f"  Case study {disease_key}: {len(seen_mets)} unique metabolites")
    
    return pd.DataFrame(records)


def enrich_with_metabolite_info(df, mpidb_lookup, mei_lookup):
    """Fill in SMILES and improve metabolite names from MPIDB/MEI lookups."""
    for idx, row in df.iterrows():
        hmdb = row["HMDB_ID"]
        if hmdb:
            info = mpidb_lookup.get(hmdb) or mei_lookup.get(hmdb)
            if info:
                if not row["SMILES"] and info.get("SMILES"):
                    df.at[idx, "SMILES"] = info["SMILES"]
                if not row["Metabolite_Name"] and info.get("Metabolite_Name"):
                    df.at[idx, "Metabolite_Name"] = info["Metabolite_Name"]
    return df


def deduplicate(df):
    """Remove duplicate metabolite-disease pairs, keeping the best evidence."""
    if df.empty:
        return df
    
    # Priority: Literature-curated > Case-study predicted
    evidence_priority = {
        "Literature-curated": 0,
        "Predicted": 1,
    }
    
    df["_priority"] = df["Evidence_Level"].map(evidence_priority).fillna(2)
    df = df.sort_values("_priority").drop_duplicates(
        subset=["HMDB_ID", "Disease_ID"], keep="first"
    )
    df = df.drop(columns=["_priority"])
    
    return df.reset_index(drop=True)


def main():
    """Main pipeline to build the MDI database."""
    print("=" * 60)
    print("Building MDI Database (Phase 2a — Seed)")
    print("=" * 60)
    
    # Step 1: Load metabolite info lookups
    print("\n[1/6] Loading metabolite info from MPIDB...")
    mpidb_lookup = load_mpidb_metabolite_info()
    print(f"  MPIDB: {len(mpidb_lookup)} metabolites with HMDB IDs")
    
    print("\n[2/6] Loading metabolite info from MEI database...")
    mei_lookup = load_mei_metabolite_info()
    print(f"  MEI: {len(mei_lookup)} metabolites with HMDB IDs")
    
    # Step 2: Extract disease-metabolite pairs from disease_mpi
    print("\n[3/6] Extracting disease-metabolite pairs from disease_mpi/...")
    df_disease = extract_disease_metabolites()
    print(f"  Total: {len(df_disease)} metabolite-disease pairs")
    
    # Step 3: Extract additional pairs from case studies
    print("\n[4/6] Extracting from case studies...")
    df_cases = extract_case_study_metabolites()
    print(f"  Total: {len(df_cases)} additional pairs from case studies")
    
    # Step 4: Combine and deduplicate
    print("\n[5/6] Combining and deduplicating...")
    df_all = pd.concat([df_disease, df_cases], ignore_index=True)
    print(f"  Combined: {len(df_all)} pairs before dedup")
    df_all = deduplicate(df_all)
    print(f"  After dedup: {len(df_all)} unique pairs")
    
    # Step 5: Enrich with metabolite info
    print("\n[6/6] Enriching with metabolite info (names, SMILES)...")
    df_all = enrich_with_metabolite_info(df_all, mpidb_lookup, mei_lookup)
    
    # Fill any remaining empty SMILES/names
    n_with_smiles = (df_all["SMILES"].str.len() > 0).sum()
    n_with_name = (df_all["Metabolite_Name"].str.len() > 0).sum()
    print(f"  With SMILES: {n_with_smiles}/{len(df_all)} ({100*n_with_smiles/len(df_all):.1f}%)")
    print(f"  With name: {n_with_name}/{len(df_all)} ({100*n_with_name/len(df_all):.1f}%)")
    
    # Export
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df_all.to_csv(OUTPUT_PATH, index=False)
    
    print(f"\n{'=' * 60}")
    print(f"MDI DATABASE COMPLETE")
    print(f"{'=' * 60}")
    print(f"  Output: {OUTPUT_PATH}")
    print(f"  Total records: {len(df_all):,}")
    print(f"  Unique metabolites: {df_all['HMDB_ID'].nunique()}")
    print(f"  Unique diseases: {df_all['Disease_Name'].nunique()}")
    print(f"  Disease categories: {df_all['Category'].nunique()}")
    print(f"\n  Breakdown by disease:")
    for disease, group in df_all.groupby("Disease_Name"):
        print(f"    {disease}: {len(group)} metabolites")
    print(f"\n  Breakdown by evidence level:")
    for evidence, group in df_all.groupby("Evidence_Level"):
        print(f"    {evidence}: {len(group)} pairs")
    print(f"\n  Breakdown by association type:")
    for assoc, group in df_all.groupby("Association_Type"):
        print(f"    {assoc}: {len(group)} pairs")


if __name__ == "__main__":
    main()
