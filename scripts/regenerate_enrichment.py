#!/usr/bin/env python3
"""
Regenerate pathway_enrichment.csv for all 40 diseases using Fisher's exact test
against the MPI database KEGG pathway annotations.

Strategy:
1. Load MPIDB_v2.csv — extract all (HMDB_ID → Pathway) associations for human entries
2. For each disease, load predictions.csv (scored interactions above threshold)
3. Build the set of metabolites in the disease network
4. For each KEGG pathway, perform Fisher's exact test:
   - a = disease metabolites in pathway
   - b = disease metabolites NOT in pathway  
   - c = background metabolites in pathway (not in disease set)
   - d = background metabolites NOT in pathway
5. Apply BH FDR correction
6. Compute fold enrichment = (a/(a+b)) / ((a+c)/(a+b+c+d))
7. Save with columns: Pathway_ID, Pathway_Name, Fold_Enrichment, P_value,
   Metabolite_Count, Background_Count, FDR, Significant
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DISEASE_DIR = PROJECT_ROOT / "data" / "mpidatabase" / "disease_mpi"
MPIDB_PATH = PROJECT_ROOT / "data" / "mpidatabase" / "MPIDB_v2.csv"


def load_pathway_annotations():
    """Build metabolite → pathway mapping from MPI database."""
    mpi = pd.read_csv(MPIDB_PATH, dtype=str)
    human = mpi[mpi["Species"] == "Homo sapiens"].copy()

    # Build HMDB_ID → set of (Pathway_ID, Pathway_Name)
    met_to_pathways = {}
    pathway_names = {}
    all_metabolites = set()

    for _, row in human.iterrows():
        hmdb = str(row.get("HMDB ID", "")).strip()
        pw_id = str(row.get("Pathway_ID", "")).strip()
        pw_name = str(row.get("Pathway_Name", "")).strip()

        if not hmdb or hmdb == "nan":
            continue
        all_metabolites.add(hmdb)

        if not pw_id or pw_id == "nan":
            continue

        # Handle multi-pathway entries (semicolon-separated)
        pw_ids = [p.strip() for p in pw_id.split(";") if p.strip() and p.strip() != "nan"]
        pw_names_list = [p.strip() for p in pw_name.split(";") if p.strip() and p.strip() != "nan"]

        for i, pid in enumerate(pw_ids):
            pname = pw_names_list[i] if i < len(pw_names_list) else pid
            met_to_pathways.setdefault(hmdb, set()).add(pid)
            pathway_names[pid] = pname

    # Build pathway → set of metabolites (background)
    pathway_to_mets = {}
    for hmdb, pws in met_to_pathways.items():
        for pw in pws:
            pathway_to_mets.setdefault(pw, set()).add(hmdb)

    logger.info(f"Pathway annotations: {len(all_metabolites)} metabolites, "
                f"{len(pathway_to_mets)} pathways")
    return met_to_pathways, pathway_to_mets, pathway_names, all_metabolites


def fishers_exact_test(a, b, c, d):
    """Fisher's exact test (one-sided, greater)."""
    try:
        from scipy.stats import fisher_exact
        table = [[a, b], [c, d]]
        _, pval = fisher_exact(table, alternative="greater")
        return pval
    except ImportError:
        # Fallback: simple hypergeometric p-value approximation
        total = a + b + c + d
        if total == 0:
            return 1.0
        expected = (a + b) * (a + c) / total
        if expected == 0:
            return 1.0
        # Chi-square approximation
        chi2 = (a - expected) ** 2 / expected
        # Very rough p-value from chi2(1)
        import math
        pval = math.exp(-chi2 / 2)
        return min(1.0, pval)


def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    fdr = [0.0] * n
    cummin = 1.0
    for rank_idx in range(n - 1, -1, -1):
        orig_idx, pval = indexed[rank_idx]
        rank = rank_idx + 1
        adjusted = pval * n / rank
        cummin = min(cummin, adjusted)
        fdr[orig_idx] = min(1.0, cummin)
    return fdr


def compute_enrichment(disease_dir, met_to_pathways, pathway_to_mets,
                       pathway_names, all_metabolites, score_threshold=0.3):
    """Compute pathway enrichment for one disease."""
    pred_path = disease_dir / "predictions.csv"
    met_path = disease_dir / "metabolites.csv"

    # Get disease metabolites (HMDB IDs)
    disease_hmdb_ids = set()

    if pred_path.exists():
        pred_df = pd.read_csv(pred_path)
        if not pred_df.empty and "Prediction Score" in pred_df.columns:
            scores = pd.to_numeric(pred_df["Prediction Score"], errors="coerce")
            above = pred_df[scores >= score_threshold]
            if "HMDB_ID" in above.columns:
                disease_hmdb_ids = set(above["HMDB_ID"].dropna().astype(str).unique())

    # Also include all metabolites from metabolites.csv
    if met_path.exists():
        met_df = pd.read_csv(met_path)
        if "hmdb_id" in met_df.columns:
            disease_hmdb_ids.update(met_df["hmdb_id"].dropna().astype(str).unique())

    if not disease_hmdb_ids:
        return pd.DataFrame()

    # Find which disease metabolites are in the annotation background
    disease_mets_annotated = disease_hmdb_ids & all_metabolites
    background_size = len(all_metabolites)
    disease_size = len(disease_mets_annotated)

    if disease_size == 0:
        return pd.DataFrame()

    rows = []
    for pw_id, pw_mets in pathway_to_mets.items():
        # a = disease mets in this pathway
        a = len(disease_mets_annotated & pw_mets)
        if a == 0:
            continue
        b = disease_size - a
        c = len(pw_mets) - a
        d = background_size - a - b - c

        pval = fishers_exact_test(a, b, c, d)

        # Fold enrichment
        observed_ratio = a / disease_size if disease_size > 0 else 0
        expected_ratio = len(pw_mets) / background_size if background_size > 0 else 0
        fold = observed_ratio / expected_ratio if expected_ratio > 0 else 0

        rows.append({
            "Pathway_ID": pw_id,
            "Pathway_Name": pathway_names.get(pw_id, pw_id),
            "Fold_Enrichment": round(fold, 4),
            "P_value": pval,
            "Metabolite_Count": a,
            "Background_Count": len(pw_mets),
        })

    if not rows:
        return pd.DataFrame(columns=[
            "Pathway_ID", "Pathway_Name", "Fold_Enrichment", "P_value",
            "Metabolite_Count", "Background_Count", "FDR", "Significant"
        ])

    df = pd.DataFrame(rows)
    df = df.sort_values("P_value").reset_index(drop=True)

    # FDR correction
    df["FDR"] = bh_fdr(df["P_value"].tolist())
    df["Significant"] = df["FDR"].apply(lambda x: "Yes" if x <= 0.25 else "No")

    return df


def main():
    met_to_pathways, pathway_to_mets, pathway_names, all_metabolites = load_pathway_annotations()

    diseases_processed = 0
    for disease_dir in sorted(DISEASE_DIR.iterdir()):
        if not disease_dir.is_dir():
            continue

        logger.info(f"Processing: {disease_dir.name}")
        enrich_df = compute_enrichment(
            disease_dir, met_to_pathways, pathway_to_mets,
            pathway_names, all_metabolites
        )

        enrich_df.to_csv(disease_dir / "pathway_enrichment.csv", index=False)
        n_sig = len(enrich_df[enrich_df.get("Significant", pd.Series()) == "Yes"]) if not enrich_df.empty else 0
        logger.info(f"  {len(enrich_df)} pathways ({n_sig} significant at FDR<=0.25)")
        diseases_processed += 1

    logger.info(f"\nDone! Processed {diseases_processed} diseases.")


if __name__ == "__main__":
    main()
