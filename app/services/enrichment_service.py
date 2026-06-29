"""
Pathway enrichment service — Fisher's exact test per KEGG pathway, BH FDR correction.
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent.parent
CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "kegg"


def load_pathway_annotations():
    """Load enzyme→pathway and pathway names from KEGG cache."""
    ec_pathways = {}
    pw_file = CACHE_DIR / "enzyme_pathway_links.json"
    if pw_file.exists():
        with open(pw_file) as f:
            ec_pathways = json.load(f)

    pw_names = {}
    name_file = CACHE_DIR / "pathway_names.json"
    if name_file.exists():
        with open(name_file) as f:
            pw_names = json.load(f)

    return ec_pathways, pw_names


def load_uniprot_to_ec():
    """Build UniProt→EC mapping from KEGG cache."""
    mapping_file = CACHE_DIR / "enzyme_uniprot_mapping.json"
    if not mapping_file.exists():
        # Try building from per-organism caches
        result = {}
        for org in ['hsa', 'mmu', 'rno', 'eco', 'bta', 'pae', 'ath', 'sce', 'dme', 'cel']:
            gene_file = CACHE_DIR / f"enzyme_genes_{org}.json"
            uniprot_file = CACHE_DIR / f"gene_uniprot_{org}.json"
            if gene_file.exists() and uniprot_file.exists():
                with open(gene_file) as f:
                    ec_genes = json.load(f)
                with open(uniprot_file) as f:
                    gene_uniprot = json.load(f)
                for ec_id, genes in ec_genes.items():
                    for gene in genes:
                        uid = gene_uniprot.get(gene, '')
                        if uid:
                            result.setdefault(uid, set()).add(ec_id)
        return {k: list(v) for k, v in result.items()}

    with open(mapping_file) as f:
        ec_to_uniprots = json.load(f)
    result = {}
    for ec_id, uniprots in ec_to_uniprots.items():
        for uid in uniprots:
            result.setdefault(uid, []).append(ec_id)
    return result


def run_enrichment(predictions_df, organism="All", fdr_threshold=0.05, score_threshold=0.5):
    """
    Run pathway enrichment analysis on prediction results.

    Args:
        predictions_df: DataFrame with at least 'Uniprot ID' or 'Protein' column
        organism: Target organism (unused currently, for future filtering)
        fdr_threshold: FDR cutoff for significance
        score_threshold: Minimum prediction score to include

    Returns:
        DataFrame with columns: Pathway, Pathway_Name, Fold_Enrichment, P_value, FDR, Protein_Count, Metabolite_Count
    """
    ec_pathways, pw_names = load_pathway_annotations()
    uniprot_to_ec = load_uniprot_to_ec()

    if not ec_pathways or not uniprot_to_ec:
        logger.warning("No pathway annotation data available")
        return pd.DataFrame()

    # Get proteins from predictions
    # Check multiple possible column names for UniProt IDs
    prot_col = None
    for candidate in ['Uniprot ID', 'Uniprot_ID', 'UniProt_ID', 'uniprot_id', 'Protein']:
        if candidate in predictions_df.columns:
            prot_col = candidate
            break
    if prot_col is None:
        logger.warning("No protein ID column found in predictions")
        return pd.DataFrame()

    # Filter by score if available
    if 'Prediction Score' in predictions_df.columns:
        scores = pd.to_numeric(predictions_df['Prediction Score'], errors='coerce')
        query_proteins = set(predictions_df.loc[scores >= score_threshold, prot_col].dropna())
    else:
        query_proteins = set(predictions_df[prot_col].dropna())

    if not query_proteins:
        return pd.DataFrame()

    # Helper: normalize pathway IDs — ec* and map* represent the same pathway,
    # but pathway_names.json only has map* keys, so merge them.
    def _normalize_pw(pw_id):
        if pw_id.startswith('ec'):
            return 'map' + pw_id[2:]
        return pw_id

    # Build pathway→protein sets for query (normalized IDs)
    query_pathway_proteins = {}
    for uid in query_proteins:
        ecs = uniprot_to_ec.get(uid, [])
        for ec in ecs:
            pws = ec_pathways.get(ec, [])
            for pw in pws:
                npw = _normalize_pw(pw)
                query_pathway_proteins.setdefault(npw, set()).add(uid)

    # Build background: all proteins with pathway annotations (normalized IDs)
    all_annotated_proteins = set()
    pathway_background = {}
    for uid, ecs in uniprot_to_ec.items():
        for ec in ecs:
            pws = ec_pathways.get(ec, [])
            if pws:
                all_annotated_proteins.add(uid)
                for pw in pws:
                    npw = _normalize_pw(pw)
                    pathway_background.setdefault(npw, set()).add(uid)

    N = len(all_annotated_proteins)  # total background
    n = len(query_proteins & all_annotated_proteins)  # query with annotations

    if N == 0 or n == 0:
        return pd.DataFrame()

    # Fisher's exact test per pathway
    results = []
    for pw_id, bg_proteins in pathway_background.items():
        K = len(bg_proteins)  # pathway size in background
        k = len(query_pathway_proteins.get(pw_id, set()))  # overlap

        if k == 0:
            continue

        # 2x2 contingency table
        a = k                    # query & pathway
        b = n - k                # query & not pathway
        c = K - k                # not query & pathway
        d = N - n - K + k        # not query & not pathway

        if d < 0:
            d = 0

        _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

        # Fold enrichment
        expected = (n * K) / N if N > 0 else 0
        fold = k / expected if expected > 0 else 0

        pw_name = pw_names.get(pw_id, pw_id)

        results.append({
            'Pathway_ID': pw_id,
            'Pathway_Name': pw_name,
            'Fold_Enrichment': round(fold, 2),
            'P_value': p_value,
            'Protein_Count': k,
            'Background_Count': K,
        })

    if not results:
        return pd.DataFrame()

    df_results = pd.DataFrame(results)

    # BH FDR correction
    reject, fdr_values, _, _ = multipletests(df_results['P_value'], method='fdr_bh')
    df_results['FDR'] = fdr_values
    df_results['Significant'] = reject

    # Sort by FDR
    df_results = df_results.sort_values('FDR').reset_index(drop=True)

    return df_results


def run_disease_enrichment(metabolite_ids, fdr_threshold=0.05):
    """
    Run disease enrichment analysis on a list of metabolites.

    Uses Fisher's exact test with BH FDR correction against the MDI database.

    Args:
        metabolite_ids: List of HMDB IDs or metabolite names
        fdr_threshold: FDR cutoff for significance

    Returns:
        DataFrame: Disease_Name, Disease_ID, Category, Fold_Enrichment, P_value,
                   Metabolite_Count, Background_Count, FDR, Significant
    """
    from app.services.mdi_service import get_mdi_db

    mdi_df = get_mdi_db()
    if mdi_df.empty:
        logger.warning("MDI database is empty — cannot run disease enrichment")
        return pd.DataFrame()

    # Build disease → metabolite sets from MDI database (background)
    disease_metabolites = {}
    disease_meta = {}
    all_mdi_metabolites = set()

    for _, row in mdi_df.iterrows():
        disease = row.get("Disease_Name", "")
        hmdb = row.get("HMDB_ID", "")
        met_name = row.get("Metabolite_Name", "")
        if not disease:
            continue

        # Use HMDB ID as primary identifier, fallback to name
        met_key = hmdb if hmdb else met_name
        if met_key:
            disease_metabolites.setdefault(disease, set()).add(met_key)
            all_mdi_metabolites.add(met_key)
            if disease not in disease_meta:
                disease_meta[disease] = {
                    "Disease_ID": row.get("Disease_ID", ""),
                    "Category": row.get("Category", ""),
                }

    N = len(all_mdi_metabolites)  # total unique metabolites in MDI background
    if N == 0:
        return pd.DataFrame()

    # Normalize query IDs — could be HMDB IDs, metabolite names, or mixed
    query_ids = set()
    for mid in metabolite_ids:
        mid_str = str(mid).strip()
        if mid_str:
            query_ids.add(mid_str)
            # Also check if it's a name that maps to an HMDB in the MDI db
            name_match = mdi_df.loc[mdi_df["Metabolite_Name"].str.lower() == mid_str.lower(), "HMDB_ID"]
            for hmdb in name_match:
                if hmdb:
                    query_ids.add(hmdb)
            # And reverse: HMDB → name
            hmdb_match = mdi_df.loc[mdi_df["HMDB_ID"] == mid_str, "Metabolite_Name"]
            for name in hmdb_match:
                if name:
                    query_ids.add(name)

    # Intersect with background
    query_annotated = query_ids & all_mdi_metabolites
    n = len(query_annotated)

    if n == 0:
        logger.info(f"No query metabolites found in MDI background (query={len(query_ids)})")
        return pd.DataFrame()

    # Fisher's exact test per disease
    results = []
    for disease, bg_mets in disease_metabolites.items():
        K = len(bg_mets)  # disease size in background
        k = len(query_annotated & bg_mets)  # overlap

        if k == 0:
            continue

        # 2x2 contingency table
        a = k               # query & disease
        b = n - k            # query & not disease
        c = K - k            # not query & disease
        d = N - n - K + k    # not query & not disease
        if d < 0:
            d = 0

        _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

        # Fold enrichment
        expected = (n * K) / N if N > 0 else 0
        fold = k / expected if expected > 0 else 0

        meta = disease_meta.get(disease, {})
        results.append({
            'Disease_Name': disease,
            'Disease_ID': meta.get('Disease_ID', ''),
            'Category': meta.get('Category', ''),
            'Fold_Enrichment': round(fold, 2),
            'P_value': p_value,
            'Metabolite_Count': k,
            'Background_Count': K,
        })

    if not results:
        return pd.DataFrame()

    df_results = pd.DataFrame(results)

    # BH FDR correction
    reject, fdr_values, _, _ = multipletests(df_results['P_value'], method='fdr_bh')
    df_results['FDR'] = fdr_values
    df_results['Significant'] = reject

    # Sort by FDR
    df_results = df_results.sort_values('FDR').reset_index(drop=True)

    return df_results


def run_microbe_enrichment(metabolite_ids, fdr_threshold=0.05):
    """
    Run microbiome enrichment analysis on a list of metabolites.

    Uses Fisher's exact test with BH FDR correction against the MMI database.

    Args:
        metabolite_ids: List of HMDB IDs or metabolite names
        fdr_threshold: FDR cutoff for significance

    Returns:
        DataFrame: Microbe_Name, Taxonomy_ID, Organism, Fold_Enrichment, P_value,
                   Metabolite_Count, Background_Count, FDR, Significant
    """
    from app.services.mmi_service import get_mmi_db

    mmi_df = get_mmi_db()
    if mmi_df.empty:
        logger.warning("MMI database is empty — cannot run microbe enrichment")
        return pd.DataFrame()

    # Build microbe → metabolite sets from MMI database (background)
    microbe_metabolites = {}
    microbe_meta = {}
    all_mmi_metabolites = set()

    for _, row in mmi_df.iterrows():
        microbe = row.get("Microbe_Name", "")
        hmdb = row.get("HMDB_ID", "")
        met_name = row.get("Metabolite_Name", "")
        if not microbe:
            continue

        met_key = hmdb if hmdb else met_name
        if met_key:
            microbe_metabolites.setdefault(microbe, set()).add(met_key)
            all_mmi_metabolites.add(met_key)
            if microbe not in microbe_meta:
                microbe_meta[microbe] = {
                    "Taxonomy_ID": row.get("Taxonomy_ID", ""),
                    "Organism": row.get("Organism", ""),
                    "Rank": row.get("Rank", ""),
                }

    N = len(all_mmi_metabolites)
    if N == 0:
        return pd.DataFrame()

    # Normalize query IDs
    query_ids = set()
    for mid in metabolite_ids:
        mid_str = str(mid).strip()
        if mid_str:
            query_ids.add(mid_str)
            name_match = mmi_df.loc[mmi_df["Metabolite_Name"].str.lower() == mid_str.lower(), "HMDB_ID"]
            for hmdb in name_match:
                if hmdb:
                    query_ids.add(hmdb)
            hmdb_match = mmi_df.loc[mmi_df["HMDB_ID"] == mid_str, "Metabolite_Name"]
            for name in hmdb_match:
                if name:
                    query_ids.add(name)

    query_annotated = query_ids & all_mmi_metabolites
    n = len(query_annotated)

    if n == 0:
        logger.info(f"No query metabolites found in MMI background (query={len(query_ids)})")
        return pd.DataFrame()

    # Fisher's exact test per microbe
    results = []
    for microbe, bg_mets in microbe_metabolites.items():
        K = len(bg_mets)
        k = len(query_annotated & bg_mets)

        if k == 0:
            continue

        a = k
        b = n - k
        c = K - k
        d = N - n - K + k
        if d < 0:
            d = 0

        _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

        expected = (n * K) / N if N > 0 else 0
        fold = k / expected if expected > 0 else 0

        meta = microbe_meta.get(microbe, {})
        results.append({
            'Microbe_Name': microbe,
            'Taxonomy_ID': meta.get('Taxonomy_ID', ''),
            'Organism': meta.get('Organism', ''),
            'Rank': meta.get('Rank', ''),
            'Fold_Enrichment': round(fold, 2),
            'P_value': p_value,
            'Metabolite_Count': k,
            'Background_Count': K,
        })

    if not results:
        return pd.DataFrame()

    df_results = pd.DataFrame(results)

    reject, fdr_values, _, _ = multipletests(df_results['P_value'], method='fdr_bh')
    df_results['FDR'] = fdr_values
    df_results['Significant'] = reject

    df_results = df_results.sort_values('FDR').reset_index(drop=True)

    return df_results


def run_drug_enrichment(metabolite_ids, fdr_threshold=0.05):
    """
    Run drug enrichment analysis on a list of metabolites.

    Uses Fisher's exact test with BH FDR correction against the MDrI database.

    Args:
        metabolite_ids: List of HMDB IDs or metabolite names
        fdr_threshold: FDR cutoff for significance

    Returns:
        DataFrame: Drug_Name, DrugBank_ID, Interaction_Type, Fold_Enrichment,
                   P_value, Metabolite_Count, Background_Count, FDR, Significant
    """
    from app.services.mdri_service import get_mdri_db

    mdri_df = get_mdri_db()
    if mdri_df.empty:
        logger.warning("MDrI database is empty — cannot run drug enrichment")
        return pd.DataFrame()

    # Build drug → metabolite sets from MDrI database (background)
    drug_metabolites = {}
    drug_meta = {}
    all_mdri_metabolites = set()

    for _, row in mdri_df.iterrows():
        drug = row.get("Drug_Name", "")
        hmdb = row.get("HMDB_ID", "")
        met_name = row.get("Metabolite_Name", "")
        if not drug:
            continue

        met_key = hmdb if hmdb else met_name
        if met_key:
            drug_metabolites.setdefault(drug, set()).add(met_key)
            all_mdri_metabolites.add(met_key)
            if drug not in drug_meta:
                drug_meta[drug] = {
                    "DrugBank_ID": row.get("DrugBank_ID", ""),
                    "Interaction_Type": row.get("Interaction_Type", ""),
                }

    N = len(all_mdri_metabolites)
    if N == 0:
        return pd.DataFrame()

    # Normalize query IDs
    query_ids = set()
    for mid in metabolite_ids:
        mid_str = str(mid).strip()
        if mid_str:
            query_ids.add(mid_str)
            name_match = mdri_df.loc[
                mdri_df["Metabolite_Name"].str.lower() == mid_str.lower(), "HMDB_ID"
            ]
            for hmdb in name_match:
                if hmdb:
                    query_ids.add(hmdb)
            hmdb_match = mdri_df.loc[mdri_df["HMDB_ID"] == mid_str, "Metabolite_Name"]
            for name in hmdb_match:
                if name:
                    query_ids.add(name)

    query_annotated = query_ids & all_mdri_metabolites
    n = len(query_annotated)

    if n == 0:
        logger.info(f"No query metabolites found in MDrI background (query={len(query_ids)})")
        return pd.DataFrame()

    # Fisher's exact test per drug
    results = []
    for drug, bg_mets in drug_metabolites.items():
        K = len(bg_mets)
        k = len(query_annotated & bg_mets)

        if k == 0:
            continue

        a = k
        b = n - k
        c = K - k
        d = N - n - K + k
        if d < 0:
            d = 0

        _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

        expected = (n * K) / N if N > 0 else 0
        fold = k / expected if expected > 0 else 0

        meta = drug_meta.get(drug, {})
        results.append({
            'Drug_Name': drug,
            'DrugBank_ID': meta.get('DrugBank_ID', ''),
            'Interaction_Type': meta.get('Interaction_Type', ''),
            'Fold_Enrichment': round(fold, 2),
            'P_value': p_value,
            'Metabolite_Count': k,
            'Background_Count': K,
        })

    if not results:
        return pd.DataFrame()

    df_results = pd.DataFrame(results)

    reject, fdr_values, _, _ = multipletests(df_results['P_value'], method='fdr_bh')
    df_results['FDR'] = fdr_values
    df_results['Significant'] = reject

    df_results = df_results.sort_values('FDR').reset_index(drop=True)

    return df_results


def run_gene_enrichment(metabolite_ids, fdr_threshold=0.05):
    """
    Run gene enrichment analysis on a list of metabolites.

    Uses Fisher's exact test with BH FDR correction against the MGI database.

    Args:
        metabolite_ids: List of HMDB IDs or metabolite names
        fdr_threshold: FDR cutoff for significance

    Returns:
        DataFrame: Gene_Symbol, Gene_ID, Organism, Interaction_Type,
                   Fold_Enrichment, P_value, Metabolite_Count,
                   Background_Count, FDR, Significant
    """
    from app.services.mgi_service import get_mgi_db

    mgi_df = get_mgi_db()
    if mgi_df.empty:
        logger.warning("MGI database is empty — cannot run gene enrichment")
        return pd.DataFrame()

    # Build gene → metabolite sets from MGI database (background)
    gene_metabolites = {}
    gene_meta = {}
    all_mgi_metabolites = set()

    for _, row in mgi_df.iterrows():
        gene = row.get("Gene_Symbol", "")
        hmdb = row.get("HMDB_ID", "")
        met_name = row.get("Metabolite_Name", "")
        if not gene:
            continue

        met_key = hmdb if hmdb else met_name
        if met_key:
            gene_metabolites.setdefault(gene, set()).add(met_key)
            all_mgi_metabolites.add(met_key)
            if gene not in gene_meta:
                gene_meta[gene] = {
                    "Gene_ID": row.get("Gene_ID", ""),
                    "Organism": row.get("Organism", ""),
                    "Interaction_Type": row.get("Interaction_Type", ""),
                }

    N = len(all_mgi_metabolites)
    if N == 0:
        return pd.DataFrame()

    # Normalize query IDs
    query_ids = set()
    for mid in metabolite_ids:
        mid_str = str(mid).strip()
        if mid_str:
            query_ids.add(mid_str)
            name_match = mgi_df.loc[
                mgi_df["Metabolite_Name"].str.lower() == mid_str.lower(), "HMDB_ID"
            ]
            for hmdb in name_match:
                if hmdb:
                    query_ids.add(hmdb)
            hmdb_match = mgi_df.loc[mgi_df["HMDB_ID"] == mid_str, "Metabolite_Name"]
            for name in hmdb_match:
                if name:
                    query_ids.add(name)

    query_annotated = query_ids & all_mgi_metabolites
    n = len(query_annotated)

    if n == 0:
        logger.info(f"No query metabolites found in MGI background (query={len(query_ids)})")
        return pd.DataFrame()

    # Fisher's exact test per gene
    results = []
    for gene, bg_mets in gene_metabolites.items():
        K = len(bg_mets)
        k = len(query_annotated & bg_mets)

        if k == 0:
            continue

        a = k
        b = n - k
        c = K - k
        d = N - n - K + k
        if d < 0:
            d = 0

        _, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')

        expected = (n * K) / N if N > 0 else 0
        fold = k / expected if expected > 0 else 0

        meta = gene_meta.get(gene, {})
        results.append({
            'Gene_Symbol': gene,
            'Gene_ID': meta.get('Gene_ID', ''),
            'Organism': meta.get('Organism', ''),
            'Interaction_Type': meta.get('Interaction_Type', ''),
            'Fold_Enrichment': round(fold, 2),
            'P_value': p_value,
            'Metabolite_Count': k,
            'Background_Count': K,
        })

    if not results:
        return pd.DataFrame()

    df_results = pd.DataFrame(results)

    reject, fdr_values, _, _ = multipletests(df_results['P_value'], method='fdr_bh')
    df_results['FDR'] = fdr_values
    df_results['Significant'] = reject

    df_results = df_results.sort_values('FDR').reset_index(drop=True)

    return df_results
