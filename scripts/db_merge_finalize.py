#!/usr/bin/env python3
"""
Database Merge + Pathway Annotation Script
=============================================
Merges MPIDB_May2024.csv (original), MPIDB_v2_kegg.csv, MPIDB_v2_rhea.csv
(and BRENDA/HMDB if they exist). Adds Evidence_Source, Pathway_ID, and
Pathway_Name columns. Deduplicates and outputs MPIDB_v2.csv.

Output: data/mpidatabase/MPIDB_v2.csv

Usage:
    python scripts/db_merge_finalize.py
"""

import json
import logging
import time
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "mpidatabase"
CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "kegg"

SOURCE_FILES = {
    'original': DATA_DIR / "MPIDB_May2024.csv",
    'KEGG': DATA_DIR / "MPIDB_v2_kegg.csv",
    'Rhea': DATA_DIR / "MPIDB_v2_rhea.csv",
    'BRENDA': DATA_DIR / "MPIDB_v2_brenda.csv",
    'HMDB': DATA_DIR / "MPIDB_v2_hmdb.csv",
}

OUTPUT_PATH = DATA_DIR / "MPIDB_v2.csv"

EXPECTED_COLUMNS = [
    'Species', 'Metabolite Name', 'HMDB ID', 'SMILES',
    'Uniprot ID', 'Protein Name', 'Gene Name',
]


def load_enzyme_pathway_cache():
    """Load enzyme→pathway mapping from KEGG cache."""
    cache_file = CACHE_DIR / "enzyme_pathway_links.json"
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)
    return {}


def fetch_pathway_names(pathway_ids, max_retries=3):
    """Fetch pathway names from KEGG for a set of pathway IDs."""
    cache_file = CACHE_DIR / "pathway_names.json"

    cache = {}
    if cache_file.exists():
        with open(cache_file) as f:
            cache = json.load(f)

    to_fetch = [pid for pid in pathway_ids if pid not in cache]
    logger.info(f"Pathway names: {len(cache)} cached, {len(to_fetch)} to fetch")

    # Batch fetch: use KEGG list pathway
    if to_fetch:
        try:
            time.sleep(0.5)
            resp = requests.get("https://rest.kegg.jp/list/pathway", timeout=60)
            if resp.status_code == 200:
                for line in resp.text.strip().split('\n'):
                    if '\t' in line:
                        pid, name = line.split('\t', 1)
                        pid = pid.replace("path:", "")
                        cache[pid] = name.strip()
        except Exception as e:
            logger.warning(f"Failed to fetch pathway list: {e}")

        # Also fetch organism-specific pathways
        for org_code in ['hsa', 'mmu', 'rno', 'eco', 'bta', 'pae', 'ath', 'sce', 'dme', 'cel']:
            try:
                time.sleep(0.5)
                resp = requests.get(f"https://rest.kegg.jp/list/pathway/{org_code}", timeout=60)
                if resp.status_code == 200:
                    for line in resp.text.strip().split('\n'):
                        if '\t' in line:
                            pid, name = line.split('\t', 1)
                            pid = pid.replace("path:", "")
                            cache[pid] = name.strip()
            except Exception as e:
                logger.warning(f"Failed to fetch pathways for {org_code}: {e}")

        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        with open(cache_file, 'w') as f:
            json.dump(cache, f)

    return cache


def build_uniprot_to_ec():
    """Build UniProt→EC mapping from KEGG per-organism caches."""
    result = {}

    # Try per-organism caches (enzyme_genes_*.json + gene_uniprot_*.json)
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

    if result:
        logger.info(f"Built UniProt→EC mapping: {len(result)} proteins with EC annotations")
    else:
        logger.warning("No per-organism KEGG caches found; pathway annotation will be limited")

    return {k: list(v) for k, v in result.items()}


def main():
    logger.info("=" * 60)
    logger.info("Database Merge + Pathway Annotation")
    logger.info("=" * 60)

    # Step 1: Load and tag each source
    dfs = []
    for source_name, path in SOURCE_FILES.items():
        if path.exists():
            df = pd.read_csv(path)
            # Ensure expected columns exist
            for col in EXPECTED_COLUMNS:
                if col not in df.columns:
                    df[col] = ''
            df['Evidence_Source'] = source_name
            dfs.append(df)
            logger.info(f"Loaded {source_name}: {len(df)} rows")
        else:
            logger.info(f"Skipping {source_name}: {path} not found")

    if not dfs:
        logger.error("No source files found!")
        return

    # Step 2: Concatenate
    df_all = pd.concat(dfs, ignore_index=True)
    logger.info(f"Total before dedup: {len(df_all)}")

    # Step 3: Remove exact duplicates (same HMDB ID + UniProt ID + Species)
    before_dedup = len(df_all)
    df_all = df_all.drop_duplicates(subset=['HMDB ID', 'Uniprot ID', 'Species'], keep='first')
    duplicates_removed = before_dedup - len(df_all)
    logger.info(f"Duplicates removed: {duplicates_removed}")
    logger.info(f"After dedup: {len(df_all)}")

    # Step 4: Add pathway annotations
    logger.info("Adding pathway annotations...")

    # Build UniProt→EC→Pathway chain
    uniprot_to_ecs = build_uniprot_to_ec()
    ec_to_pathways = load_enzyme_pathway_cache()

    # Get all pathway IDs for name lookup
    all_pathway_ids = set()
    for ec_id, pathways in ec_to_pathways.items():
        all_pathway_ids.update(pathways)

    pathway_names = fetch_pathway_names(all_pathway_ids)

    # Annotate each row
    pathway_id_col = []
    pathway_name_col = []

    for _, row in df_all.iterrows():
        uid = str(row['Uniprot ID'])
        ec_numbers = uniprot_to_ecs.get(uid, [])

        pw_ids = set()
        for ec in ec_numbers:
            pws = ec_to_pathways.get(ec, [])
            pw_ids.update(pws)

        if pw_ids:
            pw_id_str = '; '.join(sorted(pw_ids)[:5])  # Top 5
            pw_name_str = '; '.join(
                pathway_names.get(pid, pid) for pid in sorted(pw_ids)[:5]
            )
        else:
            pw_id_str = ''
            pw_name_str = ''

        pathway_id_col.append(pw_id_str)
        pathway_name_col.append(pw_name_str)

    df_all['Pathway_ID'] = pathway_id_col
    df_all['Pathway_Name'] = pathway_name_col

    annotated = sum(1 for x in pathway_id_col if x)
    logger.info(f"Pathway annotations added: {annotated}/{len(df_all)} rows have pathways")

    # Step 5: Reorder columns
    final_columns = [
        'Species', 'Metabolite Name', 'HMDB ID', 'SMILES',
        'Uniprot ID', 'Protein Name', 'Gene Name',
        'Pathway_ID', 'Pathway_Name', 'Evidence_Source',
    ]
    df_all = df_all[final_columns]

    # Step 6: Save
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df_all.to_csv(OUTPUT_PATH, index=False)
    logger.info(f"Saved merged database to {OUTPUT_PATH}")

    # Step 7: Quality report
    logger.info("\n" + "=" * 60)
    logger.info("QUALITY REPORT")
    logger.info("=" * 60)
    logger.info(f"Total MPIs: {len(df_all)}")
    logger.info(f"Duplicates removed: {duplicates_removed}")
    logger.info(f"Unique metabolites (HMDB IDs): {df_all['HMDB ID'].nunique()}")
    logger.info(f"Unique proteins (UniProt IDs): {df_all['Uniprot ID'].nunique()}")
    logger.info(f"Rows with pathway annotation: {annotated}")

    logger.info(f"\nPer-organism counts:")
    for species, count in df_all['Species'].value_counts().items():
        logger.info(f"  {species}: {count}")

    logger.info(f"\nPer-source counts:")
    for source, count in df_all['Evidence_Source'].value_counts().items():
        logger.info(f"  {source}: {count}")


if __name__ == "__main__":
    main()
