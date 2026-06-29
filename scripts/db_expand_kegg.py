#!/usr/bin/env python3
"""
KEGG Database Expansion Script
================================
Fetches enzyme-compound links per organism, maps to UniProt via
organism-specific KEGG gene→UniProt conversion, converts compound
structures to SMILES, and deduplicates against existing MPIDB.

Output: data/mpidatabase/MPIDB_v2_kegg.csv

Usage:
    python scripts/db_expand_kegg.py
"""

import os
import sys
import time
import json
import logging
from pathlib import Path

import pandas as pd
import requests
from rdkit import Chem

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "mpidatabase"
EXISTING_DB_PATH = DATA_DIR / "MPIDB_May2024.csv"
OUTPUT_PATH = DATA_DIR / "MPIDB_v2_kegg.csv"
CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "kegg"

# Supported organisms: KEGG code → species name
SUPPORTED_ORGANISMS = {
    'hsa': 'Homo sapiens',
    'mmu': 'Mus musculus',
    'rno': 'Rattus norvegicus',
    'eco': 'Escherichia coli',
    'bta': 'Bos taurus',
    'pae': 'Pseudomonas aeruginosa',
    'ath': 'Arabidopsis thaliana',
    'sce': 'Saccharomyces cerevisiae',
    'dme': 'Drosophila melanogaster',
    'cel': 'Caenorhabditis elegans',
}

KEGG_DELAY = 0.5
UNIPROT_DELAY = 0.2


def kegg_get(endpoint, max_retries=3):
    """Fetch from KEGG REST API with rate limiting."""
    url = f"https://rest.kegg.jp/{endpoint}"
    for attempt in range(max_retries):
        try:
            time.sleep(KEGG_DELAY)
            resp = requests.get(url, timeout=60)
            if resp.status_code == 200:
                return resp.text
            elif resp.status_code == 403:
                logger.warning(f"Rate limited on {endpoint}, waiting 10s...")
                time.sleep(10)
            else:
                logger.warning(f"KEGG {endpoint} returned {resp.status_code} (attempt {attempt+1})")
                return None
        except requests.RequestException as e:
            logger.warning(f"KEGG request failed: {e} (attempt {attempt+1})")
            time.sleep(5)
    return None


def parse_kegg_link(text):
    """Parse KEGG link/conv API output (tab-separated pairs)."""
    pairs = []
    if not text:
        return pairs
    for line in text.strip().split('\n'):
        if '\t' in line:
            a, b = line.strip().split('\t')
            pairs.append((a.strip(), b.strip()))
    return pairs


def load_cached_or_fetch(cache_name, fetch_func):
    """Load from cache or fetch and cache."""
    cache_file = CACHE_DIR / f"{cache_name}.json"
    if cache_file.exists():
        logger.info(f"Loading cached {cache_name}...")
        with open(cache_file) as f:
            return json.load(f)
    logger.info(f"Fetching {cache_name}...")
    data = fetch_func()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_file, 'w') as f:
        json.dump(data, f)
    return data


def mol_to_smiles(mol_text):
    """Convert KEGG MOL block to SMILES."""
    try:
        mol = Chem.MolFromMolBlock(mol_text)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception:
        pass
    return None


def fetch_enzyme_compound_links():
    """Fetch all enzyme→compound links."""
    text = kegg_get("link/compound/enzyme")
    pairs = parse_kegg_link(text)
    result = {}
    for ec, cpd in pairs:
        ec_id = ec.replace("ec:", "")
        cpd_id = cpd.replace("cpd:", "")
        result.setdefault(ec_id, []).append(cpd_id)
    logger.info(f"Fetched {len(pairs)} enzyme-compound links ({len(result)} enzymes)")
    return result


def fetch_enzyme_pathway_links():
    """Fetch enzyme→pathway links."""
    text = kegg_get("link/pathway/enzyme")
    pairs = parse_kegg_link(text)
    result = {}
    for ec, pw in pairs:
        ec_id = ec.replace("ec:", "")
        pw_id = pw.replace("path:", "")
        result.setdefault(ec_id, []).append(pw_id)
    logger.info(f"Fetched {len(pairs)} enzyme-pathway links")
    return result


def fetch_org_enzyme_genes(org_code):
    """Fetch enzyme→gene mapping for a specific organism.
    Uses: link/{org}/enzyme  →  ec:X.X.X.X  {org}:geneID
    """
    text = kegg_get(f"link/{org_code}/enzyme")
    pairs = parse_kegg_link(text)
    result = {}
    for ec, gene in pairs:
        ec_id = ec.replace("ec:", "")
        result.setdefault(ec_id, []).append(gene)
    return result


def fetch_org_gene_uniprot(org_code):
    """Fetch gene→UniProt mapping for a specific organism.
    Uses: conv/uniprot/{org}  →  {org}:geneID  up:UniProtID
    """
    text = kegg_get(f"conv/uniprot/{org_code}")
    pairs = parse_kegg_link(text)
    result = {}
    for gene, up in pairs:
        up_id = up.replace("up:", "")
        result[gene] = up_id
    return result


def fetch_uniprot_info(uniprot_id, max_retries=3):
    """Fetch protein name and gene name from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    for attempt in range(max_retries):
        try:
            time.sleep(UNIPROT_DELAY)
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                protein_name = ""
                if 'proteinDescription' in data:
                    pd_data = data['proteinDescription']
                    if 'recommendedName' in pd_data:
                        protein_name = pd_data['recommendedName'].get('fullName', {}).get('value', '')
                    elif 'submissionNames' in pd_data and pd_data['submissionNames']:
                        protein_name = pd_data['submissionNames'][0].get('fullName', {}).get('value', '')
                gene_name = ""
                if 'genes' in data and data['genes']:
                    gd = data['genes'][0]
                    gene_name = gd.get('geneName', {}).get('value', '') or \
                                (gd.get('orderedLocusNames', [{}])[0].get('value', '') if gd.get('orderedLocusNames') else '')
                return protein_name, gene_name
            elif resp.status_code == 404:
                return None, None
        except Exception:
            time.sleep(2)
    return None, None


def batch_fetch_uniprot_info(uniprot_ids):
    """Fetch UniProt info with persistent cache."""
    cache_file = CACHE_DIR / "uniprot_info.json"
    cache = {}
    if cache_file.exists():
        with open(cache_file) as f:
            cache = json.load(f)

    to_fetch = [uid for uid in uniprot_ids if uid not in cache]
    logger.info(f"UniProt info: {len(cache)} cached, {len(to_fetch)} to fetch")

    for i, uid in enumerate(to_fetch):
        if i % 100 == 0 and i > 0:
            logger.info(f"  Fetched {i}/{len(to_fetch)} UniProt entries...")
            with open(cache_file, 'w') as f:
                json.dump(cache, f)
        pname, gname = fetch_uniprot_info(uid)
        cache[uid] = {'protein_name': pname or '', 'gene_name': gname or ''}

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_file, 'w') as f:
        json.dump(cache, f)
    return cache


def fetch_compound_info_batch(compound_ids):
    """Fetch name + SMILES for compounds from KEGG."""
    cache_file = CACHE_DIR / "compound_info.json"
    cache = {}
    if cache_file.exists():
        with open(cache_file) as f:
            cache = json.load(f)

    to_fetch = [cid for cid in compound_ids if cid not in cache]
    logger.info(f"Compound info: {len(cache)} cached, {len(to_fetch)} to fetch")

    for i, cid in enumerate(to_fetch):
        if i % 50 == 0 and i > 0:
            logger.info(f"  Fetched {i}/{len(to_fetch)} compounds...")
            with open(cache_file, 'w') as f:
                json.dump(cache, f)

        # Fetch MOL for SMILES
        mol_text = kegg_get(f"get/{cid}/mol")
        smiles = mol_to_smiles(mol_text) if mol_text else ''

        # Fetch name
        text = kegg_get(f"get/{cid}")
        name = ''
        if text:
            for line in text.split('\n'):
                if line.startswith("NAME"):
                    name = line.replace("NAME", "").strip().rstrip(';').strip()
                    break

        # Fetch HMDB cross-reference from compound entry
        hmdb_id = ''
        if text:
            in_dblinks = False
            for line in text.split('\n'):
                if line.startswith("DBLINKS"):
                    in_dblinks = True
                    line = line.replace("DBLINKS", "").strip()
                if in_dblinks:
                    if line.strip().startswith("HMDB:"):
                        hmdb_id = line.strip().replace("HMDB:", "").strip().split()[0]
                        # Normalize to 7-digit
                        if hmdb_id.startswith("HMDB") and len(hmdb_id) < 11:
                            hmdb_id = "HMDB" + hmdb_id[4:].zfill(7)
                        break
                    elif not line.startswith(" ") and not line.startswith("DBLINKS"):
                        in_dblinks = False

        cache[cid] = {'smiles': smiles, 'name': name, 'hmdb_id': hmdb_id}

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_file, 'w') as f:
        json.dump(cache, f)
    return cache


def main():
    logger.info("=" * 60)
    logger.info("KEGG Database Expansion Pipeline")
    logger.info("=" * 60)

    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Fetch all enzyme→compound links (global)
    ec_compounds = load_cached_or_fetch("enzyme_compound_links", fetch_enzyme_compound_links)

    # Step 2: Fetch enzyme→pathway links (for later use)
    ec_pathways = load_cached_or_fetch("enzyme_pathway_links", fetch_enzyme_pathway_links)

    # Step 3: For each organism, get enzyme→gene→UniProt mapping
    all_mpi_pairs = []

    for org_code, species_name in SUPPORTED_ORGANISMS.items():
        logger.info(f"\nProcessing {species_name} ({org_code})...")

        # Get enzyme→genes for this organism
        cache_name = f"enzyme_genes_{org_code}"
        ec_genes = load_cached_or_fetch(cache_name, lambda oc=org_code: fetch_org_enzyme_genes(oc))
        logger.info(f"  {len(ec_genes)} enzymes have genes in {org_code}")

        # Get gene→UniProt for this organism
        cache_name = f"gene_uniprot_{org_code}"
        gene_uniprot = load_cached_or_fetch(cache_name, lambda oc=org_code: fetch_org_gene_uniprot(oc))
        logger.info(f"  {len(gene_uniprot)} gene→UniProt mappings")

        # Build MPI pairs: enzyme → compound × uniprot
        org_pairs = 0
        for ec_id, compounds in ec_compounds.items():
            genes = ec_genes.get(ec_id, [])
            for gene in genes:
                uniprot_id = gene_uniprot.get(gene, '')
                if not uniprot_id:
                    continue
                for cpd_id in compounds:
                    all_mpi_pairs.append({
                        'species': species_name,
                        'ec_number': ec_id,
                        'compound_id': cpd_id,
                        'uniprot_id': uniprot_id,
                    })
                    org_pairs += 1

        logger.info(f"  {org_pairs} MPI pairs for {species_name}")

    logger.info(f"\nTotal raw MPI pairs: {len(all_mpi_pairs)}")

    if len(all_mpi_pairs) == 0:
        logger.error("No MPI pairs found! Check KEGG API connectivity.")
        return

    # Step 4: Deduplicate pairs
    df_pairs = pd.DataFrame(all_mpi_pairs)
    df_pairs = df_pairs.drop_duplicates(subset=['compound_id', 'uniprot_id', 'species'])
    logger.info(f"After pair dedup: {len(df_pairs)}")

    # Step 5: Fetch compound info (SMILES, names, HMDB cross-refs)
    all_compounds = df_pairs['compound_id'].unique().tolist()
    compound_cache = fetch_compound_info_batch(all_compounds)

    # Step 6: Fetch UniProt info (protein names, gene names)
    all_uniprots = df_pairs['uniprot_id'].unique().tolist()
    uniprot_cache = batch_fetch_uniprot_info(all_uniprots)

    # Step 7: Build final DataFrame
    rows = []
    for _, pair in df_pairs.iterrows():
        cid = pair['compound_id']
        uid = pair['uniprot_id']
        cinfo = compound_cache.get(cid, {})
        uinfo = uniprot_cache.get(uid, {})

        hmdb_id = cinfo.get('hmdb_id', '')
        smiles = cinfo.get('smiles', '')

        # Skip entries without either HMDB ID or SMILES
        if not hmdb_id and not smiles:
            continue

        rows.append({
            'Species': pair['species'],
            'Metabolite Name': cinfo.get('name', ''),
            'HMDB ID': hmdb_id,
            'SMILES': smiles,
            'Uniprot ID': uid,
            'Protein Name': uinfo.get('protein_name', ''),
            'Gene Name': uinfo.get('gene_name', ''),
        })

    df_new = pd.DataFrame(rows)
    logger.info(f"Entries with HMDB or SMILES: {len(df_new)}")

    if df_new.empty:
        logger.warning("No valid entries found. Saving empty file.")
        df_new = pd.DataFrame(columns=['Species', 'Metabolite Name', 'HMDB ID', 'SMILES', 'Uniprot ID', 'Protein Name', 'Gene Name'])
        df_new.to_csv(OUTPUT_PATH, index=False)
        return

    # Step 8: Deduplicate within new data
    df_new = df_new.drop_duplicates(subset=['HMDB ID', 'Uniprot ID', 'Species'])
    logger.info(f"After internal dedup: {len(df_new)}")

    # Step 9: Deduplicate against existing database
    if EXISTING_DB_PATH.exists():
        df_existing = pd.read_csv(EXISTING_DB_PATH)
        existing_keys = set(
            zip(df_existing['HMDB ID'], df_existing['Uniprot ID'], df_existing['Species'])
        )
        df_new['_key'] = list(zip(df_new['HMDB ID'], df_new['Uniprot ID'], df_new['Species']))
        df_new = df_new[~df_new['_key'].isin(existing_keys)]
        df_new = df_new.drop(columns=['_key'])
        logger.info(f"After dedup against existing DB: {len(df_new)} new MPIs")

    # Step 10: Save
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df_new.to_csv(OUTPUT_PATH, index=False)
    logger.info(f"Saved to {OUTPUT_PATH}")

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total new MPIs added: {len(df_new)}")
    if len(df_new) > 0:
        logger.info(f"\nPer-organism counts:")
        for species, count in df_new['Species'].value_counts().items():
            logger.info(f"  {species}: {count}")
        logger.info(f"\nUnique metabolites (HMDB): {df_new['HMDB ID'].nunique()}")
        logger.info(f"Unique proteins (UniProt): {df_new['Uniprot ID'].nunique()}")


if __name__ == "__main__":
    main()
