#!/usr/bin/env python3
"""
Rhea Database Expansion Script
================================
Uses Rhea→UniProt (high-quality curated enzyme annotations) combined with
ChEBI→HMDB mapping to find additional metabolite-protein interactions not
covered by KEGG alone.

Approach:
  1. rhea2uniprot.tsv → Rhea reactions with UniProt IDs (391K entries)
  2. ChEBI flat files → ChEBI compound_id → HMDB ID mapping (19K entries)
  3. Use Rhea rhea-directions.tsv to link master reactions to directional IDs
  4. Download the Rhea biopax.owl which has reaction→participant→ChEBI links
  5. Join: match ChEBI/HMDB metabolites with UniProt proteins per reaction

Output: data/mpidatabase/MPIDB_v2_rhea.csv

Usage:
    python scripts/db_expand_rhea.py
"""

import gzip
import json
import logging
import re
import time
import xml.etree.ElementTree as ET
from io import BytesIO
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "mpidatabase"
EXISTING_DB_PATH = DATA_DIR / "MPIDB_May2024.csv"
OUTPUT_PATH = DATA_DIR / "MPIDB_v2_rhea.csv"
CACHE_DIR = PROJECT_ROOT / "data" / "cache" / "rhea"

SUPPORTED_ORGANISMS = {
    'Homo sapiens', 'Mus musculus', 'Rattus norvegicus',
    'Escherichia coli', 'Bos taurus', 'Pseudomonas aeruginosa',
    'Arabidopsis thaliana', 'Saccharomyces cerevisiae',
    'Drosophila melanogaster', 'Caenorhabditis elegans',
}

UNIPROT_DELAY = 0.2


def download_file(url, dest_path, description="file"):
    if dest_path.exists():
        logger.info(f"Using cached {description}")
        return
    logger.info(f"Downloading {description}...")
    resp = requests.get(url, timeout=300, stream=True)
    resp.raise_for_status()
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, 'wb') as f:
        for chunk in resp.iter_content(chunk_size=8192):
            f.write(chunk)
    logger.info(f"Downloaded ({dest_path.stat().st_size / 1024:.1f} KB)")


def fetch_uniprot_info(uniprot_id, max_retries=3):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    for attempt in range(max_retries):
        try:
            time.sleep(UNIPROT_DELAY)
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                pn = ""
                if 'proteinDescription' in data:
                    pd_data = data['proteinDescription']
                    if 'recommendedName' in pd_data:
                        pn = pd_data['recommendedName'].get('fullName', {}).get('value', '')
                    elif 'submissionNames' in pd_data and pd_data['submissionNames']:
                        pn = pd_data['submissionNames'][0].get('fullName', {}).get('value', '')
                gn = ""
                if 'genes' in data and data['genes']:
                    gd = data['genes'][0]
                    gn = gd.get('geneName', {}).get('value', '') or \
                         (gd.get('orderedLocusNames', [{}])[0].get('value', '') if gd.get('orderedLocusNames') else '')
                org = data.get('organism', {}).get('scientificName', '')
                return pn, gn, org
            elif resp.status_code == 404:
                return None, None, None
        except Exception:
            time.sleep(2)
    return None, None, None


def batch_fetch_uniprot_info(uniprot_ids):
    cache_file = CACHE_DIR / "uniprot_info_rhea.json"
    cache = {}
    if cache_file.exists():
        with open(cache_file) as f:
            cache = json.load(f)
    kegg_cache = PROJECT_ROOT / "data" / "cache" / "kegg" / "uniprot_info.json"
    if kegg_cache.exists():
        with open(kegg_cache) as f:
            for uid, info in json.load(f).items():
                if uid not in cache:
                    cache[uid] = info
    to_fetch = [uid for uid in uniprot_ids if uid not in cache]
    logger.info(f"UniProt info: {len(cache)} cached, {len(to_fetch)} to fetch")
    for i, uid in enumerate(to_fetch):
        if i % 100 == 0 and i > 0:
            logger.info(f"  Fetched {i}/{len(to_fetch)}...")
            with open(cache_file, 'w') as f:
                json.dump(cache, f)
        pn, gn, org = fetch_uniprot_info(uid)
        cache[uid] = {'protein_name': pn or '', 'gene_name': gn or '', 'organism': org or ''}
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_file, 'w') as f:
        json.dump(cache, f)
    return cache


def load_existing_metabolite_info():
    info = {}
    for p in [EXISTING_DB_PATH, DATA_DIR / "MPIDB_v2_kegg.csv"]:
        if p.exists():
            df = pd.read_csv(p)
            for _, row in df.iterrows():
                h = str(row.get('HMDB ID', ''))
                if h and h != 'nan' and h not in info:
                    info[h] = {'name': str(row.get('Metabolite Name', '')), 'smiles': str(row.get('SMILES', ''))}
    logger.info(f"Metabolite info: {len(info)} HMDB IDs")
    return info


def parse_chebi_to_hmdb():
    """Parse ChEBI→HMDB from database_accession.tsv.gz (source_id=35 = HMDB)."""
    chebi_path = CACHE_DIR / "database_accession.tsv.gz"
    download_file("https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files/database_accession.tsv.gz",
                  chebi_path, "ChEBI accessions")
    mapping = {}
    with gzip.open(chebi_path, 'rt') as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6 and parts[5] == '35' and parts[2].startswith('HMDB'):
                hmdb = parts[2].strip()
                if len(hmdb) < 11:
                    hmdb = "HMDB" + hmdb[4:].zfill(7)
                mapping[parts[1]] = hmdb  # compound_id → hmdb
    logger.info(f"ChEBI→HMDB: {len(mapping)} mappings")
    return mapping


def parse_rhea_chebi_from_biopax():
    """Parse Rhea BioPAX OWL to get reaction→ChEBI participant mappings.
    The BioPAX format explicitly links reactions to their small molecule participants.
    """
    cache_file = CACHE_DIR / "rhea_reaction_chebi.json"
    if cache_file.exists():
        with open(cache_file) as f:
            data = json.load(f)
            if data:  # Only use if non-empty
                logger.info(f"Loaded cached reaction→ChEBI: {len(data)} reactions")
                return data

    biopax_path = CACHE_DIR / "rhea-biopax.owl.gz"
    # Try to download BioPAX
    try:
        download_file("https://ftp.expasy.org/databases/rhea/biopax/rhea-biopax.owl.gz",
                      biopax_path, "Rhea BioPAX OWL")
    except Exception as e:
        logger.warning(f"Could not download BioPAX: {e}")
        return {}

    logger.info("Parsing Rhea BioPAX for reaction→ChEBI...")
    # BioPAX uses bp:left and bp:right for reaction participants
    # SmallMolecule entities reference ChEBI via bp:entityReference
    reaction_chebis = {}

    # Parse the OWL/RDF XML
    chebi_pattern = re.compile(r'CHEBI[_:](\d+)')

    # Read line by line to handle large files
    current_reaction_id = None
    current_block_type = None
    current_chebis = set()

    with gzip.open(biopax_path, 'rt', errors='replace') as f:
        for line in f:
            # Detect BiochemicalReaction blocks
            m = re.search(r'<bp:BiochemicalReaction rdf:about=".*?/(\d+)"', line)
            if m:
                if current_reaction_id and current_chebis:
                    reaction_chebis[current_reaction_id] = list(current_chebis)
                current_reaction_id = m.group(1)
                current_chebis = set()
                current_block_type = 'reaction'
                continue

            # Detect SmallMolecule blocks (not reactions)
            if '<bp:SmallMolecule ' in line or '<bp:Protein ' in line:
                if current_reaction_id and current_chebis:
                    reaction_chebis[current_reaction_id] = list(current_chebis)
                current_reaction_id = None
                current_block_type = 'molecule'

            # Collect ChEBI refs within reaction blocks
            if current_reaction_id:
                for cm in chebi_pattern.finditer(line):
                    current_chebis.add(cm.group(1))

        # Final block
        if current_reaction_id and current_chebis:
            reaction_chebis[current_reaction_id] = list(current_chebis)

    logger.info(f"Parsed {len(reaction_chebis)} reactions with ChEBI participants")

    with open(cache_file, 'w') as f:
        json.dump(reaction_chebis, f)
    return reaction_chebis


def main():
    logger.info("=" * 60)
    logger.info("Rhea Database Expansion Pipeline")
    logger.info("=" * 60)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # Step 1: Rhea→UniProt
    r2u_path = CACHE_DIR / "rhea2uniprot.tsv"
    download_file("https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv", r2u_path, "rhea2uniprot")
    df_r2u = pd.read_csv(r2u_path, sep='\t')
    logger.info(f"rhea2uniprot: {len(df_r2u)} rows")

    rhea_to_uniprots = {}
    for _, row in df_r2u.iterrows():
        mid = str(row['MASTER_ID'])
        rhea_to_uniprots.setdefault(mid, set()).add(str(row['ID']))
    # Also index by RHEA_ID (directional IDs)
    for _, row in df_r2u.iterrows():
        rid = str(row['RHEA_ID'])
        rhea_to_uniprots.setdefault(rid, set()).add(str(row['ID']))
    logger.info(f"Rhea→UniProt: {len(rhea_to_uniprots)} reaction IDs")

    # Step 2: ChEBI→HMDB
    chebi_to_hmdb = parse_chebi_to_hmdb()

    # Step 3: Rhea→ChEBI (from BioPAX)
    reaction_chebis = parse_rhea_chebi_from_biopax()

    # Step 4: Build MPI pairs
    mpi_pairs = []
    matched = 0
    for rhea_id, chebis in reaction_chebis.items():
        uniprots = rhea_to_uniprots.get(rhea_id, set())
        if not uniprots:
            continue
        hmdb_ids = [chebi_to_hmdb[cid] for cid in chebis if cid in chebi_to_hmdb]
        if not hmdb_ids:
            continue
        matched += 1
        for hmdb_id in set(hmdb_ids):
            for uid in uniprots:
                mpi_pairs.append({'hmdb_id': hmdb_id, 'uniprot_id': uid})

    logger.info(f"Matched {matched} reactions → {len(mpi_pairs)} raw MPI pairs")

    if not mpi_pairs:
        logger.warning("No MPI pairs found. Saving empty file.")
        pd.DataFrame(columns=['Species', 'Metabolite Name', 'HMDB ID', 'SMILES',
                              'Uniprot ID', 'Protein Name', 'Gene Name']).to_csv(OUTPUT_PATH, index=False)
        return

    df_pairs = pd.DataFrame(mpi_pairs).drop_duplicates()
    logger.info(f"Unique pairs: {len(df_pairs)}")

    # Step 5: UniProt info for organism filtering
    all_uniprots = df_pairs['uniprot_id'].unique().tolist()
    uniprot_cache = batch_fetch_uniprot_info(all_uniprots)

    met_info = load_existing_metabolite_info()

    # Step 6: Filter and build
    rows = []
    for _, pair in df_pairs.iterrows():
        uid = pair['uniprot_id']
        info = uniprot_cache.get(uid, {})
        org = info.get('organism', '')
        if org not in SUPPORTED_ORGANISMS:
            continue
        hmdb_id = pair['hmdb_id']
        mi = met_info.get(hmdb_id, {})
        rows.append({
            'Species': org,
            'Metabolite Name': mi.get('name', ''),
            'HMDB ID': hmdb_id,
            'SMILES': mi.get('smiles', ''),
            'Uniprot ID': uid,
            'Protein Name': info.get('protein_name', ''),
            'Gene Name': info.get('gene_name', ''),
        })

    df_new = pd.DataFrame(rows)
    logger.info(f"After organism filter: {len(df_new)}")
    df_new = df_new.drop_duplicates(subset=['HMDB ID', 'Uniprot ID', 'Species'])
    logger.info(f"After internal dedup: {len(df_new)}")

    # Dedup against existing
    existing_keys = set()
    for p in [EXISTING_DB_PATH, DATA_DIR / "MPIDB_v2_kegg.csv"]:
        if p.exists():
            df_ex = pd.read_csv(p)
            existing_keys.update(zip(df_ex['HMDB ID'], df_ex['Uniprot ID'], df_ex['Species']))
    if existing_keys and len(df_new) > 0:
        df_new['_k'] = list(zip(df_new['HMDB ID'], df_new['Uniprot ID'], df_new['Species']))
        df_new = df_new[~df_new['_k'].isin(existing_keys)].drop(columns=['_k'])
        logger.info(f"After dedup vs existing: {len(df_new)}")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df_new.to_csv(OUTPUT_PATH, index=False)
    logger.info(f"Saved to {OUTPUT_PATH}")

    logger.info("\n" + "=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total new MPIs: {len(df_new)}")
    if len(df_new) > 0:
        for sp, cnt in df_new['Species'].value_counts().items():
            logger.info(f"  {sp}: {cnt}")
        logger.info(f"Unique metabolites: {df_new['HMDB ID'].nunique()}")
        logger.info(f"Unique proteins: {df_new['Uniprot ID'].nunique()}")


if __name__ == "__main__":
    main()
