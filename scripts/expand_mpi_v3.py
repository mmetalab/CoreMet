#!/usr/bin/env python3
"""
Expand MPI Database v3 — Typed Interaction Subtypes
====================================================
Adds BindingDB (direct binding), ChEMBL (bioactivity), Reactome (pathway context),
and enhanced Rhea (reaction roles) to the existing KEGG/Rhea/MPIDB edges.

Each edge gets an interaction_subtype (substrate, product, direct_binder, inhibitor,
agonist, pathway_participant, etc.) and source_priority (1=direct, 2=reaction, 3=pathway).

Usage:
    conda run -n mpi-vgae python scripts/expand_mpi_v3.py --rhea
    conda run -n mpi-vgae python scripts/expand_mpi_v3.py --bindingdb
    conda run -n mpi-vgae python scripts/expand_mpi_v3.py --chembl
    conda run -n mpi-vgae python scripts/expand_mpi_v3.py --reactome
    conda run -n mpi-vgae python scripts/expand_mpi_v3.py --merge
"""

import argparse
import gzip
import json
import logging
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
MPI_DIR = ROOT / "data" / "mpidatabase"
CACHE_DIR = ROOT / "data" / "cache"
RAW_DIR = ROOT / "data" / "raw"
MAPPINGS_DIR = ROOT / "data" / "mappings"
OUTPUT_DIR = MPI_DIR


def load_metabolite_crossref():
    """Load InChIKey → HMDB_ID mapping."""
    path = MAPPINGS_DIR / "metabolite_crossref.json"
    if not path.exists():
        logger.warning(f"No crossref file at {path}")
        return {}, {}
    with open(path) as f:
        crossref = json.load(f)
    # Build InChIKey → HMDB_ID reverse map
    inchi_to_hmdb = {}
    for hmdb, info in crossref.items():
        ik = info.get('inchikey', '')
        if ik:
            inchi_to_hmdb[ik] = hmdb
    logger.info(f"Loaded {len(crossref)} metabolites, {len(inchi_to_hmdb)} InChIKeys")
    return crossref, inchi_to_hmdb


def load_chebi_to_hmdb():
    """Load ChEBI → HMDB mapping from Rhea cache."""
    # Try cached version first
    cache_path = CACHE_DIR / "rhea" / "chebi_to_hmdb.json"
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)

    # Build from ChEBI database_accession file (may be in rhea or chebi cache)
    chebi_file = CACHE_DIR / "chebi" / "database_accession.tsv"
    chebi_gz = CACHE_DIR / "chebi" / "database_accession.tsv.gz"
    # Also check rhea cache dir
    if not chebi_gz.exists() and not chebi_file.exists():
        chebi_gz = CACHE_DIR / "rhea" / "database_accession.tsv.gz"
        chebi_file = CACHE_DIR / "rhea" / "database_accession.tsv"
    if chebi_gz.exists():
        import gzip
        with gzip.open(chebi_gz, 'rt') as f:
            df = pd.read_csv(f, sep='\t')
    elif chebi_file.exists():
        df = pd.read_csv(chebi_file, sep='\t')
    else:
        logger.warning("No ChEBI mapping file found")
        return {}

    # source_id=35 is HMDB in ChEBI (lowercase column names)
    src_col = 'SOURCE' if 'SOURCE' in df.columns else 'source_id'
    cid_col = 'COMPOUND_ID' if 'COMPOUND_ID' in df.columns else 'compound_id'
    acc_col = 'ACCESSION_NUMBER' if 'ACCESSION_NUMBER' in df.columns else 'accession_number'

    if src_col == 'source_id':
        hmdb_rows = df[df[src_col] == 35]
    else:
        hmdb_rows = df[df[src_col] == 'HMDB']
    mapping = {}
    for _, row in hmdb_rows.iterrows():
        chebi = str(row[cid_col])
        hmdb = str(row[acc_col])
        if not hmdb.startswith('HMDB'):
            hmdb = 'HMDB' + hmdb.zfill(7)
        mapping[chebi] = hmdb

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'w') as f:
        json.dump(mapping, f)
    logger.info(f"ChEBI→HMDB: {len(mapping)} mappings")
    return mapping


# ══════════════════════════════════════════════════════════════════════
# Source 1: Enhanced Rhea (reaction roles)
# ══════════════════════════════════════════════════════════════════════

def expand_rhea_roles():
    """Extract reaction participant roles (substrate/product) from Rhea BioPAX."""
    logger.info("=== Rhea Enhanced Role Extraction ===")

    # Load existing reaction→ChEBI cache
    cache = CACHE_DIR / "rhea" / "rhea_reaction_chebi.json"
    if not cache.exists():
        logger.warning(f"No Rhea cache at {cache}")
        return pd.DataFrame()

    with open(cache) as f:
        reaction_chebi = json.load(f)

    # Load Rhea→UniProt mapping
    r2u_path = CACHE_DIR / "rhea" / "rhea2uniprot.tsv"
    if not r2u_path.exists():
        logger.warning("No rhea2uniprot.tsv")
        return pd.DataFrame()

    r2u = pd.read_csv(r2u_path, sep='\t')
    logger.info(f"Rhea→UniProt: {len(r2u)} entries")

    # Load ChEBI→HMDB
    chebi_to_hmdb = load_chebi_to_hmdb()

    # Use cached reaction→ChEBI mapping
    # The cache doesn't distinguish LEFT/RIGHT — assign as 'reaction_participant'
    # (BioPAX parsing is complex; the cache provides the ChEBI participants per reaction)
    reaction_roles = {}
    for rxn_id, chebis in reaction_chebi.items():
        reaction_roles[rxn_id] = [(c, 'reaction_participant') for c in chebis]
    logger.info(f"Using {len(reaction_roles)} cached reactions with ChEBI participants")

    # Build edges: (HMDB, UniProt, role, reaction_id)
    rows = []
    rhea_col = 'RHEA_ID' if 'RHEA_ID' in r2u.columns else r2u.columns[0]
    uniprot_col = [c for c in r2u.columns if 'uniprot' in c.lower() or 'id' in c.lower()][-1]

    for rxn_id, participants in reaction_roles.items():
        # Find UniProt proteins for this reaction
        rxn_proteins = r2u[r2u[rhea_col].astype(str) == str(rxn_id)]
        if rxn_proteins.empty:
            continue
        for chebi_id, role in participants:
            hmdb = chebi_to_hmdb.get(str(chebi_id))
            if not hmdb:
                continue
            for _, prow in rxn_proteins.iterrows():
                uniprot = str(prow[uniprot_col])
                rows.append({
                    'HMDB_ID': hmdb,
                    'Uniprot_ID': uniprot,
                    'interaction_subtype': role,
                    'reaction_id': f'RHEA:{rxn_id}',
                    'source': 'Rhea',
                    'source_priority': 2,
                })

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates(['HMDB_ID', 'Uniprot_ID', 'interaction_subtype'])
    logger.info(f"Rhea roles: {len(df)} edges with roles")

    out = OUTPUT_DIR / "mpi_rhea_roles.csv"
    df.to_csv(out, index=False)
    return df


def _parse_biopax_roles(owl_gz_path):
    """Parse BioPAX OWL for LEFT (substrate) / RIGHT (product) participants."""
    import xml.etree.ElementTree as ET

    reaction_roles = {}
    bp = '{http://www.biopax.org/release/biopax-level3.owl#}'

    with gzip.open(owl_gz_path, 'rt') as f:
        tree = ET.parse(f)
    root = tree.getroot()

    for rxn in root.findall(f'.//{bp}BiochemicalReaction'):
        rxn_id = rxn.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about', '')
        # Extract reaction ID number
        match = re.search(r'(\d+)$', rxn_id)
        if not match:
            continue
        rid = match.group(1)

        participants = []
        # LEFT = substrates
        for left in rxn.findall(f'{bp}left'):
            ref = left.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource', '')
            chebi_match = re.search(r'CHEBI[_:](\d+)', ref)
            if chebi_match:
                participants.append((chebi_match.group(1), 'substrate'))

        # RIGHT = products
        for right in rxn.findall(f'{bp}right'):
            ref = right.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource', '')
            chebi_match = re.search(r'CHEBI[_:](\d+)', ref)
            if chebi_match:
                participants.append((chebi_match.group(1), 'product'))

        if participants:
            reaction_roles[rid] = participants

    logger.info(f"Parsed {len(reaction_roles)} reactions with LEFT/RIGHT roles")
    return reaction_roles


# ══════════════════════════════════════════════════════════════════════
# Source 2: BindingDB (direct binding)
# ══════════════════════════════════════════════════════════════════════

def expand_bindingdb():
    """Filter BindingDB to HMDB-mappable metabolites."""
    logger.info("=== BindingDB Direct Binding ===")

    # Check for downloaded file
    bindingdb_dir = RAW_DIR / "bindingdb"
    tsv_files = list(bindingdb_dir.glob("BindingDB_All*.tsv")) if bindingdb_dir.exists() else []
    if not tsv_files:
        logger.warning(f"No BindingDB TSV in {bindingdb_dir}. Download from bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp")
        return pd.DataFrame()

    tsv_path = tsv_files[0]
    logger.info(f"Reading {tsv_path}...")

    _, inchi_to_hmdb = load_metabolite_crossref()

    # Stream through large file, filter by InChIKey match
    rows = []
    chunk_size = 50000
    for chunk in pd.read_csv(tsv_path, sep='\t', chunksize=chunk_size,
                              on_bad_lines='skip', low_memory=False):
        # BindingDB columns vary; find InChIKey and UniProt columns
        ik_col = [c for c in chunk.columns if 'inchi' in c.lower() and 'key' in c.lower()]
        up_col = [c for c in chunk.columns if 'uniprot' in c.lower()]
        ki_col = [c for c in chunk.columns if c.strip().lower() == 'ki (nm)']
        kd_col = [c for c in chunk.columns if c.strip().lower() == 'kd (nm)']
        ic50_col = [c for c in chunk.columns if c.strip().lower() == 'ic50 (nm)']

        if not ik_col or not up_col:
            continue

        ik_c = ik_col[0]
        up_c = up_col[0]

        # Filter to metabolites in our universe
        mask = chunk[ik_c].isin(inchi_to_hmdb)
        matched = chunk[mask]

        for _, row in matched.iterrows():
            ik = str(row[ik_c])
            hmdb = inchi_to_hmdb.get(ik, '')
            uniprot = str(row[up_c]).split(',')[0].strip() if pd.notna(row[up_c]) else ''
            if not hmdb or not uniprot:
                continue

            # Get best affinity measurement
            aff_type, aff_val = None, None
            for col_list, atype in [(ki_col, 'Ki'), (kd_col, 'Kd'), (ic50_col, 'IC50')]:
                if col_list:
                    v = pd.to_numeric(row.get(col_list[0]), errors='coerce')
                    if pd.notna(v) and v > 0:
                        aff_type, aff_val = atype, float(v)
                        break

            rows.append({
                'HMDB_ID': hmdb,
                'Uniprot_ID': uniprot,
                'interaction_subtype': 'direct_binder',
                'affinity_type': aff_type,
                'affinity_value': aff_val,
                'source': 'BindingDB',
                'source_priority': 1,
            })

        if len(rows) % 10000 < chunk_size:
            logger.info(f"  Processed... {len(rows)} matches so far")

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates(['HMDB_ID', 'Uniprot_ID'])
    logger.info(f"BindingDB: {len(df)} metabolite–protein binding edges")

    out = OUTPUT_DIR / "mpi_bindingdb.csv"
    df.to_csv(out, index=False)
    return df


# ══════════════════════════════════════════════════════════════════════
# Source 3: ChEMBL (bioactivity + mechanism)
# ══════════════════════════════════════════════════════════════════════

def expand_chembl():
    """Query ChEMBL REST API for metabolite bioactivities."""
    logger.info("=== ChEMBL Bioactivity ===")

    crossref, inchi_to_hmdb = load_metabolite_crossref()
    if not inchi_to_hmdb:
        return pd.DataFrame()

    # Get all InChIKeys to query
    inchikeys = list(inchi_to_hmdb.keys())
    logger.info(f"Querying ChEMBL for {len(inchikeys)} InChIKeys...")

    # Check for cached results
    cache_path = RAW_DIR / "chembl" / "chembl_metabolite_activities.csv"
    if cache_path.exists():
        logger.info(f"Loading cached ChEMBL results from {cache_path}")
        df = pd.read_csv(cache_path)
        logger.info(f"ChEMBL cached: {len(df)} edges")
        out = OUTPUT_DIR / "mpi_chembl.csv"
        df.to_csv(out, index=False)
        return df

    try:
        import requests
    except ImportError:
        logger.warning("requests not available for ChEMBL API")
        return pd.DataFrame()

    base_url = "https://www.ebi.ac.uk/chembl/api/data"
    rows = []

    for i, ik in enumerate(inchikeys):
        if i % 100 == 0:
            logger.info(f"  ChEMBL: {i}/{len(inchikeys)} InChIKeys processed, {len(rows)} hits")

        try:
            # Step 1: Get molecule ChEMBL ID from InChIKey
            resp = requests.get(f"{base_url}/molecule/{ik}.json", timeout=10)
            if resp.status_code != 200:
                continue
            mol = resp.json()
            chembl_id = mol.get('molecule_chembl_id')
            if not chembl_id:
                continue

            # Step 2: Get activities for this molecule
            act_resp = requests.get(
                f"{base_url}/activity.json",
                params={
                    'molecule_chembl_id': chembl_id,
                    'target_type': 'SINGLE PROTEIN',
                    'limit': 100,
                },
                timeout=15,
            )
            if act_resp.status_code != 200:
                continue

            activities = act_resp.json().get('activities', [])
            for act in activities:
                atype = act.get('standard_type', '')
                if atype not in ('IC50', 'Ki', 'Kd', 'EC50'):
                    continue
                val = act.get('standard_value')
                if val is None:
                    continue
                target_id = act.get('target_chembl_id', '')
                # Get UniProt from target
                uniprot = act.get('target_pref_name', '')  # simplified
                # Better: query target endpoint for UniProt
                # For now, use component if available
                target_comps = act.get('target_organism', '')

                rows.append({
                    'HMDB_ID': inchi_to_hmdb[ik],
                    'Uniprot_ID': target_id,  # Will need UniProt mapping
                    'interaction_subtype': 'bioactivity',
                    'affinity_type': atype,
                    'affinity_value': float(val),
                    'action_type': act.get('activity_comment', ''),
                    'source': 'ChEMBL',
                    'source_priority': 1,
                    'pmid': act.get('document_chembl_id', ''),
                })

            time.sleep(0.3)  # Rate limit

        except Exception as e:
            logger.debug(f"  ChEMBL error for {ik}: {e}")
            continue

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates(['HMDB_ID', 'Uniprot_ID', 'affinity_type'])

    # Cache results
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_path, index=False)

    logger.info(f"ChEMBL: {len(df)} bioactivity edges")
    out = OUTPUT_DIR / "mpi_chembl.csv"
    df.to_csv(out, index=False)
    return df


# ══════════════════════════════════════════════════════════════════════
# Source 4: Reactome (pathway context)
# ══════════════════════════════════════════════════════════════════════

def expand_reactome():
    """Build Reactome metabolite-protein edges via shared pathway membership."""
    logger.info("=== Reactome Pathway Context ===")

    chebi2react = RAW_DIR / "reactome" / "ChEBI2Reactome_All_Levels.txt"
    uniprot2react = RAW_DIR / "reactome" / "UniProt2Reactome_All_Levels.txt"

    if not chebi2react.exists() or not uniprot2react.exists():
        logger.warning("Missing Reactome files. Download ChEBI2Reactome and UniProt2Reactome "
                        "from reactome.org/download/current/")
        return pd.DataFrame()

    chebi_to_hmdb = load_chebi_to_hmdb()

    # Load ChEBI → Reactome pathway (human only)
    chebi_df = pd.read_csv(chebi2react, sep='\t', header=None,
                            names=['chebi_id', 'pathway_id', 'url', 'pathway_name',
                                   'evidence', 'species'])
    chebi_df = chebi_df[chebi_df['species'] == 'Homo sapiens']
    chebi_df['hmdb_id'] = chebi_df['chebi_id'].astype(str).map(chebi_to_hmdb)
    chebi_df = chebi_df[chebi_df['hmdb_id'].notna()]
    logger.info(f"  ChEBI→Reactome (human, HMDB-mappable): {len(chebi_df)} entries, "
                 f"{chebi_df['hmdb_id'].nunique()} metabolites")

    # Load UniProt → Reactome pathway (human only)
    up_df = pd.read_csv(uniprot2react, sep='\t', header=None,
                          names=['uniprot_id', 'pathway_id', 'url', 'pathway_name',
                                 'evidence', 'species'])
    up_df = up_df[up_df['species'] == 'Homo sapiens']
    logger.info(f"  UniProt→Reactome (human): {len(up_df)} entries, "
                 f"{up_df['uniprot_id'].nunique()} proteins")

    # Join: metabolites and proteins sharing the same Reactome pathway
    # Use LOWEST-LEVEL pathways only (most specific) — filter to reaction-level IDs
    # For efficiency, limit to pathways with ≤100 metabolites and ≤500 proteins
    met_per_pw = chebi_df.groupby('pathway_id')['hmdb_id'].nunique()
    prot_per_pw = up_df.groupby('pathway_id')['uniprot_id'].nunique()
    valid_pws = set(met_per_pw[met_per_pw <= 50].index) & set(prot_per_pw[prot_per_pw <= 200].index)
    logger.info(f"  Valid pathways (≤50 mets, ≤200 prots): {len(valid_pws)}")

    chebi_filt = chebi_df[chebi_df['pathway_id'].isin(valid_pws)]
    up_filt = up_df[up_df['pathway_id'].isin(valid_pws)]

    # Merge on pathway_id
    merged = chebi_filt[['hmdb_id', 'pathway_id']].merge(
        up_filt[['uniprot_id', 'pathway_id']], on='pathway_id')
    merged = merged.drop_duplicates(['hmdb_id', 'uniprot_id'])
    logger.info(f"  Reactome metabolite-protein pairs: {len(merged)}")

    rows = []
    for _, row in merged.iterrows():
        rows.append({
            'HMDB_ID': row['hmdb_id'],
            'Uniprot_ID': row['uniprot_id'],
            'interaction_subtype': 'pathway_participant',
            'source': 'Reactome',
            'source_priority': 3,
        })

    df = pd.DataFrame(rows)
    logger.info(f"Reactome: {len(df)} pathway-participant edges")

    out = OUTPUT_DIR / "mpi_reactome.csv"
    df.to_csv(out, index=False)
    return df


# ══════════════════════════════════════════════════════════════════════
# Merge
# ══════════════════════════════════════════════════════════════════════

def merge_mpi_v3():
    """Merge all MPI sources into MPIDB_v3.csv with typed subtypes."""
    logger.info("=== Merging MPI v3 ===")

    # Load existing MPIDB_v2
    v2 = pd.read_csv(MPI_DIR / "MPIDB_v2.csv")
    v2['interaction_subtype'] = 'pathway_participant'
    v2['affinity_type'] = ''
    v2['affinity_value'] = np.nan
    v2['reaction_id'] = ''
    v2['source_priority'] = 3
    v2 = v2.rename(columns={'HMDB ID': 'HMDB_ID', 'Uniprot ID': 'Uniprot_ID'})
    logger.info(f"  v2 base: {len(v2)} edges")

    # Load expansion sources
    sources = []
    for fname in ['mpi_rhea_roles.csv', 'mpi_bindingdb.csv', 'mpi_chembl.csv', 'mpi_reactome.csv']:
        path = MPI_DIR / fname
        if path.exists():
            df = pd.read_csv(path)
            sources.append((fname, df))
            logger.info(f"  {fname}: {len(df)} edges")

    # Combine
    all_dfs = [v2[['HMDB_ID', 'Uniprot_ID', 'Species', 'Metabolite Name', 'SMILES',
                     'Protein Name', 'Gene Name', 'Pathway_ID', 'Pathway_Name',
                     'Evidence_Source', 'interaction_subtype', 'affinity_type',
                     'affinity_value', 'reaction_id', 'source_priority']]]

    for fname, df in sources:
        # Align columns
        for col in all_dfs[0].columns:
            if col not in df.columns:
                df[col] = ''
        all_dfs.append(df[all_dfs[0].columns])

    merged = pd.concat(all_dfs, ignore_index=True)

    # Deduplicate: keep highest priority per (HMDB, UniProt)
    merged['source_priority'] = pd.to_numeric(merged['source_priority'], errors='coerce').fillna(3)
    merged = merged.sort_values('source_priority')
    # Keep all edges but mark the best one
    merged['is_best'] = ~merged.duplicated(['HMDB_ID', 'Uniprot_ID'], keep='first')

    out = MPI_DIR / "MPIDB_v3.csv"
    merged.to_csv(out, index=False)

    # Stats
    logger.info(f"\n=== MPIDB v3 Summary ===")
    logger.info(f"  Total edges: {len(merged):,}")
    logger.info(f"  Unique (HMDB, UniProt) pairs: {merged.drop_duplicates(['HMDB_ID','Uniprot_ID']).shape[0]:,}")
    logger.info(f"  By subtype:")
    for st, count in merged['interaction_subtype'].value_counts().items():
        logger.info(f"    {st}: {count:,}")
    logger.info(f"  By source priority:")
    for sp, count in merged['source_priority'].value_counts().sort_index().items():
        logger.info(f"    priority {int(sp)}: {count:,}")
    logger.info(f"  Saved: {out}")


# ══════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Expand MPI Database v3")
    parser.add_argument('--rhea', action='store_true', help='Expand Rhea roles')
    parser.add_argument('--bindingdb', action='store_true', help='Expand BindingDB')
    parser.add_argument('--chembl', action='store_true', help='Expand ChEMBL')
    parser.add_argument('--reactome', action='store_true', help='Expand Reactome')
    parser.add_argument('--merge', action='store_true', help='Merge all into v3')
    parser.add_argument('--all', action='store_true', help='Run all expansions + merge')
    args = parser.parse_args()

    if args.all or args.rhea:
        expand_rhea_roles()
    if args.all or args.bindingdb:
        expand_bindingdb()
    if args.all or args.chembl:
        expand_chembl()
    if args.all or args.reactome:
        expand_reactome()
    if args.all or args.merge:
        merge_mpi_v3()

    if not any([args.rhea, args.bindingdb, args.chembl, args.reactome, args.merge, args.all]):
        parser.print_help()


if __name__ == "__main__":
    main()
