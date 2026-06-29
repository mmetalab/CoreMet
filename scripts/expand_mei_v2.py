#!/usr/bin/env python3
"""
Expand MEI Database v2 — Reaction Roles + Kinetics
===================================================
Adds enzyme_role (substrate/product/cofactor) from Rhea and kinetic parameters
from SABIO-RK to the existing KEGG-based MEI.

Usage:
    conda run -n mpi-vgae python scripts/expand_mei_v2.py --rhea-roles
    conda run -n mpi-vgae python scripts/expand_mei_v2.py --sabiork
    conda run -n mpi-vgae python scripts/expand_mei_v2.py --merge
"""

import argparse
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"
CACHE_DIR = ROOT / "data" / "cache"
SABIORK_CACHE = CACHE_DIR / "sabiork"
SABIORK_CACHE.mkdir(parents=True, exist_ok=True)


def load_chebi_to_hmdb():
    """Load ChEBI → HMDB mapping."""
    cache_path = CACHE_DIR / "rhea" / "chebi_to_hmdb.json"
    if cache_path.exists():
        with open(cache_path) as f:
            return json.load(f)
    # Build from database_accession
    import gzip
    gz = CACHE_DIR / "rhea" / "database_accession.tsv.gz"
    if not gz.exists():
        return {}
    with gzip.open(gz, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    hmdb_rows = df[df['source_id'] == 35]
    mapping = {}
    for _, row in hmdb_rows.iterrows():
        chebi = str(row['compound_id'])
        hmdb = str(row['accession_number'])
        if not hmdb.startswith('HMDB'):
            hmdb = 'HMDB' + hmdb.zfill(7)
        mapping[chebi] = hmdb
    with open(cache_path, 'w') as f:
        json.dump(mapping, f)
    return mapping


# ══════════════════════════════════════════════════════════════════════
# Rhea Reaction Roles
# ══════════════════════════════════════════════════════════════════════

def expand_rhea_mei_roles():
    """Assign substrate/product roles to MEI edges via Rhea reaction→ChEBI→EC mapping."""
    logger.info("=== Rhea MEI Role Extraction ===")

    # Load Rhea reaction→ChEBI cache
    rc_path = CACHE_DIR / "rhea" / "rhea_reaction_chebi.json"
    with open(rc_path) as f:
        reaction_chebi = json.load(f)  # {reaction_id: [chebi_id, ...]}

    # Load Rhea→EC mapping from xrefs
    xrefs = pd.read_csv(CACHE_DIR / "rhea" / "rhea2xrefs.tsv", sep='\t')
    ec_map = xrefs[xrefs['DB'] == 'EC'][['RHEA_ID', 'MASTER_ID', 'ID', 'DIRECTION']]
    # Build master_id → EC mapping
    master_to_ec = {}
    for _, row in ec_map.iterrows():
        mid = str(row['MASTER_ID'])
        master_to_ec.setdefault(mid, set()).add(row['ID'])
    logger.info(f"Rhea→EC: {len(master_to_ec)} master reactions with EC numbers")

    # Load Rhea→UniProt
    r2u = pd.read_csv(CACHE_DIR / "rhea" / "rhea2uniprot.tsv", sep='\t')
    rhea_col = r2u.columns[0]  # RHEA_ID or MASTER_ID
    up_col = r2u.columns[-1]  # ID (UniProt)
    master_to_uniprot = {}
    for _, row in r2u.iterrows():
        mid = str(row[rhea_col])
        master_to_uniprot.setdefault(mid, set()).add(str(row[up_col]))

    # Load ChEBI→HMDB
    chebi_to_hmdb = load_chebi_to_hmdb()
    logger.info(f"ChEBI→HMDB: {len(chebi_to_hmdb)} mappings")

    # Build role assignments
    # reaction_chebi uses reaction IDs that may be MASTER or directional
    rows = []
    for rxn_id, chebi_list in reaction_chebi.items():
        # Get EC numbers for this reaction
        ecs = master_to_ec.get(rxn_id, set())
        # Get UniProt proteins
        uniprots = master_to_uniprot.get(rxn_id, set())

        if not ecs or not uniprots:
            continue

        for chebi_id in chebi_list:
            hmdb = chebi_to_hmdb.get(str(chebi_id))
            if not hmdb:
                continue
            for ec in ecs:
                for up in uniprots:
                    rows.append({
                        'HMDB_ID': hmdb,
                        'EC_Number': ec,
                        'Uniprot_ID': up,
                        'enzyme_role': 'reaction_participant',  # No LEFT/RIGHT from cache
                        'reaction_id': f'RHEA:{rxn_id}',
                        'source': 'Rhea',
                    })

    df = pd.DataFrame(rows)
    df = df.drop_duplicates(['HMDB_ID', 'EC_Number', 'Uniprot_ID'])
    logger.info(f"Rhea MEI roles: {len(df)} role-annotated edges")

    out = DB_DIR / "mei_rhea_roles.csv"
    df.to_csv(out, index=False)
    return df


# ══════════════════════════════════════════════════════════════════════
# SABIO-RK Kinetics
# ══════════════════════════════════════════════════════════════════════

def expand_sabiork():
    """Query SABIO-RK REST API for kinetic parameters per EC number."""
    logger.info("=== SABIO-RK Kinetic Parameters ===")

    # Get unique EC numbers from MEI
    mei = pd.read_csv(DB_DIR / "mei_database.csv", usecols=['EC_Number'])
    ec_numbers = mei['EC_Number'].dropna().unique()
    logger.info(f"Unique EC numbers to query: {len(ec_numbers)}")

    # Check cache
    cache_file = SABIORK_CACHE / "sabiork_kinetics.csv"
    if cache_file.exists():
        logger.info(f"Loading cached SABIO-RK results from {cache_file}")
        return pd.read_csv(cache_file)

    try:
        import requests
    except ImportError:
        logger.warning("requests not available")
        return pd.DataFrame()

    base_url = "http://sabiork.h-its.org/sabioRestWebServices"
    rows = []

    for i, ec in enumerate(ec_numbers[:500]):  # Limit to first 500 for speed
        if i % 50 == 0:
            logger.info(f"  SABIO-RK: {i}/{min(len(ec_numbers), 500)} ECs, {len(rows)} kinetic entries")

        try:
            # Query kinetic laws by EC number
            resp = requests.get(
                f"{base_url}/kineticLaws/organisms",
                params={'ecnumber': ec, 'format': 'json'},
                timeout=10,
                headers={'Accept': 'application/json'},
            )
            if resp.status_code != 200:
                continue

            data = resp.json() if resp.headers.get('content-type', '').startswith('application/json') else []
            if not isinstance(data, list):
                continue

            for entry in data[:20]:  # Cap per EC
                rows.append({
                    'EC_Number': ec,
                    'param_type': entry.get('type', ''),
                    'param_value': entry.get('startValue', ''),
                    'param_unit': entry.get('unit', ''),
                    'organism': entry.get('organism', ''),
                    'substrate': entry.get('substrate', ''),
                    'source': 'SABIO-RK',
                })

            time.sleep(0.5)  # Rate limit

        except Exception as e:
            logger.debug(f"  SABIO-RK error for {ec}: {e}")
            continue

    df = pd.DataFrame(rows)
    if len(df):
        df.to_csv(cache_file, index=False)
    logger.info(f"SABIO-RK: {len(df)} kinetic entries")
    return df


# ══════════════════════════════════════════════════════════════════════
# Merge
# ══════════════════════════════════════════════════════════════════════

def merge_mei_v2():
    """Merge MEI base + Rhea roles + SABIO-RK kinetics."""
    logger.info("=== Merging MEI v2 ===")

    mei = pd.read_csv(DB_DIR / "mei_database.csv", low_memory=False)
    mei['enzyme_role'] = ''
    mei['reaction_id'] = ''
    mei['param_type'] = ''
    mei['param_value'] = np.nan
    mei['param_unit'] = ''
    logger.info(f"MEI base: {len(mei):,} edges")

    # Add Rhea role edges as new entries (they complement, not overlap, KEGG MEI)
    rhea_path = DB_DIR / "mei_rhea_roles.csv"
    if rhea_path.exists():
        rhea = pd.read_csv(rhea_path)
        logger.info(f"Rhea roles: {len(rhea):,} new entries to append")
        # Align columns
        for col in mei.columns:
            if col not in rhea.columns:
                rhea[col] = ''
        rhea['Evidence_Source'] = 'Rhea'
        mei = pd.concat([mei, rhea[mei.columns]], ignore_index=True)
        mei = mei.drop_duplicates(['HMDB_ID', 'EC_Number', 'Uniprot_ID'], keep='last')

    # Summary
    role_filled = (mei['enzyme_role'] != '').sum()
    logger.info(f"\nMEI v2 Summary:")
    logger.info(f"  Total: {len(mei):,}")
    logger.info(f"  With enzyme_role: {role_filled:,} ({100*role_filled/len(mei):.1f}%)")
    if role_filled:
        logger.info(f"  Roles: {dict(mei[mei['enzyme_role']!='']['enzyme_role'].value_counts())}")

    out = DB_DIR / "mei_database_v2_enriched.csv"
    mei.to_csv(out, index=False)
    logger.info(f"  Saved: {out}")


def main():
    parser = argparse.ArgumentParser(description="Expand MEI Database v2")
    parser.add_argument('--rhea-roles', action='store_true')
    parser.add_argument('--sabiork', action='store_true')
    parser.add_argument('--merge', action='store_true')
    parser.add_argument('--all', action='store_true')
    args = parser.parse_args()

    if args.all or args.rhea_roles:
        expand_rhea_mei_roles()
    if args.all or args.sabiork:
        expand_sabiork()
    if args.all or args.merge:
        merge_mei_v2()

    if not any([args.rhea_roles, args.sabiork, args.merge, args.all]):
        parser.print_help()


if __name__ == "__main__":
    main()
