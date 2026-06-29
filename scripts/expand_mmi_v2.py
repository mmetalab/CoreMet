#!/usr/bin/env python3
"""
Expand MMI Database — Role-Aware Microbiome Metabolism
======================================================
Adds relation_subtype (producer/consumer/bidirectional/causal/correlative) and
evidence_origin (experimental/computational/database_curated) to MMI edges.
Attempts to add MiMeDB 2.0 for broader coverage.

Usage:
    conda run -n mpi-vgae python scripts/expand_mmi_v2.py --refine
    conda run -n mpi-vgae python scripts/expand_mmi_v2.py --mimedb
    conda run -n mpi-vgae python scripts/expand_mmi_v2.py --merge
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"
RAW_DIR = ROOT / "data" / "raw"


def refine_mmi_subtypes():
    """Add relation_subtype and evidence_origin to existing MMI edges."""
    logger.info("=== Refining MMI Subtypes ===")

    mmi = pd.read_csv(DB_DIR / "mmi_database.csv", low_memory=False)
    logger.info(f"MMI base: {len(mmi):,} edges")

    # Map Relationship_Type → relation_subtype
    mmi['relation_subtype'] = mmi['Relationship_Type'].map({
        'exchange': 'bidirectional',  # AGORA2: all exchange reactions have lb=-1000, ub=1000
        'causal': 'causal',
        'correlative': 'correlative',
    }).fillna('associated')

    # Map Source → evidence_origin
    mmi['evidence_origin'] = mmi['Source'].map({
        'gutMGene': 'experimental',
        'AGORA2': 'computational',
    }).fillna('database_curated')

    # Note: AGORA2 exchange reactions are bidirectional by default (SBML bounds ±1000).
    # True producer/consumer classification requires FBA simulation, not static parsing.
    # We mark them as 'bidirectional' to be honest about the evidence level.

    logger.info(f"Subtypes: {dict(mmi['relation_subtype'].value_counts())}")
    logger.info(f"Evidence: {dict(mmi['evidence_origin'].value_counts())}")

    out = DB_DIR / "mmi_refined.csv"
    mmi.to_csv(out, index=False)
    logger.info(f"Saved: {out}")
    return mmi


def expand_mimedb():
    """Parse MiMeDB download for new microbe-metabolite edges."""
    logger.info("=== MiMeDB Expansion ===")

    mimedb_dir = RAW_DIR / "mimedb"
    files = list(mimedb_dir.glob("*.csv")) + list(mimedb_dir.glob("*.tsv")) \
        if mimedb_dir.exists() else []

    if not files:
        logger.warning(f"No MiMeDB files in {mimedb_dir}. "
                        "Download from https://mimedb.org/downloads")
        # Try wget
        try:
            import subprocess
            mimedb_dir.mkdir(parents=True, exist_ok=True)
            # MiMeDB compounds download
            logger.info("Attempting MiMeDB compounds download...")
            result = subprocess.run(
                ['curl', '-sL', 'https://mimedb.org/system/downloads/current/compounds.csv.zip',
                 '-o', str(mimedb_dir / 'compounds.csv.zip')],
                timeout=60, capture_output=True)
            if (mimedb_dir / 'compounds.csv.zip').exists():
                import zipfile
                with zipfile.ZipFile(mimedb_dir / 'compounds.csv.zip', 'r') as z:
                    z.extractall(mimedb_dir)
                files = list(mimedb_dir.glob("*.csv"))
                logger.info(f"Downloaded and extracted: {[f.name for f in files]}")
            else:
                logger.warning("Download failed — MiMeDB may require manual download")
                return pd.DataFrame()
        except Exception as e:
            logger.warning(f"Download attempt failed: {e}")
            return pd.DataFrame()

    # Load HMDB synonym table for name matching
    syn_path = ROOT / "data" / "mappings" / "hmdb_synonyms.csv"
    name_to_hmdb = {}
    if syn_path.exists():
        syn = pd.read_csv(syn_path, usecols=['hmdb_id', 'name'])
        for _, row in syn.drop_duplicates('name').iterrows():
            name_to_hmdb[str(row['name']).lower().strip()] = str(row['hmdb_id'])
        logger.info(f"Name→HMDB mappings: {len(name_to_hmdb):,}")

    rows = []
    for fpath in files:
        logger.info(f"  Parsing {fpath.name}...")
        sep = '\t' if fpath.suffix == '.tsv' else ','
        try:
            df = pd.read_csv(fpath, sep=sep, on_bad_lines='skip', low_memory=False)
        except Exception as e:
            logger.warning(f"  Error reading {fpath}: {e}")
            continue

        logger.info(f"    Columns: {list(df.columns)[:15]}")
        logger.info(f"    Rows: {len(df):,}")

        # Look for metabolite and microbe columns
        met_cols = [c for c in df.columns if any(x in c.lower() for x in ['compound', 'metabolite', 'chemical'])]
        micro_cols = [c for c in df.columns if any(x in c.lower() for x in ['microbe', 'organism', 'species', 'taxon'])]
        hmdb_cols = [c for c in df.columns if 'hmdb' in c.lower()]

        if not met_cols:
            logger.info(f"    No metabolite column found, skipping")
            continue

        for _, row in df.iterrows():
            met_name = str(row[met_cols[0]]).strip() if met_cols else ''
            microbe = str(row[micro_cols[0]]).strip() if micro_cols else ''
            hmdb = str(row[hmdb_cols[0]]).strip() if hmdb_cols else ''

            if met_name == 'nan' or not met_name:
                continue

            # Map to HMDB if not already present
            if not hmdb or hmdb == 'nan':
                hmdb = name_to_hmdb.get(met_name.lower(), '')

            if hmdb:
                rows.append({
                    'HMDB_ID': hmdb,
                    'Metabolite_Name': met_name,
                    'Microbe_Name': microbe,
                    'relation_subtype': 'associated',
                    'evidence_origin': 'database_curated',
                    'Source': 'MiMeDB',
                })

    df = pd.DataFrame(rows)
    if len(df):
        df = df.drop_duplicates(['HMDB_ID', 'Microbe_Name'])
        logger.info(f"MiMeDB: {len(df):,} unique edges")

    out = DB_DIR / "mmi_mimedb.csv"
    df.to_csv(out, index=False)
    return df


def merge_mmi_v2():
    """Merge refined MMI + MiMeDB into v3."""
    logger.info("=== Merging MMI v3 ===")

    # Load refined base
    refined_path = DB_DIR / "mmi_refined.csv"
    if refined_path.exists():
        base = pd.read_csv(refined_path, low_memory=False)
    else:
        base = pd.read_csv(DB_DIR / "mmi_database.csv", low_memory=False)
        base['relation_subtype'] = base['Relationship_Type'].map({
            'exchange': 'bidirectional', 'causal': 'causal', 'correlative': 'correlative'
        }).fillna('associated')
        base['evidence_origin'] = base['Source'].map({
            'gutMGene': 'experimental', 'AGORA2': 'computational'
        }).fillna('database_curated')
    logger.info(f"Base: {len(base):,}")

    # Add MiMeDB if available
    mimedb_path = DB_DIR / "mmi_mimedb.csv"
    if mimedb_path.exists():
        mimedb = pd.read_csv(mimedb_path)
        logger.info(f"MiMeDB: {len(mimedb):,}")
        # Align columns
        for col in base.columns:
            if col not in mimedb.columns:
                mimedb[col] = ''
        base = pd.concat([base, mimedb[base.columns]], ignore_index=True)
        base = base.drop_duplicates(['HMDB_ID', 'Microbe_Name'], keep='first')

    out = DB_DIR / "mmi_database_v3.csv"
    base.to_csv(out, index=False)

    logger.info(f"\n=== MMI v3 Summary ===")
    logger.info(f"  Total: {len(base):,}")
    logger.info(f"  By subtype: {dict(base['relation_subtype'].value_counts())}")
    logger.info(f"  By evidence: {dict(base['evidence_origin'].value_counts())}")
    logger.info(f"  Unique metabolites: {base['HMDB_ID'].nunique():,}")
    logger.info(f"  Unique microbes: {base['Microbe_Name'].nunique():,}")
    logger.info(f"  Saved: {out}")


def main():
    parser = argparse.ArgumentParser(description="Expand MMI Database")
    parser.add_argument('--refine', action='store_true', help='Refine AGORA2 subtypes')
    parser.add_argument('--mimedb', action='store_true', help='Add MiMeDB')
    parser.add_argument('--merge', action='store_true', help='Merge into v3')
    parser.add_argument('--all', action='store_true')
    args = parser.parse_args()

    if args.all or args.refine:
        refine_mmi_subtypes()
    if args.all or args.mimedb:
        expand_mimedb()
    if args.all or args.merge:
        merge_mmi_v2()

    if not any([args.refine, args.mimedb, args.merge, args.all]):
        parser.print_help()


if __name__ == "__main__":
    main()
