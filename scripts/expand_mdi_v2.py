#!/usr/bin/env python3
"""
Expand MDI Database — Typed Disease Subtypes
=============================================
Adds HMDB disease annotations (clinical_annotation) and MarkerDB biomarkers
(clinical_biomarker) to the existing CTD-based MDI.

Usage:
    conda run -n mpi-vgae python scripts/expand_mdi_v2.py --hmdb
    conda run -n mpi-vgae python scripts/expand_mdi_v2.py --markerdb
    conda run -n mpi-vgae python scripts/expand_mdi_v2.py --merge
"""

import argparse
import logging
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"
RAW_DIR = ROOT / "data" / "raw"


# ══════════════════════════════════════════════════════════════════════
# HMDB Disease Annotations
# ══════════════════════════════════════════════════════════════════════

def expand_hmdb_diseases():
    """Stream-parse HMDB XML for disease associations."""
    logger.info("=== HMDB Disease Extraction ===")

    hmdb_xml = RAW_DIR / "hmdb" / "hmdb_metabolites.xml"
    if not hmdb_xml.exists():
        logger.warning(f"HMDB XML not found at {hmdb_xml}")
        return pd.DataFrame()

    logger.info(f"Streaming {hmdb_xml} ({hmdb_xml.stat().st_size / 1e9:.1f} GB)...")

    ns = '{http://www.hmdb.ca}'
    rows = []
    met_count = 0

    # Use iterparse to handle 6 GB file without loading into memory
    context = ET.iterparse(str(hmdb_xml), events=('end',))
    for event, elem in context:
        if elem.tag == f'{ns}metabolite':
            met_count += 1
            if met_count % 10000 == 0:
                logger.info(f"  Processed {met_count:,} metabolites, {len(rows):,} disease edges")

            # Extract HMDB ID
            accession = elem.findtext(f'{ns}accession', '')
            if not accession:
                elem.clear()
                continue

            # Extract metabolite name and SMILES
            met_name = elem.findtext(f'{ns}name', '')
            smiles = elem.findtext(f'{ns}smiles', '')

            # Extract diseases
            diseases = elem.find(f'{ns}diseases')
            if diseases is not None:
                for disease in diseases.findall(f'{ns}disease'):
                    disease_name = disease.findtext(f'{ns}name', '').strip()
                    if not disease_name:
                        continue

                    # Extract references
                    refs = disease.find(f'{ns}references')
                    pmids = []
                    if refs is not None:
                        for ref in refs.findall(f'{ns}reference'):
                            pmid = ref.findtext(f'{ns}pubmed_id', '').strip()
                            if pmid:
                                pmids.append(pmid)

                    rows.append({
                        'HMDB_ID': accession,
                        'Metabolite_Name': met_name,
                        'SMILES': smiles,
                        'Disease_Name': disease_name,
                        'Disease_ID': '',  # HMDB doesn't provide MeSH/DOID
                        'association_subtype': 'clinical_annotation',
                        'pmid': ';'.join(pmids[:5]),  # Cap at 5 PMIDs
                        'source': 'HMDB',
                        'evidence_type': 'curated',
                        'confidence': 0.7,
                    })

            # Free memory
            elem.clear()

    df = pd.DataFrame(rows)
    logger.info(f"HMDB: {len(df):,} disease edges from {met_count:,} metabolites")

    if len(df):
        df = df.drop_duplicates(['HMDB_ID', 'Disease_Name'])
        logger.info(f"  After dedup: {len(df):,} unique (metabolite, disease) pairs")
        logger.info(f"  Unique metabolites: {df['HMDB_ID'].nunique():,}")
        logger.info(f"  Unique diseases: {df['Disease_Name'].nunique():,}")
        logger.info(f"  With PMID: {(df['pmid'] != '').sum():,}")

    out = DB_DIR / "mdi_hmdb.csv"
    df.to_csv(out, index=False)
    logger.info(f"  Saved: {out}")
    return df


# ══════════════════════════════════════════════════════════════════════
# MarkerDB
# ══════════════════════════════════════════════════════════════════════

def expand_markerdb():
    """Parse MarkerDB TSV for clinical biomarkers."""
    logger.info("=== MarkerDB Clinical Biomarkers ===")

    markerdb_dir = RAW_DIR / "markerdb"
    tsv_files = list(markerdb_dir.glob("*.tsv")) + list(markerdb_dir.glob("*.csv")) \
        if markerdb_dir.exists() else []

    if not tsv_files:
        logger.warning(f"No MarkerDB files in {markerdb_dir}. "
                        "Download from https://markerdb.ca/downloads")
        return pd.DataFrame()

    # Load HMDB synonym table for name matching
    syn_path = ROOT / "data" / "mappings" / "hmdb_synonyms.csv"
    name_to_hmdb = {}
    if syn_path.exists():
        syn = pd.read_csv(syn_path, usecols=['hmdb_id', 'name'])
        for _, row in syn.drop_duplicates('name').iterrows():
            name_to_hmdb[str(row['name']).lower().strip()] = str(row['hmdb_id'])
        logger.info(f"Loaded {len(name_to_hmdb):,} name→HMDB mappings")

    rows = []
    for fpath in tsv_files:
        logger.info(f"  Parsing {fpath.name}...")
        sep = '\t' if fpath.suffix == '.tsv' else ','
        df = pd.read_csv(fpath, sep=sep, on_bad_lines='skip', low_memory=False)
        logger.info(f"    Columns: {list(df.columns)[:10]}")
        logger.info(f"    Rows: {len(df):,}")

        # Find relevant columns (names vary)
        name_cols = [c for c in df.columns if 'name' in c.lower() and 'bio' not in c.lower()]
        disease_cols = [c for c in df.columns if 'condition' in c.lower() or 'disease' in c.lower()]
        cat_cols = [c for c in df.columns if 'category' in c.lower() or 'type' in c.lower()]
        specimen_cols = [c for c in df.columns if 'specimen' in c.lower() or 'biofluid' in c.lower()]

        if not name_cols or not disease_cols:
            logger.info(f"    Skipping — can't find name/disease columns")
            continue

        for _, row in df.iterrows():
            met_name = str(row[name_cols[0]]).strip() if name_cols else ''
            disease = str(row[disease_cols[0]]).strip() if disease_cols else ''
            category = str(row[cat_cols[0]]).strip() if cat_cols else ''
            specimen = str(row[specimen_cols[0]]).strip() if specimen_cols else ''

            if not met_name or not disease or met_name == 'nan' or disease == 'nan':
                continue

            # Map to HMDB
            hmdb = name_to_hmdb.get(met_name.lower(), '')

            rows.append({
                'HMDB_ID': hmdb,
                'Metabolite_Name': met_name,
                'Disease_Name': disease,
                'Disease_ID': '',
                'biomarker_category': category,
                'specimen': specimen,
                'association_subtype': 'clinical_biomarker',
                'source': 'MarkerDB',
                'evidence_type': 'curated',
                'confidence': 0.85,
            })

    df = pd.DataFrame(rows)
    if len(df):
        mapped = (df['HMDB_ID'] != '').sum()
        logger.info(f"MarkerDB: {len(df):,} total, {mapped:,} HMDB-mapped ({100*mapped/len(df):.1f}%)")
        # Keep only HMDB-mapped entries
        df = df[df['HMDB_ID'] != '']
        df = df.drop_duplicates(['HMDB_ID', 'Disease_Name'])
        logger.info(f"  Final: {len(df):,} unique mapped edges")

    out = DB_DIR / "mdi_markerdb.csv"
    df.to_csv(out, index=False)
    return df


# ══════════════════════════════════════════════════════════════════════
# Merge
# ══════════════════════════════════════════════════════════════════════

def merge_mdi_v2():
    """Merge CTD + HMDB + MarkerDB into MDI v3."""
    logger.info("=== Merging MDI v3 ===")

    # Load CTD base
    ctd = pd.read_csv(DB_DIR / "mdi_database.csv", low_memory=False)
    ctd['association_subtype'] = ctd['Association_Type'].map({
        'Biomarker': 'mechanistic_curated',
        'Therapeutic': 'therapeutic',
        'Disease-associated biomarker': 'mechanistic_curated',
    }).fillna('mechanistic_curated')
    ctd['biomarker_category'] = ''
    ctd['specimen'] = ''
    logger.info(f"CTD base: {len(ctd):,}")

    # Standard columns for merge
    std_cols = ['HMDB_ID', 'Metabolite_Name', 'SMILES', 'Disease_Name', 'Disease_ID',
                'association_subtype', 'biomarker_category', 'specimen',
                'source', 'evidence_type', 'confidence', 'pmid']

    for col in std_cols:
        if col not in ctd.columns:
            ctd[col] = ''
    ctd['source'] = ctd.get('Source', 'CTD').fillna('CTD')
    ctd['evidence_type'] = 'curated'
    ctd['confidence'] = 0.7
    ctd['pmid'] = ''

    all_dfs = [ctd[std_cols]]

    # Add HMDB diseases
    hmdb_path = DB_DIR / "mdi_hmdb.csv"
    if hmdb_path.exists():
        hmdb = pd.read_csv(hmdb_path)
        for col in std_cols:
            if col not in hmdb.columns:
                hmdb[col] = ''
        all_dfs.append(hmdb[std_cols])
        logger.info(f"HMDB diseases: {len(hmdb):,}")

    # Add MarkerDB
    marker_path = DB_DIR / "mdi_markerdb.csv"
    if marker_path.exists():
        marker = pd.read_csv(marker_path)
        for col in std_cols:
            if col not in marker.columns:
                marker[col] = ''
        all_dfs.append(marker[std_cols])
        logger.info(f"MarkerDB: {len(marker):,}")

    merged = pd.concat(all_dfs, ignore_index=True)
    merged = merged.drop_duplicates(['HMDB_ID', 'Disease_Name'], keep='first')

    out = DB_DIR / "mdi_database_v3.csv"
    merged.to_csv(out, index=False)

    logger.info(f"\n=== MDI v3 Summary ===")
    logger.info(f"  Total: {len(merged):,}")
    logger.info(f"  By subtype:")
    for st, count in merged['association_subtype'].value_counts().items():
        logger.info(f"    {st}: {count:,}")
    logger.info(f"  Unique metabolites: {merged['HMDB_ID'].nunique():,}")
    logger.info(f"  Unique diseases: {merged['Disease_Name'].nunique():,}")
    logger.info(f"  Saved: {out}")


def main():
    parser = argparse.ArgumentParser(description="Expand MDI Database")
    parser.add_argument('--hmdb', action='store_true')
    parser.add_argument('--markerdb', action='store_true')
    parser.add_argument('--merge', action='store_true')
    parser.add_argument('--all', action='store_true')
    args = parser.parse_args()

    if args.all or args.hmdb:
        expand_hmdb_diseases()
    if args.all or args.markerdb:
        expand_markerdb()
    if args.all or args.merge:
        merge_mdi_v2()

    if not any([args.hmdb, args.markerdb, args.merge, args.all]):
        parser.print_help()


if __name__ == "__main__":
    main()
