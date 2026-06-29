#!/usr/bin/env python3
"""
Standardize Evidence Model for CoreMet
==========================================
Adds uniform provenance fields to all 7 database CSVs:
  source_db, evidence_type, species, tissue, pmid, curation_type, confidence

Usage:
    conda run -n mpi-vgae python scripts/standardize_evidence.py
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"
MPI_PATH = ROOT / "data" / "mpidatabase" / "MPIDB_v2.csv"
OUT_DIR = ROOT / "data" / "databases_v2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Standard evidence fields to add
STANDARD_FIELDS = [
    'source_db', 'evidence_type', 'curation_type', 'confidence',
    'species_id', 'tissue', 'pmid',
]


def _norm_confidence(series, method='minmax'):
    """Normalize a numeric series to [0, 1]."""
    vals = pd.to_numeric(series, errors='coerce')
    if vals.isna().all():
        return pd.Series(np.nan, index=series.index)
    vmin, vmax = vals.min(), vals.max()
    if vmax == vmin:
        return pd.Series(0.5, index=series.index)
    return (vals - vmin) / (vmax - vmin)


def standardize_mpi():
    """MPI: KEGG/Rhea/MPIDB evidence."""
    logger.info("Standardizing MPI...")
    df = pd.read_csv(MPI_PATH)

    df['source_db'] = df['Evidence_Source'].fillna('MPIDB')
    df['evidence_type'] = df['Evidence_Source'].map({
        'KEGG': 'curated', 'Rhea': 'curated', 'Original': 'curated',
    }).fillna('curated')
    df['curation_type'] = 'imported'
    df['confidence'] = np.where(df['Pathway_Name'].notna(), 0.9, 0.7)
    df['species_id'] = ''  # Would need NCBI TaxID mapping
    df['tissue'] = ''
    df['pmid'] = ''

    out = OUT_DIR / "mpi_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  MPI: {len(df):,} rows → {out}")
    return df


def standardize_mei():
    """MEI: KEGG enzyme annotations."""
    logger.info("Standardizing MEI...")
    df = pd.read_csv(DB_DIR / "mei_database.csv")

    df['source_db'] = df.get('Evidence_Source', 'KEGG').fillna('KEGG')
    df['evidence_type'] = 'curated'
    df['curation_type'] = 'imported'
    df['confidence'] = 0.9
    df['species_id'] = ''
    df['tissue'] = ''
    df['pmid'] = ''

    out = OUT_DIR / "mei_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  MEI: {len(df):,} rows → {out}")
    return df


def standardize_mdi():
    """MDI: CTD disease associations."""
    logger.info("Standardizing MDI...")
    df = pd.read_csv(DB_DIR / "mdi_database.csv")

    df['source_db'] = df.get('Source', 'CTD').fillna('CTD')
    df['evidence_type'] = df['Association_Type'].map({
        'marker/mechanism': 'experimental',
        'therapeutic': 'experimental',
        'marker': 'experimental',
    }).fillna('inferred')
    # Override with Evidence_Level if available
    if 'Evidence_Level' in df.columns:
        mask = df['Evidence_Level'].str.lower() == 'direct'
        df.loc[mask, 'evidence_type'] = 'experimental'
    df['curation_type'] = 'imported'
    df['confidence'] = _norm_confidence(df.get('Avg_Network_Score', pd.Series()))
    df['confidence'] = df['confidence'].fillna(0.5)
    df['species_id'] = '9606'  # CTD is human-centric
    df['tissue'] = ''
    df['pmid'] = ''

    out = OUT_DIR / "mdi_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  MDI: {len(df):,} rows → {out}")
    return df


def standardize_mmi():
    """MMI: gutMGene + AGORA2 microbe interactions."""
    logger.info("Standardizing MMI...")
    df = pd.read_csv(DB_DIR / "mmi_database.csv")

    df['source_db'] = df.get('Source', '').fillna('gutMGene')
    df['evidence_type'] = df.get('Evidence_Level', '').map({
        'direct': 'experimental',
        'experimental': 'experimental',
        'computational': 'computational',
    }).fillna('curated')
    df['curation_type'] = np.where(df['source_db'] == 'AGORA2', 'computational', 'imported')
    df['confidence'] = np.where(df['evidence_type'] == 'experimental', 0.9, 0.6)
    df['species_id'] = df.get('Taxonomy_ID', '').astype(str)
    df['tissue'] = df.get('Tissue', '').fillna('')
    df['pmid'] = df.get('PMID', '').astype(str).replace('nan', '')

    out = OUT_DIR / "mmi_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  MMI: {len(df):,} rows → {out}")
    return df


def standardize_mdri():
    """MDrI: curated + DrugBank drug interactions."""
    logger.info("Standardizing MDrI...")
    df = pd.read_csv(DB_DIR / "mdri_database.csv")

    df['source_db'] = df.get('Source', '').fillna('DrugBank')
    df['evidence_type'] = df.get('Evidence_Level', '').map({
        'direct': 'experimental',
        'experimental': 'experimental',
    }).fillna('curated')
    df['curation_type'] = np.where(df['source_db'] == 'curated', 'manual', 'imported')
    df['confidence'] = np.where(df['evidence_type'] == 'experimental', 0.9, 0.7)
    df['species_id'] = '9606'
    df['tissue'] = df.get('Tissue', '').fillna('')
    df['pmid'] = df.get('PMID', '').astype(str).replace('nan', '')

    out = OUT_DIR / "mdri_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  MDrI: {len(df):,} rows → {out}")
    return df


def standardize_mgi():
    """MGI: CTD chemical-gene interactions."""
    logger.info("Standardizing MGI...")
    df = pd.read_csv(DB_DIR / "mgi_database.csv")

    df['source_db'] = df.get('Source', 'CTD').fillna('CTD')
    # Map CTD interaction actions to evidence types
    actions = df.get('Interaction_Actions', '').astype(str).str.lower()
    df['evidence_type'] = np.where(
        actions.str.contains('expression|activity|phosphorylation', na=False),
        'experimental', 'curated')
    df['curation_type'] = 'imported'
    df['confidence'] = 0.8  # CTD has curated interactions
    df['species_id'] = ''  # Would need organism → TaxID mapping
    df['tissue'] = ''
    df['pmid'] = df.get('PMID', '').astype(str).replace('nan', '')

    out = OUT_DIR / "mgi_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  MGI: {len(df):,} rows → {out}")
    return df


def standardize_mgwas():
    """mGWAS: GWAS Catalog + curated supplements."""
    logger.info("Standardizing mGWAS...")
    df = pd.read_csv(DB_DIR / "mgwas_database.csv")

    df['source_db'] = df.get('Source', 'GWAS_Catalog').fillna('GWAS_Catalog')
    df['evidence_type'] = 'experimental'  # All GWAS hits are experimental
    df['curation_type'] = 'imported'
    # Confidence from p-value: -log10(p) normalized to [0, 1]
    pvals = pd.to_numeric(df.get('P_Value', pd.Series()), errors='coerce')
    log_p = -np.log10(pvals.clip(lower=1e-300))
    df['confidence'] = (log_p / log_p.max()).clip(0, 1).fillna(0.5)
    df['species_id'] = '9606'
    df['tissue'] = ''
    df['pmid'] = df.get('PMID', '').astype(str).replace('nan', '')

    out = OUT_DIR / "mgwas_database_v2.csv"
    df.to_csv(out, index=False)
    logger.info(f"  mGWAS: {len(df):,} rows → {out}")
    return df


def main():
    logger.info("CoreMet Evidence Standardization")
    logger.info(f"Output: {OUT_DIR}")

    results = {}
    results['MPI'] = standardize_mpi()
    results['MEI'] = standardize_mei()
    results['MDI'] = standardize_mdi()
    results['MMI'] = standardize_mmi()
    results['MDrI'] = standardize_mdri()
    results['MGI'] = standardize_mgi()
    results['mGWAS'] = standardize_mgwas()

    # Summary
    logger.info("\n=== Summary ===")
    total = 0
    for name, df in results.items():
        n = len(df)
        total += n
        has_pmid = (df['pmid'].astype(str) != '') & (df['pmid'].astype(str) != 'nan')
        has_conf = df['confidence'].notna()
        logger.info(f"  {name:8s}: {n:>10,} rows | PMID: {has_pmid.sum():>8,} ({100*has_pmid.mean():.1f}%) | Confidence: {has_conf.sum():>10,} ({100*has_conf.mean():.1f}%)")
    logger.info(f"  {'TOTAL':8s}: {total:>10,} rows")


if __name__ == "__main__":
    main()
