#!/usr/bin/env python3
"""
Expand mGWAS Database — Study-Aware Metadata
=============================================
Adds trait_category, effect_direction, significance_tier, genome_build,
and locus_type to existing mGWAS edges.

Usage:
    conda run -n mpi-vgae python scripts/expand_mgwas_v2.py
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"

# Trait → category keyword mapping
TRAIT_KEYWORDS = {
    'lipid': ['cholesterol', 'triglyceride', 'lipoprotein', 'hdl', 'ldl', 'vldl',
              'fatty acid', 'phospholipid', 'sphingomyelin', 'ceramide', 'lipid',
              'apolipoprotein', 'phosphatidylcholine', 'lysophosphatidylcholine'],
    'amino_acid': ['amino acid', 'alanine', 'glycine', 'valine', 'leucine', 'isoleucine',
                   'proline', 'serine', 'threonine', 'phenylalanine', 'tyrosine',
                   'tryptophan', 'glutamine', 'glutamate', 'asparagine', 'aspartate',
                   'histidine', 'lysine', 'arginine', 'methionine', 'cysteine',
                   'homocysteine', 'creatinine', 'creatine'],
    'sugar': ['glucose', 'fructose', 'galactose', 'lactose', 'glycated', 'hba1c'],
    'steroid': ['testosterone', 'estradiol', 'cortisol', 'steroid', 'androgen',
                'dehydroepiandrosterone', 'dhea', 'sex hormone'],
    'nucleotide': ['urate', 'uric acid', 'purine', 'xanthine', 'hypoxanthine',
                   'inosine', 'adenosine', 'guanosine'],
    'bile_acid': ['bile acid', 'cholic', 'deoxycholic', 'chenodeoxycholic',
                  'ursodeoxycholic', 'lithocholic'],
    'vitamin': ['vitamin', 'folate', 'retinol', 'tocopherol', 'ascorbic'],
    'xenobiotic': ['caffeine', 'cotinine', 'nicotine', 'drug', 'xenobiotic'],
    'kidney': ['creatinine', 'cystatin', 'egfr', 'kidney', 'renal'],
    'liver': ['bilirubin', 'alt', 'ast', 'ggt', 'alkaline phosphatase', 'liver'],
}


def classify_trait(trait):
    """Map trait text to metabolite category."""
    trait_lower = str(trait).lower()
    for category, keywords in TRAIT_KEYWORDS.items():
        for kw in keywords:
            if kw in trait_lower:
                return category
    return 'other'


def main():
    logger.info("=== mGWAS Study-Aware Metadata ===")

    df = pd.read_csv(DB_DIR / "mgwas_database.csv", low_memory=False)
    logger.info(f"mGWAS base: {len(df):,} edges")

    # 1. Trait category
    df['trait_category'] = df['Trait'].apply(classify_trait)

    # 2. Effect direction from Beta
    beta = pd.to_numeric(df['Beta'], errors='coerce')
    df['effect_direction'] = np.where(beta > 0, 'positive',
                              np.where(beta < 0, 'negative', 'unknown'))
    df.loc[beta.isna(), 'effect_direction'] = 'unknown'

    # 3. Significance tier
    pval = pd.to_numeric(df['P_Value'], errors='coerce')
    df['significance_tier'] = np.where(pval < 5e-8, 'genome_wide',
                               np.where(pval < 1e-5, 'suggestive', 'nominal'))

    # 4. Genome build
    df['genome_build'] = 'GRCh38'

    # 5. Locus type (sentinel = lowest p per gene region)
    df['locus_type'] = 'secondary'
    if 'Mapped_Gene' in df.columns:
        sentinel_idx = df.groupby(['Mapped_Gene', 'HMDB_ID'])['P_Value'].idxmin()
        df.loc[sentinel_idx.dropna(), 'locus_type'] = 'sentinel'

    # 6. Evidence origin
    df['evidence_origin'] = 'experimental'  # All GWAS hits are experimental

    # Summary
    logger.info(f"\n=== mGWAS v3 Summary ===")
    logger.info(f"  Total: {len(df):,}")
    logger.info(f"  By trait_category:")
    for tc, count in df['trait_category'].value_counts().items():
        logger.info(f"    {tc}: {count:,}")
    logger.info(f"  By effect_direction:")
    for ed, count in df['effect_direction'].value_counts().items():
        logger.info(f"    {ed}: {count:,}")
    logger.info(f"  By significance_tier:")
    for st, count in df['significance_tier'].value_counts().items():
        logger.info(f"    {st}: {count:,}")
    logger.info(f"  By locus_type:")
    for lt, count in df['locus_type'].value_counts().items():
        logger.info(f"    {lt}: {count:,}")

    out = DB_DIR / "mgwas_database_v3.csv"
    df.to_csv(out, index=False)
    logger.info(f"  Saved: {out}")


if __name__ == "__main__":
    main()
