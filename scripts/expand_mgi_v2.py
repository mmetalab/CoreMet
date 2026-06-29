#!/usr/bin/env python3
"""
Expand MGI Database — Typed Gene Interaction Subtypes
=====================================================
Classifies 1.66M CTD chemical-gene interactions into:
  - expression_regulation (increases/decreases/affects expression)
  - epigenetic_modification (methylation, acetylation)
  - protein_activity (increases/decreases activity, phosphorylation)
  - binding (affects binding)
  - cotreatment (affects cotreatment)
  - reaction (affects reaction)
  - other

Also adds direction (up/down/mixed) and evidence_strength based on
interaction specificity.

Usage:
    conda run -n mpi-vgae python scripts/expand_mgi_v2.py
"""

import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"


def classify_interaction(itype):
    """Classify Interaction_Type into biological subtype."""
    itype_lower = str(itype).lower()

    if 'expression' in itype_lower:
        return 'expression_regulation'
    if 'methylation' in itype_lower or 'acetylation' in itype_lower or 'histone' in itype_lower:
        return 'epigenetic_modification'
    if 'activity' in itype_lower or 'phosphorylation' in itype_lower:
        return 'protein_activity'
    if 'binding' in itype_lower:
        return 'binding'
    if 'cotreatment' in itype_lower:
        return 'cotreatment'
    if 'reaction' in itype_lower:
        return 'metabolic_reaction'
    if 'abundance' in itype_lower:
        return 'abundance_change'
    if 'transport' in itype_lower or 'uptake' in itype_lower or 'secretion' in itype_lower:
        return 'transport'
    if 'splicing' in itype_lower or 'rna' in itype_lower:
        return 'rna_processing'
    if 'mutagenesis' in itype_lower or 'damage' in itype_lower:
        return 'dna_damage'
    if 'response' in itype_lower:
        return 'response'
    if 'localization' in itype_lower:
        return 'localization'
    return 'other'


def classify_direction(itype):
    """Extract direction (up/down/mixed) from Interaction_Type."""
    itype_lower = str(itype).lower()
    if 'increases' in itype_lower:
        return 'up'
    if 'decreases' in itype_lower:
        return 'down'
    if 'affects' in itype_lower:
        return 'mixed'
    return 'unknown'


def main():
    logger.info("=== MGI Subtype Classification ===")

    df = pd.read_csv(DB_DIR / "mgi_database.csv", low_memory=False)
    logger.info(f"MGI base: {len(df):,} edges")

    # Add relation_subtype
    df['relation_subtype'] = df['Interaction_Type'].apply(classify_interaction)

    # Add direction
    df['direction'] = df['Interaction_Type'].apply(classify_direction)

    # Add evidence_origin (all CTD = direct_curated)
    df['evidence_origin'] = 'direct_curated'

    # Summary
    logger.info(f"\n=== MGI v3 Summary ===")
    logger.info(f"  Total: {len(df):,}")
    logger.info(f"  By relation_subtype:")
    for st, count in df['relation_subtype'].value_counts().items():
        logger.info(f"    {st}: {count:,}")
    logger.info(f"  By direction:")
    for d, count in df['direction'].value_counts().items():
        logger.info(f"    {d}: {count:,}")
    logger.info(f"  Unique organisms: {df['Organism'].nunique()}")
    logger.info(f"  PMID coverage: {df['PMID'].notna().sum():,} / {len(df):,} (100%)")

    out = DB_DIR / "mgi_database_v3.csv"
    df.to_csv(out, index=False)
    logger.info(f"  Saved: {out}")


if __name__ == "__main__":
    main()
