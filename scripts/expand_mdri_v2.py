#!/usr/bin/env python3
"""
Expand MDrI Database — Typed Drug-Metabolite Subtypes
=====================================================
Derives relation_subtype from existing Description/Interaction_Type fields.
Subtypes: biotransformation, metabolic_inhibition, transport, enzyme_induction,
          abundance_change, pathway_context.

Usage:
    conda run -n mpi-vgae python scripts/expand_mdri_v2.py
"""

import logging
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_DIR = ROOT / "data" / "databases"


def classify_relation(row):
    """Derive relation_subtype from Description and Interaction_Type."""
    desc = str(row.get('Description', '')).lower()
    itype = str(row.get('Interaction_Type', '')).lower()

    # Priority order: most specific first
    if 'transport' in desc:
        return 'transport_interference'
    if 'inhibit' in desc and ('enzyme' in desc or 'metaboliz' in desc or 'cyp' in desc):
        return 'metabolic_inhibition'
    if 'induc' in desc and ('enzyme' in desc or 'cyp' in desc):
        return 'enzyme_induction'
    if 'substrate' in desc or 'metaboliz' in desc:
        return 'biotransformation'
    if 'produc' in desc or 'generat' in desc:
        return 'biotransformation'
    if 'alter' in desc or 'change' in desc or 'affect' in desc:
        return 'abundance_change'
    if itype == 'pharmacokinetic':
        return 'biotransformation'
    if itype == 'pharmacodynamic':
        return 'abundance_change'
    return 'associated'


def classify_mechanism(row):
    """Extract mechanism category from Description."""
    desc = str(row.get('Description', '')).lower()

    if 'cyp' in desc:
        # Extract CYP enzyme
        import re
        cyp = re.search(r'cyp\d+[a-z]*\d*', desc)
        return cyp.group(0).upper() if cyp else 'CYP-mediated'
    if 'ugt' in desc:
        return 'UGT-mediated'
    if 'p-glycoprotein' in desc or 'p-gp' in desc or 'abcb1' in desc:
        return 'P-gp transport'
    if 'oat' in desc or 'oct' in desc:
        return 'transporter-mediated'
    return ''


def main():
    logger.info("=== MDrI Subtype Classification ===")

    df = pd.read_csv(DB_DIR / "mdri_database.csv")
    logger.info(f"MDrI base: {len(df):,} edges")

    # Add relation_subtype
    df['relation_subtype'] = df.apply(classify_relation, axis=1)

    # Add mechanism_category
    df['mechanism_category'] = df.apply(classify_mechanism, axis=1)

    # Add evidence_origin
    df['evidence_origin'] = df['Source'].map({
        'CoreMet_curated': 'manual',
        'DrugBank_cross_ref': 'imported',
        'DrugBank_enzyme_bridge': 'inferred',
    }).fillna('imported')

    # Summary
    logger.info(f"\n=== MDrI v3 Summary ===")
    logger.info(f"  Total: {len(df):,}")
    logger.info(f"  By relation_subtype:")
    for st, count in df['relation_subtype'].value_counts().items():
        logger.info(f"    {st}: {count:,}")
    logger.info(f"  By mechanism_category (non-empty):")
    mech = df[df['mechanism_category'] != '']
    for mc, count in mech['mechanism_category'].value_counts().head(10).items():
        logger.info(f"    {mc}: {count:,}")
    logger.info(f"  By evidence_origin:")
    for eo, count in df['evidence_origin'].value_counts().items():
        logger.info(f"    {eo}: {count:,}")

    out = DB_DIR / "mdri_database_v3.csv"
    df.to_csv(out, index=False)
    logger.info(f"  Saved: {out}")


if __name__ == "__main__":
    main()
