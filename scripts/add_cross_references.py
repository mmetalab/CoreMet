#!/usr/bin/env python3
"""
Add Cross-Database Identifiers to CoreMet
=============================================
Enriches metabolite entities with: InChIKey, PubChem CID, ChEBI ID, KEGG Compound ID.
Uses SMILES→InChIKey via RDKit and existing MMI cross-refs as seed.

Usage:
    conda run -n mpi-vgae python scripts/add_cross_references.py
"""

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
DB_V2 = ROOT / "data" / "databases_v2"
MAP_DIR = ROOT / "data" / "mappings"
MAP_DIR.mkdir(parents=True, exist_ok=True)


def build_metabolite_crossref():
    """Build a master HMDB_ID → cross-references mapping."""
    logger.info("Building metabolite cross-reference table...")

    crossref = {}  # HMDB_ID → {inchikey, pubchem_cid, chebi_id, kegg_id, smiles}

    # 1. Collect SMILES from all databases
    smiles_map = {}
    for csv_name in ['mpi_database_v2.csv', 'mdi_database_v2.csv', 'mmi_database_v2.csv',
                      'mdri_database_v2.csv', 'mgi_database_v2.csv', 'mgwas_database_v2.csv']:
        path = DB_V2 / csv_name
        if not path.exists():
            continue
        df = pd.read_csv(path, usecols=lambda c: c in ['HMDB_ID', 'HMDB ID', 'SMILES'],
                          low_memory=False)
        hmdb_col = 'HMDB_ID' if 'HMDB_ID' in df.columns else 'HMDB ID'
        for _, row in df.drop_duplicates(hmdb_col).iterrows():
            hmdb = str(row[hmdb_col])
            smi = row.get('SMILES', '')
            if hmdb and pd.notna(smi) and str(smi).strip():
                smiles_map[hmdb] = str(smi).strip()

    logger.info(f"  SMILES collected: {len(smiles_map)} metabolites")

    # 2. Generate InChIKey from SMILES via RDKit
    try:
        from rdkit import Chem
        from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
        inchikey_map = {}
        for hmdb, smi in smiles_map.items():
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    inchi = MolToInchi(mol)
                    if inchi:
                        inchikey_map[hmdb] = InchiToInchiKey(inchi)
            except Exception:
                pass
        logger.info(f"  InChIKey generated: {len(inchikey_map)} metabolites")
    except ImportError:
        logger.warning("  RDKit not available — skipping InChIKey generation")
        inchikey_map = {}

    # 3. Collect PubChem CID, ChEBI ID, KEGG ID from MMI (richest cross-refs)
    mmi_path = DB_V2 / "mmi_database_v2.csv"
    pubchem_map = {}
    chebi_map = {}
    kegg_map = {}
    if mmi_path.exists():
        mmi = pd.read_csv(mmi_path, usecols=lambda c: c in [
            'HMDB_ID', 'PubChem_CID', 'ChEBI_ID', 'KEGG_ID'], low_memory=False)
        for _, row in mmi.drop_duplicates('HMDB_ID').iterrows():
            hmdb = str(row['HMDB_ID'])
            if pd.notna(row.get('PubChem_CID')) and str(row['PubChem_CID']).strip():
                pubchem_map[hmdb] = str(row['PubChem_CID']).strip()
            if pd.notna(row.get('ChEBI_ID')) and str(row['ChEBI_ID']).strip():
                chebi_map[hmdb] = str(row['ChEBI_ID']).strip()
            if pd.notna(row.get('KEGG_ID')) and str(row['KEGG_ID']).strip():
                kegg_map[hmdb] = str(row['KEGG_ID']).strip()
        logger.info(f"  From MMI: PubChem={len(pubchem_map)}, ChEBI={len(chebi_map)}, KEGG={len(kegg_map)}")

    # 4. Merge into crossref dict
    all_hmdb = set(smiles_map.keys())
    for hmdb in all_hmdb:
        crossref[hmdb] = {
            'smiles': smiles_map.get(hmdb, ''),
            'inchikey': inchikey_map.get(hmdb, ''),
            'pubchem_cid': pubchem_map.get(hmdb, ''),
            'chebi_id': chebi_map.get(hmdb, ''),
            'kegg_id': kegg_map.get(hmdb, ''),
        }

    # 5. Stats
    n_total = len(crossref)
    n_smiles = sum(1 for v in crossref.values() if v['smiles'])
    n_inchi = sum(1 for v in crossref.values() if v['inchikey'])
    n_pubchem = sum(1 for v in crossref.values() if v['pubchem_cid'])
    n_chebi = sum(1 for v in crossref.values() if v['chebi_id'])
    n_kegg = sum(1 for v in crossref.values() if v['kegg_id'])

    logger.info(f"\n  Cross-reference coverage ({n_total} metabolites):")
    logger.info(f"    SMILES:      {n_smiles:>6,} ({100*n_smiles/n_total:.1f}%)")
    logger.info(f"    InChIKey:    {n_inchi:>6,} ({100*n_inchi/n_total:.1f}%)")
    logger.info(f"    PubChem CID: {n_pubchem:>6,} ({100*n_pubchem/n_total:.1f}%)")
    logger.info(f"    ChEBI ID:    {n_chebi:>6,} ({100*n_chebi/n_total:.1f}%)")
    logger.info(f"    KEGG ID:     {n_kegg:>6,} ({100*n_kegg/n_total:.1f}%)")

    # 6. Save
    out_json = MAP_DIR / "metabolite_crossref.json"
    with open(out_json, 'w') as f:
        json.dump(crossref, f, indent=2)
    logger.info(f"  Saved: {out_json}")

    # Also save as CSV for easy inspection
    out_csv = MAP_DIR / "metabolite_crossref.csv"
    rows = [{'HMDB_ID': k, **v} for k, v in crossref.items()]
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    logger.info(f"  Saved: {out_csv}")

    return crossref


def enrich_databases(crossref):
    """Add cross-reference columns to all v2 database CSVs."""
    logger.info("\nEnriching databases with cross-references...")

    for csv_name in DB_V2.glob("*_database_v2.csv"):
        df = pd.read_csv(csv_name, low_memory=False)
        hmdb_col = 'HMDB_ID' if 'HMDB_ID' in df.columns else 'HMDB ID'

        if hmdb_col not in df.columns:
            logger.info(f"  {csv_name.name}: no HMDB column, skipping")
            continue

        # Add cross-ref columns if not already present
        for field in ['inchikey', 'pubchem_cid', 'chebi_id', 'kegg_id']:
            col_name = field.upper() if field == 'inchikey' else field.replace('_', ' ').title().replace(' ', '_')
            if col_name not in df.columns:
                df[col_name] = df[hmdb_col].map(
                    lambda h: crossref.get(str(h), {}).get(field, ''))

        df.to_csv(csv_name, index=False)
        logger.info(f"  {csv_name.name}: enriched ({len(df):,} rows)")


def main():
    logger.info("CoreMet Cross-Reference Enrichment")
    crossref = build_metabolite_crossref()
    enrich_databases(crossref)
    logger.info("\nDone!")


if __name__ == "__main__":
    main()
