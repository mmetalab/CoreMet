#!/usr/bin/env python3
"""
expand_mei_hmdb.py
==================
Improve MEI HMDB_ID coverage by cross-referencing KEGG Compound IDs
against the full HMDB XML (6 GB).  Uses iterparse to stream through
the XML without loading it all into memory.

Maps KEGG Compound → HMDB ID, then backfills HMDB_ID in mei_database.csv.

Usage:
    python scripts/expand_mei_hmdb.py
"""

import os, re, csv, time, sys
from collections import defaultdict
from lxml import etree

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.join(os.path.dirname(__file__), '..')
HMDB_XML = os.path.join(BASE, 'data', 'raw', 'hmdb', 'hmdb_metabolites.xml')
MEI_DB   = os.path.join(BASE, 'data', 'databases', 'mei_database.csv')
MEI_BAK  = MEI_DB + '.bak'

NS = '{http://www.hmdb.ca}'


def build_kegg_to_hmdb_map():
    """Stream-parse HMDB XML and build KEGG Compound → HMDB ID mapping.
    
    Also builds a metabolite name → HMDB ID mapping for fuzzy matching.
    """
    print("Phase 1: Building KEGG → HMDB mapping from HMDB XML...")
    print(f"  Reading: {HMDB_XML}")
    
    kegg_to_hmdb = {}       # KEGG_ID -> HMDB_ID (primary)
    name_to_hmdb = {}       # lowercase metabolite name -> HMDB_ID
    chebi_to_hmdb = {}      # ChEBI ID -> HMDB_ID
    
    count = 0
    t0 = time.time()
    
    context = etree.iterparse(HMDB_XML, events=('end',), tag=f'{NS}metabolite')
    
    for event, elem in context:
        count += 1
        
        # Extract HMDB accession
        acc_el = elem.find(f'{NS}accession')
        if acc_el is None or not acc_el.text:
            elem.clear()
            continue
        hmdb_id = acc_el.text.strip()
        
        # Normalize to 7-digit HMDB format
        if hmdb_id.startswith('HMDB') and len(hmdb_id) < 11:
            digits = hmdb_id[4:]
            hmdb_id = 'HMDB' + digits.zfill(7)
        
        # Extract name
        name_el = elem.find(f'{NS}name')
        if name_el is not None and name_el.text:
            met_name = name_el.text.strip()
            name_lower = met_name.lower()
            if name_lower not in name_to_hmdb:
                name_to_hmdb[name_lower] = hmdb_id
        
        # Extract synonyms
        synonyms_el = elem.find(f'{NS}synonyms')
        if synonyms_el is not None:
            for syn_el in synonyms_el.findall(f'{NS}synonym'):
                if syn_el.text:
                    syn_lower = syn_el.text.strip().lower()
                    if syn_lower not in name_to_hmdb:
                        name_to_hmdb[syn_lower] = hmdb_id
        
        # Extract KEGG ID from <kegg_id> element
        kegg_el = elem.find(f'{NS}kegg_id')
        if kegg_el is not None and kegg_el.text and kegg_el.text.strip():
            kegg_id = kegg_el.text.strip()
            if kegg_id not in kegg_to_hmdb:
                kegg_to_hmdb[kegg_id] = hmdb_id
        
        # Extract ChEBI ID
        chebi_el = elem.find(f'{NS}chebi_id')
        if chebi_el is not None and chebi_el.text and chebi_el.text.strip():
            chebi_id = chebi_el.text.strip()
            if chebi_id not in chebi_to_hmdb:
                chebi_to_hmdb[chebi_id] = hmdb_id
        
        # Free memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
        
        if count % 50000 == 0:
            elapsed = time.time() - t0
            print(f"  Processed {count:,} metabolites ({elapsed:.1f}s)")
    
    elapsed = time.time() - t0
    print(f"  Done: {count:,} metabolites in {elapsed:.1f}s")
    print(f"  KEGG → HMDB mappings: {len(kegg_to_hmdb):,}")
    print(f"  Name → HMDB mappings: {len(name_to_hmdb):,}")
    print(f"  ChEBI → HMDB mappings: {len(chebi_to_hmdb):,}")
    
    return kegg_to_hmdb, name_to_hmdb, chebi_to_hmdb


def update_mei_database(kegg_to_hmdb, name_to_hmdb, chebi_to_hmdb):
    """Read MEI database and backfill HMDB_ID for records that are missing it."""
    print("\nPhase 2: Backfilling HMDB_ID in MEI database...")
    
    import shutil
    shutil.copy2(MEI_DB, MEI_BAK)
    print(f"  Backup saved to {MEI_BAK}")
    
    rows = []
    total = 0
    already_has = 0
    filled_kegg = 0
    filled_name = 0
    still_missing = 0
    
    with open(MEI_BAK, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        
        for row in reader:
            total += 1
            hmdb_id = (row.get('HMDB_ID') or '').strip()
            
            if hmdb_id:
                already_has += 1
                rows.append(row)
                continue
            
            # Try KEGG mapping first (most reliable)
            kegg_id = (row.get('KEGG_Compound') or '').strip()
            if kegg_id and kegg_id in kegg_to_hmdb:
                row['HMDB_ID'] = kegg_to_hmdb[kegg_id]
                filled_kegg += 1
                rows.append(row)
                continue
            
            # Try name matching
            met_name = (row.get('Metabolite_Name') or '').strip().lower()
            if met_name and met_name in name_to_hmdb:
                row['HMDB_ID'] = name_to_hmdb[met_name]
                filled_name += 1
                rows.append(row)
                continue
            
            still_missing += 1
            rows.append(row)
    
    # Write updated database
    with open(MEI_DB, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    
    filled_total = filled_kegg + filled_name
    new_coverage = (already_has + filled_total) / total * 100
    
    print(f"\n  Total MEI records: {total:,}")
    print(f"  Already had HMDB_ID: {already_has:,}")
    print(f"  Filled via KEGG: {filled_kegg:,}")
    print(f"  Filled via Name: {filled_name:,}")
    print(f"  Still missing: {still_missing:,}")
    print(f"\n  HMDB coverage: {already_has/total*100:.1f}% → {new_coverage:.1f}%")
    print(f"  Records with HMDB_ID: {already_has + filled_total:,} / {total:,}")
    
    return {
        'total': total,
        'already_had': already_has,
        'filled_kegg': filled_kegg,
        'filled_name': filled_name,
        'still_missing': still_missing,
        'old_coverage': already_has / total * 100,
        'new_coverage': new_coverage,
    }


def main():
    print("=" * 70)
    print("MEI HMDB Coverage Expansion via HMDB XML")
    print("=" * 70)
    
    if not os.path.isfile(HMDB_XML):
        print(f"ERROR: HMDB XML not found: {HMDB_XML}")
        sys.exit(1)
    
    if not os.path.isfile(MEI_DB):
        print(f"ERROR: MEI database not found: {MEI_DB}")
        sys.exit(1)
    
    kegg_to_hmdb, name_to_hmdb, chebi_to_hmdb = build_kegg_to_hmdb_map()
    stats = update_mei_database(kegg_to_hmdb, name_to_hmdb, chebi_to_hmdb)
    
    print(f"\n{'='*70}")
    print(f"  HMDB COVERAGE: {stats['old_coverage']:.1f}% → {stats['new_coverage']:.1f}%")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
