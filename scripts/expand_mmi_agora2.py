#!/usr/bin/env python3
"""
expand_mmi_agora2.py
====================
Expand the MMI (Metabolite–Microbe Interaction) database by extracting exchange
reactions from 818 AGORA2 genome-scale metabolic models (SBML format).

Each AGORA2 model represents one human gut microbe.  Exchange reactions (R_EX_*)
represent metabolites that the microbe can import/export.  We extract these
metabolites, map them to HMDB / KEGG / PubChem / ChEBI via SBML annotations,
and merge with the existing curated MMI database.

Usage:
    python scripts/expand_mmi_agora2.py
"""

import os, re, sys, csv, time
from collections import defaultdict
from lxml import etree

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.path.join(os.path.dirname(__file__), '..')
SBML_DIR = os.path.join(BASE, 'data', 'raw', 'agora2', 'reconstructions', 'sbml')
MMI_DB   = os.path.join(BASE, 'data', 'databases', 'mmi_database.csv')
MMI_OUT  = MMI_DB  # overwrite in-place (backup first)
MMI_BAK  = MMI_DB + '.bak'

# SBML namespaces
NS_SBML  = 'http://www.sbml.org/sbml/level3/version1/core'
NS_FBC   = 'http://www.sbml.org/sbml/level3/version1/fbc/version2'
NS_RDF   = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
NS_BQBIOL = 'http://biomodels.net/biology-qualifiers/'
NS_HTML  = 'http://www.w3.org/1999/xhtml'

# identifiers.org patterns
RE_HMDB    = re.compile(r'identifiers\.org/hmdb/(HMDB\d+)', re.I)
RE_KEGG    = re.compile(r'identifiers\.org/kegg\.compound/(C\d+)', re.I)
RE_PUBCHEM = re.compile(r'identifiers\.org/pubchem\.compound/(\d+)', re.I)
RE_CHEBI   = re.compile(r'identifiers\.org/chebi/CHEBI:(\d+)', re.I)
RE_TAXONOMY = re.compile(r'identifiers\.org/taxonomy/(\d+)', re.I)

# SBML encoded chars
def decode_sbml_id(s):
    """Decode SBML-encoded ID: __40__ → (, __41__ → ), __91__ → [, __93__ → ]"""
    s = s.replace('__40__', '(').replace('__41__', ')')
    s = s.replace('__91__', '[').replace('__93__', ']')
    s = s.replace('__45__', '-').replace('__46__', '.')
    return s

# ---------------------------------------------------------------------------
# Parse one SBML file
# ---------------------------------------------------------------------------
def parse_agora2_sbml(fpath):
    """
    Parse one AGORA2 SBML file and return:
      - organism_name (str)
      - taxonomy_id (str or '')
      - list of dicts: {metabolite_name, hmdb_id, kegg_id, pubchem_cid, chebi_id}
    for each exchange reaction metabolite.
    """
    tree = etree.parse(fpath)
    root = tree.getroot()
    ns = root.nsmap.get(None, NS_SBML)

    model = root.find(f'{{{ns}}}model')
    if model is None:
        return None, None, []

    organism_name = (model.get('name') or model.get('id') or '').replace('_', ' ')

    # Try to find taxonomy ID from model annotation
    taxonomy_id = ''
    model_annot = model.find(f'{{{ns}}}annotation')
    if model_annot is not None:
        annot_str = etree.tostring(model_annot, encoding='unicode')
        m = RE_TAXONOMY.search(annot_str)
        if m:
            taxonomy_id = m.group(1)

    # Build species lookup: species_id -> annotation info
    species_map = {}
    for sp in model.findall(f'.//{{{ns}}}species'):
        sid = sp.get('id', '')
        sname = sp.get('name', '')
        comp = sp.get('compartment', '')

        # Only consider extracellular species
        if comp != 'e':
            continue

        hmdb_id = ''
        kegg_id = ''
        pubchem_cid = ''
        chebi_id = ''

        # Extract from annotation (RDF)
        annotation = sp.find(f'{{{ns}}}annotation')
        if annotation is not None:
            annot_str = etree.tostring(annotation, encoding='unicode')
            m = RE_HMDB.search(annot_str)
            if m:
                hmdb_id = m.group(1)
                # Normalize old 5-digit HMDB IDs to 7-digit
                if len(hmdb_id) == 9:  # HMDB00112 → already 9 chars → OK (HMDB + 5 digits)
                    hmdb_id = 'HMDB00' + hmdb_id[4:]  # HMDB00112 → HMDB0000112
                elif len(hmdb_id) < 11 and hmdb_id.startswith('HMDB'):
                    # Pad to 11 chars total (HMDB + 7 digits)
                    digits = hmdb_id[4:]
                    hmdb_id = 'HMDB' + digits.zfill(7)
            m = RE_KEGG.search(annot_str)
            if m:
                kegg_id = m.group(1)
            m = RE_PUBCHEM.search(annot_str)
            if m:
                pubchem_cid = m.group(1)
            m = RE_CHEBI.search(annot_str)
            if m:
                chebi_id = 'CHEBI:' + m.group(1)

        # Also try extracting ChEBI from notes
        if not chebi_id:
            notes = sp.find(f'{{{ns}}}notes')
            if notes is not None:
                notes_str = etree.tostring(notes, encoding='unicode')
                m_chebi = re.search(r'ChEBIID:(\d+)', notes_str)
                if m_chebi:
                    chebi_id = 'CHEBI:' + m_chebi.group(1)

        species_map[sid] = {
            'metabolite_name': sname,
            'hmdb_id': hmdb_id,
            'kegg_id': kegg_id,
            'pubchem_cid': pubchem_cid,
            'chebi_id': chebi_id,
        }

    # Find exchange reactions (R_EX_*)
    exchange_metabolites = []
    seen_species = set()

    for rxn in model.findall(f'.//{{{ns}}}reaction'):
        rid = rxn.get('id', '')
        if not rid.startswith('R_EX_'):
            continue

        # Get the single reactant species
        lr = rxn.find(f'{{{ns}}}listOfReactants')
        if lr is not None:
            for sr in lr.findall(f'{{{ns}}}speciesReference'):
                sp_id = sr.get('species', '')
                if sp_id in species_map and sp_id not in seen_species:
                    seen_species.add(sp_id)
                    exchange_metabolites.append(species_map[sp_id])

        # Some exchange reactions might only have products
        lp = rxn.find(f'{{{ns}}}listOfProducts')
        if lp is not None:
            for sr in lp.findall(f'{{{ns}}}speciesReference'):
                sp_id = sr.get('species', '')
                if sp_id in species_map and sp_id not in seen_species:
                    seen_species.add(sp_id)
                    exchange_metabolites.append(species_map[sp_id])

    return organism_name, taxonomy_id, exchange_metabolites


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("AGORA2 → MMI Expansion")
    print("=" * 70)

    if not os.path.isdir(SBML_DIR):
        print(f"ERROR: SBML directory not found: {SBML_DIR}")
        sys.exit(1)

    # List SBML files
    sbml_files = sorted(f for f in os.listdir(SBML_DIR) if f.endswith('.xml'))
    print(f"Found {len(sbml_files)} AGORA2 SBML files")

    # -----------------------------------------------------------------------
    # Parse all SBML files
    # -----------------------------------------------------------------------
    new_rows = []
    organisms_parsed = 0
    total_exchanges = 0
    annotated_exchanges = 0

    t0 = time.time()
    for i, fname in enumerate(sbml_files):
        fpath = os.path.join(SBML_DIR, fname)
        org_name, tax_id, metabolites = parse_agora2_sbml(fpath)

        if org_name is None:
            continue

        organisms_parsed += 1

        for met in metabolites:
            total_exchanges += 1
            # Only include metabolites with at least one identifier that can
            # cross-reference to our DB (HMDB, KEGG, PubChem, or ChEBI)
            has_id = bool(met['hmdb_id'] or met['kegg_id'] or met['pubchem_cid'] or met['chebi_id'])
            if not has_id:
                continue
            annotated_exchanges += 1

            row = {
                'Metabolite_Name': met['metabolite_name'],
                'HMDB_ID': met['hmdb_id'],
                'KEGG_ID': met['kegg_id'],
                'PubChem_CID': met['pubchem_cid'],
                'ChEBI_ID': met['chebi_id'],
                'SMILES': '',
                'Microbe_Name': org_name,
                'Taxonomy_ID': tax_id,
                'Rank': 'strain',
                'Substrate': '',
                'Substrate_PubChem_CID': '',
                'Relationship_Type': 'exchange',
                'Tissue': 'gut',
                'Organism': 'human',
                'Evidence_Level': 'computational',
                'Experimental_Method': 'genome-scale metabolic model',
                'PMID': '30543047',  # AGORA2 paper PMID
                'Source': 'AGORA2',
            }
            new_rows.append(row)

        if (i + 1) % 100 == 0:
            elapsed = time.time() - t0
            print(f"  Parsed {i+1}/{len(sbml_files)} files ({elapsed:.1f}s)")

    elapsed = time.time() - t0
    print(f"\nParsing complete in {elapsed:.1f}s")
    print(f"  Organisms parsed: {organisms_parsed}")
    print(f"  Total exchange metabolites: {total_exchanges}")
    print(f"  With identifiers: {annotated_exchanges}")

    # -----------------------------------------------------------------------
    # Deduplicate new rows (same microbe + metabolite HMDB/KEGG/ChEBI)
    # -----------------------------------------------------------------------
    seen = set()
    deduped = []
    for row in new_rows:
        # Create composite key
        key_parts = [row['Microbe_Name']]
        if row['HMDB_ID']:
            key_parts.append(row['HMDB_ID'])
        elif row['KEGG_ID']:
            key_parts.append(row['KEGG_ID'])
        elif row['ChEBI_ID']:
            key_parts.append(row['ChEBI_ID'])
        elif row['PubChem_CID']:
            key_parts.append(row['PubChem_CID'])
        else:
            key_parts.append(row['Metabolite_Name'])
        key = tuple(key_parts)
        if key not in seen:
            seen.add(key)
            deduped.append(row)

    print(f"  After deduplication: {len(deduped)} new rows")

    # -----------------------------------------------------------------------
    # Load existing database
    # -----------------------------------------------------------------------
    existing = []
    existing_keys = set()
    if os.path.isfile(MMI_DB):
        with open(MMI_DB, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                existing.append(row)
                # Build existing key for dedup against new rows
                key_parts = [row.get('Microbe_Name', '')]
                if row.get('HMDB_ID'):
                    key_parts.append(row['HMDB_ID'])
                elif row.get('KEGG_ID'):
                    key_parts.append(row['KEGG_ID'])
                elif row.get('ChEBI_ID'):
                    key_parts.append(row['ChEBI_ID'])
                elif row.get('PubChem_CID'):
                    key_parts.append(row['PubChem_CID'])
                else:
                    key_parts.append(row.get('Metabolite_Name', ''))
                existing_keys.add(tuple(key_parts))
        print(f"  Existing MMI rows: {len(existing)}")

    # -----------------------------------------------------------------------
    # Merge: only add rows not already in existing DB
    # -----------------------------------------------------------------------
    truly_new = []
    for row in deduped:
        key_parts = [row['Microbe_Name']]
        if row['HMDB_ID']:
            key_parts.append(row['HMDB_ID'])
        elif row['KEGG_ID']:
            key_parts.append(row['KEGG_ID'])
        elif row['ChEBI_ID']:
            key_parts.append(row['ChEBI_ID'])
        elif row['PubChem_CID']:
            key_parts.append(row['PubChem_CID'])
        else:
            key_parts.append(row['Metabolite_Name'])
        if tuple(key_parts) not in existing_keys:
            truly_new.append(row)

    print(f"  Truly new (not in existing): {len(truly_new)}")

    # -----------------------------------------------------------------------
    # Stats summary
    # -----------------------------------------------------------------------
    unique_microbes = len(set(r['Microbe_Name'] for r in truly_new))
    unique_mets_hmdb = len(set(r['HMDB_ID'] for r in truly_new if r['HMDB_ID']))
    unique_mets_kegg = len(set(r['KEGG_ID'] for r in truly_new if r['KEGG_ID']))
    print(f"\n  New unique microbes: {unique_microbes}")
    print(f"  New unique metabolites (HMDB): {unique_mets_hmdb}")
    print(f"  New unique metabolites (KEGG): {unique_mets_kegg}")

    # -----------------------------------------------------------------------
    # Write merged database
    # -----------------------------------------------------------------------
    # Backup existing
    if os.path.isfile(MMI_DB):
        import shutil
        shutil.copy2(MMI_DB, MMI_BAK)
        print(f"\n  Backup saved to {MMI_BAK}")

    columns = [
        'Metabolite_Name', 'HMDB_ID', 'KEGG_ID', 'PubChem_CID', 'ChEBI_ID',
        'SMILES', 'Microbe_Name', 'Taxonomy_ID', 'Rank', 'Substrate',
        'Substrate_PubChem_CID', 'Relationship_Type', 'Tissue', 'Organism',
        'Evidence_Level', 'Experimental_Method', 'PMID', 'Source'
    ]

    merged = existing + truly_new
    with open(MMI_OUT, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(merged)

    total = len(merged)
    print(f"\n{'='*70}")
    print(f"  FINAL MMI DATABASE: {total} rows")
    print(f"    Existing curated: {len(existing)}")
    print(f"    New from AGORA2:  {len(truly_new)}")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
