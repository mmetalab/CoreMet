#!/usr/bin/env python3
"""Extract drug SMILES from DrugBank XML for MDrI model training."""
import os, csv, time
from lxml import etree

BASE = os.path.join(os.path.dirname(__file__), '..')
DRUGBANK_XML = os.path.join(BASE, 'data', 'raw', 'drugbank', 'drug_bank_database.xml')
OUTPUT = os.path.join(BASE, 'data', 'processed', 'drug_smiles.csv')
NS = '{http://www.drugbank.ca}'

def main():
    print("Extracting drug SMILES from DrugBank XML...")
    t0 = time.time()
    
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
    
    results = []
    count = 0
    context = etree.iterparse(DRUGBANK_XML, events=('end',), tag=f'{NS}drug')
    
    for event, elem in context:
        # Only top-level drug elements (check for drugbank-id child)
        dbids = elem.findall(f'{NS}drugbank-id')
        if not dbids:
            elem.clear()
            continue
            
        primary_id = None
        for dbid in dbids:
            if dbid.get('primary') == 'true':
                primary_id = dbid.text
                break
        if not primary_id:
            primary_id = dbids[0].text if dbids[0].text else None
        
        if not primary_id:
            elem.clear()
            continue
            
        count += 1
        
        # Get name
        name = elem.findtext(f'{NS}name', '')
        
        # Get SMILES from calculated-properties or experimental-properties
        smiles = ''
        for prop_group in [f'{NS}calculated-properties', f'{NS}experimental-properties']:
            props = elem.find(prop_group)
            if props is not None:
                for prop in props.findall(f'{NS}property'):
                    kind = prop.findtext(f'{NS}kind', '')
                    if kind == 'SMILES':
                        smiles = prop.findtext(f'{NS}value', '')
                        if smiles:
                            break
            if smiles:
                break
        
        if smiles:
            results.append({
                'DrugBank_ID': primary_id,
                'Drug_Name': name,
                'SMILES': smiles,
            })
        
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    
    del context
    
    # Write CSV
    with open(OUTPUT, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['DrugBank_ID', 'Drug_Name', 'SMILES'])
        writer.writeheader()
        writer.writerows(results)
    
    elapsed = time.time() - t0
    print(f"Parsed {count:,} drugs in {elapsed:.1f}s")
    print(f"Extracted SMILES for {len(results):,} drugs")
    print(f"Saved to {OUTPUT}")

if __name__ == '__main__':
    main()
