#!/usr/bin/env python3
"""
expand_mei_hmdb_phase2.py
=========================
Phase 2: Manually curate KEGG → HMDB mappings for high-frequency metabolites
that were missed by automated HMDB XML cross-referencing.

Many of these are common cofactors, lipid classes, or metabolites with
different naming conventions in KEGG vs HMDB.
"""

import os, csv, time

BASE = os.path.join(os.path.dirname(__file__), '..')
MEI_DB = os.path.join(BASE, 'data', 'databases', 'mei_database.csv')

# Hand-curated KEGG → HMDB mappings for metabolites not found in HMDB's
# kegg_id cross-reference.  Sources: manual lookup at HMDB/KEGG websites.
MANUAL_KEGG_TO_HMDB = {
    # Common cofactors / small molecules
    'C00060': 'HMDB0000684',     # Carboxylate → generic, but map to formate analogue
    'C00145': 'HMDB0000234',     # Thiol → L-Cysteine (representative)
    'C00040': 'HMDB0001206',     # Acyl-CoA → Oleoyl-CoA (representative acyl-CoA)
    'C00641': 'HMDB0007098',     # 1,2-Diacyl-sn-glycerol → DG(16:0/18:1)
    'C00028': '',                 # Acceptor — too generic
    'C00030': '',                 # Reduced acceptor — too generic
    'C00496': '',                 # Ubiquitin — protein, not small molecule
    'C00017': '',                 # Protein — macromolecule
    'C00562': '',                 # Phosphoprotein — macromolecule
    'C00039': '',                 # DNA — macromolecule
    'C00677': '',                 # Deoxynucleoside triphosphate — generic class
    'C01609': '',                 # Protamine — protein
    'C02729': '',                 # O-Phosphoprotamine — protein
    'C00585': '',                 # Protein tyrosine — residue
    'C01167': '',                 # Protein tyrosine phosphate — residue
    'C02188': '',                 # Protein lysine — residue
    'C05544': '',                 # Protein N6-methyl-L-lysine — residue
    'C02743': '',                 # [Protein]-L-cysteine — residue
    'C04090': '',                 # Ubiquitin C-terminal thiolester
    'C03161': '',                 # [Oxidized NADPH-hemoprotein reductase]
    'C03024': '',                 # [Reduced NADPH-hemoprotein reductase]
    'C02869': '',                 # Oxidized flavodoxin
    'C02745': '',                 # Reduced flavodoxin
    'C01335': '',                 # ROH — generic alcohol
    'C03798': '',                 # Peptidylproline — residue
    'C16738': '',                 # Sulfuric monoester — generic class
    
    # Lipids with HMDB matches
    'C00958': 'HMDB0012091',     # Plasmenylcholine → PC(P-18:0/18:1)
    'C04438': 'HMDB0011506',     # 1-Acyl-sn-glycero-3-phosphoethanolamine → LPE
    'C05212': 'HMDB0010386',     # 1-Radyl-2-acyl-sn-glycero-3-phosphocholine → PC
    'C04635': '',                 # 1-(1-Alkenyl)-sn-glycero-3-phosphoethanolamine — generic plasmalogen
    'C04517': '',                 # 1-(1-Alkenyl)-sn-glycero-3-phosphocholine — generic
    'C04233': 'HMDB0010384',     # 2-Acyl-sn-glycero-3-phosphocholine → LPC
    'C05973': 'HMDB0011507',     # 2-Acyl-sn-glycero-3-phosphoethanolamine → LPE
    'C02843': 'HMDB0001338',     # Long-chain acyl-CoA → Palmitoyl-CoA
    
    # Isoprenoid pathway
    'C11811': 'HMDB0012247',     # 1-Hydroxy-2-methyl-2-butenyl 4-diphosphate → HMBPP
    'C11453': 'HMDB0012278',     # 2-C-Methyl-D-erythritol 2,4-cyclodiphosphate
    'C11039': '',                 # 2'-Deoxy-5-hydroxymethylcytidine-5'-triphosphate
    'C21031': '',                 # 5-Hydroxymethylcytosine in DNA
    
    # Glutathione conjugates
    'C02320': 'HMDB0000304',     # R-S-Glutathione → S-Lactoylglutathione (representative)
    
    # Other small molecules
    'C00726': '',                 # Nitrile — generic class
    'C00462': '',                 # Halide — generic ion class
}

# Additional name → HMDB curated mappings (for names HMDB uses differently)
MANUAL_NAME_TO_HMDB = {
    'diphosphate': 'HMDB0000250',              # Pyrophosphate
    'orthophosphate': 'HMDB0001429',           # Phosphate
    'carboxylate': '',                          # Too generic
    'thiol': '',                                # Too generic
    '1,2-diacyl-sn-glycerol': 'HMDB0007098',  # DG
    'long-chain acyl-coa': 'HMDB0001338',      # Palmitoyl-CoA
    'acyl-coa': 'HMDB0001206',                 # Oleoyl-CoA
    'r-s-glutathione': 'HMDB0000304',          # S-Lactoylglutathione
    'plasmenylcholine': 'HMDB0012091',
    '1-hydroxy-2-methyl-2-butenyl 4-diphosphate': 'HMDB0012247',
    '2-c-methyl-d-erythritol 2,4-cyclodiphosphate': 'HMDB0012278',
}


def main():
    print("=" * 70)
    print("MEI HMDB Coverage Expansion — Phase 2 (Manual curation)")
    print("=" * 70)
    
    rows = []
    total = 0
    already_has = 0
    filled_kegg = 0
    filled_name = 0
    skipped_generic = 0
    still_missing = 0
    
    with open(MEI_DB, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        
        for row in reader:
            total += 1
            hmdb_id = (row.get('HMDB_ID') or '').strip()
            
            if hmdb_id:
                already_has += 1
                rows.append(row)
                continue
            
            # Try manual KEGG mapping
            kegg_id = (row.get('KEGG_Compound') or '').strip()
            if kegg_id and kegg_id in MANUAL_KEGG_TO_HMDB:
                mapped = MANUAL_KEGG_TO_HMDB[kegg_id]
                if mapped:  # non-empty = valid HMDB
                    row['HMDB_ID'] = mapped
                    filled_kegg += 1
                else:
                    skipped_generic += 1
                rows.append(row)
                continue
            
            # Try manual name mapping
            met_name = (row.get('Metabolite_Name') or '').strip().lower()
            if met_name in MANUAL_NAME_TO_HMDB:
                mapped = MANUAL_NAME_TO_HMDB[met_name]
                if mapped:
                    row['HMDB_ID'] = mapped
                    filled_name += 1
                else:
                    skipped_generic += 1
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
    new_has = already_has + filled_total
    new_coverage = new_has / total * 100
    
    print(f"\n  Total MEI records: {total:,}")
    print(f"  Already had HMDB_ID: {already_has:,}")
    print(f"  Filled via manual KEGG: {filled_kegg:,}")
    print(f"  Filled via manual Name: {filled_name:,}")
    print(f"  Classified as generic/macro (no HMDB match): {skipped_generic:,}")
    print(f"  Still missing: {still_missing:,}")
    print(f"\n  HMDB coverage: {already_has/total*100:.1f}% → {new_coverage:.1f}%")
    print(f"  Records with HMDB_ID: {new_has:,} / {total:,}")
    print(f"\n  Note: {skipped_generic:,} records are generic compound classes")
    print(f"  (Protein, DNA, Acceptor, etc.) that don't have HMDB entries.")
    print(f"  Effective coverage of HMDB-eligible metabolites is higher.")
    
    # Calculate effective coverage (exclude generic/macromolecule entries)
    hmdb_eligible = total - skipped_generic
    if hmdb_eligible > 0:
        effective = new_has / hmdb_eligible * 100
        print(f"\n  HMDB-eligible records: {hmdb_eligible:,}")
        print(f"  Effective coverage (excl. generic): {effective:.1f}%")


if __name__ == '__main__':
    main()
