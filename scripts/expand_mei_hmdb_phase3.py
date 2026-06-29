#!/usr/bin/env python3
"""
expand_mei_hmdb_phase3.py
=========================
Phase 3: Comprehensive generic/macromolecule classification + extended
HMDB mapping via fuzzy name matching and HMDB XML synonym re-search.

Strategy:
1. Pattern-based generic classification (regex) for protein residues,
   macromolecules, generic compound classes, carrier proteins, etc.
2. Extended manual KEGG→HMDB curation for remaining specific compounds
3. Fuzzy/partial name matching against full HMDB synonym dictionary
4. Report accurate effective coverage
"""

import os, csv, re, time
from lxml import etree

BASE = os.path.join(os.path.dirname(__file__), '..')
MEI_DB = os.path.join(BASE, 'data', 'databases', 'mei_database.csv')
HMDB_XML = os.path.join(BASE, 'data', 'raw', 'hmdb', 'hmdb_metabolites.xml')
NS = '{http://www.hmdb.ca}'

# ─── Pattern-based generic classification ────────────────────────────────
# These patterns identify metabolites that are generic compound classes,
# macromolecules, protein residues, or enzyme complexes that inherently
# do not have HMDB small-molecule entries.
GENERIC_PATTERNS = [
    # Protein-related
    r'^Protein\b',                # Protein, Protein tyrosine, etc.
    r'^\[Protein\]',             # [Protein]-L-cysteine, etc.
    r'^Phosphoprotein',          # Phosphoprotein
    r'^Protamine',               # Protamine
    r'^O-Phosphoprotamine',      # Modified protamine
    r'^S-Palmitoylprotein',      # Lipidated protein
    r'^Histone',                 # Histone residues
    r'^Ubiquitin',               # Ubiquitin, Ubiquitin C-terminal thiolester
    r'carrier protein',          # Acyl-carrier protein, Acyl-[acyl-carrier protein]
    r'\[acp\]',                  # Hexadecanoyl-[acp], etc.
    
    # DNA/RNA modifications
    r'\bDNA\b',                  # DNA, 5-Hydroxymethylcytosine in DNA
    r'^RNA$',                    # RNA
    r'^mRNA$',                   # mRNA
    r'^tRNA',                    # tRNA variants
    r'^rRNA',                    # rRNA
    
    # Enzyme complexes and cofactor forms
    r'^\[Oxidized ',             # [Oxidized NADPH-hemoprotein reductase]
    r'^\[Reduced ',              # [Reduced NADPH-hemoprotein reductase]
    r'^Oxidized flavodoxin',     # Oxidized flavodoxin
    r'^Reduced flavodoxin',      # Reduced flavodoxin
    r'^Thioredoxin',             # Thioredoxin/Thioredoxin disulfide (protein)
    r'^Ferrocytochrome',         # Ferrocytochrome b5
    r'^Ferricytochrome',         # Ferricytochrome b5
    r'^Cytochrome',              # Cytochrome c variants
    r'^Glutaredoxin',            # Glutaredoxin (protein)
    r'^Reduced ferredoxin',      # Ferredoxin (protein)
    r'^Oxidized ferredoxin',     # Ferredoxin (protein)
    
    # Generic placeholders
    r'^Acceptor$',               # Generic electron acceptor
    r'^Reduced acceptor$',       # Generic reduced acceptor
    r'^Donor$',                  # Generic donor
    r'^ROH$',                    # Generic alcohol
    r'^RX$',                     # Generic halide
    r'^R\'H$',                   # Generic substrate
    r'^RCO-',                    # Generic acyl groups
    r'^NDP$',                    # Generic nucleoside diphosphate
    r'^NTP$',                    # Generic nucleoside triphosphate
    
    # Generic compound classes (not specific enough for HMDB)
    r'^Nitrile$',                # Generic nitrile class
    r'^Halide$',                 # Generic halide class
    r'^Alkane$',                 # Generic alkane class
    r'^Aldehyde$',               # Generic aldehyde class
    r'^Primary alcohol$',        # Generic alcohol class
    r'^Secondary alcohol$',      # Generic
    r'^Long-chain fatty acid$',  # Generic fatty acid class
    r'^Long-chain aldehyde$',    # Generic
    r'^Long-chain alcohol$',     # Generic
    r'^Medium-chain fatty acid$', # Generic
    r'^alpha-Hydroxy fatty acid$', # Generic fatty acid class
    r'^Sulfuric monoester$',     # Generic sulfate ester
    r'^Hepatotoxins$',           # Generic class
    r'^Nucleoside triphosphate$', # Generic NTP
    r'^Deoxynucleoside triphosphate$', # Generic dNTP
    r'^2-Oxo acid$',             # Generic 2-oxo acid class
    r'^3-Oxoacyl-CoA$',          # Generic oxoacyl-CoA class
    r'^2,3-Dehydroacyl-CoA$',    # Generic dehydroacyl-CoA
    r'^trans-2,3-Dehydroacyl-CoA$', # Generic
    r'^Acyl-CoA$',               # Generic acyl-CoA
    r'^Long-chain acyl-CoA$',    # Generic acyl-CoA
    r'^Short-chain acyl-CoA$',   # Generic acyl-CoA
    r'^3-Oxo acid$',             # Generic
    r'^Omega-hydroxy fatty acid$', # Generic
    r'^Primary amine$',          # Generic
    r'^Secondary amine$',        # Generic
    r'^Sterol$',                 # Generic
    r'^Phenol$',                 # Wait, phenol is HMDB0000228 - NOT generic
    r'^D-Hexose$',               # Generic hexose
    r'^D-Hexose 6-phosphate$',   # Generic
    r'^Glycolipid$',             # Generic
    r'^Sphingolipid$',           # Generic class
    r'^Peptidylproline',         # Peptide residues
    r'^N-Acylsphingosine$',      # Generic ceramide
    r'^1-Acyl-sn-glycerol$',     # Generic MAG
    r'^1,2-Diacyl-sn-glycerol$', # Generic DAG
    r'^1-Acyl-sn-glycerol 3-phosphate$', # Generic LPA
    r'^1-Acylglycerol$',         # Generic MAG
    r'^1-\(1-Alkenyl\)-sn-glycero-3-phospho', # Generic plasmalogen classes
    r'^1-Organyl-2-acyl-sn-glycero-3-phospho', # Generic phospholipid
    r'^1-Radyl-2-acyl-sn-glycero', # Generic phospholipid
    r'^2-Acyl-1-alkyl-sn-glycero', # Generic ether lipid
    r'^Amide$',                  # Generic
    r'^Carboxylate$',            # Generic
    r'^Thiol$',                  # Generic
    r'^An alcohol$',             # Generic
    r'^A ketone$',               # Generic
    r'^Carboxylic acid$',        # Too generic (but specific HMDB name matches exist)
    r'^Heme$',                   # HMDB has heme actually... but as a cofactor class? Let me check
]

# Remove patterns that match real HMDB metabolites
# Phenol IS in HMDB: HMDB0000228, don't classify as generic
GENERIC_PATTERNS = [p for p in GENERIC_PATTERNS if 'Phenol' not in p and 'Heme' not in p and 'Carboxylic acid' not in p]

# Compile all patterns
_generic_re = [re.compile(p, re.IGNORECASE) for p in GENERIC_PATTERNS]

def is_generic(name):
    """Check if a metabolite name matches any generic pattern."""
    for pat in _generic_re:
        if pat.search(name):
            return True
    return False


# ─── Extended manual KEGG→HMDB curation ──────────────────────────────────
# These are specific KEGG compounds that have known HMDB matches but were
# not found by automated cross-referencing (different naming, synonyms).
MANUAL_KEGG_TO_HMDB = {
    # Bicarbonate / Carbon dioxide related
    'C00288': 'HMDB0000595',     # HCO3- = Bicarbonate
    
    # Specific nucleotides & derivatives
    'C00201': '',                 # Nucleoside triphosphate — generic
    
    # CoA derivatives (specific)
    'C16169': 'HMDB0013622',     # (4Z,7Z,10Z,13Z,16Z,19Z)-Docosahexaenoyl-CoA = DHA-CoA
    'C00264': '',                 # 3-Oxoacyl-CoA — generic
    'C00605': '',                 # 2,3-Dehydroacyl-CoA — generic
    'C05764': 'HMDB0001338',     # Hexadecanoyl-[acp] → Palmitoyl-CoA (equivalent)
    'C00173': '',                 # Acyl-[acyl-carrier protein] — generic
    
    # Mevalonate pathway
    'C00772': 'HMDB0001387',     # Mevaldate = Mevaldehyde
    
    # Aflatoxins
    'C19595': 'HMDB0060553',     # Aflatoxin B1-endo-8,9-epoxide
    'C19594': 'HMDB0060581',     # Aflatoxin-M1-8,9-epoxide
    
    # NNK/NNAL metabolites
    'C19577': 'HMDB0041823',     # 1-(Methylnitrosoamino)-4-(3-pyridinyl)-1,4-butanediol
    
    # Common metabolites with naming differences
    'C00161': '',                 # 2-Oxo acid — generic
    'C00638': '',                 # Long-chain fatty acid — generic
    'C00226': '',                 # Primary alcohol — generic
    'C01371': '',                 # Alkane — generic
    'C00681': '',                 # 1-Acyl-sn-glycerol 3-phosphate — generic LPA
    'C01885': '',                 # 1-Acylglycerol — generic MAG
    'C05102': '',                 # alpha-Hydroxy fatty acid — generic
    'C21758': '',                 # Hepatotoxins — generic
    'C01322': '',                 # RX — generic halide
    
    # Additional specific compounds
    'C05546': '',                 # Protein N6,N6,N6-trimethyl-L-lysine — residue
    'C02415': '',                 # Histone-L-lysine — residue
    'C01997': '',                 # Histone N6-acetyl-L-lysine — residue
    'C02189': '',                 # [Protein]-L-serine — residue
    'C00613': '',                 # Protein-L-arginine — residue
    'C00343': '',                 # Thioredoxin disulfide — protein
    'C00342': '',                 # Thioredoxin — protein
    'C00999': '',                 # Ferrocytochrome b5 — protein
    'C00996': '',                 # Ferricytochrome b5 — protein
    'C00229': '',                 # Acyl-carrier protein — protein
    'C17202': '',                 # S-Palmitoylprotein — protein
    'C00454': '',                 # NDP — generic
    'C00046': '',                 # RNA — macro

    # Sphingolipids that have representative HMDB entries
    'C00195': 'HMDB0004947',     # N-Acylsphingosine → Ceramide (d18:1/16:0)
    
    # More lipids
    'C00550': 'HMDB0000037',     # Sphingomyelin → SM
    'C01120': 'HMDB0010400',     # Sphinganine 1-phosphate → Sphinganine 1-phosphate
    
    # Amino acid derivatives
    'C00049': 'HMDB0000191',     # L-Aspartate
    'C04640': 'HMDB0004017',     # 2-Acyl-sn-glycero-3-phosphocholine → LPC(16:0)
}


def build_hmdb_maps():
    """Parse HMDB XML to build comprehensive lookup dictionaries."""
    print("  Parsing HMDB XML for comprehensive synonym matching...")
    t0 = time.time()
    
    name_to_hmdb = {}
    kegg_to_hmdb = {}
    inchikey_to_hmdb = {}
    count = 0
    
    context = etree.iterparse(HMDB_XML, events=('end',), tag=f'{NS}metabolite')
    
    for event, elem in context:
        accession = elem.findtext(f'{NS}accession')
        if not accession:
            elem.clear()
            continue
        
        count += 1
        
        # Primary name
        pname = elem.findtext(f'{NS}name', '')
        if pname:
            key = pname.strip().lower()
            if key not in name_to_hmdb:
                name_to_hmdb[key] = accession
        
        # IUPAC name
        iupac = elem.findtext(f'{NS}iupac_name', '')
        if iupac:
            key = iupac.strip().lower()
            if key not in name_to_hmdb:
                name_to_hmdb[key] = accession
        
        # Traditional name
        trad = elem.findtext(f'{NS}traditional_iupac', '')
        if trad:
            key = trad.strip().lower()
            if key not in name_to_hmdb:
                name_to_hmdb[key] = accession
        
        # Synonyms
        syns_el = elem.find(f'{NS}synonyms')
        if syns_el is not None:
            for syn in syns_el.findall(f'{NS}synonym'):
                if syn.text:
                    key = syn.text.strip().lower()
                    if key not in name_to_hmdb:
                        name_to_hmdb[key] = accession
        
        # KEGG ID
        kegg_el = elem.findtext(f'{NS}kegg_id', '')
        if kegg_el and kegg_el.strip():
            kegg_to_hmdb[kegg_el.strip()] = accession
        
        # InChIKey
        inchikey = elem.findtext(f'{NS}inchikey', '')
        if inchikey and inchikey.strip():
            inchikey_to_hmdb[inchikey.strip()] = accession
        
        # Memory cleanup
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
    
    del context
    
    elapsed = time.time() - t0
    print(f"    Parsed {count:,} metabolites in {elapsed:.1f}s")
    print(f"    Name/synonym→HMDB: {len(name_to_hmdb):,}")
    print(f"    KEGG→HMDB: {len(kegg_to_hmdb):,}")
    print(f"    InChIKey→HMDB: {len(inchikey_to_hmdb):,}")
    
    return name_to_hmdb, kegg_to_hmdb, inchikey_to_hmdb


def normalize_name(name):
    """
    Generate normalized name variants for fuzzy matching.
    Returns list of candidate strings to match against HMDB synonyms.
    """
    variants = []
    n = name.strip().lower()
    variants.append(n)
    
    # Remove stereochemistry prefixes
    for prefix in ['(s)-', '(r)-', '(+)-', '(-)-', 'l-', 'd-', 'trans-', 'cis-',
                    'alpha-', 'beta-', 'gamma-', 'delta-', 'omega-']:
        if n.startswith(prefix):
            variants.append(n[len(prefix):])
    
    # Remove trailing charge indicators
    for suffix in ['(2-)', '(1-)', '(3-)', '(2+)', '(1+)', '(3+)']:
        if n.endswith(suffix):
            variants.append(n[:-len(suffix)].strip())
    
    # Remove text in brackets at start
    m = re.match(r'^\[.*?\]\s*', n)
    if m:
        variants.append(n[m.end():])
    
    # Replace common naming differences
    n2 = n.replace('-', ' ').replace('_', ' ')
    variants.append(n2)
    
    # Handle "ate" vs "ic acid" endings
    if n.endswith('ate'):
        variants.append(n[:-3] + 'ic acid')
    if n.endswith('ic acid'):
        variants.append(n[:-7] + 'ate')
    
    return list(set(variants))


def main():
    print("=" * 70)
    print("MEI HMDB Coverage Expansion — Phase 3 (Comprehensive)")
    print("=" * 70)
    
    # Build HMDB lookup maps
    name_to_hmdb, kegg_to_hmdb, inchikey_to_hmdb = build_hmdb_maps()
    
    # Merge manual KEGG mappings
    for k, v in MANUAL_KEGG_TO_HMDB.items():
        if v:
            kegg_to_hmdb[k] = v
    
    print("\n  Processing MEI database...")
    
    rows = []
    total = 0
    already_has = 0
    filled_kegg = 0
    filled_name_exact = 0
    filled_name_fuzzy = 0
    classified_generic = 0
    still_missing = 0
    
    # Track what we fill for reporting
    fill_examples = {'kegg': [], 'name_exact': [], 'name_fuzzy': []}
    generic_examples = []
    missing_examples = []
    
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
            
            met_name = (row.get('Metabolite_Name') or '').strip()
            kegg_id = (row.get('KEGG_Compound') or '').strip()
            
            # Step 1: Check if generic/macromolecule
            if is_generic(met_name):
                classified_generic += 1
                rows.append(row)
                continue
            
            # Step 2: Check manual KEGG mapping (empty = classified generic already)
            if kegg_id and kegg_id in MANUAL_KEGG_TO_HMDB:
                mapped = MANUAL_KEGG_TO_HMDB[kegg_id]
                if mapped:
                    row['HMDB_ID'] = mapped
                    filled_kegg += 1
                    if len(fill_examples['kegg']) < 5:
                        fill_examples['kegg'].append(f"{met_name} ({kegg_id}) → {mapped}")
                else:
                    classified_generic += 1
                rows.append(row)
                continue
            
            # Step 3: Try KEGG→HMDB from HMDB XML
            if kegg_id and kegg_id in kegg_to_hmdb:
                row['HMDB_ID'] = kegg_to_hmdb[kegg_id]
                filled_kegg += 1
                rows.append(row)
                continue
            
            # Step 4: Exact name match
            name_lower = met_name.lower()
            if name_lower in name_to_hmdb:
                row['HMDB_ID'] = name_to_hmdb[name_lower]
                filled_name_exact += 1
                rows.append(row)
                continue
            
            # Step 5: Fuzzy name matching
            matched = False
            for variant in normalize_name(met_name):
                if variant in name_to_hmdb:
                    row['HMDB_ID'] = name_to_hmdb[variant]
                    filled_name_fuzzy += 1
                    if len(fill_examples['name_fuzzy']) < 10:
                        fill_examples['name_fuzzy'].append(
                            f"{met_name} → '{variant}' → {name_to_hmdb[variant]}")
                    matched = True
                    break
            
            if not matched:
                still_missing += 1
                if len(missing_examples) < 30:
                    missing_examples.append(f"{met_name} (KEGG={kegg_id})")
            
            rows.append(row)
    
    # Write updated database
    print("  Writing updated database...")
    with open(MEI_DB, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    
    # Report
    filled_total = filled_kegg + filled_name_exact + filled_name_fuzzy
    new_has = already_has + filled_total
    new_coverage = new_has / total * 100
    hmdb_eligible = total - classified_generic
    effective = new_has / hmdb_eligible * 100 if hmdb_eligible > 0 else 0
    
    print(f"\n{'='*70}")
    print(f"  RESULTS")
    print(f"{'='*70}")
    print(f"  Total MEI records:        {total:,}")
    print(f"  Already had HMDB_ID:      {already_has:,}")
    print(f"  Filled via KEGG mapping:  {filled_kegg:,}")
    print(f"  Filled via exact name:    {filled_name_exact:,}")
    print(f"  Filled via fuzzy name:    {filled_name_fuzzy:,}")
    print(f"  Classified as generic:    {classified_generic:,}")
    print(f"  Still missing:            {still_missing:,}")
    print(f"\n  Overall HMDB coverage:    {new_coverage:.1f}% ({new_has:,} / {total:,})")
    print(f"  HMDB-eligible records:    {hmdb_eligible:,}")
    print(f"  Effective coverage:       {effective:.1f}% ({new_has:,} / {hmdb_eligible:,})")
    
    if fill_examples['kegg']:
        print(f"\n  KEGG fill examples:")
        for ex in fill_examples['kegg']:
            print(f"    {ex}")
    
    if fill_examples['name_fuzzy']:
        print(f"\n  Fuzzy name match examples:")
        for ex in fill_examples['name_fuzzy']:
            print(f"    {ex}")
    
    if missing_examples:
        print(f"\n  Still missing examples ({still_missing:,} total):")
        for ex in missing_examples[:20]:
            print(f"    {ex}")


if __name__ == '__main__':
    main()
