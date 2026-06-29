#!/usr/bin/env python3
"""
build_mgwas_database.py — Build Metabolite–SNP (mGWAS) Database
================================================================

Sources (each row in final CSV carries a 'Source' provenance tag):

  1. GWAS Catalog (EBI)       — genome-wide significant metabolite QTL associations
     URL: https://www.ebi.ac.uk/gwas/api/search/downloads/full
     Tag: "GWAS_Catalog"

  2. Shin et al. 2014 (Nat Genet 46:543)   — 486 metabolites, 145 loci
     PMID: 24816252
     Data: Supplementary Table 1 (significant loci)
     Tag: "Shin2014"

  3. Long et al. 2017 (Nat Genet 49:568)   — 147 metabolites, multi-ethnic
     PMID: 28263315
     Data: Supplementary Table (significant loci)
     Tag: "Long2017"

  4. Chen et al. 2023 (Nat Genet 55:1903)  — 1,091 metabolites, 248 loci
     PMID: 36635386
     Data: Supplementary Data
     Tag: "Chen2023"

  5. Yin et al. 2024 (Nature 628:130)      — 233 NMR traits, 443 loci
     PMID: 38448586
     Data: Supplementary Table 5 (all lead SNP–trait associations)
     Tag: "Yin2024"

Pipeline:
  Step 1: Build metabolite synonym table from HMDB XML + existing xref mapping
  Step 2: Download & filter GWAS Catalog for metabolite-related traits
  Step 3: Parse curated mGWAS supplementary tables (if available locally)
  Step 4: Map metabolite names → HMDB IDs via synonym table
  Step 5: Deduplicate on (HMDB_ID, rsID), keep best p-value, merge sources
  Step 6: Attach SMILES from xref mapping
  Step 7: Output mgwas_database.csv

Usage:
    python scripts/build_mgwas_database.py                    # Full pipeline
    python scripts/build_mgwas_database.py --skip-download    # Use cached GWAS file
    python scripts/build_mgwas_database.py --synonyms-only    # Only build synonym table
"""

import os
import re
import sys
import gzip
import json
import time
import logging
import argparse
from pathlib import Path
from collections import defaultdict

import pandas as pd
import numpy as np

# Optional: for GWAS Catalog download
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
)
logger = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache" / "mgwas"
MAPPINGS_DIR = DATA_DIR / "mappings"
DB_DIR = DATA_DIR / "databases"

HMDB_XML = RAW_DIR / "hmdb" / "hmdb_metabolites.xml"
XREF_CSV = MAPPINGS_DIR / "hmdb_xref_mapping.csv"
SYNONYM_CSV = MAPPINGS_DIR / "hmdb_synonyms.csv"
GWAS_CACHE = CACHE_DIR / "gwas_catalog_associations.tsv"
OUTPUT_CSV = DB_DIR / "mgwas_database.csv"


# ══════════════════════════════════════════════════════════════════════
# Step 1: Build metabolite synonym table from HMDB XML
# ══════════════════════════════════════════════════════════════════════

def build_synonym_table(force=False):
    """Extract metabolite names + synonyms from HMDB XML.

    Uses streaming XML parsing to handle the 6 GB file efficiently.
    Output: hmdb_synonyms.csv with columns [hmdb_id, name, synonym_type]
    where synonym_type ∈ {primary, synonym, iupac, traditional, cas}.
    """
    if SYNONYM_CSV.exists() and not force:
        logger.info(f"Synonym table exists: {SYNONYM_CSV} ({_count_lines(SYNONYM_CSV)-1:,} entries)")
        return pd.read_csv(SYNONYM_CSV, dtype=str)

    if not HMDB_XML.exists():
        logger.warning(f"HMDB XML not found at {HMDB_XML}; building from xref mapping only")
        return _build_synonym_from_xref()

    logger.info(f"Parsing HMDB XML ({HMDB_XML}, ~6 GB) for synonyms — streaming mode...")
    from lxml import etree

    NS = '{http://www.hmdb.ca}'
    rows = []
    count = 0

    # Streaming parse
    context = etree.iterparse(str(HMDB_XML), events=('end',), tag=f'{NS}metabolite')
    for event, elem in context:
        # Extract HMDB accession
        acc_el = elem.find(f'{NS}accession')
        if acc_el is None or not acc_el.text:
            elem.clear()
            continue
        hmdb_id = acc_el.text.strip()

        # Primary name
        name_el = elem.find(f'{NS}name')
        if name_el is not None and name_el.text:
            rows.append((hmdb_id, name_el.text.strip().lower(), 'primary'))

        # IUPAC name
        iupac_el = elem.find(f'{NS}iupac_name')
        if iupac_el is not None and iupac_el.text:
            rows.append((hmdb_id, iupac_el.text.strip().lower(), 'iupac'))

        # Traditional IUPAC
        trad_el = elem.find(f'{NS}traditional_iupac')
        if trad_el is not None and trad_el.text:
            rows.append((hmdb_id, trad_el.text.strip().lower(), 'traditional'))

        # CAS number
        cas_el = elem.find(f'{NS}cas_registry_number')
        if cas_el is not None and cas_el.text:
            rows.append((hmdb_id, cas_el.text.strip(), 'cas'))

        # Synonyms
        syn_parent = elem.find(f'{NS}synonyms')
        if syn_parent is not None:
            for syn_el in syn_parent.findall(f'{NS}synonym'):
                if syn_el.text:
                    rows.append((hmdb_id, syn_el.text.strip().lower(), 'synonym'))

        # Chemical formula — useful for matching
        formula_el = elem.find(f'{NS}chemical_formula')
        if formula_el is not None and formula_el.text:
            rows.append((hmdb_id, formula_el.text.strip(), 'formula'))

        count += 1
        if count % 10000 == 0:
            logger.info(f"  Parsed {count:,} metabolites...")

        # Free memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    df = pd.DataFrame(rows, columns=['hmdb_id', 'name', 'synonym_type'])
    df = df.drop_duplicates(subset=['hmdb_id', 'name'])

    MAPPINGS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(SYNONYM_CSV, index=False)
    logger.info(f"Synonym table saved: {SYNONYM_CSV} ({len(df):,} entries for {df['hmdb_id'].nunique():,} metabolites)")
    return df


def _build_synonym_from_xref():
    """Fallback: build a minimal synonym table from xref mapping."""
    xref = pd.read_csv(XREF_CSV, dtype=str)
    rows = []
    for _, r in xref.iterrows():
        hmdb_id = r['hmdb_id']
        name = r.get('metabolite_name', '')
        if pd.notna(name) and name.strip():
            rows.append((hmdb_id, name.strip().lower(), 'primary'))
    df = pd.DataFrame(rows, columns=['hmdb_id', 'name', 'synonym_type'])
    MAPPINGS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(SYNONYM_CSV, index=False)
    logger.info(f"Fallback synonym table: {len(df):,} entries")
    return df


# ══════════════════════════════════════════════════════════════════════
# Step 2: Download & filter GWAS Catalog
# ══════════════════════════════════════════════════════════════════════

# Keywords that indicate metabolite-related GWAS traits
METABOLITE_TRAIT_KEYWORDS = [
    # General metabolomics
    'metabolite', 'metabolic', 'metabolomics', 'metabotype', 'metabolome',
    # Lipids
    'lipid', 'cholesterol', 'triglyceride', 'hdl', 'ldl', 'vldl',
    'sphingolipid', 'sphingomyelin', 'ceramide', 'glycerophospholipid',
    'phosphatidylcholine', 'lysophosphatidylcholine', 'phospholipid',
    'fatty acid', 'acylcarnitine', 'lipoprotein',
    # Amino acids
    'amino acid', 'alanine', 'arginine', 'asparagine', 'aspartate',
    'cysteine', 'glutamate', 'glutamine', 'glycine', 'histidine',
    'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine',
    'proline', 'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine',
    'homocysteine', 'creatinine', 'creatine',
    # Sugars / energy
    'glucose', 'fructose', 'lactate', 'pyruvate', 'citrate', 'succinate',
    'glycerol',
    # Nucleotides / purines
    'urate', 'uric acid', 'xanthine', 'hypoxanthine',
    # Vitamins / cofactors
    'vitamin', 'folate', 'betaine', 'choline',
    # Bile acids
    'bile acid', 'taurine', 'glycocholate',
    # Hormones (steroid metabolites)
    'cortisol', 'testosterone', 'estradiol', 'progesterone', 'aldosterone',
    # Other
    'carnitine', 'acetylcarnitine', 'bilirubin', 'albumin',
    'lipoprotein(a)', 'apolipoprotein',
    # Specific measurement terms
    'blood levels', 'serum levels', 'plasma levels', 'urinary levels',
    'concentration measurement',
]

# EFO trait URIs that map to metabolite measurements
METABOLITE_EFO_PREFIXES = [
    'EFO_0004529',   # metabolite measurement
    'EFO_0004530',   # lipid measurement
    'EFO_0004469',   # amino acid measurement
    'EFO_0009789',   # metabolomics measurement
]


def download_gwas_catalog(force=False):
    """Download GWAS Catalog full associations TSV."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    if GWAS_CACHE.exists() and not force:
        size_mb = GWAS_CACHE.stat().st_size / (1024 * 1024)
        logger.info(f"GWAS Catalog cache exists: {GWAS_CACHE} ({size_mb:.0f} MB)")
        return GWAS_CACHE

    if not HAS_REQUESTS:
        logger.error("requests package not available; cannot download GWAS Catalog")
        sys.exit(1)

    import zipfile

    url = ("https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
           "gwas-catalog-associations_ontology-annotated-full.zip")
    zip_path = CACHE_DIR / "gwas-catalog-associations_ontology-annotated-full.zip"

    logger.info(f"Downloading GWAS Catalog from EBI FTP ...")
    logger.info(f"  URL: {url}")
    logger.info("  (This may take several minutes)")

    resp = requests.get(url, stream=True, timeout=600)
    resp.raise_for_status()

    total = int(resp.headers.get('content-length', 0))
    downloaded = 0
    with open(zip_path, 'wb') as f:
        for chunk in resp.iter_content(chunk_size=8192 * 16):
            f.write(chunk)
            downloaded += len(chunk)
            if total > 0 and downloaded % (50 * 1024 * 1024) == 0:
                pct = downloaded / total * 100
                logger.info(f"  Downloaded {downloaded / 1e6:.0f} / {total / 1e6:.0f} MB ({pct:.1f}%)")

    size_mb = zip_path.stat().st_size / (1024 * 1024)
    logger.info(f"  Download complete: {size_mb:.0f} MB ZIP → {zip_path}")

    # Extract TSV from ZIP
    logger.info("  Extracting ZIP...")
    with zipfile.ZipFile(zip_path, 'r') as zf:
        tsv_names = [n for n in zf.namelist() if n.endswith('.tsv')]
        if not tsv_names:
            logger.error(f"No TSV files found in ZIP: {zf.namelist()}")
            sys.exit(1)
        tsv_name = tsv_names[0]
        logger.info(f"  Extracting {tsv_name}...")
        zf.extract(tsv_name, CACHE_DIR)
        extracted = CACHE_DIR / tsv_name
        if extracted != GWAS_CACHE:
            extracted.rename(GWAS_CACHE)

    size_mb = GWAS_CACHE.stat().st_size / (1024 * 1024)
    logger.info(f"  Extracted: {size_mb:.0f} MB → {GWAS_CACHE}")

    # Clean up ZIP to save disk space
    zip_path.unlink(missing_ok=True)

    return GWAS_CACHE


def filter_gwas_catalog(gwas_path):
    """Filter GWAS Catalog for metabolite-related traits.

    Returns DataFrame with columns: rsID, Chromosome, Position,
    Mapped_Gene, P_Value, Beta, Trait, PMID, EFO_Trait_URI.
    """
    logger.info(f"Reading GWAS Catalog: {gwas_path}")

    # Read in chunks for memory efficiency
    usecols = [
        'SNPS', 'CHR_ID', 'CHR_POS', 'MAPPED_GENE',
        'P-VALUE', 'OR or BETA', 'MAPPED_TRAIT',
        'MAPPED_TRAIT_URI', 'PUBMEDID', 'STUDY',
        'DISEASE/TRAIT', 'STRONGEST SNP-RISK ALLELE',
    ]

    try:
        df = pd.read_csv(gwas_path, sep='\t', low_memory=False,
                         usecols=lambda c: c in usecols, dtype=str)
    except Exception as e:
        logger.error(f"Failed to read GWAS Catalog: {e}")
        return pd.DataFrame()

    logger.info(f"  Total GWAS associations: {len(df):,}")

    # --- Filter by trait keywords ---
    trait_col = 'MAPPED_TRAIT' if 'MAPPED_TRAIT' in df.columns else 'DISEASE/TRAIT'
    trait_lower = df[trait_col].fillna('').str.lower()

    # Build regex pattern for efficient matching
    pattern = '|'.join(re.escape(kw) for kw in METABOLITE_TRAIT_KEYWORDS)
    mask_keyword = trait_lower.str.contains(pattern, regex=True, na=False)

    # Also filter by EFO URI if available
    mask_efo = pd.Series(False, index=df.index)
    if 'MAPPED_TRAIT_URI' in df.columns:
        uri_str = df['MAPPED_TRAIT_URI'].fillna('')
        for efo in METABOLITE_EFO_PREFIXES:
            mask_efo |= uri_str.str.contains(efo, na=False)

    mask = mask_keyword | mask_efo
    filtered = df[mask].copy()
    logger.info(f"  Metabolite-related: {len(filtered):,} associations "
                f"({mask_keyword.sum():,} keyword + {mask_efo.sum():,} EFO)")

    # --- Clean up ---
    # Extract rsID from SNPS column (may contain "rs12345-A" or "rs12345 x rs67890")
    if 'SNPS' in filtered.columns:
        filtered['rsID'] = filtered['SNPS'].str.extract(r'(rs\d+)', expand=False)
    elif 'STRONGEST SNP-RISK ALLELE' in filtered.columns:
        filtered['rsID'] = filtered['STRONGEST SNP-RISK ALLELE'].str.extract(r'(rs\d+)', expand=False)
    else:
        filtered['rsID'] = ''

    # Drop rows without valid rsID
    filtered = filtered[filtered['rsID'].notna() & (filtered['rsID'] != '')]

    # Rename columns
    result = pd.DataFrame({
        'rsID': filtered['rsID'],
        'Chromosome': filtered.get('CHR_ID', pd.Series('', index=filtered.index)),
        'Position': filtered.get('CHR_POS', pd.Series('', index=filtered.index)),
        'Mapped_Gene': filtered.get('MAPPED_GENE', pd.Series('', index=filtered.index)),
        'P_Value': filtered.get('P-VALUE', pd.Series('', index=filtered.index)),
        'Beta': filtered.get('OR or BETA', pd.Series('', index=filtered.index)),
        'Trait': filtered[trait_col],
        'PMID': filtered.get('PUBMEDID', pd.Series('', index=filtered.index)),
        'Source': 'GWAS_Catalog',
    })

    # Clean p-value
    result['P_Value'] = pd.to_numeric(result['P_Value'], errors='coerce')

    # Filter genome-wide significance (p < 5e-8)
    sig_mask = result['P_Value'] < 5e-8
    result_sig = result[sig_mask].copy()
    logger.info(f"  Genome-wide significant (p < 5e-8): {len(result_sig):,}")

    return result_sig


# ══════════════════════════════════════════════════════════════════════
# Step 3: Parse curated mGWAS publications
# ══════════════════════════════════════════════════════════════════════

def load_curated_mgwas():
    """Load curated mGWAS associations from local supplementary tables.

    Checks for both TSV and Excel (.xlsx) variants of each file.
    Each returns standardized columns: rsID, Chromosome, Position,
    Mapped_Gene, P_Value, Beta, Trait (metabolite name), PMID, Source.
    """
    curated_dir = RAW_DIR / "mgwas"
    all_curated = []

    def _find_file(names):
        """Return first existing path from a list of candidate filenames."""
        for name in names:
            p = curated_dir / name
            if p.exists():
                return p
        return None

    # --- Shin et al. 2014 ---
    shin_path = _find_file([
        'shin2014_supplement.xlsx', 'shin2014_table_s1.tsv',
        'shin2014_table_s1.xlsx', 'shin2014_table_s1.csv',
    ])
    if shin_path:
        logger.info(f"Loading Shin 2014: {shin_path}")
        shin = _parse_shin2014(shin_path)
        all_curated.append(shin)
        logger.info(f"  Shin 2014: {len(shin):,} associations")
    else:
        logger.info("Shin 2014 not found — skipping "
                     "(download from PMID:24816252 Supplementary Table)")

    # --- Long et al. 2017 ---
    long_path = _find_file([
        'long2017_table.tsv', 'long2017_table.xlsx',
        'long2017_table.csv', 'long2017_supplement.xlsx',
    ])
    if long_path:
        logger.info(f"Loading Long 2017: {long_path}")
        long_ = _parse_long2017(long_path)
        all_curated.append(long_)
        logger.info(f"  Long 2017: {len(long_):,} associations")
    else:
        logger.info("Long 2017 not found — skipping "
                     "(requires institutional access to PMID:28263315)")

    # --- Chen et al. 2023 ---
    chen_path = _find_file([
        'chen2023_supplement.xlsx', 'chen2023_supp_data.tsv',
        'chen2023_supp_data.xlsx', 'chen2023_supp_data.csv',
    ])
    if chen_path:
        logger.info(f"Loading Chen 2023: {chen_path}")
        chen = _parse_chen2023(chen_path)
        all_curated.append(chen)
        logger.info(f"  Chen 2023: {len(chen):,} associations")
    else:
        logger.info("Chen 2023 not found — skipping "
                     "(download from PMID:36635386 Supplementary Data)")

    # --- Yin et al. 2024 ---
    yin_path = _find_file([
        'yin2024_supp_data6.xlsx', 'yin2024_supplement.xlsx',
        'yin2024_table.tsv', 'yin2024_table.csv',
    ])
    if yin_path:
        logger.info(f"Loading Yin 2024: {yin_path}")
        yin = _parse_yin2024(yin_path)
        all_curated.append(yin)
        logger.info(f"  Yin 2024: {len(yin):,} associations")
    else:
        logger.info("Yin 2024 not found — skipping "
                     "(download from PMID:38448586 Supplementary Table 5)")

    if all_curated:
        return pd.concat(all_curated, ignore_index=True)
    return pd.DataFrame()


def _parse_shin2014(path):
    """Parse Shin et al. 2014 supplementary table.

    Supports two formats:
    - Excel (.xlsx): Table S4 sheet with header row at index 2
    - TSV/CSV: flat file with columns SNP, Metabolite, P-value, etc.
    """
    path = Path(path)
    if path.suffix in ('.xlsx', '.xls'):
        # Downloaded PMC Excel: Table S4 has significant loci
        df = _read_flexible(path, sheet_name='Table S4', header=2)
    else:
        df = _read_flexible(path)
    col_map = _find_columns(df, {
        'rsID': ['snp', 'rsid', 'rs_id', 'variant', 'marker'],
        'Trait': ['biochemical', 'metabolite', 'trait', 'compound', 'analyte'],
        'Chromosome': ['chr', 'chromosome', 'chrom'],
        'Position': ['pos', 'position', 'bp', 'base_pair'],
        'P_Value': ['p', 'p_value', 'pvalue', 'p-value', 'pval'],
        'Beta': ['beta', 'effect', 'effect_size', 'b'],
        'Mapped_Gene': ['gene', 'mapped_gene', 'nearest_gene', 'locus'],
    })
    result = _standardize(df, col_map, source='Shin2014', pmid='24816252')
    return result


def _parse_long2017(path):
    """Parse Long et al. 2017 supplementary table."""
    df = _read_flexible(path)
    col_map = _find_columns(df, {
        'rsID': ['snp', 'rsid', 'rs_id', 'variant', 'marker'],
        'Trait': ['metabolite', 'trait', 'compound', 'analyte'],
        'Chromosome': ['chr', 'chromosome', 'chrom'],
        'Position': ['pos', 'position', 'bp', 'base_pair'],
        'P_Value': ['p', 'p_value', 'pvalue', 'p-value', 'pval'],
        'Beta': ['beta', 'effect', 'effect_size', 'b'],
        'Mapped_Gene': ['gene', 'mapped_gene', 'nearest_gene', 'locus'],
    })
    return _standardize(df, col_map, source='Long2017', pmid='28263315')


def _parse_chen2023(path):
    """Parse Chen et al. 2023 supplementary data."""
    df = _read_flexible(path)
    col_map = _find_columns(df, {
        'rsID': ['snp', 'rsid', 'rs_id', 'variant', 'marker', 'lead_snp'],
        'Trait': ['metabolite', 'trait', 'compound', 'analyte', 'metabolite_name'],
        'Chromosome': ['chr', 'chromosome', 'chrom'],
        'Position': ['pos', 'position', 'bp', 'base_pair'],
        'P_Value': ['p', 'p_value', 'pvalue', 'p-value', 'pval'],
        'Beta': ['beta', 'effect', 'effect_size', 'b'],
        'Mapped_Gene': ['gene', 'mapped_gene', 'nearest_gene', 'locus'],
    })
    return _standardize(df, col_map, source='Chen2023', pmid='36635386')


def _parse_yin2024(path):
    """Parse Yin et al. 2024 (Nature 628:130) supplementary table.

    TableS5 has 8,795 lead SNP–metabolic trait associations.
    Columns: Region_ID, Chromosome, Metabolite, Rsid, P-value, Effect, etc.
    NMR abbreviations (Ala, Crea, Gln, etc.) need expansion for HMDB matching.
    """
    path = Path(path)
    if path.suffix in ('.xlsx', '.xls'):
        df = _read_flexible(path, sheet_name='TableS5_all_lead_SNP_assoc', header=1)
    else:
        df = _read_flexible(path)

    # Expand NMR abbreviations to full metabolite names for HMDB matching
    _NMR_ABBREV = {
        'Ala': 'Alanine', 'Gln': 'Glutamine', 'Gly': 'Glycine', 'His': 'Histidine',
        'Ile': 'Isoleucine', 'Leu': 'Leucine', 'Val': 'Valine', 'Phe': 'Phenylalanine',
        'Tyr': 'Tyrosine', 'Crea': 'Creatinine', 'Alb': 'Albumin', 'ApoA1': 'Apolipoprotein A-I',
        'ApoB': 'Apolipoprotein B', 'GlycA': 'GlycA', 'Lac': 'Lactate', 'Glc': 'Glucose',
        'Ace': 'Acetate', 'AcAce': 'Acetoacetate', 'bOHBut': '3-Hydroxybutyric acid',
        'Cit': 'Citrate', 'Glol': 'Glycerol', 'Urea': 'Urea', 'Pyr': 'Pyruvate',
        'DHA': 'Docosahexaenoic acid', 'LA': 'Linoleic acid', 'FAw3': 'Omega-3 fatty acids',
        'FAw6': 'Omega-6 fatty acids', 'MUFA': 'Monounsaturated fatty acids',
        'PUFA': 'Polyunsaturated fatty acids', 'SFA': 'Saturated fatty acids',
        'UnSat': 'Unsaturated fatty acids', 'TotFA': 'Total fatty acids',
        'Acetone': 'Acetone', 'EtOH': 'Ethanol',
        'TotCho': 'Total cholesterol', 'EstC': 'Cholesteryl esters',
        'FreeC': 'Free cholesterol', 'remnant_C': 'Remnant cholesterol',
        'Clinical_LDL_C': 'LDL cholesterol', 'TotTG': 'Triglycerides',
    }

    col_map = _find_columns(df, {
        'rsID': ['rsid', 'snp', 'rs_id', 'variant'],
        'Trait': ['metabolite', 'trait'],
        'Chromosome': ['chromosome', 'chr', 'chrom'],
        'Position': ['lead variant position', 'pos', 'position'],
        'P_Value': ['p-value', 'p', 'pvalue'],
        'Beta': ['effect', 'beta', 'effect_size'],
        'Mapped_Gene': ['gene', 'nearest_gene', 'candidate_gene'],
    })

    result = _standardize(df, col_map, source='Yin2024', pmid='38448586')

    # Expand NMR abbreviations in Trait column
    if not result.empty:
        result['Trait'] = result['Trait'].map(
            lambda x: _NMR_ABBREV.get(x, x) if pd.notna(x) else x
        )

    return result


def _read_flexible(path, sheet_name=0, header=0):
    """Read a TSV/CSV/Excel file with flexible delimiter detection.

    For Excel files, sheet_name selects the sheet (name or index).
    header selects the header row (0-indexed).
    """
    path = Path(path)
    if path.suffix in ('.xlsx', '.xls'):
        return pd.read_excel(path, sheet_name=sheet_name, header=header, dtype=str)
    # Try TSV first, then CSV
    try:
        df = pd.read_csv(path, sep='\t', dtype=str)
        if len(df.columns) > 2:
            return df
    except Exception:
        pass
    return pd.read_csv(path, dtype=str)


def _find_columns(df, target_map):
    """Fuzzy-match DataFrame columns to target column names.

    target_map: {target_name: [possible_column_names_lowercase]}
    Returns: {target_name: actual_column_name_or_None}
    """
    cols_lower = {c.lower().strip(): c for c in df.columns}
    result = {}
    for target, candidates in target_map.items():
        found = None
        for cand in candidates:
            if cand in cols_lower:
                found = cols_lower[cand]
                break
            # Partial match
            for cl, orig in cols_lower.items():
                if cand in cl:
                    found = orig
                    break
            if found:
                break
        result[target] = found
    return result


def _standardize(df, col_map, source, pmid):
    """Standardize a curated mGWAS table to common schema."""
    rows = []
    for _, row in df.iterrows():
        rsid = str(row.get(col_map.get('rsID', ''), '')).strip() if col_map.get('rsID') else ''
        if not rsid or not rsid.startswith('rs'):
            # Try to extract rsID
            for val in row.values:
                m = re.search(r'(rs\d+)', str(val))
                if m:
                    rsid = m.group(1)
                    break
        if not rsid.startswith('rs'):
            continue

        trait = str(row.get(col_map.get('Trait', ''), '')).strip() if col_map.get('Trait') else ''
        chrom = str(row.get(col_map.get('Chromosome', ''), '')).strip() if col_map.get('Chromosome') else ''
        pos = str(row.get(col_map.get('Position', ''), '')).strip() if col_map.get('Position') else ''
        pval = str(row.get(col_map.get('P_Value', ''), '')).strip() if col_map.get('P_Value') else ''
        beta = str(row.get(col_map.get('Beta', ''), '')).strip() if col_map.get('Beta') else ''
        gene = str(row.get(col_map.get('Mapped_Gene', ''), '')).strip() if col_map.get('Mapped_Gene') else ''

        rows.append({
            'rsID': rsid,
            'Chromosome': chrom,
            'Position': pos,
            'Mapped_Gene': gene,
            'P_Value': pval,
            'Beta': beta,
            'Trait': trait,
            'PMID': pmid,
            'Source': source,
        })

    result = pd.DataFrame(rows)
    if not result.empty:
        result['P_Value'] = pd.to_numeric(result['P_Value'], errors='coerce')
    return result


# ══════════════════════════════════════════════════════════════════════
# Step 4: Map metabolite trait names → HMDB IDs
# ══════════════════════════════════════════════════════════════════════

def build_name_to_hmdb(synonym_df):
    """Build a lowercase name → HMDB_ID lookup dictionary.

    Priority: primary > synonym > iupac > traditional > cas > formula.
    For duplicates, prefer the HMDB ID with more synonyms (more data).
    """
    # Priority ordering
    priority = {'primary': 0, 'synonym': 1, 'iupac': 2, 'traditional': 3, 'cas': 4, 'formula': 5}

    name_to_hmdb = {}
    name_to_priority = {}

    for _, row in synonym_df.iterrows():
        name = str(row['name']).strip().lower()
        hmdb = str(row['hmdb_id']).strip()
        stype = str(row.get('synonym_type', 'synonym'))
        prio = priority.get(stype, 5)

        if name not in name_to_hmdb or prio < name_to_priority.get(name, 99):
            name_to_hmdb[name] = hmdb
            name_to_priority[name] = prio

    logger.info(f"  Name→HMDB lookup: {len(name_to_hmdb):,} unique names")
    return name_to_hmdb


def _clean_trait_name(trait):
    """Clean GWAS trait string to extract metabolite name.

    Examples:
      "Valine levels" → "valine"
      "L-Carnitine measurement" → "l-carnitine"
      "Serum metabolite levels (X-12345)" → "x-12345"
      "Glycerophosphocholine (16:0) levels" → "glycerophosphocholine"
    """
    if not isinstance(trait, str):
        return ''
    t = trait.lower().strip()

    # Remove common suffixes
    for suffix in [' levels', ' level', ' measurement', ' measurements',
                   ' concentration', ' concentrations',
                   ' in blood', ' in serum', ' in plasma', ' in urine',
                   ' (metabolomics)', ' metabolite', ' amount',
                   ' ratio', ' percentage', ' to total fatty acids',
                   ' to total lipids']:
        t = t.replace(suffix, '')

    # Remove parenthetical lipid chain annotations like (16:0), (18:1/18:2)
    # but keep the base compound name
    t = re.sub(r'\s*\(\d+:\d+[^)]*\)\s*', '', t)

    # Remove leading "serum ", "plasma ", "urinary "
    for prefix in ['serum ', 'plasma ', 'urinary ', 'blood ', 'circulating ']:
        if t.startswith(prefix):
            t = t[len(prefix):]

    return t.strip()


def _convert_lipid_species(cleaned):
    """Convert GWAS lipid species names to HMDB-style shorthand.

    Examples:
      "level of phosphatidylcholine (16:0_18:1)" → "pc(16:0/18:1)"
      "phosphatidylcholine diacyl c36:2"         → "pc(36:2)"
      "phosphatidylcholine acyl-alkyl c34:2"     → "pc(o-34:2)"
      "sphingomyelin c18:1"                      → "sm(d18:1)"
      "lysophosphatidylcholine a c18:0"           → "lpc(18:0)"
    """
    _LIPID_CLASS_MAP = {
        'phosphatidylcholine': 'pc',
        'lysophosphatidylcholine': 'lpc',
        'phosphatidylethanolamine': 'pe',
        'sphingomyelin': 'sm',
        'diacylglycerol': 'dg',
        'triacylglycerol': 'tg',
        'ceramide': 'cer',
    }

    t = cleaned.lower().strip()

    # Remove "level of " prefix and " in blood serum/plasma" suffix
    t = re.sub(r'^level\s+of\s+', '', t)
    t = re.sub(r'\s+in\s+(blood\s+)?(serum|plasma)$', '', t)

    # Pattern 1: "phosphatidylcholine (16:0_18:1)" → "pc(16:0/18:1)"
    m = re.match(r'(\w+)\s*\((\d+:\d+)[_/](\d+:\d+)\)', t)
    if m:
        cls_name, chain1, chain2 = m.groups()
        abbrev = _LIPID_CLASS_MAP.get(cls_name)
        if abbrev:
            return f"{abbrev}({chain1}/{chain2})"

    # Pattern 2: "phosphatidylcholine (O-16:0_18:1)" → "pc(o-16:0/18:1)"
    m = re.match(r'(\w+)\s*\(o-(\d+:\d+)[_/](\d+:\d+)\)', t)
    if m:
        cls_name, chain1, chain2 = m.groups()
        abbrev = _LIPID_CLASS_MAP.get(cls_name)
        if abbrev:
            return f"{abbrev}(o-{chain1}/{chain2})"

    # Pattern 3: "phosphatidylcholine diacyl c36:2" → "pc(36:2)"
    m = re.match(r'(\w+)\s+diacyl\s+c?(\d+:\d+)', t)
    if m:
        cls_name, chain = m.groups()
        abbrev = _LIPID_CLASS_MAP.get(cls_name)
        if abbrev:
            return f"{abbrev}({chain})"

    # Pattern 4: "phosphatidylcholine acyl-alkyl c34:2" → "pc(o-34:2)"
    m = re.match(r'(\w+)\s+acyl-alkyl\s+c?(\d+:\d+)', t)
    if m:
        cls_name, chain = m.groups()
        abbrev = _LIPID_CLASS_MAP.get(cls_name)
        if abbrev:
            return f"{abbrev}(o-{chain})"

    # Pattern 5: "lysophosphatidylcholine a c18:0" → "lpc(18:0)"
    m = re.match(r'(\w+)\s+a\s+c?(\d+:\d+)', t)
    if m:
        cls_name, chain = m.groups()
        abbrev = _LIPID_CLASS_MAP.get(cls_name)
        if abbrev:
            return f"{abbrev}({chain})"

    # Pattern 6: "sphingomyelin c18:1" → "sm(d18:1)"
    m = re.match(r'sphingomyelin\s+c?(\d+:\d+)', t)
    if m:
        return f"sm(d{m.group(1)})"

    return None


def map_traits_to_hmdb(df, name_to_hmdb):
    """Map Trait column to HMDB_ID using fuzzy matching.

    Strategy:
    1. Exact match of cleaned trait name
    2. Exact match of original trait (lowercased)
    3. Try removing common prefixes (L-, D-, (+)-, (-)-, alpha-, beta-)
    4. Mark unmapped traits for manual curation

    Returns: DataFrame with new HMDB_ID and Metabolite_Name columns.
    """
    mapped = []
    unmapped_traits = defaultdict(int)

    for _, row in df.iterrows():
        trait = str(row.get('Trait', ''))
        cleaned = _clean_trait_name(trait)

        hmdb_id = None
        met_name = ''

        # Strategy 1: exact match on cleaned name
        if cleaned and cleaned in name_to_hmdb:
            hmdb_id = name_to_hmdb[cleaned]
            met_name = cleaned

        # Strategy 2: exact match on original trait (lowercased)
        if not hmdb_id:
            t_lower = trait.lower().strip()
            if t_lower in name_to_hmdb:
                hmdb_id = name_to_hmdb[t_lower]
                met_name = t_lower

        # Strategy 3: try removing stereochemistry prefixes
        if not hmdb_id and cleaned:
            for prefix in ['l-', 'd-', '(+)-', '(-)-', 'alpha-', 'beta-',
                           'gamma-', 'delta-', 'n-', 'o-', 's-']:
                if cleaned.startswith(prefix):
                    stripped = cleaned[len(prefix):]
                    if stripped in name_to_hmdb:
                        hmdb_id = name_to_hmdb[stripped]
                        met_name = stripped
                        break

        # Strategy 4: try each word (for compound traits like "valine/leucine")
        # Blocklist: generic terms that match HMDB entries but are too ambiguous
        _GENERIC_BLOCKLIST = {
            'acid', 'base', 'salt', 'ester', 'amine', 'amide', 'sugar',
            'lipid', 'fatty', 'amino', 'ratio', 'total', 'free', 'level',
            'amount', 'percentage', 'rate', 'index', 'score', 'body',
            'mass', 'weight', 'size', 'count', 'number', 'type',
            'protein', 'enzyme', 'receptor', 'factor', 'gene',
            'omega-3', 'omega-6', 'polyunsaturated', 'monounsaturated',
            'saturated', 'unsaturated',
        }
        if not hmdb_id and cleaned:
            parts = re.split(r'[/,;\s]+', cleaned)
            for part in parts:
                part = part.strip()
                if len(part) >= 4 and part not in _GENERIC_BLOCKLIST and part in name_to_hmdb:
                    hmdb_id = name_to_hmdb[part]
                    met_name = part
                    break

        # Strategy 5: lipid species name conversion
        # "level of phosphatidylcholine (16:0_18:1) in blood serum" → "pc(16:0/18:1)"
        # "phosphatidylcholine diacyl c36:2" → "pc(36:2)"
        if not hmdb_id and cleaned:
            lipid_name = _convert_lipid_species(cleaned)
            if lipid_name and lipid_name in name_to_hmdb:
                hmdb_id = name_to_hmdb[lipid_name]
                met_name = lipid_name

        if hmdb_id:
            new_row = row.to_dict()
            new_row['HMDB_ID'] = hmdb_id
            new_row['Metabolite_Name'] = met_name
            mapped.append(new_row)
        else:
            unmapped_traits[trait] += 1

    result = pd.DataFrame(mapped)

    # Report unmapped
    n_unmapped = sum(unmapped_traits.values())
    total = len(df)
    pct = len(result) / total * 100 if total > 0 else 0
    logger.info(f"  Trait→HMDB mapping: {len(result):,} / {total:,} mapped ({pct:.1f}%)")

    if unmapped_traits:
        # Save unmapped for manual curation
        unmapped_path = CACHE_DIR / "unmapped_traits.tsv"
        unmapped_df = pd.DataFrame([
            {'Trait': t, 'Count': c} for t, c in
            sorted(unmapped_traits.items(), key=lambda x: -x[1])
        ])
        unmapped_df.to_csv(unmapped_path, sep='\t', index=False)
        logger.info(f"  Unmapped traits ({n_unmapped:,} associations, "
                    f"{len(unmapped_traits):,} unique) saved to {unmapped_path}")

    return result


# ══════════════════════════════════════════════════════════════════════
# Step 5: Deduplicate & finalize
# ══════════════════════════════════════════════════════════════════════

def deduplicate_and_finalize(df, xref_df):
    """Deduplicate on (HMDB_ID, rsID), keep best p-value, merge sources.

    Attach SMILES from xref mapping.
    """
    if df.empty:
        return df

    # Ensure P_Value is numeric
    df['P_Value'] = pd.to_numeric(df['P_Value'], errors='coerce')

    # Group by (HMDB_ID, rsID) — keep row with smallest p-value,
    # but merge Source and PMID from all rows
    groups = df.groupby(['HMDB_ID', 'rsID'], sort=False)

    deduped_rows = []
    for (hmdb, rsid), group in groups:
        # Sort by p-value (smallest first)
        group_sorted = group.sort_values('P_Value', ascending=True)
        best = group_sorted.iloc[0].to_dict()

        # Merge sources
        sources = sorted(group['Source'].dropna().unique())
        best['Source'] = '; '.join(sources)

        # Merge PMIDs
        pmids = sorted(group['PMID'].dropna().astype(str).unique())
        best['PMID'] = '; '.join(pmids)

        deduped_rows.append(best)

    result = pd.DataFrame(deduped_rows)

    # Attach SMILES from xref mapping
    if xref_df is not None and not xref_df.empty:
        smiles_map = dict(zip(xref_df['hmdb_id'], xref_df['smiles']))
        name_map = dict(zip(xref_df['hmdb_id'], xref_df['metabolite_name']))
        result['SMILES'] = result['HMDB_ID'].map(smiles_map).fillna('')
        # Use canonical metabolite name from xref if available
        result['Metabolite_Name'] = result['HMDB_ID'].map(name_map).fillna(result['Metabolite_Name'])

    # Clean up column order
    col_order = [
        'HMDB_ID', 'Metabolite_Name', 'SMILES',
        'rsID', 'Chromosome', 'Position', 'Mapped_Gene',
        'P_Value', 'Beta', 'Trait', 'Source', 'PMID',
    ]
    for c in col_order:
        if c not in result.columns:
            result[c] = ''
    result = result[col_order]

    # Sort by HMDB_ID, then p-value
    result = result.sort_values(['HMDB_ID', 'P_Value'], ascending=[True, True])
    result = result.reset_index(drop=True)

    logger.info(f"  Deduplicated: {len(result):,} unique (HMDB_ID, rsID) pairs")
    logger.info(f"  Metabolites: {result['HMDB_ID'].nunique():,}, SNPs: {result['rsID'].nunique():,}")

    return result


# ══════════════════════════════════════════════════════════════════════
# Main pipeline
# ══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Build mGWAS database')
    parser.add_argument('--skip-download', action='store_true',
                        help='Use cached GWAS Catalog file')
    parser.add_argument('--synonyms-only', action='store_true',
                        help='Only build the synonym table')
    parser.add_argument('--force-synonyms', action='store_true',
                        help='Rebuild synonym table even if it exists')
    args = parser.parse_args()

    logger.info("=" * 70)
    logger.info("Building mGWAS (Metabolite–SNP) Database")
    logger.info("=" * 70)

    # --- Step 1: Build synonym table ---
    logger.info("\n--- Step 1: Building metabolite synonym table ---")
    synonym_df = build_synonym_table(force=args.force_synonyms)

    if args.synonyms_only:
        logger.info("Synonyms-only mode — done.")
        return

    # Build lookup
    name_to_hmdb = build_name_to_hmdb(synonym_df)

    # Load xref mapping for SMILES
    xref_df = pd.read_csv(XREF_CSV, dtype=str) if XREF_CSV.exists() else pd.DataFrame()

    all_dfs = []

    # --- Step 2: GWAS Catalog ---
    logger.info("\n--- Step 2: GWAS Catalog ---")
    if args.skip_download and not GWAS_CACHE.exists():
        logger.warning("--skip-download but no cache found; will download")
    gwas_path = download_gwas_catalog(force=not args.skip_download and not GWAS_CACHE.exists())
    gwas_df = filter_gwas_catalog(gwas_path)
    if not gwas_df.empty:
        logger.info("  Mapping GWAS Catalog traits to HMDB IDs...")
        gwas_mapped = map_traits_to_hmdb(gwas_df, name_to_hmdb)
        if not gwas_mapped.empty:
            all_dfs.append(gwas_mapped)
            logger.info(f"  GWAS Catalog contributed: {len(gwas_mapped):,} mapped associations")

    # --- Step 3: Curated mGWAS publications ---
    logger.info("\n--- Step 3: Curated mGWAS publications ---")
    curated_df = load_curated_mgwas()
    if not curated_df.empty:
        logger.info("  Mapping curated traits to HMDB IDs...")
        curated_mapped = map_traits_to_hmdb(curated_df, name_to_hmdb)
        if not curated_mapped.empty:
            all_dfs.append(curated_mapped)
            logger.info(f"  Curated sources contributed: {len(curated_mapped):,} mapped associations")

    # --- Step 4: Merge & deduplicate ---
    logger.info("\n--- Step 4: Merge & deduplicate ---")
    if not all_dfs:
        logger.error("No data collected from any source — check data availability")
        sys.exit(1)

    merged = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"  Total before dedup: {len(merged):,}")

    final = deduplicate_and_finalize(merged, xref_df)

    # --- Step 5: Save ---
    logger.info("\n--- Step 5: Saving output ---")
    DB_DIR.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUTPUT_CSV, index=False)
    logger.info(f"  Saved: {OUTPUT_CSV}")

    # --- Summary ---
    logger.info("\n" + "=" * 70)
    logger.info("mGWAS Database Summary")
    logger.info("=" * 70)
    logger.info(f"  Total associations:  {len(final):,}")
    logger.info(f"  Unique metabolites:  {final['HMDB_ID'].nunique():,}")
    logger.info(f"  Unique SNPs:         {final['rsID'].nunique():,}")
    logger.info(f"  Unique genes:        {final['Mapped_Gene'].nunique():,}")
    if 'Source' in final.columns:
        for src in sorted(final['Source'].unique()):
            n = (final['Source'] == src).sum()
            # Handle merged sources like "GWAS_Catalog; Shin2014"
            logger.info(f"    {src}: {n:,}")
    logger.info(f"  With SMILES:         {(final['SMILES'] != '').sum():,} / {len(final):,}")
    logger.info(f"  Chromosomes:         {sorted(final['Chromosome'].dropna().unique())}")
    logger.info(f"  Output:              {OUTPUT_CSV}")


def _count_lines(path):
    """Count lines in a file."""
    with open(path) as f:
        return sum(1 for _ in f)


if __name__ == '__main__':
    main()
