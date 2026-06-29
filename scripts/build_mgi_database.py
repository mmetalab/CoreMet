#!/usr/bin/env python3
"""
Build MGI (Metabolite–Gene Interaction) database.

Data sources:
  1. CTD Chemical–Gene Interactions (583 MB TSV, ~3M rows)
     — Stream-filtered to HMDB metabolites via CAS + name matching.
  2. KEGG Enzyme–Gene + Enzyme–Compound links
     — Compound→HMDB mapping via xref table, gene→symbol via KEGG cache.

Output: data/databases/mgi_database.csv
Schema: HMDB_ID, Metabolite_Name, SMILES, Gene_ID, Gene_Symbol,
        Organism, Interaction_Type, Interaction_Actions, Source, PMID
"""

import csv
import json
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
from lxml import etree

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

BASE = Path(__file__).resolve().parent.parent
CTD_FILE = BASE / "data" / "raw" / "ctd" / "CTD_chem_gene_ixns.tsv"
HMDB_XML = BASE / "data" / "raw" / "hmdb" / "hmdb_metabolites.xml"
SYNONYM_CSV = BASE / "data" / "mappings" / "hmdb_synonyms.csv"
XREF_CSV = BASE / "data" / "mappings" / "hmdb_xref_mapping.csv"
KEGG_CACHE = BASE / "data" / "cache" / "kegg"
OUTPUT = BASE / "data" / "databases" / "mgi_database.csv"

HMDB_NS = "{http://www.hmdb.ca}"

# Organisms we care about (CTD has many)
ORGANISM_MAP = {
    "9606": "Homo sapiens",
    "10090": "Mus musculus",
    "10116": "Rattus norvegicus",
    "9913": "Bos taurus",
    "7955": "Danio rerio",
    "7227": "Drosophila melanogaster",
    "6239": "Caenorhabditis elegans",
    "4932": "Saccharomyces cerevisiae",
    "3702": "Arabidopsis thaliana",
    "511145": "Escherichia coli",
    "287": "Pseudomonas aeruginosa",
}

# KEGG organism codes → scientific names
KEGG_ORG_MAP = {
    "hsa": "Homo sapiens",
    "mmu": "Mus musculus",
    "rno": "Rattus norvegicus",
    "bta": "Bos taurus",
    "dme": "Drosophila melanogaster",
    "cel": "Caenorhabditis elegans",
    "sce": "Saccharomyces cerevisiae",
    "ath": "Arabidopsis thaliana",
    "eco": "Escherichia coli",
    "pae": "Pseudomonas aeruginosa",
}

# Simplify CTD interaction actions
ACTION_MAP = {
    "increases^expression": "increases expression",
    "decreases^expression": "decreases expression",
    "increases^activity": "increases activity",
    "decreases^activity": "decreases activity",
    "affects^expression": "affects expression",
    "affects^activity": "affects activity",
    "increases^phosphorylation": "increases phosphorylation",
    "decreases^phosphorylation": "decreases phosphorylation",
    "affects^binding": "affects binding",
    "increases^secretion": "increases secretion",
    "decreases^secretion": "decreases secretion",
    "increases^metabolic processing": "increases metabolism",
    "decreases^metabolic processing": "decreases metabolism",
    "affects^metabolic processing": "affects metabolism",
    "increases^transport": "increases transport",
    "decreases^transport": "decreases transport",
    "affects^reaction": "affects reaction",
    "increases^reaction": "affects reaction",
    "decreases^reaction": "affects reaction",
    "affects^folding": "affects folding",
    "affects^localization": "affects localization",
    "increases^abundance": "increases abundance",
    "decreases^abundance": "decreases abundance",
    "affects^abundance": "affects abundance",
    "increases^oxidation": "increases oxidation",
    "decreases^oxidation": "decreases oxidation",
    "increases^reduction": "increases reduction",
    "affects^splicing": "affects splicing",
    "increases^uptake": "increases uptake",
    "decreases^uptake": "decreases uptake",
    "affects^mutagenesis": "affects mutagenesis",
    "increases^stability": "increases stability",
    "decreases^stability": "decreases stability",
    "increases^methylation": "increases methylation",
    "decreases^methylation": "decreases methylation",
}


# ── Step 1: Build metabolite lookup tables ────────────────────────────

def build_cas_lookup() -> dict:
    """Parse HMDB XML → CAS → (HMDB_ID, name, SMILES)."""
    if not HMDB_XML.exists():
        log.warning("HMDB XML not found — CAS matching disabled")
        return {}
    log.info("Parsing HMDB XML for CAS → metabolite mapping …")
    t0 = time.time()
    cas_map = {}
    n = 0
    context = etree.iterparse(str(HMDB_XML), events=("end",), tag=f"{HMDB_NS}metabolite")
    for _, mel in context:
        parent = mel.getparent()
        if parent is None or parent.tag != f"{HMDB_NS}hmdb":
            mel.clear()
            continue
        n += 1
        hmdb_id = mel.findtext(f"{HMDB_NS}accession", "")
        name = mel.findtext(f"{HMDB_NS}name", "")
        cas = mel.findtext(f"{HMDB_NS}cas_registry_number", "") or ""
        smiles = mel.findtext(f"{HMDB_NS}smiles", "") or ""
        if cas and cas != "0":
            cas_map[cas] = (hmdb_id, name, smiles)
        mel.clear()
        if n % 50000 == 0:
            log.info(f"  … {n:,} metabolites")
    log.info(f"  Done: {n:,} metabolites, {len(cas_map):,} with CAS ({time.time()-t0:.1f}s)")
    return cas_map


def build_name_lookup() -> dict:
    """Build lowered-name → (HMDB_ID, metabolite_name) from synonyms table."""
    if not SYNONYM_CSV.exists():
        log.warning("HMDB synonyms CSV not found — name matching disabled")
        return {}
    log.info("Loading HMDB synonym table for name matching …")
    t0 = time.time()
    name_map = {}  # lowercase name → (hmdb_id, canonical_name)

    # First load xref for SMILES
    smiles_map = {}
    if XREF_CSV.exists():
        xref = pd.read_csv(XREF_CSV, dtype=str).fillna("")
        for _, row in xref.iterrows():
            if row["smiles"]:
                smiles_map[row["hmdb_id"]] = row["smiles"]

    # Load primary names first (highest priority)
    syn_df = pd.read_csv(SYNONYM_CSV, dtype=str).fillna("")
    primary = syn_df[syn_df["synonym_type"] == "primary"]
    canonical = {}  # hmdb_id → canonical name
    for _, row in primary.iterrows():
        canonical[row["hmdb_id"]] = row["name"]
        name_map[row["name"].lower()] = (row["hmdb_id"], row["name"])

    # Then add all synonyms (don't overwrite primary)
    for _, row in syn_df.iterrows():
        key = row["name"].lower().strip()
        if key and key not in name_map and len(key) > 2:
            canon = canonical.get(row["hmdb_id"], row["name"])
            name_map[key] = (row["hmdb_id"], canon)

    log.info(f"  Name lookup: {len(name_map):,} entries ({time.time()-t0:.1f}s)")
    return name_map


def load_smiles_map() -> dict:
    """HMDB_ID → SMILES from xref CSV + HMDB XML fallback."""
    smiles = {}
    if XREF_CSV.exists():
        xref = pd.read_csv(XREF_CSV, dtype=str).fillna("")
        for _, row in xref.iterrows():
            if row.get("smiles"):
                smiles[row["hmdb_id"]] = row["smiles"]
    return smiles


# ── Step 2: Stream CTD Chemical–Gene Interactions ─────────────────────

def stream_ctd_genes(cas_map: dict, name_map: dict) -> list:
    """
    Stream CTD_chem_gene_ixns.tsv and match chemicals to HMDB.

    Matching strategy:
      1. CAS number match (most reliable)
      2. Chemical name exact match against HMDB synonyms (fallback)

    Returns list of dicts.
    """
    if not CTD_FILE.exists():
        log.error(f"CTD file not found: {CTD_FILE}")
        return []

    log.info("Streaming CTD Chemical–Gene Interactions …")
    t0 = time.time()
    rows = []
    matched_cas = 0
    matched_name = 0
    skipped_no_match = 0
    skipped_no_organism = 0
    total = 0

    with open(CTD_FILE, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue

            chem_name = parts[0]
            chem_id = parts[1]     # MeSH ID (C-prefixed)
            cas_rn = parts[2]      # CAS Registry Number (may be empty)
            gene_symbol = parts[3]
            gene_id = parts[4]     # Entrez Gene ID
            gene_forms = parts[5]  # protein, mRNA, gene, etc.
            organism = parts[6]
            org_id = parts[7]
            interaction = parts[8]
            actions = parts[9]     # pipe-separated action codes
            pmids = parts[10] if len(parts) > 10 else ""

            # Skip rows without organism (in vitro / unspecified)
            if not org_id:
                skipped_no_organism += 1
                continue

            # Only keep organisms we care about
            if org_id not in ORGANISM_MAP:
                skipped_no_organism += 1
                continue

            # Match to HMDB
            hmdb_id = None
            met_name = None
            match_method = None

            # Strategy 1: CAS match
            if cas_rn and cas_rn in cas_map:
                hmdb_id, met_name, _ = cas_map[cas_rn]
                match_method = "CAS"
                matched_cas += 1
            # Strategy 2: Name match
            elif chem_name.lower() in name_map:
                hmdb_id, met_name = name_map[chem_name.lower()]
                match_method = "name"
                matched_name += 1
            else:
                skipped_no_match += 1
                continue

            # Simplify interaction actions
            action_parts = actions.split("|") if actions else []
            simplified = []
            for a in action_parts:
                simplified.append(ACTION_MAP.get(a, a.replace("^", " ")))
            # Deduplicate and pick the primary action
            simplified = list(dict.fromkeys(simplified))  # preserve order, deduplicate
            interaction_type = simplified[0] if simplified else "interacts"
            interaction_actions = "; ".join(simplified)

            rows.append({
                "HMDB_ID": hmdb_id,
                "Metabolite_Name": met_name,
                "Gene_ID": gene_id,
                "Gene_Symbol": gene_symbol,
                "Organism": ORGANISM_MAP[org_id],
                "Interaction_Type": interaction_type,
                "Interaction_Actions": interaction_actions,
                "Source": "CTD",
                "PMID": pmids.replace("|", "; "),
            })

            if total % 500000 == 0:
                log.info(f"  … {total:,} rows processed, {len(rows):,} matched")

    log.info(
        f"  CTD done: {total:,} rows → {len(rows):,} matched "
        f"(CAS: {matched_cas:,}, name: {matched_name:,}, "
        f"no-match: {skipped_no_match:,}, no-org: {skipped_no_organism:,}) "
        f"({time.time()-t0:.1f}s)"
    )
    return rows


# ── Step 3: KEGG Enzyme–Compound–Gene links ──────────────────────────

def load_kegg_mgi() -> list:
    """
    Build metabolite–gene links from KEGG cache:
      enzyme_compound_links.json  →  enzyme EC → [KEGG compound IDs]
      enzyme_genes_<org>.json     →  enzyme EC → [org:geneID]
      compound_info.json          →  KEGG compound → {name, smiles, hmdb_id}

    Link path: compound → enzyme → gene
    """
    ec_compound_path = KEGG_CACHE / "enzyme_compound_links.json"
    compound_info_path = KEGG_CACHE / "compound_info.json"

    if not ec_compound_path.exists() or not compound_info_path.exists():
        log.warning("KEGG cache files not found — KEGG MGI disabled")
        return []

    log.info("Loading KEGG enzyme–compound–gene links …")
    t0 = time.time()

    with open(ec_compound_path) as f:
        ec_compounds = json.load(f)  # EC → [compound IDs]
    with open(compound_info_path) as f:
        compound_info = json.load(f)  # compound → {name, smiles, hmdb_id}

    # Load xref to map KEGG compounds → HMDB
    kegg_to_hmdb = {}  # kegg_id → hmdb_id
    if XREF_CSV.exists():
        xref = pd.read_csv(XREF_CSV, dtype=str).fillna("")
        # xref has hmdb_id, kegg_id columns
        if "kegg_id" in xref.columns:
            for _, row in xref.iterrows():
                if row.get("kegg_id"):
                    kegg_to_hmdb[row["kegg_id"]] = row["hmdb_id"]

    # Also use compound_info hmdb_id field
    for cpd_id, info in compound_info.items():
        if info.get("hmdb_id") and cpd_id not in kegg_to_hmdb:
            kegg_to_hmdb[cpd_id] = info["hmdb_id"]

    # Build compound → EC mapping (reverse of ec_compounds)
    compound_to_ec = defaultdict(set)
    for ec, cpds in ec_compounds.items():
        for cpd in cpds:
            compound_to_ec[cpd].add(ec)

    # Load all organism enzyme→gene files
    ec_genes = defaultdict(list)  # EC → [(gene_id, gene_symbol_prefix, organism)]
    for org_code, org_name in KEGG_ORG_MAP.items():
        gene_file = KEGG_CACHE / f"enzyme_genes_{org_code}.json"
        if not gene_file.exists():
            continue
        with open(gene_file) as f:
            org_ec_genes = json.load(f)  # EC → ["hsa:12345", ...]
        for ec, genes in org_ec_genes.items():
            for g in genes:
                # g is like "hsa:12345"
                parts = g.split(":", 1)
                if len(parts) == 2:
                    ec_genes[ec].append((parts[1], org_code, org_name))

    # Now build MGI rows: for each compound with HMDB mapping,
    # find enzymes, then find genes
    rows = []
    xref_names = {}
    if XREF_CSV.exists():
        xref = pd.read_csv(XREF_CSV, dtype=str).fillna("")
        for _, row in xref.iterrows():
            xref_names[row["hmdb_id"]] = row.get("metabolite_name", "")

    for cpd_id, ecs in compound_to_ec.items():
        hmdb_id = kegg_to_hmdb.get(cpd_id)
        if not hmdb_id:
            continue
        cpd_info = compound_info.get(cpd_id, {})
        met_name = xref_names.get(hmdb_id, cpd_info.get("name", ""))

        for ec in ecs:
            gene_list = ec_genes.get(ec, [])
            for gene_id_str, org_code, org_name in gene_list:
                rows.append({
                    "HMDB_ID": hmdb_id,
                    "Metabolite_Name": met_name,
                    "Gene_ID": gene_id_str,
                    "Gene_Symbol": "",  # Will be filled later if possible
                    "Organism": org_name,
                    "Interaction_Type": "enzymatic reaction",
                    "Interaction_Actions": f"enzyme {ec}",
                    "Source": "KEGG",
                    "PMID": "",
                })

    log.info(
        f"  KEGG MGI: {len(rows):,} links from "
        f"{len(kegg_to_hmdb):,} compound→HMDB mappings "
        f"({time.time()-t0:.1f}s)"
    )
    return rows


# ── Step 4: Resolve KEGG gene IDs to symbols ─────────────────────────

def resolve_kegg_gene_symbols(rows: list) -> list:
    """
    For KEGG-sourced rows, try to resolve gene IDs to symbols
    using the gene_uniprot cache files (which map org:geneID → uniprot).
    We use the gene ID as the symbol if no better mapping exists.
    """
    # KEGG gene lists have format "hsa:12345" but we store just "12345"
    # For gene symbols, we'll use the Entrez gene ID as-is since KEGG
    # doesn't directly provide gene symbols in its cache
    for row in rows:
        if row["Source"] == "KEGG" and not row["Gene_Symbol"]:
            row["Gene_Symbol"] = row["Gene_ID"]
    return rows


# ── Step 5: Deduplicate and finalize ──────────────────────────────────

def deduplicate_and_finalize(ctd_rows: list, kegg_rows: list, smiles_map: dict) -> pd.DataFrame:
    """
    Merge CTD + KEGG rows, deduplicate by (HMDB_ID, Gene_Symbol, Organism),
    keeping unique interaction types and merging PMIDs.
    """
    log.info("Deduplicating and finalizing …")
    all_rows = ctd_rows + kegg_rows
    if not all_rows:
        log.warning("No rows to finalize!")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)
    log.info(f"  Raw rows: {len(df):,}")

    # Group by (HMDB_ID, Gene_Symbol, Organism, Interaction_Type) and merge
    # Keep unique PMIDs and sources
    grouped = df.groupby(
        ["HMDB_ID", "Gene_Symbol", "Organism", "Interaction_Type"],
        as_index=False
    ).agg({
        "Metabolite_Name": "first",
        "Gene_ID": "first",
        "Interaction_Actions": lambda x: "; ".join(sorted(set(
            a.strip() for val in x for a in str(val).split(";") if a.strip()
        ))),
        "Source": lambda x: "; ".join(sorted(set(x))),
        "PMID": lambda x: "; ".join(sorted(set(
            p.strip() for val in x for p in str(val).split(";") if p.strip()
        ))[:5]),  # Cap at 5 PMIDs
    })

    # Add SMILES
    grouped["SMILES"] = grouped["HMDB_ID"].map(smiles_map).fillna("")

    # Reorder columns
    col_order = [
        "HMDB_ID", "Metabolite_Name", "SMILES",
        "Gene_ID", "Gene_Symbol", "Organism",
        "Interaction_Type", "Interaction_Actions",
        "Source", "PMID",
    ]
    for c in col_order:
        if c not in grouped.columns:
            grouped[c] = ""
    grouped = grouped[col_order].sort_values(
        ["HMDB_ID", "Gene_Symbol", "Organism"]
    ).reset_index(drop=True)

    log.info(f"  Final: {len(grouped):,} unique interactions")
    log.info(f"  Metabolites: {grouped['HMDB_ID'].nunique():,}")
    log.info(f"  Genes: {grouped['Gene_Symbol'].nunique():,}")
    log.info(f"  Organisms: {grouped['Organism'].nunique():,}")
    log.info(f"  Interaction types: {grouped['Interaction_Type'].nunique():,}")
    sources = grouped["Source"].value_counts()
    for src, cnt in sources.items():
        log.info(f"    {src}: {cnt:,}")

    return grouped


# ── Main ──────────────────────────────────────────────────────────────

def main():
    log.info("=" * 60)
    log.info("Building MGI (Metabolite–Gene Interaction) database")
    log.info("=" * 60)
    t_start = time.time()

    # Step 1: Build lookup tables
    cas_map = build_cas_lookup()
    name_map = build_name_lookup()
    smiles_map = load_smiles_map()

    # Step 2: Stream CTD
    ctd_rows = stream_ctd_genes(cas_map, name_map)

    # Step 3: KEGG links
    kegg_rows = load_kegg_mgi()
    kegg_rows = resolve_kegg_gene_symbols(kegg_rows)

    # Step 4: Deduplicate and finalize
    df = deduplicate_and_finalize(ctd_rows, kegg_rows, smiles_map)

    if df.empty:
        log.error("No data produced — aborting.")
        sys.exit(1)

    # Step 5: Write output
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT, index=False)
    log.info(f"Wrote {len(df):,} rows → {OUTPUT}")
    log.info(f"Total time: {time.time()-t_start:.1f}s")

    # Summary
    print("\n" + "=" * 60)
    print("MGI Database Summary")
    print("=" * 60)
    print(f"  Total interactions: {len(df):,}")
    print(f"  Unique metabolites: {df['HMDB_ID'].nunique():,}")
    print(f"  Unique genes:       {df['Gene_Symbol'].nunique():,}")
    print(f"  Unique organisms:   {df['Organism'].nunique():,}")
    print(f"  Interaction types:  {df['Interaction_Type'].nunique():,}")
    print(f"  Sources: {dict(df['Source'].value_counts())}")
    print(f"  Output: {OUTPUT}")


if __name__ == "__main__":
    main()
