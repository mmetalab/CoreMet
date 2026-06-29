"""
Build the Metabolite–Enzyme Interaction (MEI) database from KEGG + Rhea caches.

Pipeline:
  1.  EC number  →  KEGG compound IDs  (enzyme_compound_links.json)
  2.  EC number  →  gene IDs per organism  (enzyme_genes_{org}.json)
  3.  Gene IDs   →  UniProt accessions  (gene_uniprot_{org}.json)
  4.  Compound   →  name, SMILES  (compound_info.json)
  5.  UniProt    →  protein_name, gene_name  (uniprot_info.json)
  6.  HMDB ID mapping via MPIDB_v2 reverse-lookup
  7.  Rhea EC bridge  →  extra UniProt enzymes
  8.  Pathway annotations  (enzyme_pathway_links, pathway_names)
  9.  Deduplicate, validate, export
"""

import json
import csv
import os
from pathlib import Path
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).parent.parent
CACHE_KEGG = BASE_DIR / "data" / "cache" / "kegg"
CACHE_RHEA = BASE_DIR / "data" / "cache" / "rhea"
MPI_DB     = BASE_DIR / "data" / "mpidatabase" / "MPIDB_v2.csv"
OUT_DIR    = BASE_DIR / "data" / "databases"
OUT_FILE   = OUT_DIR / "mei_database.csv"

ORGANISMS = {
    "hsa": "Homo sapiens",
    "mmu": "Mus musculus",
    "rno": "Rattus norvegicus",
    "eco": "Escherichia coli",
    "bta": "Bos taurus",
    "pae": "Pseudomonas aeruginosa",
    "ath": "Arabidopsis thaliana",
    "sce": "Saccharomyces cerevisiae",
    "dme": "Drosophila melanogaster",
    "cel": "Caenorhabditis elegans",
}


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def load_tsv_lines(path: Path, max_lines=None):
    """Yield dicts from a TSV file."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            if max_lines and i >= max_lines:
                break
            yield row


# ── Step 1–2: Load KEGG caches ────────────────────────────────────────────────

def load_kegg_caches():
    """Return all KEGG data structures."""
    print("[1/9] Loading KEGG caches ...")

    ec_compounds = load_json(CACHE_KEGG / "enzyme_compound_links.json")
    compound_info = load_json(CACHE_KEGG / "compound_info.json")
    uniprot_info = load_json(CACHE_KEGG / "uniprot_info.json")

    # Per-organism: EC→genes and gene→uniprot
    ec_genes_by_org = {}
    gene_uniprot_by_org = {}
    for code in ORGANISMS:
        gf = CACHE_KEGG / f"enzyme_genes_{code}.json"
        uf = CACHE_KEGG / f"gene_uniprot_{code}.json"
        if gf.exists() and uf.exists():
            ec_genes_by_org[code] = load_json(gf)
            gene_uniprot_by_org[code] = load_json(uf)

    # Pathways
    ec_pathways = {}
    pw_names = {}
    pw_file = CACHE_KEGG / "enzyme_pathway_links.json"
    nm_file = CACHE_KEGG / "pathway_names.json"
    if pw_file.exists():
        ec_pathways = load_json(pw_file)
    if nm_file.exists():
        pw_names = load_json(nm_file)

    print(f"     {len(ec_compounds):,} ECs with compound links")
    print(f"     {len(compound_info):,} compounds with metadata")
    print(f"     {len(uniprot_info):,} UniProt entries with metadata")
    print(f"     {len(ec_genes_by_org)} organisms with gene mappings")

    return (ec_compounds, compound_info, uniprot_info,
            ec_genes_by_org, gene_uniprot_by_org,
            ec_pathways, pw_names)


# ── Step 3: Build HMDB reverse-lookup from MPIDB ──────────────────────────────

def build_hmdb_lookup():
    """Build metabolite name → HMDB ID and SMILES → HMDB ID from MPI DB."""
    print("[2/9] Building HMDB ID lookup from MPIDB_v2.csv ...")
    name_to_hmdb = {}
    smiles_to_hmdb = {}
    name_to_smiles = {}

    if not MPI_DB.exists():
        print("     WARNING: MPIDB_v2.csv not found, HMDB mapping will be empty")
        return name_to_hmdb, smiles_to_hmdb, name_to_smiles

    with open(MPI_DB) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get("Metabolite Name") or "").strip()
            hmdb = (row.get("HMDB ID") or "").strip()
            smiles = (row.get("SMILES") or "").strip()
            if name and hmdb:
                name_to_hmdb[name.lower()] = hmdb
            if smiles and hmdb:
                smiles_to_hmdb[smiles] = hmdb
            if name and smiles:
                name_to_smiles[name.lower()] = smiles

    print(f"     {len(name_to_hmdb):,} name→HMDB mappings")
    print(f"     {len(smiles_to_hmdb):,} SMILES→HMDB mappings")
    return name_to_hmdb, smiles_to_hmdb, name_to_smiles


# ── Step 4: Build Rhea → EC bridge ────────────────────────────────────────────

def build_rhea_ec_bridge():
    """Build Rhea reaction → EC number and Rhea → UniProt mappings."""
    print("[3/9] Building Rhea → EC bridge ...")
    rhea_to_ec = defaultdict(set)
    xref_file = CACHE_RHEA / "rhea2xrefs.tsv"
    if xref_file.exists():
        for row in load_tsv_lines(xref_file):
            if row.get("DB") == "EC":
                master = row.get("MASTER_ID", row.get("RHEA_ID", ""))
                ec = row.get("ID", "")
                if master and ec:
                    rhea_to_ec[master].add(ec)

    rhea_to_uniprot = defaultdict(set)
    uni_file = CACHE_RHEA / "rhea2uniprot.tsv"
    if uni_file.exists():
        for row in load_tsv_lines(uni_file):
            master = row.get("MASTER_ID", row.get("RHEA_ID", ""))
            uid = row.get("ID", "")
            if master and uid:
                rhea_to_uniprot[master].add(uid)

    print(f"     {len(rhea_to_ec):,} Rhea reactions with EC links")
    print(f"     {len(rhea_to_uniprot):,} Rhea reactions with UniProt links")
    return rhea_to_ec, rhea_to_uniprot


# ── Step 5: Build MEI records ─────────────────────────────────────────────────

def build_mei_records(ec_compounds, compound_info, uniprot_info,
                      ec_genes_by_org, gene_uniprot_by_org,
                      ec_pathways, pw_names,
                      name_to_hmdb, smiles_to_hmdb, name_to_smiles,
                      rhea_to_ec, rhea_to_uniprot):
    """Build the full MEI record list."""

    # ── 5a: For each EC, gather all UniProt IDs per organism (from KEGG) ──
    print("[4/9] Building EC → UniProt mapping from KEGG ...")
    ec_uniprots = defaultdict(lambda: defaultdict(set))  # ec → org_code → {uniprot_ids}
    for org_code, ec_genes in ec_genes_by_org.items():
        gene_uni = gene_uniprot_by_org.get(org_code, {})
        for ec, genes in ec_genes.items():
            for gene in genes:
                uid = gene_uni.get(gene, "")
                if uid:
                    ec_uniprots[ec][org_code].add(uid)

    # ── 5b: Add Rhea UniProt IDs via EC bridge ──
    print("[5/9] Adding Rhea enzyme links via EC bridge ...")
    rhea_additions = 0
    for rhea_id, ecs in rhea_to_ec.items():
        unis = rhea_to_uniprot.get(rhea_id, set())
        if not unis:
            continue
        for ec in ecs:
            # We add Rhea UniProts but mark them as organism="unknown" for now
            # unless we can determine organism from uniprot_info_rhea
            for uid in unis:
                ec_uniprots[ec]["rhea"].add(uid)
                rhea_additions += 1
    print(f"     {rhea_additions:,} Rhea enzyme additions")

    # Load Rhea UniProt info for organism detection
    rhea_uni_info = {}
    rhea_info_file = CACHE_RHEA / "uniprot_info_rhea.json"
    if rhea_info_file.exists():
        rhea_uni_info = load_json(rhea_info_file)

    # ── 5c: Build MEI rows ──
    print("[6/9] Generating MEI records ...")
    records = []
    seen = set()  # (hmdb_or_name, ec, uniprot) dedup key

    for ec, compounds in ec_compounds.items():
        # Get pathway info for this EC
        pathways = ec_pathways.get(ec, [])
        # Filter to map* pathways only
        map_pws = [p for p in pathways if p.startswith("map")]
        pw_id = map_pws[0] if map_pws else ""
        pw_name = pw_names.get(pw_id, "")

        for cpd_id in compounds:
            cpd = compound_info.get(cpd_id, {})
            met_name = cpd.get("name", "")
            smiles = cpd.get("smiles", "")
            if not met_name:
                continue

            # Resolve HMDB ID
            hmdb_id = ""
            name_lower = met_name.lower()
            if name_lower in name_to_hmdb:
                hmdb_id = name_to_hmdb[name_lower]
            elif smiles and smiles in smiles_to_hmdb:
                hmdb_id = smiles_to_hmdb[smiles]

            # If no SMILES from KEGG but found in our lookup, use that
            if not smiles and name_lower in name_to_smiles:
                smiles = name_to_smiles[name_lower]

            # For each organism with this EC
            for org_code, uniprots in ec_uniprots.get(ec, {}).items():
                org_name = ORGANISMS.get(org_code, "")

                for uid in uniprots:
                    dedup_key = (met_name.lower(), ec, uid)
                    if dedup_key in seen:
                        continue
                    seen.add(dedup_key)

                    # Get protein info
                    info = uniprot_info.get(uid, {})
                    if not info and uid in rhea_uni_info:
                        info = rhea_uni_info[uid]
                        if not org_name:
                            org_name = info.get("organism", "")

                    prot_name = info.get("protein_name", "")
                    gene_name = info.get("gene_name", "")

                    # Skip if we don't have organism info for Rhea-only entries
                    # unless they have one in rhea_uni_info
                    if org_code == "rhea" and not org_name:
                        continue

                    source = "KEGG" if org_code != "rhea" else "Rhea"

                    records.append({
                        "Species": org_name,
                        "Metabolite_Name": met_name,
                        "HMDB_ID": hmdb_id,
                        "SMILES": smiles,
                        "EC_Number": ec,
                        "Enzyme_Name": prot_name,
                        "Uniprot_ID": uid,
                        "Gene_Name": gene_name,
                        "KEGG_Compound": cpd_id,
                        "Pathway_ID": pw_id,
                        "Pathway_Name": pw_name,
                        "Evidence_Source": source,
                    })

    print(f"     {len(records):,} raw MEI records")
    return records


# ── Step 6: Deduplicate against MPIDB ─────────────────────────────────────────

def load_mpi_keys():
    """Load (metabolite_name_lower, uniprot_id) from MPIDB for dedup."""
    keys = set()
    if not MPI_DB.exists():
        return keys
    with open(MPI_DB) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = (row.get("Metabolite Name") or "").strip().lower()
            uid = (row.get("Uniprot ID") or "").strip()
            if name and uid:
                keys.add((name, uid))
    return keys


def deduplicate(records):
    """Remove records that already exist in MPIDB (same metabolite + UniProt)."""
    print("[7/9] Deduplicating against MPIDB_v2.csv ...")
    mpi_keys = load_mpi_keys()
    before = len(records)
    filtered = []
    for r in records:
        key = (r["Metabolite_Name"].lower(), r["Uniprot_ID"])
        if key not in mpi_keys:
            filtered.append(r)
    print(f"     {before:,} → {len(filtered):,} after dedup ({before - len(filtered):,} overlap)")
    return filtered


# ── Step 7: Statistics ─────────────────────────────────────────────────────────

def print_stats(records):
    """Print summary statistics."""
    print("[8/9] MEI Database Statistics:")
    print(f"     Total records:     {len(records):,}")
    ecs = set(r["EC_Number"] for r in records)
    mets = set(r["Metabolite_Name"].lower() for r in records)
    unis = set(r["Uniprot_ID"] for r in records)
    orgs = set(r["Species"] for r in records if r["Species"])
    hmdb_count = sum(1 for r in records if r["HMDB_ID"])
    smiles_count = sum(1 for r in records if r["SMILES"])
    kegg_count = sum(1 for r in records if r["Evidence_Source"] == "KEGG")
    rhea_count = sum(1 for r in records if r["Evidence_Source"] == "Rhea")

    print(f"     Unique ECs:        {len(ecs):,}")
    print(f"     Unique metabolites:{len(mets):,}")
    print(f"     Unique enzymes:    {len(unis):,}")
    print(f"     Organisms:         {len(orgs):,}: {sorted(orgs)}")
    print(f"     With HMDB ID:      {hmdb_count:,} ({100*hmdb_count/max(len(records),1):.1f}%)")
    print(f"     With SMILES:       {smiles_count:,} ({100*smiles_count/max(len(records),1):.1f}%)")
    print(f"     Source KEGG:       {kegg_count:,}")
    print(f"     Source Rhea:       {rhea_count:,}")


# ── Step 8: Export ─────────────────────────────────────────────────────────────

def export_csv(records):
    """Write MEI records to CSV."""
    print(f"[9/9] Writing {len(records):,} records to {OUT_FILE} ...")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "Species", "Metabolite_Name", "HMDB_ID", "SMILES",
        "EC_Number", "Enzyme_Name", "Uniprot_ID", "Gene_Name",
        "KEGG_Compound", "Pathway_ID", "Pathway_Name", "Evidence_Source",
    ]
    with open(OUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print(f"     Done! → {OUT_FILE}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Building CoreMet MEI Database")
    print("=" * 60)

    (ec_compounds, compound_info, uniprot_info,
     ec_genes_by_org, gene_uniprot_by_org,
     ec_pathways, pw_names) = load_kegg_caches()

    name_to_hmdb, smiles_to_hmdb, name_to_smiles = build_hmdb_lookup()
    rhea_to_ec, rhea_to_uniprot = build_rhea_ec_bridge()

    records = build_mei_records(
        ec_compounds, compound_info, uniprot_info,
        ec_genes_by_org, gene_uniprot_by_org,
        ec_pathways, pw_names,
        name_to_hmdb, smiles_to_hmdb, name_to_smiles,
        rhea_to_ec, rhea_to_uniprot,
    )

    records = deduplicate(records)
    print_stats(records)
    export_csv(records)

    print("\n✅ MEI database build complete!")


if __name__ == "__main__":
    main()
