#!/usr/bin/env python3
"""
Expand MDrI (Metabolite-Drug Interactions) database from 97 → 1,000+.

Strategy:
  1. Parse DrugBank XML → drug metadata, pathway metabolites, enzymes/targets
  2. Parse HMDB XML → metabolite metadata, DrugBank cross-refs, protein associations
  3. Cross-reference:
     a) HMDB metabolites that ARE drugs (have DrugBank IDs) → extract their enzyme/target interactions
     b) DrugBank pathway metabolites → direct drug-metabolite links
     c) Shared enzyme bridge: drug targets enzyme E ∧ metabolite is processed by enzyme E
  4. Merge with existing curated MDrI + deduplicate.

Output: data/databases/mdri_database.csv  (preserves existing schema)
"""
import os, sys, csv, time, json
from collections import defaultdict
from pathlib import Path
from lxml import etree

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = Path(__file__).resolve().parent.parent
DRUGBANK_XML = BASE / "data" / "raw" / "drugbank" / "drug_bank_database.xml"
HMDB_XML     = BASE / "data" / "raw" / "hmdb" / "hmdb_metabolites.xml"
EXISTING_DB  = BASE / "data" / "databases" / "mdri_database.csv"
OUTPUT_DB    = BASE / "data" / "databases" / "mdri_database.csv"
BACKUP_DB    = BASE / "data" / "databases" / "mdri_database.csv.bak"

DB_NS   = "{http://www.drugbank.ca}"
HMDB_NS = "{http://www.hmdb.ca}"

# ── Step 1: Parse DrugBank ───────────────────────────────────────────────────
def parse_drugbank():
    """
    Extract per-drug:
      - Basic info: DrugBank_ID, name, CAS, type (small molecule), SMILES, groups
      - Pathway metabolites: list of metabolite names/HMDB IDs tied to each drug
      - Enzymes / targets / transporters / carriers: UniProt IDs and actions
    Returns dict keyed by DrugBank_ID.
    """
    print("[1/3] Parsing DrugBank XML …")
    t0 = time.time()
    drugs = {}      # {DB_ID: {...}}
    pathway_links = []  # (drug_dbid, metabolite_name, metabolite_hmdb)

    context = etree.iterparse(str(DRUGBANK_XML), events=("end",), tag=f"{DB_NS}drug")
    n = 0

    for _ev, drug_el in context:
        # Only top-level <drug> elements
        parent = drug_el.getparent()
        if parent is None or parent.tag != f"{DB_NS}drugbank":
            drug_el.clear()
            continue

        n += 1
        dtype = drug_el.get("type", "")
        dbid  = drug_el.findtext(f"{DB_NS}drugbank-id[@primary='true']", "")
        name  = drug_el.findtext(f"{DB_NS}name", "")
        cas   = drug_el.findtext(f"{DB_NS}cas-number", "") or ""

        # Groups (approved, investigational, …)
        groups_el = drug_el.find(f"{DB_NS}groups")
        groups = []
        if groups_el is not None:
            groups = [g.text for g in groups_el.findall(f"{DB_NS}group") if g.text]

        # SMILES from calculated-properties
        smiles = ""
        calc = drug_el.find(f"{DB_NS}calculated-properties")
        if calc is not None:
            for prop in calc.findall(f"{DB_NS}property"):
                if prop.findtext(f"{DB_NS}kind") == "SMILES":
                    smiles = prop.findtext(f"{DB_NS}value", "")
                    break

        # Categories
        cats_el = drug_el.find(f"{DB_NS}categories")
        categories = []
        if cats_el is not None:
            categories = [c.findtext(f"{DB_NS}category") for c in cats_el.findall(f"{DB_NS}category")]
            categories = [c for c in categories if c]

        # External IDs (look for HMDB, PubChem, KEGG)
        ext_map = {}
        ext_el = drug_el.find(f"{DB_NS}external-identifiers")
        if ext_el is not None:
            for eid in ext_el.findall(f"{DB_NS}external-identifier"):
                res = eid.findtext(f"{DB_NS}resource", "")
                val = eid.findtext(f"{DB_NS}identifier", "")
                ext_map[res] = val

        # Pathway metabolites
        pathways_el = drug_el.find(f"{DB_NS}pathways")
        if pathways_el is not None:
            for pw in pathways_el.findall(f"{DB_NS}pathway"):
                metabs_el = pw.find(f"{DB_NS}metabolites") if pw is not None else None
                if metabs_el is not None:
                    for m in metabs_el:
                        mname = m.findtext(f"{DB_NS}name", "")
                        mhmdb = m.findtext(f"{DB_NS}hmdb-id", "")
                        if mname or mhmdb:
                            pathway_links.append((dbid, mname, mhmdb))

        # Enzymes / targets / transporters / carriers
        protein_links = []
        for section_name, child_tag in [
            ("targets", "target"),
            ("enzymes", "enzyme"),
            ("transporters", "transporter"),
            ("carriers", "carrier"),
        ]:
            sec_el = drug_el.find(f"{DB_NS}{section_name}")
            if sec_el is None:
                continue
            for item in sec_el.findall(f"{DB_NS}{child_tag}"):
                pname = item.findtext(f"{DB_NS}name", "")
                organism = item.findtext(f"{DB_NS}organism", "")
                # Skip non-human
                if organism and "human" not in organism.lower() and "homo" not in organism.lower():
                    continue
                known = item.findtext(f"{DB_NS}known-action", "")
                actions_el = item.find(f"{DB_NS}actions")
                actions = []
                if actions_el is not None:
                    actions = [a.text for a in actions_el.findall(f"{DB_NS}action") if a.text]
                # UniProt IDs
                poly_el = item.find(f"{DB_NS}polypeptide")
                uniprot = ""
                if poly_el is not None:
                    uniprot = poly_el.get("id", "")
                if pname:
                    protein_links.append({
                        "role": section_name, "name": pname,
                        "actions": actions, "known_action": known,
                        "uniprot": uniprot,
                    })

        drugs[dbid] = {
            "name": name, "cas": cas, "type": dtype,
            "smiles": smiles, "groups": groups,
            "categories": categories, "ext": ext_map,
            "proteins": protein_links,
        }

        drug_el.clear()
        if n % 3000 == 0:
            print(f"  … {n} drugs parsed")

    print(f"  Done: {n} drugs, {len(pathway_links)} pathway-metabolite links  ({time.time()-t0:.1f}s)")
    return drugs, pathway_links


# ── Step 2: Parse HMDB ───────────────────────────────────────────────────────
def parse_hmdb():
    """
    Extract per-metabolite:
      - HMDB_ID, name, CAS, SMILES, DrugBank_ID
      - Protein associations (name, uniprot_id)
    Returns dict keyed by HMDB_ID, and a name→HMDB_ID index.
    """
    print("[2/3] Parsing HMDB XML …")
    t0 = time.time()
    metabs = {}      # {HMDB_ID: {...}}
    name_idx = {}    # {lowercase_name: HMDB_ID}

    context = etree.iterparse(str(HMDB_XML), events=("end",), tag=f"{HMDB_NS}metabolite")
    n = 0

    for _ev, mel in context:
        parent = mel.getparent()
        if parent is None or parent.tag != f"{HMDB_NS}hmdb":
            mel.clear()
            continue

        n += 1
        hmdb_id = mel.findtext(f"{HMDB_NS}accession", "")
        name    = mel.findtext(f"{HMDB_NS}name", "")
        cas     = mel.findtext(f"{HMDB_NS}cas_registry_number", "") or ""
        smiles  = mel.findtext(f"{HMDB_NS}smiles", "") or ""
        db_id   = mel.findtext(f"{HMDB_NS}drugbank_id", "") or ""

        # Protein associations
        prot_assoc = []
        pa_el = mel.find(f"{HMDB_NS}protein_associations")
        if pa_el is not None:
            for p in pa_el.findall(f"{HMDB_NS}protein"):
                pname = p.findtext(f"{HMDB_NS}name", "")
                uniprot = p.findtext(f"{HMDB_NS}uniprot_id", "") or ""
                if pname:
                    prot_assoc.append({"name": pname, "uniprot": uniprot})

        # Synonyms for better matching
        syn_el = mel.find(f"{HMDB_NS}synonyms")
        synonyms = []
        if syn_el is not None:
            synonyms = [s.text for s in syn_el.findall(f"{HMDB_NS}synonym") if s.text]

        metabs[hmdb_id] = {
            "name": name, "cas": cas, "smiles": smiles,
            "drugbank_id": db_id, "proteins": prot_assoc,
            "synonyms": synonyms,
        }
        name_idx[name.lower()] = hmdb_id
        for syn in synonyms:
            name_idx[syn.lower()] = hmdb_id

        mel.clear()
        if n % 40000 == 0:
            print(f"  … {n} metabolites parsed")

    print(f"  Done: {n} metabolites, {sum(1 for v in metabs.values() if v['drugbank_id'])} with DrugBank ID  ({time.time()-t0:.1f}s)")
    return metabs, name_idx


# ── Step 3: Build expanded MDrI ──────────────────────────────────────────────
def build_mdri(drugs, pathway_links, metabs, name_idx):
    """
    Link strategy:
      A) Pathway metabolites: DrugBank pathways directly list metabolite names/HMDB
      B) Metabolite-IS-drug: HMDB metabolites with DrugBank IDs → take drug's enzyme/target actions
      C) Shared enzyme bridge: Drug X acts on enzyme E; metabolite Y is associated with enzyme E
    """
    print("[3/3] Building MDrI interactions …")
    rows = []   # list of dicts (matches CSV schema)
    seen = set()  # (HMDB_ID, DrugBank_ID, interaction_type)  dedup key

    # Helper
    def add_row(hmdb_id, metab_name, smiles, drug_name, drug_dbid,
                interaction_type, tissue, evidence_level, description, pmid, source):
        key = (hmdb_id, drug_dbid, interaction_type, metab_name.lower())
        if key in seen:
            return
        seen.add(key)
        rows.append({
            "Metabolite_Name": metab_name,
            "HMDB_ID": hmdb_id,
            "SMILES": smiles,
            "Drug_Name": drug_name,
            "DrugBank_ID": drug_dbid,
            "Interaction_Type": interaction_type,
            "Tissue": tissue,
            "Cell_Location": "",
            "Biospecimen": "",
            "Evidence_Level": evidence_level,
            "Description": description,
            "PMID": pmid,
            "Source": source,
        })

    # ── A) Pathway metabolites ──────────────────────────────────────────────
    print("  Strategy A: pathway-metabolite links …")
    for drug_dbid, mname, mhmdb in pathway_links:
        drug_info = drugs.get(drug_dbid)
        if not drug_info:
            continue
        # Resolve HMDB
        hmdb_id = mhmdb
        if not hmdb_id:
            hmdb_id = name_idx.get(mname.lower(), "")
        metab_info = metabs.get(hmdb_id, {})
        mname_resolved = metab_info.get("name", mname) if metab_info else mname
        smiles = metab_info.get("smiles", "") if metab_info else ""
        if hmdb_id or mname:
            add_row(
                hmdb_id=hmdb_id,
                metab_name=mname_resolved,
                smiles=smiles,
                drug_name=drug_info["name"],
                drug_dbid=drug_dbid,
                interaction_type="Pharmacokinetic",
                tissue="",
                evidence_level="Literature" if "approved" in drug_info["groups"] else "Predicted",
                description=f"Metabolite in {drug_info['name']} pharmacokinetic pathway",
                pmid="",
                source="DrugBank_pathway",
            )
    print(f"    → {len(rows)} after pathway links")

    # ── B) Metabolite-IS-drug: enzyme/target interactions ────────────────────
    print("  Strategy B: HMDB metabolites that are DrugBank drugs …")
    before = len(rows)
    for hmdb_id, minfo in metabs.items():
        db_id = minfo.get("drugbank_id", "")
        if not db_id or db_id not in drugs:
            continue
        drug_info = drugs[db_id]
        # The metabolite IS a drug — record its enzyme/target interactions
        for pl in drug_info.get("proteins", []):
            role = pl["role"]
            actions = pl.get("actions", [])
            action_str = "/".join(actions) if actions else "unknown"
            if role == "targets":
                itype = "Pharmacodynamic"
                desc = f"{minfo['name']} (as drug {drug_info['name']}) {action_str} target {pl['name']}"
            elif role == "enzymes":
                itype = "Pharmacokinetic"
                desc = f"{minfo['name']} metabolized by enzyme {pl['name']} ({action_str})"
            elif role == "transporters":
                itype = "Pharmacokinetic"
                desc = f"{minfo['name']} transported by {pl['name']} ({action_str})"
            elif role == "carriers":
                itype = "Pharmacokinetic"
                desc = f"{minfo['name']} carried by {pl['name']} ({action_str})"
            else:
                continue
            add_row(
                hmdb_id=hmdb_id,
                metab_name=minfo["name"],
                smiles=minfo.get("smiles", ""),
                drug_name=drug_info["name"],
                drug_dbid=db_id,
                interaction_type=itype,
                tissue="",
                evidence_level="Literature",
                description=desc,
                pmid="",
                source="DrugBank_cross_ref",
            )
    print(f"    → +{len(rows) - before} new ({len(rows)} total)")

    # ── C) Shared enzyme bridge (stringent) ────────────────────────────────
    print("  Strategy C: shared-enzyme bridge (stringent) …")
    before = len(rows)

    # Build enzyme→metabolites index from HMDB protein associations
    enzyme_to_metabs = defaultdict(list)  # uniprot_id → [(hmdb_id, name, smiles)]
    for hmdb_id, minfo in metabs.items():
        for pa in minfo.get("proteins", []):
            uid = pa.get("uniprot", "")
            if uid:
                enzyme_to_metabs[uid].append((hmdb_id, minfo["name"], minfo.get("smiles", "")))

    # Build enzyme→drugs index to measure specificity
    enzyme_to_drugs = defaultdict(int)
    for dbid, dinfo in drugs.items():
        if "approved" not in dinfo.get("groups", []):
            continue
        for pl in dinfo.get("proteins", []):
            uid = pl.get("uniprot", "")
            if uid:
                enzyme_to_drugs[uid] += 1

    # Stringent filters to avoid combinatorial explosion
    MAX_METABS_PER_ENZYME = 50   # skip overly promiscuous enzymes
    MAX_DRUGS_PER_ENZYME  = 50   # skip very common drug targets
    REQUIRE_KNOWN_ACTION  = True  # drug must have known-action=yes

    skipped_promiscuous = 0

    for dbid, dinfo in drugs.items():
        if "approved" not in dinfo.get("groups", []):
            continue
        for pl in dinfo.get("proteins", []):
            uid = pl.get("uniprot", "")
            if not uid or uid not in enzyme_to_metabs:
                continue
            # Stringent: require known action
            if REQUIRE_KNOWN_ACTION and pl.get("known_action", "").lower() != "yes":
                continue
            # Only enzymes and transporters (most relevant to metabolite pharmacokinetics)
            role = pl["role"]
            if role not in ("enzymes", "transporters"):
                continue
            # Skip promiscuous enzymes
            n_metabs = len(enzyme_to_metabs[uid])
            n_drugs = enzyme_to_drugs.get(uid, 0)
            if n_metabs > MAX_METABS_PER_ENZYME or n_drugs > MAX_DRUGS_PER_ENZYME:
                skipped_promiscuous += 1
                continue

            actions = pl.get("actions", [])
            action_str = "/".join(actions) if actions else "modulates"
            for hmdb_id, mname, msmiles in enzyme_to_metabs[uid]:
                mdb = metabs.get(hmdb_id, {}).get("drugbank_id", "")
                if mdb == dbid:
                    continue
                desc = f"{dinfo['name']} {action_str} {pl['name']} ({role}); " \
                       f"{mname} also associated with this protein"
                add_row(
                    hmdb_id=hmdb_id,
                    metab_name=mname,
                    smiles=msmiles,
                    drug_name=dinfo["name"],
                    drug_dbid=dbid,
                    interaction_type="Pharmacokinetic",
                    tissue="",
                    evidence_level="Predicted",
                    description=desc,
                    pmid="",
                    source="DrugBank_enzyme_bridge",
                )
    print(f"    Skipped {skipped_promiscuous} promiscuous enzyme pairs")
    print(f"    → +{len(rows) - before} new ({len(rows)} total)")

    return rows


# ── Step 4: Merge with existing curated DB ───────────────────────────────────
def merge_and_save(new_rows):
    """Merge new rows with existing curated MDrI, preserving curated entries."""
    print("Merging with existing curated database …")

    # Read existing
    existing = []
    if EXISTING_DB.exists():
        with open(EXISTING_DB) as f:
            reader = csv.DictReader(f)
            existing = list(reader)
        print(f"  Existing curated entries: {len(existing)}")
        # Backup
        import shutil
        shutil.copy2(EXISTING_DB, BACKUP_DB)
        print(f"  Backup saved to {BACKUP_DB}")

    # Build dedup set from existing
    existing_keys = set()
    for r in existing:
        key = (r.get("HMDB_ID",""), r.get("DrugBank_ID",""),
               r.get("Interaction_Type",""), r.get("Metabolite_Name","").lower())
        existing_keys.add(key)

    # Filter new rows to avoid duplicates with curated ones
    added = 0
    for r in new_rows:
        key = (r["HMDB_ID"], r["DrugBank_ID"], r["Interaction_Type"], r["Metabolite_Name"].lower())
        if key not in existing_keys:
            existing.append(r)
            existing_keys.add(key)
            added += 1

    # Sort
    existing.sort(key=lambda r: (r.get("Metabolite_Name","").lower(), r.get("Drug_Name","").lower()))

    # Write
    fields = [
        "Metabolite_Name","HMDB_ID","SMILES","Drug_Name","DrugBank_ID",
        "Interaction_Type","Tissue","Cell_Location","Biospecimen",
        "Evidence_Level","Description","PMID","Source",
    ]
    with open(OUTPUT_DB, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(existing)

    total = len(existing)
    print(f"\n{'='*60}")
    print(f"MDrI database expanded: {total - added} curated + {added} new = {total} total")
    print(f"Saved to {OUTPUT_DB}")

    # Print source breakdown
    from collections import Counter
    src_counts = Counter(r.get("Source","") for r in existing)
    for src, cnt in src_counts.most_common():
        print(f"  {src}: {cnt}")

    return total


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    t_start = time.time()

    # Validate paths
    for path, label in [(DRUGBANK_XML, "DrugBank"), (HMDB_XML, "HMDB")]:
        if not path.exists():
            print(f"ERROR: {label} XML not found at {path}")
            sys.exit(1)

    drugs, pw_links = parse_drugbank()
    metabs, name_idx = parse_hmdb()
    new_rows = build_mdri(drugs, pw_links, metabs, name_idx)
    total = merge_and_save(new_rows)

    print(f"\nTotal time: {time.time() - t_start:.1f}s")
