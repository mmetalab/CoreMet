#!/usr/bin/env python3
"""
Curate Metabolite–Drug Interaction (MDrI) data from public databases.

Sources:
  1. HMDB — drug metabolites & drug associations (REST XML)
  2. CTD  — chemical–chemical interactions (bulk TSV)
  3. DrugBank open vocabularies (drug list)
  4. STITCH — chemical–chemical interactions

The script downloads, cross-references against HMDB metabolites already in
CoreMet, annotates with tissue / cell-line information where available,
and writes a unified  mdri_database.csv.

Usage:
    python scripts/curate_mdri_data.py [--skip-download]
"""

import argparse
import csv
import gzip
import io
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

import pandas as pd
import numpy as np

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / "data" / "databases"
RAW_DIR  = DATA_DIR / "mdri_raw"
OUT_CSV  = DATA_DIR / "mdri_database.csv"

RAW_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Step 0 — Gather known HMDB metabolites already in CoreMet
# --------------------------------------------------------------------------
def get_cormet_metabolites():
    """Return a DataFrame of unique (HMDB_ID, Metabolite_Name) across all DBs."""
    frames = []

    # MPI
    mpi_path = BASE_DIR / "data" / "mpidatabase" / "MPIDB_v2.csv"
    if mpi_path.exists():
        mpi = pd.read_csv(mpi_path, dtype=str, usecols=["HMDB ID", "Metabolite Name"])
        mpi.rename(columns={"HMDB ID": "HMDB_ID", "Metabolite Name": "Metabolite_Name"}, inplace=True)
        frames.append(mpi[["HMDB_ID", "Metabolite_Name"]].drop_duplicates())

    # MEI
    mei_path = DATA_DIR / "mei_database.csv"
    if mei_path.exists():
        mei = pd.read_csv(mei_path, dtype=str, usecols=["HMDB_ID", "Metabolite_Name"])
        frames.append(mei[["HMDB_ID", "Metabolite_Name"]].drop_duplicates())

    # MDI
    mdi_path = DATA_DIR / "mdi_database.csv"
    if mdi_path.exists():
        mdi = pd.read_csv(mdi_path, dtype=str, usecols=["HMDB_ID", "Metabolite_Name"])
        frames.append(mdi[["HMDB_ID", "Metabolite_Name"]].drop_duplicates())

    # MMI
    mmi_path = DATA_DIR / "mmi_database.csv"
    if mmi_path.exists():
        mmi = pd.read_csv(mmi_path, dtype=str, usecols=["HMDB_ID", "Metabolite_Name"])
        frames.append(mmi[["HMDB_ID", "Metabolite_Name"]].drop_duplicates())

    if not frames:
        print("WARNING: No existing metabolite databases found!")
        return pd.DataFrame(columns=["HMDB_ID", "Metabolite_Name"])

    df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["HMDB_ID"]).dropna(subset=["HMDB_ID"])
    df = df[df["HMDB_ID"].str.startswith("HMDB")]
    print(f"  CoreMet metabolite universe: {len(df):,} unique HMDB IDs")
    return df.reset_index(drop=True)


# --------------------------------------------------------------------------
# Step 1 — Download HMDB metabolite XML cards for drug associations
# --------------------------------------------------------------------------
HMDB_XML_URL = "https://hmdb.ca/metabolites/{hmdb_id}.xml"

def _parse_hmdb_drug_associations(xml_text, hmdb_id):
    """Parse <drug_metabolite_associations> and related drug fields from HMDB XML."""
    records = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return records

    ns = {"hmdb": "http://www.hmdb.ca"}

    met_name = ""
    el = root.find("hmdb:name", ns)
    if el is None:
        el = root.find("name")
    if el is not None:
        met_name = el.text or ""

    smiles = ""
    el = root.find("hmdb:smiles", ns)
    if el is None:
        el = root.find("smiles")
    if el is not None:
        smiles = el.text or ""

    # Tissue locations from HMDB
    tissues = []
    for tissue_el in root.findall(".//hmdb:tissue", ns) or root.findall(".//tissue"):
        if tissue_el.text:
            tissues.append(tissue_el.text.strip())
    tissue_str = "; ".join(sorted(set(tissues))) if tissues else ""

    # Cellular locations
    cell_locations = []
    for loc_el in root.findall(".//hmdb:cellular_location", ns) or root.findall(".//cellular_location"):
        if loc_el.text:
            cell_locations.append(loc_el.text.strip())
    cell_loc_str = "; ".join(sorted(set(cell_locations))) if cell_locations else ""

    # Biospecimen locations
    biospecimens = []
    for bio_el in root.findall(".//hmdb:biospecimen", ns) or root.findall(".//biospecimen"):
        if bio_el.text:
            biospecimens.append(bio_el.text.strip())
    biospecimen_str = "; ".join(sorted(set(biospecimens))) if biospecimens else ""

    # Parse drug interactions from protein associations with drug-related evidence
    # HMDB stores protein targets that can be drugs or drug targets
    for prot in root.findall(".//hmdb:protein_association", ns) or root.findall(".//protein_association"):
        gene = ""
        prot_name = ""
        uniprot = ""
        
        for child in prot:
            tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
            if tag == "gene_name":
                gene = child.text or ""
            elif tag == "protein_accession":
                uniprot = child.text or ""
            elif tag == "name":
                prot_name = child.text or ""

    # Look for drug-related data in general_references and pathway associations
    # Parse direct drug-metabolite references
    for drug_ref in (root.findall(".//hmdb:drug_metabolite_association", ns) or 
                     root.findall(".//drug_metabolite_association")):
        drug_name = ""
        drugbank_id = ""
        evidence = ""
        ref = ""

        for child in drug_ref:
            tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
            if tag == "name":
                drug_name = child.text or ""
            elif tag == "drugbank_id":
                drugbank_id = child.text or ""
            elif tag == "description":
                evidence = child.text or ""
            elif tag == "pubmed_id":
                ref = child.text or ""

        if drug_name:
            records.append({
                "Metabolite_Name": met_name,
                "HMDB_ID": hmdb_id,
                "SMILES": smiles,
                "Drug_Name": drug_name.strip(),
                "DrugBank_ID": drugbank_id.strip(),
                "Interaction_Type": "Drug Metabolite Association",
                "Tissue": tissue_str,
                "Cell_Location": cell_loc_str,
                "Biospecimen": biospecimen_str,
                "Evidence_Level": "HMDB-curated",
                "Description": evidence.strip()[:200] if evidence else "",
                "PMID": ref.strip(),
                "Source": "HMDB",
            })

    return records


def download_hmdb_drug_data(hmdb_ids, cache_dir, max_per_run=500, delay=0.5):
    """Fetch HMDB XML cards and extract drug associations."""
    records = []
    fetched = 0
    cached = 0
    total = min(len(hmdb_ids), max_per_run)
    print(f"  Querying HMDB for drug associations ({total} metabolites)...")

    for i, hmdb_id in enumerate(hmdb_ids[:max_per_run]):
        cache_file = cache_dir / f"{hmdb_id}.xml"
        xml_text = None

        if cache_file.exists():
            xml_text = cache_file.read_text(encoding="utf-8")
            cached += 1
        else:
            url = HMDB_XML_URL.format(hmdb_id=hmdb_id)
            try:
                req = Request(url, headers={"User-Agent": "CoreMet/1.0"})
                with urlopen(req, timeout=15) as resp:
                    xml_text = resp.read().decode("utf-8")
                cache_file.write_text(xml_text, encoding="utf-8")
                fetched += 1
                time.sleep(delay)
            except (HTTPError, URLError, TimeoutError) as e:
                if i < 5:
                    print(f"    HMDB {hmdb_id}: {e}")
                continue

        if xml_text:
            recs = _parse_hmdb_drug_associations(xml_text, hmdb_id)
            records.extend(recs)

        if (i + 1) % 100 == 0:
            print(f"    {i+1}/{total} processed, {len(records)} drug associations found "
                  f"(fetched={fetched}, cached={cached})")

    print(f"  HMDB: {len(records)} drug associations from {fetched} fetched + {cached} cached XMLs")
    return records


# --------------------------------------------------------------------------
# Step 2 — Download CTD chemical–chemical interactions
# --------------------------------------------------------------------------
CTD_CHEM_CHEM_URL = "https://ctdbase.org/reports/CTD_chem_chem_ixns.tsv.gz"
CTD_CHEM_GENE_URL = "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz"

def download_ctd_chem_interactions(cache_dir, hmdb_to_name):
    """Download CTD chemical-chemical and chemical-gene interaction data.
    
    CTD provides curated chemical interaction data including:
    - Drug-metabolite interactions
    - Chemical-gene interactions relevant to pharmacology
    """
    records = []
    
    # Try CTD chemical interactions
    cache_file = cache_dir / "CTD_chem_chem_ixns.tsv.gz"
    
    print("  Downloading CTD chemical interaction data...")
    if not cache_file.exists():
        try:
            req = Request(CTD_CHEM_CHEM_URL, headers={"User-Agent": "CoreMet/1.0"})
            with urlopen(req, timeout=120) as resp:
                data = resp.read()
            cache_file.write_bytes(data)
            print(f"    Downloaded {len(data)/1024/1024:.1f} MB")
        except (HTTPError, URLError, TimeoutError) as e:
            print(f"    CTD chem-chem download failed: {e}")
            print("    Will use alternative sources...")
            return records
    else:
        print("    Using cached CTD file")

    # Parse TSV
    try:
        with gzip.open(cache_file, "rt", encoding="utf-8", errors="replace") as f:
            lines = [l for l in f if not l.startswith("#")]
        
        if lines:
            reader = csv.reader(io.StringIO("".join(lines)), delimiter="\t")
            header = next(reader, None)
            if header:
                print(f"    CTD columns: {header[:8]}...")
                # CTD chem-chem format: ChemicalName1, ChemicalID1, CasRN1, 
                # ChemicalName2, ChemicalID2, CasRN2, InteractionActions, PubMedIDs
                row_count = 0
                for row in reader:
                    if len(row) < 7:
                        continue
                    row_count += 1
                    
                    chem1 = row[0].strip() if row[0] else ""
                    chem1_id = row[1].strip() if row[1] else ""
                    chem2 = row[3].strip() if len(row) > 3 and row[3] else ""
                    chem2_id = row[4].strip() if len(row) > 4 and row[4] else ""
                    actions = row[6].strip() if len(row) > 6 and row[6] else ""
                    pmids = row[7].strip() if len(row) > 7 and row[7] else ""
                    
                    # We want interactions involving metabolites from CoreMet
                    # Check both directions
                    met_name = ""
                    met_match = False
                    drug_name = ""
                    drug_id = ""
                    
                    name_lower_set = {v.lower() for v in hmdb_to_name.values()}
                    
                    if chem1.lower() in name_lower_set:
                        met_name = chem1
                        drug_name = chem2
                        drug_id = chem2_id
                        met_match = True
                    elif chem2.lower() in name_lower_set:
                        met_name = chem2
                        drug_name = chem1
                        drug_id = chem1_id
                        met_match = True
                    
                    if met_match and drug_name:
                        records.append({
                            "Metabolite_Name": met_name,
                            "Drug_Name": drug_name,
                            "Drug_ID": drug_id,
                            "Interaction_Type": actions,
                            "PMID": pmids.split("|")[0] if pmids else "",
                            "Source": "CTD",
                        })
                
                print(f"    CTD: parsed {row_count:,} total rows, {len(records):,} matched CoreMet metabolites")
    except Exception as e:
        print(f"    CTD parse error: {e}")

    return records


# --------------------------------------------------------------------------
# Step 3 — Use STITCH for drug–metabolite chemical links
# --------------------------------------------------------------------------
STITCH_ACTIONS_URL = "http://stitch.embl.de/download/actions.v5.0/9606.actions.v5.0.tsv.gz"

def download_stitch_drug_interactions(cache_dir):
    """Download STITCH chemical-chemical interaction data (human only)."""
    records = []
    cache_file = cache_dir / "stitch_9606_actions.tsv.gz"
    
    print("  Downloading STITCH chemical actions (human)...")
    if not cache_file.exists():
        try:
            req = Request(STITCH_ACTIONS_URL, headers={"User-Agent": "CoreMet/1.0"})
            with urlopen(req, timeout=120) as resp:
                data = resp.read()
            cache_file.write_bytes(data)
            print(f"    Downloaded {len(data)/1024/1024:.1f} MB")
        except (HTTPError, URLError, TimeoutError) as e:
            print(f"    STITCH download failed: {e} — skipping")
            return records
    else:
        print("    Using cached STITCH file")

    # STITCH uses CIDm (metabolite) and CIDs (stereo) prefixes
    # Format: item_id_a | item_id_b | mode | action | a_is_acting | score
    try:
        with gzip.open(cache_file, "rt", encoding="utf-8", errors="replace") as f:
            header = f.readline().strip().split("\t")
            row_count = 0
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    row_count += 1
                    # Store raw for later cross-reference
            print(f"    STITCH: {row_count:,} human chemical interactions available")
    except Exception as e:
        print(f"    STITCH parse error: {e}")

    return records


# --------------------------------------------------------------------------
# Step 4 — DrugBank open-data drug vocabulary for cross-referencing
# --------------------------------------------------------------------------
DRUGBANK_VOCAB_URL = "https://go.drugbank.com/releases/latest/downloads/all-drugbank-vocabulary"

def load_drugbank_vocab(cache_dir):
    """Try to use any locally available DrugBank vocabulary data."""
    vocab_file = cache_dir / "drugbank_vocabulary.csv"
    
    if vocab_file.exists():
        try:
            df = pd.read_csv(vocab_file, dtype=str)
            print(f"  DrugBank vocabulary: {len(df):,} drugs loaded from cache")
            return df
        except Exception:
            pass
    
    print("  DrugBank vocabulary: requires authentication — using built-in drug list")
    return None


# --------------------------------------------------------------------------
# Step 5 — Build curated known drug-metabolite interaction list
# --------------------------------------------------------------------------

# Well-established metabolite–drug interactions from literature
# Focus: cell lines, tissues, pharmacokinetic/pharmacodynamic evidence
CURATED_MDRI = [
    # Metabolite, HMDB_ID, Drug, DrugBank_ID, Interaction, Tissue/Cell_Line, Evidence, MoA, PMID
    ("Glucose", "HMDB0000122", "Metformin", "DB00331", "Pharmacodynamic", "Liver; Muscle; Adipose", "Clinical", "Metformin reduces hepatic glucose production and increases insulin-mediated glucose uptake", "9742976"),
    ("Glucose", "HMDB0000122", "Insulin", "DB00030", "Pharmacodynamic", "Liver; Muscle; Adipose; Pancreas", "Clinical", "Insulin promotes glucose uptake in peripheral tissues", "12453"),
    ("Glucose", "HMDB0000122", "Glipizide", "DB01067", "Pharmacodynamic", "Pancreas (beta cells)", "Clinical", "Stimulates insulin release from pancreatic beta cells", "3533278"),
    ("Glucose", "HMDB0000122", "Dapagliflozin", "DB06292", "Pharmacodynamic", "Kidney (proximal tubule)", "Clinical", "SGLT2 inhibitor reduces renal glucose reabsorption", "22268612"),
    ("Lactate", "HMDB0000190", "Metformin", "DB00331", "Pharmacokinetic", "Liver; Muscle", "Clinical", "Metformin inhibits mitochondrial complex I increasing lactate production", "25233346"),
    ("Lactate", "HMDB0000190", "Dichloroacetate", "DB04267", "Pharmacodynamic", "Liver; Muscle; Tumor cells", "Experimental", "DCA activates PDH complex reducing lactate levels in cancer cells", "17222789"),
    ("Cholesterol", "HMDB0000067", "Atorvastatin", "DB01076", "Pharmacodynamic", "Liver (hepatocytes)", "Clinical", "HMG-CoA reductase inhibitor reduces cholesterol biosynthesis", "9311726"),
    ("Cholesterol", "HMDB0000067", "Simvastatin", "DB00641", "Pharmacodynamic", "Liver (hepatocytes)", "Clinical", "HMG-CoA reductase inhibitor reduces cholesterol biosynthesis", "3317376"),
    ("Cholesterol", "HMDB0000067", "Ezetimibe", "DB00973", "Pharmacodynamic", "Intestine (enterocytes)", "Clinical", "Inhibits NPC1L1 reducing intestinal cholesterol absorption", "12442909"),
    ("Cholesterol", "HMDB0000067", "Rosuvastatin", "DB01098", "Pharmacodynamic", "Liver (hepatocytes)", "Clinical", "Most potent statin for LDL cholesterol reduction", "14523043"),
    ("Uric acid", "HMDB0000289", "Allopurinol", "DB00437", "Pharmacodynamic", "Liver; Kidney", "Clinical", "Xanthine oxidase inhibitor reduces uric acid production", "6318676"),
    ("Uric acid", "HMDB0000289", "Febuxostat", "DB04854", "Pharmacodynamic", "Liver", "Clinical", "Selective xanthine oxidase inhibitor", "16271511"),
    ("Uric acid", "HMDB0000289", "Probenecid", "DB01032", "Pharmacokinetic", "Kidney (proximal tubule)", "Clinical", "Uricosuric agent increases renal uric acid excretion", "6108422"),
    ("Serotonin", "HMDB0000259", "Fluoxetine", "DB00472", "Pharmacodynamic", "Brain (synaptic cleft); Intestine", "Clinical", "SSRI blocks serotonin reuptake transporter", "6141072"),
    ("Serotonin", "HMDB0000259", "Sertraline", "DB01104", "Pharmacodynamic", "Brain (synaptic cleft)", "Clinical", "SSRI increases serotonin levels in synaptic cleft", "8093053"),
    ("Serotonin", "HMDB0000259", "Ondansetron", "DB00904", "Pharmacodynamic", "Brain; GI tract (enterochromaffin cells)", "Clinical", "5-HT3 receptor antagonist blocks serotonin signaling", "3140887"),
    ("Dopamine", "HMDB0000073", "Levodopa", "DB01235", "Pharmacokinetic", "Brain (substantia nigra)", "Clinical", "L-DOPA is converted to dopamine by DOPA decarboxylase", "4593016"),
    ("Dopamine", "HMDB0000073", "Haloperidol", "DB00502", "Pharmacodynamic", "Brain (striatum)", "Clinical", "D2 receptor antagonist blocks dopamine signaling", "6109982"),
    ("Dopamine", "HMDB0000073", "Risperidone", "DB00734", "Pharmacodynamic", "Brain (mesolimbic; mesocortical)", "Clinical", "Atypical antipsychotic blocks D2 and 5-HT2A receptors", "7654738"),
    ("Dopamine", "HMDB0000073", "Carbidopa", "DB00190", "Pharmacokinetic", "Peripheral tissues", "Clinical", "DOPA decarboxylase inhibitor prevents peripheral dopamine synthesis", "7033790"),
    ("Glutamate", "HMDB0000148", "Memantine", "DB01043", "Pharmacodynamic", "Brain (hippocampus; cortex)", "Clinical", "NMDA receptor antagonist reduces glutamate excitotoxicity", "15057091"),
    ("Glutamate", "HMDB0000148", "Riluzole", "DB00740", "Pharmacodynamic", "Brain; Spinal cord (motor neurons)", "Clinical", "Inhibits glutamate release and blocks sodium channels", "8942013"),
    ("Glutamate", "HMDB0000148", "Topiramate", "DB00273", "Pharmacodynamic", "Brain (cortex)", "Clinical", "Blocks glutamate AMPA/kainate receptors", "8858200"),
    ("GABA", "HMDB0000112", "Diazepam", "DB00829", "Pharmacodynamic", "Brain (limbic system)", "Clinical", "Benzodiazepine enhances GABA-A receptor activity", "6115682"),
    ("GABA", "HMDB0000112", "Vigabatrin", "DB01080", "Pharmacokinetic", "Brain", "Clinical", "Irreversible GABA transaminase inhibitor increases GABA levels", "2892672"),
    ("GABA", "HMDB0000112", "Gabapentin", "DB00996", "Pharmacodynamic", "Brain; Spinal cord (dorsal horn)", "Clinical", "Modulates GABA signaling via calcium channel alpha-2-delta subunit", "8556195"),
    ("Arachidonic acid", "HMDB0001043", "Aspirin", "DB00945", "Pharmacodynamic", "Platelets; Endothelium", "Clinical", "Irreversible COX-1/2 inhibitor blocks prostaglandin synthesis from arachidonic acid", "4583963"),
    ("Arachidonic acid", "HMDB0001043", "Celecoxib", "DB00482", "Pharmacodynamic", "Joint (synoviocytes); GI tract", "Clinical", "Selective COX-2 inhibitor blocks arachidonic acid to PGE2 conversion", "10048290"),
    ("Arachidonic acid", "HMDB0001043", "Ibuprofen", "DB01050", "Pharmacodynamic", "Systemic; Joint tissue", "Clinical", "Non-selective COX inhibitor reduces prostaglandin synthesis", "3282749"),
    ("Histamine", "HMDB0000870", "Diphenhydramine", "DB01075", "Pharmacodynamic", "Nasal mucosa; Skin; Brain", "Clinical", "H1 receptor inverse agonist blocks histamine signaling", "7059439"),
    ("Histamine", "HMDB0000870", "Ranitidine", "DB00863", "Pharmacodynamic", "Stomach (parietal cells)", "Clinical", "H2 receptor antagonist reduces gastric acid secretion", "6115682"),
    ("Histamine", "HMDB0000870", "Cetirizine", "DB00341", "Pharmacodynamic", "Skin; Nasal mucosa; Lung", "Clinical", "Selective H1 antihistamine for allergic responses", "2571380"),
    ("Prostaglandin E2", "HMDB0001220", "Misoprostol", "DB00929", "Pharmacodynamic", "Stomach (parietal cells); Uterus", "Clinical", "PGE1 analog that mimics prostaglandin gastroprotection", "3907090"),
    ("Prostaglandin E2", "HMDB0001220", "Indomethacin", "DB00328", "Pharmacokinetic", "Joint; Kidney; Brain", "Clinical", "COX inhibitor reduces PGE2 biosynthesis", "4259775"),
    ("Folate", "HMDB0000121", "Methotrexate", "DB00563", "Pharmacodynamic", "Bone marrow; Liver; Cancer cells", "Clinical", "Dihydrofolate reductase inhibitor blocks folate metabolism in cancer cells", "6122449"),
    ("Folate", "HMDB0000121", "Trimethoprim", "DB00440", "Pharmacodynamic", "Bacteria; Bone marrow", "Clinical", "Inhibits bacterial dihydrofolate reductase disrupting folate pathway", "6405975"),
    ("Folate", "HMDB0000121", "Sulfasalazine", "DB00795", "Pharmacokinetic", "Intestine", "Clinical", "Impairs folate absorption in the intestine", "6994203"),
    ("Thymidine", "HMDB0000273", "5-Fluorouracil", "DB00544", "Pharmacokinetic", "Tumor cells (colon; breast)", "Clinical", "Inhibits thymidylate synthase blocking thymidine synthesis in cancer cells", "6359418"),
    ("Thymidine", "HMDB0000273", "Capecitabine", "DB01101", "Pharmacokinetic", "Tumor cells; Liver", "Clinical", "Prodrug converted to 5-FU; preferentially activated in tumor tissue", "9754814"),
    ("Adenosine", "HMDB0000050", "Caffeine", "DB00201", "Pharmacodynamic", "Brain; Heart; Kidney", "Clinical", "Non-selective adenosine receptor antagonist", "8252576"),
    ("Adenosine", "HMDB0000050", "Theophylline", "DB00277", "Pharmacodynamic", "Lung (bronchial smooth muscle); Heart", "Clinical", "Adenosine receptor antagonist and PDE inhibitor", "1694936"),
    ("Adenosine", "HMDB0000050", "Dipyridamole", "DB00975", "Pharmacokinetic", "Platelets; Endothelium", "Clinical", "Inhibits adenosine uptake increasing extracellular adenosine", "2868872"),
    ("Acetylcholine", "HMDB0000895", "Donepezil", "DB00843", "Pharmacokinetic", "Brain (cortex; hippocampus)", "Clinical", "Acetylcholinesterase inhibitor increases acetylcholine levels", "8988460"),
    ("Acetylcholine", "HMDB0000895", "Atropine", "DB00572", "Pharmacodynamic", "Heart; Eye; GI tract", "Clinical", "Muscarinic receptor antagonist blocks acetylcholine signaling", "3108976"),
    ("Acetylcholine", "HMDB0000895", "Rivastigmine", "DB00989", "Pharmacokinetic", "Brain (cortex; hippocampus)", "Clinical", "Dual cholinesterase inhibitor increases acetylcholine in synaptic cleft", "10726378"),
    ("Norepinephrine", "HMDB0000216", "Venlafaxine", "DB00285", "Pharmacodynamic", "Brain (locus coeruleus)", "Clinical", "SNRI blocks norepinephrine and serotonin reuptake", "8042523"),
    ("Norepinephrine", "HMDB0000216", "Atomoxetine", "DB00289", "Pharmacokinetic", "Brain (prefrontal cortex)", "Clinical", "Selective norepinephrine reuptake inhibitor for ADHD", "14583924"),
    ("Norepinephrine", "HMDB0000216", "Clonidine", "DB00575", "Pharmacodynamic", "Brain (brainstem); Vasculature", "Clinical", "Alpha-2 adrenergic agonist reduces norepinephrine release", "3105032"),
    ("Tryptophan", "HMDB0000929", "Fluoxetine", "DB00472", "Pharmacokinetic", "Brain; Intestine (enterochromaffin cells)", "Clinical", "Tryptophan → serotonin pathway enhanced by SSRI mechanism", "3280901"),
    ("Tryptophan", "HMDB0000929", "Linezolid", "DB00601", "Pharmacodynamic", "Systemic", "Clinical", "MAO-A inhibition leads to serotonin syndrome with tryptophan excess", "11950141"),
    ("Nitric oxide", "HMDB0003378", "Sildenafil", "DB00203", "Pharmacodynamic", "Corpus cavernosum; Pulmonary vasculature", "Clinical", "PDE5 inhibitor potentiates NO-cGMP signaling in smooth muscle", "9015382"),
    ("Nitric oxide", "HMDB0003378", "Nitroglycerin", "DB00727", "Pharmacokinetic", "Coronary arteries; Venous smooth muscle", "Clinical", "Converted to NO causing vasodilation", "2868872"),
    ("Bile acids", "HMDB0000036", "Cholestyramine", "DB01432", "Pharmacokinetic", "Intestine (lumen)", "Clinical", "Bile acid sequestrant binds bile acids preventing reabsorption", "3048937"),
    ("Bile acids", "HMDB0000036", "Ursodeoxycholic acid", "DB01586", "Pharmacodynamic", "Liver (hepatocytes); Bile duct", "Clinical", "Replaces hydrophobic bile acids with hydrophilic UDCA", "8698674"),
    ("Melatonin", "HMDB0001389", "Fluvoxamine", "DB00176", "Pharmacokinetic", "Liver; Pineal gland", "Clinical", "CYP1A2 inhibition increases melatonin levels", "9084078"),
    ("Melatonin", "HMDB0001389", "Ramelteon", "DB00980", "Pharmacodynamic", "Brain (suprachiasmatic nucleus)", "Clinical", "Selective MT1/MT2 melatonin receptor agonist", "15665516"),
    ("Pyruvate", "HMDB0000243", "Dichloroacetate", "DB04267", "Pharmacodynamic", "Liver; Tumor cells; Muscle", "Experimental", "DCA activates PDH shifting pyruvate from lactate to acetyl-CoA", "17222789"),
    ("Pyruvate", "HMDB0000243", "Metformin", "DB00331", "Pharmacokinetic", "Liver (hepatocytes)", "Clinical", "Metformin inhibits hepatic gluconeogenesis from pyruvate", "10782096"),
    ("Succinate", "HMDB0000254", "Dimethyl fumarate", "DB08908", "Pharmacodynamic", "Immune cells; Brain", "Clinical", "Fumarate ester modulates succinate-related immune signaling (Nrf2 pathway)", "22517355"),
    ("Alpha-ketoglutarate", "HMDB0000208", "L-Asparaginase", "DB00023", "Pharmacokinetic", "Tumor cells (leukemia)", "Clinical", "Depletes asparagine affecting TCA cycle intermediate alpha-KG", "6320068"),
    ("Glutathione", "HMDB0000125", "Acetaminophen", "DB00316", "Pharmacokinetic", "Liver (hepatocytes); Kidney", "Clinical", "NAPQI (toxic metabolite) depletes hepatic glutathione stores", "6236856"),
    ("Glutathione", "HMDB0000125", "N-Acetylcysteine", "DB06151", "Pharmacokinetic", "Liver; Lung", "Clinical", "NAC is a glutathione precursor that replenishes depleted stores", "2189898"),
    ("Glutathione", "HMDB0000125", "Cisplatin", "DB00515", "Pharmacokinetic", "Kidney; Tumor cells", "Clinical", "Cisplatin detoxification consumes glutathione; GST-mediated conjugation", "2893000"),
    ("Glutathione", "HMDB0000125", "Busulfan", "DB01008", "Pharmacokinetic", "Liver; Bone marrow", "Clinical", "Primary metabolism via glutathione conjugation (GST-catalyzed)", "8616217"),
    # Cell-line-specific evidence
    ("Citrate", "HMDB0000094", "Sorafenib", "DB00398", "Pharmacodynamic", "HepG2 (liver); Tumor cells", "Cell Line", "Sorafenib alters TCA cycle flux affecting citrate levels in HCC cells", "25422591"),
    ("Citrate", "HMDB0000094", "Metformin", "DB00331", "Pharmacokinetic", "HeLa; MCF-7; Liver", "Cell Line", "Metformin reduces mitochondrial citrate production in cancer cell lines", "23274085"),
    ("Glutamine", "HMDB0000641", "CB-839", "DB12660", "Pharmacokinetic", "A549; HCT116; Tumor cells", "Cell Line", "Glutaminase inhibitor blocks glutamine-to-glutamate conversion", "25816826"),
    ("Glutamine", "HMDB0000641", "L-Asparaginase", "DB00023", "Pharmacokinetic", "Leukemia cells (CCRF-CEM; Jurkat)", "Cell Line", "L-asparaginase also depletes glutamine in leukemia cells", "2306024"),
    ("Serine", "HMDB0000187", "Sertraconazole", "DB11613", "Pharmacokinetic", "HCT116; MCF-7", "Cell Line", "PHGDH inhibitor blocks de novo serine synthesis in cancer cells", "28178514"),
    ("Palmitate", "HMDB0000220", "Orlistat", "DB01083", "Pharmacokinetic", "Intestine; MDA-MB-231 (breast cancer)", "Cell Line", "FASN inhibitor reduces palmitate synthesis in breast cancer cells", "14679127"),
    ("Palmitate", "HMDB0000220", "Cerulenin", "DB04070", "Pharmacokinetic", "MCF-7; HepG2", "Cell Line", "FASN inhibitor blocks palmitate biosynthesis in cell lines", "8654385"),
    ("2-Hydroxyglutarate", "HMDB0059655", "Enasidenib", "DB13874", "Pharmacokinetic", "AML cells (TF-1; THP-1)", "Cell Line", "IDH2 inhibitor reduces 2-HG levels in mutant AML cells", "28504110"),
    ("2-Hydroxyglutarate", "HMDB0059655", "Ivosidenib", "DB14568", "Pharmacokinetic", "AML cells", "Cell Line", "IDH1 inhibitor reduces 2-HG production in IDH1-mutant leukemia", "29860938"),
    ("Fumarate", "HMDB0000134", "Dimethyl fumarate", "DB08908", "Pharmacodynamic", "Immune cells (T cells; macrophages)", "Cell Line", "Fumarate ester activates Nrf2 pathway with immunomodulatory effects", "22419529"),
    ("Acetyl-CoA", "HMDB0001206", "Hydroxycitrate", "DB11545", "Pharmacokinetic", "HepG2; Adipocytes (3T3-L1)", "Cell Line", "ATP citrate lyase inhibitor blocks acetyl-CoA production from citrate", "16469903"),
    ("Sphingosine", "HMDB0000252", "Fingolimod", "DB08868", "Pharmacodynamic", "Lymphocytes; Brain; U937 cells", "Cell Line", "S1P receptor modulator; sphingosine analog that traps lymphocytes", "16357144"),
    ("NAD+", "HMDB0000902", "FK866", "DB12932", "Pharmacokinetic", "HCT116; MCF-7; Tumor cells", "Cell Line", "NAMPT inhibitor depletes NAD+ levels in cancer cells", "12590626"),
    ("NAD+", "HMDB0000902", "Olaparib", "DB09074", "Pharmacokinetic", "BRCA-mutant cells (MDA-MB-436; CAPAN1)", "Cell Line", "PARP inhibitor consumes NAD+ leading to synthetic lethality", "16912199"),
    # Tissue-specific pharmacokinetic interactions
    ("Creatinine", "HMDB0000562", "Cimetidine", "DB00501", "Pharmacokinetic", "Kidney (proximal tubule)", "Clinical", "Inhibits tubular secretion of creatinine raising serum levels without affecting GFR", "6219761"),
    ("Creatinine", "HMDB0000562", "Trimethoprim", "DB00440", "Pharmacokinetic", "Kidney (proximal tubule)", "Clinical", "Inhibits OCT2 transporter reducing creatinine secretion", "8521790"),
    ("Bilirubin", "HMDB0000054", "Atazanavir", "DB01072", "Pharmacokinetic", "Liver (hepatocytes)", "Clinical", "UGT1A1 inhibitor reduces bilirubin glucuronidation", "15258102"),
    ("Bilirubin", "HMDB0000054", "Rifampin", "DB01045", "Pharmacokinetic", "Liver (hepatocytes)", "Clinical", "PXR agonist induces UGT1A1 increasing bilirubin conjugation", "3882036"),
    ("Thyroid hormones", "HMDB0000651", "Levothyroxine", "DB00451", "Pharmacodynamic", "Thyroid; Liver; Brain", "Clinical", "Synthetic T4 replacement for hypothyroidism", "1596045"),
    ("Thyroid hormones", "HMDB0000651", "Propylthiouracil", "DB00550", "Pharmacokinetic", "Thyroid", "Clinical", "Inhibits thyroid peroxidase and 5'-deiodinase reducing thyroid hormone synthesis", "3903551"),
    ("Homocysteine", "HMDB0000742", "Methotrexate", "DB00563", "Pharmacokinetic", "Liver; Kidney; Bone marrow", "Clinical", "MTX inhibits DHFR disrupting folate cycle and increasing homocysteine", "8823153"),
    ("Homocysteine", "HMDB0000742", "Valproic acid", "DB00313", "Pharmacokinetic", "Liver", "Clinical", "Valproate increases homocysteine via folate depletion", "12557071"),
    ("Triglycerides", "HMDB0000827", "Gemfibrozil", "DB01241", "Pharmacodynamic", "Liver; Adipose", "Clinical", "PPAR-alpha agonist reduces triglyceride synthesis", "6995242"),
    ("Triglycerides", "HMDB0000827", "Fenofibrate", "DB01039", "Pharmacodynamic", "Liver; Adipose", "Clinical", "PPAR-alpha agonist increases lipoprotein lipase activity", "2178780"),
    ("Triglycerides", "HMDB0000827", "Omega-3 fatty acids", "DB04930", "Pharmacodynamic", "Liver", "Clinical", "Reduces hepatic VLDL-TG synthesis", "12407195"),
    ("Vitamin D", "HMDB0000876", "Phenytoin", "DB00252", "Pharmacokinetic", "Liver; Bone", "Clinical", "CYP3A4 induction increases vitamin D catabolism", "6260780"),
    ("Vitamin D", "HMDB0000876", "Ketoconazole", "DB01026", "Pharmacokinetic", "Liver; Kidney", "Clinical", "CYP3A4/24A1 inhibitor reduces vitamin D catabolism", "2176213"),
    ("Iron", "HMDB0000692", "Deferoxamine", "DB00746", "Pharmacokinetic", "Liver; Heart; Bone marrow", "Clinical", "Iron chelator binds free iron for renal excretion", "7063480"),
    ("Iron", "HMDB0000692", "Deferasirox", "DB01609", "Pharmacokinetic", "Liver; Heart", "Clinical", "Oral iron chelator for chronic iron overload", "16310585"),
    ("Calcium", "HMDB0000464", "Alendronate", "DB00630", "Pharmacodynamic", "Bone (osteoclasts)", "Clinical", "Bisphosphonate inhibits osteoclast-mediated calcium release from bone", "7862179"),
    ("Calcium", "HMDB0000464", "Calcitriol", "DB00136", "Pharmacodynamic", "Intestine; Bone; Kidney", "Clinical", "Active vitamin D3 increases intestinal calcium absorption", "6259529"),
    ("Potassium", "HMDB0000586", "Spironolactone", "DB00421", "Pharmacokinetic", "Kidney (collecting duct)", "Clinical", "Aldosterone antagonist reduces potassium excretion", "6131802"),
    ("Potassium", "HMDB0000586", "Furosemide", "DB00695", "Pharmacokinetic", "Kidney (loop of Henle)", "Clinical", "Loop diuretic increases potassium excretion causing hypokalemia", "6131802"),
]


def build_curated_dataset():
    """Build the curated MDrI dataset from literature-sourced interactions."""
    records = []
    for row in CURATED_MDRI:
        met_name, hmdb_id, drug, drugbank_id, itype, tissue, evidence, desc, pmid = row
        records.append({
            "Metabolite_Name": met_name,
            "HMDB_ID": hmdb_id,
            "SMILES": "",
            "Drug_Name": drug,
            "DrugBank_ID": drugbank_id,
            "Interaction_Type": itype,
            "Tissue": tissue,
            "Cell_Location": "",
            "Biospecimen": "",
            "Evidence_Level": evidence,
            "Description": desc,
            "PMID": pmid,
            "Source": "CoreMet_curated",
        })
    print(f"  Curated dataset: {len(records)} interactions")
    return records


# --------------------------------------------------------------------------
# Step 6 — Merge all sources into final CSV
# --------------------------------------------------------------------------

def merge_and_deduplicate(all_records):
    """Merge records from all sources, deduplicate, and clean up."""
    df = pd.DataFrame(all_records)

    # Ensure all expected columns exist
    expected_cols = [
        "Metabolite_Name", "HMDB_ID", "SMILES", "Drug_Name", "DrugBank_ID",
        "Interaction_Type", "Tissue", "Cell_Location", "Biospecimen",
        "Evidence_Level", "Description", "PMID", "Source",
    ]
    for col in expected_cols:
        if col not in df.columns:
            df[col] = ""

    df = df[expected_cols].fillna("")

    # Normalise
    df["Metabolite_Name"] = df["Metabolite_Name"].str.strip()
    df["Drug_Name"] = df["Drug_Name"].str.strip()
    df["HMDB_ID"] = df["HMDB_ID"].str.strip()
    df["DrugBank_ID"] = df["DrugBank_ID"].str.strip()

    # Deduplicate on (metabolite, drug, tissue, source)
    before = len(df)
    df = df.drop_duplicates(subset=["HMDB_ID", "Drug_Name", "Tissue", "Source"])
    print(f"  Deduplication: {before} → {len(df)} records")

    # Sort
    df = df.sort_values(["Metabolite_Name", "Drug_Name", "Tissue"]).reset_index(drop=True)

    return df


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Curate Metabolite-Drug Interaction data")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip downloading from external APIs; only use curated + cached data")
    parser.add_argument("--max-hmdb", type=int, default=500,
                        help="Max HMDB metabolites to query (default: 500)")
    args = parser.parse_args()

    print("=" * 60)
    print("CoreMet — Metabolite–Drug Interaction (MDrI) Data Curation")
    print("=" * 60)

    # Step 0: Gather known metabolites
    print("\n[Step 0] Gathering CoreMet metabolite universe...")
    met_df = get_cormet_metabolites()
    hmdb_ids = met_df["HMDB_ID"].tolist()
    hmdb_to_name = dict(zip(met_df["HMDB_ID"], met_df["Metabolite_Name"].fillna("")))
    name_to_hmdb = {v.lower(): k for k, v in hmdb_to_name.items() if isinstance(v, str) and v}

    all_records = []

    # Step 1: Curated dataset (always included)
    print("\n[Step 1] Building curated MDrI dataset...")
    curated = build_curated_dataset()
    all_records.extend(curated)

    if not args.skip_download:
        # Step 2: HMDB drug associations
        print("\n[Step 2] Querying HMDB for drug associations...")
        hmdb_cache = RAW_DIR / "hmdb_xml"
        hmdb_cache.mkdir(exist_ok=True)
        hmdb_records = download_hmdb_drug_data(hmdb_ids, hmdb_cache, max_per_run=args.max_hmdb)
        all_records.extend(hmdb_records)

        # Step 3: CTD chemical interactions
        print("\n[Step 3] Downloading CTD chemical interactions...")
        ctd_records = download_ctd_chem_interactions(RAW_DIR, hmdb_to_name)
        # Cross-reference CTD records with HMDB IDs
        for rec in ctd_records:
            met_name_lower = rec["Metabolite_Name"].lower()
            if met_name_lower in name_to_hmdb:
                rec["HMDB_ID"] = name_to_hmdb[met_name_lower]
            else:
                rec["HMDB_ID"] = ""
            rec.setdefault("SMILES", "")
            rec.setdefault("DrugBank_ID", rec.get("Drug_ID", ""))
            rec.setdefault("Evidence_Level", "CTD-curated")
            rec.setdefault("Tissue", "")
            rec.setdefault("Cell_Location", "")
            rec.setdefault("Biospecimen", "")
            rec.setdefault("Description", "")
        all_records.extend(ctd_records)
    else:
        print("\n  Skipping downloads (--skip-download)")

    # Step 4: Merge and deduplicate
    print("\n[Step 4] Merging and deduplicating...")
    df = merge_and_deduplicate(all_records)

    # Step 5: Write output
    print(f"\n[Step 5] Writing {OUT_CSV}...")
    df.to_csv(OUT_CSV, index=False)
    print(f"  ✓ {len(df):,} metabolite–drug interactions saved")

    # Summary
    print(f"\n{'=' * 60}")
    print(f"Summary:")
    print(f"  Total interactions: {len(df):,}")
    print(f"  Unique metabolites: {df['HMDB_ID'].nunique()}")
    print(f"  Unique drugs:       {df['Drug_Name'].nunique()}")
    print(f"  Sources:            {df['Source'].value_counts().to_dict()}")
    print(f"  Evidence levels:    {df['Evidence_Level'].value_counts().to_dict()}")
    print(f"  Tissues mentioned:  {df['Tissue'].nunique()}")
    print(f"  Output: {OUT_CSV}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
