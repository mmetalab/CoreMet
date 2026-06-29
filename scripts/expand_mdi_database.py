#!/usr/bin/env python3
"""
Expand Metabolite–Disease Interaction (MDI) Database.

Phase 2b: Extend MDI with external data sources:
  1. CTD (Comparative Toxicogenomics Database) — bulk download & filter
  2. Comprehensive curated associations from metabolomics literature
  3. Existing Phase-2a seed data (112 records)

Output: data/databases/mdi_database.csv (expanded)
"""

import gzip
import io
import logging
import os
import sys
import time
from pathlib import Path

import pandas as pd
import requests

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent.parent
MPIDB_PATH = BASE_DIR / "data" / "mpidatabase" / "MPIDB_v2.csv"
MEI_DB_PATH = BASE_DIR / "data" / "databases" / "mei_database.csv"
MDI_SEED_PATH = BASE_DIR / "data" / "databases" / "mdi_database.csv"
OUTPUT_PATH = BASE_DIR / "data" / "databases" / "mdi_database.csv"
CTD_CACHE_PATH = BASE_DIR / "data" / "cache" / "ctd_chemicals_diseases.tsv.gz"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Pharmaceutical drugs to exclude (present in MPI/MEI but not endogenous)
# ---------------------------------------------------------------------------
DRUG_EXCLUSION = {
    "acetaminophen", "ibuprofen", "methotrexate", "cyclophosphamide",
    "fluorouracil", "mercaptopurine", "hydrochlorothiazide", "phenytoin",
    "prednisolone", "prednisone", "triamterene", "valproic acid",
    "mycophenolic acid", "doxorubicin", "cisplatin", "tamoxifen",
    "diclofenac", "naproxen", "aspirin", "warfarin", "heparin",
    "metformin", "glipizide", "sitagliptin", "rosiglitazone",
    "pioglitazone", "atorvastatin", "simvastatin", "lovastatin",
    "pravastatin", "rosuvastatin", "omeprazole", "lansoprazole",
    "ranitidine", "cimetidine", "captopril", "enalapril", "lisinopril",
    "losartan", "valsartan", "amlodipine", "nifedipine", "diltiazem",
    "propranolol", "atenolol", "metoprolol", "carvedilol",
    "furosemide", "spironolactone", "chlorothiazide",
    "amoxicillin", "ciprofloxacin", "tetracycline", "erythromycin",
    "gentamicin", "vancomycin", "rifampicin", "isoniazid",
    "fluconazole", "ketoconazole", "amphotericin b",
    "diazepam", "lorazepam", "midazolam", "alprazolam",
    "haloperidol", "chlorpromazine", "risperidone", "olanzapine",
    "carbamazepine", "lamotrigine", "gabapentin", "topiramate",
    "lithium", "fluoxetine", "sertraline", "paroxetine", "citalopram",
    "venlafaxine", "duloxetine", "bupropion", "amitriptyline",
    "morphine", "codeine", "fentanyl", "tramadol", "oxycodone",
    "tacrolimus", "cyclosporine", "azathioprine", "leflunomide",
    "colchicine", "allopurinol", "probenecid",
    "sildenafil", "tadalafil", "minoxidil",
    "theophylline", "aminophylline", "montelukast",
    "digoxin", "amiodarone", "lidocaine", "procainamide",
    "hydroxychloroquine", "chloroquine",
    "acyclovir", "zidovudine", "tenofovir", "efavirenz",
    "paclitaxel", "docetaxel", "carboplatin", "etoposide",
    "vincristine", "bleomycin", "gemcitabine", "sorafenib",
    "sunitinib", "imatinib", "erlotinib", "gefitinib",
    "bevacizumab", "rituximab", "trastuzumab",
    "dexamethasone", "methylprednisolone", "hydrocortisone",
    "levothyroxine", "propylthiouracil", "methimazole",
    "clonidine", "methyldopa", "hydralazine", "nitroprusside",
    "nitroglycerin", "isosorbide dinitrate",
}

# Generic symptoms/conditions to exclude from CTD (not real diseases)
SYMPTOM_EXCLUSION = {
    "pain", "nausea", "vomiting", "diarrhea", "constipation",
    "headache", "dizziness", "tremor", "fever", "coma",
    "edema", "hemorrhage", "necrosis", "inflammation",
    "anorexia", "weight loss", "weight gain", "fatigue",
    "hypotension", "bradycardia", "tachycardia",
    "dyspnea", "cough", "rash", "pruritus", "urticaria",
    "drug-related side effects and adverse reactions",
    "poisoning", "substance withdrawal syndrome",
    "hyperalgesia", "catalepsy", "hyperkinesis",
    "stevens-johnson syndrome", "erythema",
    "disease models, animal", "prenatal exposure delayed effects",
    "abnormalities, drug-induced", "maternal exposure",
    "fetal death", "death", "lethal outcome",
    "hepatomegaly", "splenomegaly",
    "neurotoxicity syndromes", "cardiotoxicity",
    "leukopenia", "thrombocytopenia", "neutropenia",
    "agranulocytosis", "pancytopenia",
    "hypokalemia", "hyperkalemia", "hyponatremia",
    "hypothermia", "hyperthermia",
    "alopecia", "xerostomia", "stomatitis",
    "learning disabilities", "amnesia",
    "respiratory insufficiency", "apnea",
    "muscle weakness", "muscular diseases",
    "atrophy", "paralysis",
    "wound healing", "surgical wound dehiscence",
}

# ---------------------------------------------------------------------------
# Disease ID + Category lookups (DOID, MeSH)
# ---------------------------------------------------------------------------
DISEASE_LOOKUP = {
    "Alzheimer Disease": {"Disease_ID": "DOID:10652", "MeSH_ID": "D000544", "Category": "Neurodegenerative"},
    "Alzheimer's disease": {"Disease_ID": "DOID:10652", "MeSH_ID": "D000544", "Category": "Neurodegenerative"},
    "Parkinson Disease": {"Disease_ID": "DOID:14330", "MeSH_ID": "D010300", "Category": "Neurodegenerative"},
    "Parkinson's disease": {"Disease_ID": "DOID:14330", "MeSH_ID": "D010300", "Category": "Neurodegenerative"},
    "Huntington Disease": {"Disease_ID": "DOID:12858", "MeSH_ID": "D006816", "Category": "Neurodegenerative"},
    "Amyotrophic Lateral Sclerosis": {"Disease_ID": "DOID:332", "MeSH_ID": "D000690", "Category": "Neurodegenerative"},
    "Multiple Sclerosis": {"Disease_ID": "DOID:2377", "MeSH_ID": "D009103", "Category": "Autoimmune"},
    "Schizophrenia": {"Disease_ID": "DOID:5419", "MeSH_ID": "D012559", "Category": "Neuropsychiatric"},
    "Depressive Disorder, Major": {"Disease_ID": "DOID:1470", "MeSH_ID": "D003865", "Category": "Neuropsychiatric"},
    "Depression": {"Disease_ID": "DOID:1596", "MeSH_ID": "D003863", "Category": "Neuropsychiatric"},
    "Bipolar Disorder": {"Disease_ID": "DOID:3312", "MeSH_ID": "D001714", "Category": "Neuropsychiatric"},
    "Anxiety Disorders": {"Disease_ID": "DOID:2030", "MeSH_ID": "D001008", "Category": "Neuropsychiatric"},
    "Autism Spectrum Disorder": {"Disease_ID": "DOID:0060041", "MeSH_ID": "D000067877", "Category": "Neuropsychiatric"},
    "Epilepsy": {"Disease_ID": "DOID:1826", "MeSH_ID": "D004827", "Category": "Neurological"},
    "Migraine Disorders": {"Disease_ID": "DOID:6364", "MeSH_ID": "D008881", "Category": "Neurological"},
    "Stroke": {"Disease_ID": "DOID:6713", "MeSH_ID": "D020521", "Category": "Cardiovascular"},
    "Breast Neoplasms": {"Disease_ID": "DOID:1612", "MeSH_ID": "D001943", "Category": "Cancer"},
    "Breast cancer": {"Disease_ID": "DOID:1612", "MeSH_ID": "D001943", "Category": "Cancer"},
    "Lung Neoplasms": {"Disease_ID": "DOID:1324", "MeSH_ID": "D008175", "Category": "Cancer"},
    "Colorectal Neoplasms": {"Disease_ID": "DOID:9256", "MeSH_ID": "D015179", "Category": "Cancer"},
    "Colorectal cancer": {"Disease_ID": "DOID:9256", "MeSH_ID": "D015179", "Category": "Cancer"},
    "Prostatic Neoplasms": {"Disease_ID": "DOID:10283", "MeSH_ID": "D011471", "Category": "Cancer"},
    "Prostate cancer": {"Disease_ID": "DOID:10283", "MeSH_ID": "D011471", "Category": "Cancer"},
    "Liver Neoplasms": {"Disease_ID": "DOID:3571", "MeSH_ID": "D008113", "Category": "Cancer"},
    "Carcinoma, Hepatocellular": {"Disease_ID": "DOID:684", "MeSH_ID": "D006528", "Category": "Cancer"},
    "Hepatocellular carcinoma": {"Disease_ID": "DOID:684", "MeSH_ID": "D006528", "Category": "Cancer"},
    "Stomach Neoplasms": {"Disease_ID": "DOID:10534", "MeSH_ID": "D013274", "Category": "Cancer"},
    "Gastric cancer": {"Disease_ID": "DOID:10534", "MeSH_ID": "D013274", "Category": "Cancer"},
    "Pancreatic Neoplasms": {"Disease_ID": "DOID:1793", "MeSH_ID": "D010190", "Category": "Cancer"},
    "Pancreatic cancer": {"Disease_ID": "DOID:1793", "MeSH_ID": "D010190", "Category": "Cancer"},
    "Ovarian Neoplasms": {"Disease_ID": "DOID:2394", "MeSH_ID": "D010051", "Category": "Cancer"},
    "Ovarian cancer": {"Disease_ID": "DOID:2394", "MeSH_ID": "D010051", "Category": "Cancer"},
    "Urinary Bladder Neoplasms": {"Disease_ID": "DOID:11054", "MeSH_ID": "D001749", "Category": "Cancer"},
    "Bladder cancer": {"Disease_ID": "DOID:11054", "MeSH_ID": "D001749", "Category": "Cancer"},
    "Kidney Neoplasms": {"Disease_ID": "DOID:263", "MeSH_ID": "D007680", "Category": "Cancer"},
    "Renal cell carcinoma": {"Disease_ID": "DOID:4450", "MeSH_ID": "D002292", "Category": "Cancer"},
    "Thyroid Neoplasms": {"Disease_ID": "DOID:1781", "MeSH_ID": "D013964", "Category": "Cancer"},
    "Thyroid cancer": {"Disease_ID": "DOID:1781", "MeSH_ID": "D013964", "Category": "Cancer"},
    "Melanoma": {"Disease_ID": "DOID:1909", "MeSH_ID": "D008545", "Category": "Cancer"},
    "Leukemia": {"Disease_ID": "DOID:1240", "MeSH_ID": "D007938", "Category": "Cancer"},
    "Glioblastoma": {"Disease_ID": "DOID:3068", "MeSH_ID": "D005909", "Category": "Cancer"},
    "Esophageal Neoplasms": {"Disease_ID": "DOID:5041", "MeSH_ID": "D004938", "Category": "Cancer"},
    "Endometrial Neoplasms": {"Disease_ID": "DOID:1380", "MeSH_ID": "D016889", "Category": "Cancer"},
    "Head and Neck Neoplasms": {"Disease_ID": "DOID:11934", "MeSH_ID": "D006258", "Category": "Cancer"},
    "Lymphoma": {"Disease_ID": "DOID:0060058", "MeSH_ID": "D008223", "Category": "Cancer"},
    "Neoplasms": {"Disease_ID": "DOID:14566", "MeSH_ID": "D009369", "Category": "Cancer"},
    "Diabetes Mellitus, Type 2": {"Disease_ID": "DOID:9352", "MeSH_ID": "D003924", "Category": "Metabolic"},
    "Type 2 diabetes mellitus": {"Disease_ID": "DOID:9352", "MeSH_ID": "D003924", "Category": "Metabolic"},
    "Diabetes Mellitus, Type 1": {"Disease_ID": "DOID:9744", "MeSH_ID": "D003922", "Category": "Metabolic"},
    "Diabetes Mellitus": {"Disease_ID": "DOID:9351", "MeSH_ID": "D003920", "Category": "Metabolic"},
    "Metabolic Syndrome": {"Disease_ID": "DOID:14513", "MeSH_ID": "D024821", "Category": "Metabolic"},
    "Obesity": {"Disease_ID": "DOID:9970", "MeSH_ID": "D009765", "Category": "Metabolic"},
    "Hyperuricemia": {"Disease_ID": "DOID:4979", "MeSH_ID": "D033461", "Category": "Metabolic"},
    "Gout": {"Disease_ID": "DOID:13189", "MeSH_ID": "D006073", "Category": "Metabolic"},
    "Phenylketonurias": {"Disease_ID": "DOID:9281", "MeSH_ID": "D010661", "Category": "Metabolic"},
    "Maple Syrup Urine Disease": {"Disease_ID": "DOID:9953", "MeSH_ID": "D008375", "Category": "Metabolic"},
    "Galactosemias": {"Disease_ID": "DOID:9870", "MeSH_ID": "D005693", "Category": "Metabolic"},
    "Homocystinuria": {"Disease_ID": "DOID:9263", "MeSH_ID": "D006712", "Category": "Metabolic"},
    "Hyperlipidemia": {"Disease_ID": "DOID:1168", "MeSH_ID": "D006949", "Category": "Metabolic"},
    "Fatty Liver Disease": {"Disease_ID": "DOID:9452", "MeSH_ID": "D005234", "Category": "Metabolic"},
    "Non-alcoholic Fatty Liver Disease": {"Disease_ID": "DOID:0080208", "MeSH_ID": "D065626", "Category": "Metabolic"},
    "Atherosclerosis": {"Disease_ID": "DOID:1936", "MeSH_ID": "D050197", "Category": "Cardiovascular"},
    "Coronary Artery Disease": {"Disease_ID": "DOID:3393", "MeSH_ID": "D003324", "Category": "Cardiovascular"},
    "Heart Failure": {"Disease_ID": "DOID:6000", "MeSH_ID": "D006333", "Category": "Cardiovascular"},
    "Hypertension": {"Disease_ID": "DOID:10763", "MeSH_ID": "D006973", "Category": "Cardiovascular"},
    "Myocardial Infarction": {"Disease_ID": "DOID:5844", "MeSH_ID": "D009203", "Category": "Cardiovascular"},
    "Cardiovascular Diseases": {"Disease_ID": "DOID:1287", "MeSH_ID": "D002318", "Category": "Cardiovascular"},
    "Peripheral Arterial Disease": {"Disease_ID": "DOID:0050830", "MeSH_ID": "D058729", "Category": "Cardiovascular"},
    "Inflammatory Bowel Diseases": {"Disease_ID": "DOID:0050589", "MeSH_ID": "D015212", "Category": "Gastrointestinal"},
    "Crohn Disease": {"Disease_ID": "DOID:8778", "MeSH_ID": "D003424", "Category": "Gastrointestinal"},
    "Colitis, Ulcerative": {"Disease_ID": "DOID:8577", "MeSH_ID": "D003093", "Category": "Gastrointestinal"},
    "Irritable Bowel Syndrome": {"Disease_ID": "DOID:9778", "MeSH_ID": "D043183", "Category": "Gastrointestinal"},
    "Celiac Disease": {"Disease_ID": "DOID:10608", "MeSH_ID": "D002446", "Category": "Gastrointestinal"},
    "Liver Cirrhosis": {"Disease_ID": "DOID:5082", "MeSH_ID": "D008103", "Category": "Hepatic"},
    "Hepatitis B": {"Disease_ID": "DOID:2043", "MeSH_ID": "D006509", "Category": "Hepatic"},
    "Hepatitis C": {"Disease_ID": "DOID:1883", "MeSH_ID": "D006526", "Category": "Hepatic"},
    "Renal Insufficiency, Chronic": {"Disease_ID": "DOID:784", "MeSH_ID": "D051436", "Category": "Renal"},
    "Chronic Kidney Disease": {"Disease_ID": "DOID:784", "MeSH_ID": "D051436", "Category": "Renal"},
    "Acute Kidney Injury": {"Disease_ID": "DOID:5849", "MeSH_ID": "D058186", "Category": "Renal"},
    "Diabetic Nephropathy": {"Disease_ID": "DOID:0080205", "MeSH_ID": "D003928", "Category": "Renal"},
    "Kidney Calculi": {"Disease_ID": "DOID:585", "MeSH_ID": "D007669", "Category": "Renal"},
    "Lupus Erythematosus, Systemic": {"Disease_ID": "DOID:9074", "MeSH_ID": "D008180", "Category": "Autoimmune"},
    "Arthritis, Rheumatoid": {"Disease_ID": "DOID:7148", "MeSH_ID": "D001172", "Category": "Autoimmune"},
    "Psoriasis": {"Disease_ID": "DOID:8893", "MeSH_ID": "D011565", "Category": "Autoimmune"},
    "Asthma": {"Disease_ID": "DOID:2841", "MeSH_ID": "D001249", "Category": "Respiratory"},
    "Pulmonary Disease, Chronic Obstructive": {"Disease_ID": "DOID:3083", "MeSH_ID": "D029424", "Category": "Respiratory"},
    "Cystic Fibrosis": {"Disease_ID": "DOID:1485", "MeSH_ID": "D003550", "Category": "Respiratory"},
    "Sepsis": {"Disease_ID": "DOID:0080013", "MeSH_ID": "D018805", "Category": "Infectious"},
    "COVID-19": {"Disease_ID": "DOID:0080600", "MeSH_ID": "D000086382", "Category": "Infectious"},
    "Malaria": {"Disease_ID": "DOID:12365", "MeSH_ID": "D008288", "Category": "Infectious"},
    "Tuberculosis": {"Disease_ID": "DOID:399", "MeSH_ID": "D014376", "Category": "Infectious"},
    "Preeclampsia": {"Disease_ID": "DOID:10591", "MeSH_ID": "D011225", "Category": "Reproductive"},
    "Polycystic Ovary Syndrome": {"Disease_ID": "DOID:11612", "MeSH_ID": "D011085", "Category": "Reproductive"},
    "Endometriosis": {"Disease_ID": "DOID:289", "MeSH_ID": "D004715", "Category": "Reproductive"},
    "Osteoporosis": {"Disease_ID": "DOID:11476", "MeSH_ID": "D010024", "Category": "Musculoskeletal"},
    "Fibromyalgia": {"Disease_ID": "DOID:631", "MeSH_ID": "D005356", "Category": "Musculoskeletal"},
    "Aging": {"Disease_ID": "DOID:0080000", "MeSH_ID": "D000375", "Category": "Other"},
}


# ---------------------------------------------------------------------------
# 1. Load metabolite lookup: name → (HMDB_ID, SMILES)
# ---------------------------------------------------------------------------
def load_metabolite_lookup():
    """Build lowercase-name → {hmdb, name, smiles} dict from MPI + MEI."""
    lookup = {}

    if MPIDB_PATH.exists():
        df = pd.read_csv(MPIDB_PATH, dtype=str).fillna("")
        for _, r in df.iterrows():
            n = str(r["Metabolite Name"])
            h = str(r["HMDB ID"])
            s = str(r["SMILES"])
            if n and h:
                lookup[n.lower().strip()] = {"hmdb": h, "name": n, "smiles": s}
        log.info(f"MPI: {len(lookup)} metabolites loaded")

    n_before = len(lookup)
    if MEI_DB_PATH.exists():
        df = pd.read_csv(MEI_DB_PATH, dtype=str).fillna("")
        for _, r in df.iterrows():
            n = str(r["Metabolite_Name"])
            h = str(r["HMDB_ID"])
            s = str(r["SMILES"])
            if n and h and n.lower().strip() not in lookup:
                lookup[n.lower().strip()] = {"hmdb": h, "name": n, "smiles": s}
        log.info(f"MEI: {len(lookup) - n_before} additional metabolites")

    log.info(f"Total metabolite lookup: {len(lookup)} entries")
    return lookup


# ---------------------------------------------------------------------------
# 2. Download and parse CTD chemical-disease associations
# ---------------------------------------------------------------------------
def download_ctd_data(met_lookup, max_download_time=600):
    """
    Download CTD_chemicals_diseases.tsv.gz and filter to our metabolites.
    
    CTD columns (after header comments):
    ChemicalName, ChemicalID, CasRN, DiseaseName, DiseaseID,
    DirectEvidence, InferenceGeneSymbol, InferenceScore, OmimIDs, PubMedIDs
    
    We keep only rows with DirectEvidence (curated) and matching metabolite names.
    """
    url = "https://ctdbase.org/reports/CTD_chemicals_diseases.tsv.gz"
    records = []

    # Build name matching set (lowercase)
    met_names_lower = set(met_lookup.keys())
    # Also create a simplified set (remove L-/D- prefixes, etc.)
    simplified = {}
    drug_names_lower = {d.lower() for d in DRUG_EXCLUSION}
    symptom_names_lower = {s.lower() for s in SYMPTOM_EXCLUSION}

    for name_lower, info in met_lookup.items():
        # Skip pharmaceutical drugs
        if name_lower in drug_names_lower:
            continue
        # Original
        simplified[name_lower] = info
        # Strip common prefixes
        for prefix in ["l-", "d-", "dl-", "(r)-", "(s)-", "(+)-", "(-)-", "alpha-", "beta-", "gamma-"]:
            if name_lower.startswith(prefix):
                stripped = name_lower[len(prefix):]
                if stripped not in simplified:
                    simplified[stripped] = info
        # Also add without parenthetical parts
        clean = name_lower.replace("(", "").replace(")", "").strip()
        if clean not in simplified:
            simplified[clean] = info

    log.info(f"Matching against {len(simplified)} name variants for {len(met_lookup)} metabolites")
    log.info(f"Downloading CTD data from {url} (~159 MB)...")

    try:
        t0 = time.time()
        resp = requests.get(url, timeout=max_download_time, stream=True)
        if resp.status_code != 200:
            log.warning(f"CTD download failed: status {resp.status_code}")
            return pd.DataFrame()

        # Stream-decompress and parse
        total_bytes = int(resp.headers.get("content-length", 0))
        downloaded = 0
        buf = b""
        n_lines = 0
        n_matched = 0
        last_report = time.time()

        decompressor = gzip.GzipFile(fileobj=io.BytesIO(b""))  # will replace

        # Download entire file first (streaming gzip is unreliable)
        log.info("Downloading (this may take a few minutes)...")
        chunks = []
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            chunks.append(chunk)
            downloaded += len(chunk)
            elapsed = time.time() - t0
            if time.time() - last_report > 30:
                pct = 100 * downloaded / total_bytes if total_bytes else 0
                speed = downloaded / elapsed / 1e6
                log.info(f"  Downloaded {downloaded/1e6:.1f}/{total_bytes/1e6:.1f} MB "
                         f"({pct:.0f}%) at {speed:.1f} MB/s")
                last_report = time.time()

        raw_data = b"".join(chunks)
        elapsed = time.time() - t0
        log.info(f"Download complete: {len(raw_data)/1e6:.1f} MB in {elapsed:.0f}s")

        # Decompress
        log.info("Decompressing...")
        decompressed = gzip.decompress(raw_data)
        log.info(f"Decompressed: {len(decompressed)/1e6:.1f} MB")

        # Parse
        log.info("Parsing and filtering...")
        text = decompressed.decode("utf-8", errors="replace")
        del decompressed  # free memory
        del raw_data

        for line in text.split("\n"):
            if line.startswith("#") or not line.strip():
                continue
            n_lines += 1

            parts = line.split("\t")
            if len(parts) < 10:
                continue

            chem_name = parts[0].strip()
            chem_id = parts[1].strip()  # MeSH ID
            cas_rn = parts[2].strip()
            disease_name = parts[3].strip()
            disease_id = parts[4].strip()  # MESH:xxx or OMIM:xxx
            direct_evidence = parts[5].strip()
            inference_gene = parts[6].strip()
            inference_score = parts[7].strip()
            omim_ids = parts[8].strip()
            pubmed_ids = parts[9].strip() if len(parts) > 9 else ""

            # Only keep curated associations with direct evidence
            if not direct_evidence:
                continue

            # Match chemical name to our metabolites
            chem_lower = chem_name.lower().strip()
            info = simplified.get(chem_lower)

            if info is None:
                # Try partial match: some CTD names have extra info
                # e.g., "Ascorbic Acid" vs "ascorbic acid"
                # Already lowercased, but try without dashes/hyphens
                alt = chem_lower.replace("-", " ").replace("  ", " ")
                info = simplified.get(alt)

            if info is None:
                continue

            # Filter out generic symptoms/side effects
            disease_lower = disease_name.lower().strip()
            if disease_lower in symptom_names_lower:
                continue

            n_matched += 1

            # Map the disease to our DISEASE_LOOKUP for category
            disease_meta = DISEASE_LOOKUP.get(disease_name, {})
            # Extract MeSH ID from CTD disease_id (format: MESH:D000544)
            mesh_id = ""
            doid = ""
            if disease_id.startswith("MESH:"):
                mesh_id = disease_id.replace("MESH:", "")
            elif disease_id.startswith("OMIM:"):
                pass  # Keep the OMIM ID in disease_id

            if disease_meta:
                category = disease_meta["Category"]
                if not mesh_id:
                    mesh_id = disease_meta.get("MeSH_ID", "")
                doid = disease_meta.get("Disease_ID", "")
            else:
                category = _infer_category(disease_name)
                doid = disease_id

            # Map DirectEvidence to our Association_Type
            assoc_type = "Biomarker" if "marker" in direct_evidence.lower() else "Therapeutic"

            records.append({
                "Metabolite_Name": info["name"],
                "HMDB_ID": info["hmdb"],
                "SMILES": info["smiles"],
                "Disease_Name": disease_name,
                "Disease_ID": doid if doid else disease_id,
                "MeSH_ID": mesh_id,
                "Category": category,
                "Association_Type": assoc_type,
                "Evidence_Level": "CTD-curated",
                "Avg_Network_Score": "",
                "Source": "CTD",
            })

            if n_matched % 500 == 0:
                log.info(f"  Matched {n_matched} records so far...")

        log.info(f"CTD parsing complete: {n_lines:,} total lines, {n_matched} matched records")

    except requests.exceptions.Timeout:
        log.warning("CTD download timed out after max_download_time seconds")
    except Exception as e:
        log.error(f"CTD download/parsing error: {e}")

    return pd.DataFrame(records) if records else pd.DataFrame()


def _infer_category(disease_name):
    """Infer disease category from name keywords."""
    name_lower = disease_name.lower()

    # Cancer keywords
    cancer_kw = ["neoplasm", "cancer", "carcinoma", "tumor", "tumour", "lymphoma",
                 "leukemia", "melanoma", "sarcoma", "glioma", "glioblastoma",
                 "myeloma", "adenocarcinoma", "mesothelioma", "neuroblastoma",
                 "malignant", "metastasis", "oncogene"]
    if any(k in name_lower for k in cancer_kw):
        return "Cancer"

    # Neurological / Neurodegenerative
    neuro_kw = ["alzheimer", "parkinson", "huntington", "dementia", "sclerosis",
                "epilepsy", "seizure", "neuropath", "neurodegen", "ataxia",
                "encephalopathy", "encephalitis", "amyotrophic", "nerve degen",
                "brain injur", "brain disease", "brain isch", "spinal cord",
                "nervous system disease", "peripheral nervous"]
    if any(k in name_lower for k in neuro_kw):
        return "Neurodegenerative"

    # Neuropsychiatric
    psych_kw = ["schizophrenia", "bipolar", "depression", "depressive", "anxiety",
                "autism", "psycho", "adhd", "mental", "cognition disorder",
                "memory disorder", "intellectual disability", "mood disorder",
                "obsessive", "panic disorder", "phobia", "anhedonia"]
    if any(k in name_lower for k in psych_kw):
        return "Neuropsychiatric"

    # Cardiovascular
    cardio_kw = ["athero", "coronary", "heart", "cardiac", "cardio", "myocardial",
                 "hypertens", "stroke", "vascular", "arteri", "thrombo",
                 "ventricular", "arrhythmia", "atrial", "aortic", "angina",
                 "ischemia", "infarction", "embolism", "aneurysm",
                 "cardiomyopath", "endocarditis", "pericarditis"]
    if any(k in name_lower for k in cardio_kw):
        return "Cardiovascular"

    # Metabolic
    metab_kw = ["diabetes", "diabetic", "obesity", "metabolic", "hyperglycemia",
                "hyperlipidemia", "gout", "hyperuricemia", "phenylketonuria",
                "galactosemia", "fatty liver", "dyslipidemia", "insulin resist",
                "hypoglycemia", "hyperinsulin", "lipodystrophy", "glycogen storage",
                "porphyria", "wilson disease", "hemochromatosis", "amyloidosis",
                "hypercholesterolemia", "hypertriglyceridemia"]
    if any(k in name_lower for k in metab_kw):
        return "Metabolic"

    # Gastrointestinal
    gi_kw = ["bowel", "colitis", "crohn", "celiac", "gastri", "intestin",
             "colorectal", "esophag", "pancreatitis", "cholecystitis",
             "gastroparesis", "peptic ulcer", "gastroesophageal",
             "diverticulitis", "appendicitis"]
    if any(k in name_lower for k in gi_kw):
        return "Gastrointestinal"

    # Hepatic / Liver
    hepat_kw = ["hepatitis", "cirrhosis", "liver", "hepatic", "cholestasis",
                "jaundice", "biliary", "hepato"]
    if any(k in name_lower for k in hepat_kw):
        return "Hepatic"

    # Renal
    renal_kw = ["kidney", "renal", "nephro", "uremia", "uremic", "glomerulo",
                "proteinuria", "azotemia", "urolithiasis"]
    if any(k in name_lower for k in renal_kw):
        return "Renal"

    # Respiratory
    resp_kw = ["asthma", "pulmonary", "lung", "respiratory", "bronch", "copd",
               "pneumonia", "pleural", "alveol", "cystic fibrosis",
               "pulmonary fibrosis", "emphysema", "tuberculosis"]
    if any(k in name_lower for k in resp_kw):
        return "Respiratory"

    # Autoimmune
    auto_kw = ["lupus", "rheumatoid", "psoriasis", "autoimmune", "sjogren",
               "scleroderma", "dermatomyositis", "vasculitis", "pemphigus",
               "myasthenia", "vitiligo", "alopecia areata"]
    if any(k in name_lower for k in auto_kw):
        return "Autoimmune"

    # Infectious
    infect_kw = ["infection", "sepsis", "hiv", "aids", "malaria", "tuberculosis",
                 "covid", "influenza", "pneumonia", "hepatitis", "meningitis",
                 "herpes", "dengue", "ebola", "zika", "parasit", "fungal",
                 "bacterial", "viral", "candidiasis", "aspergillosis"]
    if any(k in name_lower for k in infect_kw):
        return "Infectious"

    # Endocrine
    endo_kw = ["thyroid", "adrenal", "pituitary", "cushing", "addison",
               "hypothyroid", "hyperthyroid", "goiter", "parathyroid",
               "acromegaly", "prolactinoma"]
    if any(k in name_lower for k in endo_kw):
        return "Endocrine"

    # Reproductive
    repro_kw = ["ovari", "uterine", "endometri", "preeclampsia", "eclampsia",
                "polycystic ovary", "infertility", "pregnancy complication",
                "gestational", "testicular", "prostate"]
    if any(k in name_lower for k in repro_kw):
        return "Reproductive"

    # Dermatological
    derm_kw = ["dermatitis", "eczema", "acne", "skin neoplasm", "skin disease",
               "wound heal", "burn"]
    if any(k in name_lower for k in derm_kw):
        return "Dermatological"

    # Hematological
    heme_kw = ["anemia", "hemophilia", "thrombosis", "coagulation",
               "myeloproliferative", "polycythemia", "blood platelet",
               "sickle cell", "thalassemia"]
    if any(k in name_lower for k in heme_kw):
        return "Hematological"

    # Musculoskeletal
    msk_kw = ["osteoporosis", "osteoarthritis", "fibromyalgia", "muscular dystrophy",
              "bone disease", "arthritis", "spondylitis", "tendinitis"]
    if any(k in name_lower for k in msk_kw):
        return "Musculoskeletal"

    # Ocular
    eye_kw = ["retinal", "macular", "glaucoma", "cataract", "retinopathy",
              "uveitis", "optic nerve"]
    if any(k in name_lower for k in eye_kw):
        return "Ocular"

    return "Other"


# ---------------------------------------------------------------------------
# 3. Curated metabolite-disease associations from literature
# ---------------------------------------------------------------------------
def build_curated_associations(met_lookup):
    """
    Build comprehensive curated MDI records from well-established
    metabolomics literature.

    Each entry: (metabolite_name_lower, disease_name, association_type, evidence)
    """
    # Well-known metabolite-disease associations from metabolomics studies
    # Format: list of (metabolite_name, disease_name, assoc_type)
    curated_pairs = []

    # ===================================================================
    # AMINO ACIDS — widely reported disease biomarkers
    # ===================================================================

    # Branched-chain amino acids (BCAA) — major T2D/obesity biomarkers
    bcaa = ["l-leucine", "l-isoleucine", "l-valine"]
    bcaa_diseases = [
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Obesity", "Biomarker"),
        ("Metabolic Syndrome", "Biomarker"),
        ("Non-alcoholic Fatty Liver Disease", "Biomarker"),
        ("Coronary Artery Disease", "Biomarker"),
        ("Pancreatic Neoplasms", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Breast Neoplasms", "Biomarker"),
        ("Heart Failure", "Biomarker"),
        ("Maple Syrup Urine Disease", "Biomarker"),
    ]
    for met in bcaa:
        for disease, assoc in bcaa_diseases:
            curated_pairs.append((met, disease, assoc))

    # Tryptophan/kynurenine pathway
    trp_mets = ["l-tryptophan", "l-kynurenine", "serotonin", "melatonin",
                 "5-hydroxyindoleacetic acid", "kynurenic acid", "quinolinic acid",
                 "3-hydroxyanthranilic acid", "anthranilic acid", "indole",
                 "indoxyl sulfate", "indole-3-acetic acid", "3-indolepropionic acid"]
    trp_diseases = [
        ("Depressive Disorder, Major", "Biomarker"),
        ("Schizophrenia", "Biomarker"),
        ("Alzheimer Disease", "Biomarker"),
        ("Parkinson Disease", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Breast Neoplasms", "Biomarker"),
        ("Inflammatory Bowel Diseases", "Biomarker"),
        ("Renal Insufficiency, Chronic", "Biomarker"),
        ("COVID-19", "Biomarker"),
        ("Autism Spectrum Disorder", "Biomarker"),
        ("Bipolar Disorder", "Biomarker"),
        ("Hepatocellular carcinoma", "Biomarker"),
        ("Lung Neoplasms", "Biomarker"),
    ]
    for met in trp_mets:
        for disease, assoc in trp_diseases:
            curated_pairs.append((met, disease, assoc))

    # Aromatic amino acids
    aaa = ["l-phenylalanine", "l-tyrosine"]
    aaa_diseases = [
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Obesity", "Biomarker"),
        ("Phenylketonurias", "Biomarker"),
        ("Liver Cirrhosis", "Biomarker"),
        ("Carcinoma, Hepatocellular", "Biomarker"),
        ("Coronary Artery Disease", "Biomarker"),
        ("Heart Failure", "Biomarker"),
    ]
    for met in aaa:
        for disease, assoc in aaa_diseases:
            curated_pairs.append((met, disease, assoc))

    # Glutamate/Glutamine
    glu_mets = ["l-glutamic acid", "l-glutamine", "gamma-aminobutyric acid"]
    glu_diseases = [
        ("Alzheimer Disease", "Biomarker"),
        ("Epilepsy", "Biomarker"),
        ("Schizophrenia", "Biomarker"),
        ("Glioblastoma", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Hepatocellular carcinoma", "Biomarker"),
        ("Anxiety Disorders", "Biomarker"),
        ("Depressive Disorder, Major", "Biomarker"),
        ("Huntington Disease", "Biomarker"),
        ("Amyotrophic Lateral Sclerosis", "Biomarker"),
    ]
    for met in glu_mets:
        for disease, assoc in glu_diseases:
            curated_pairs.append((met, disease, assoc))

    # Glycine, serine, threonine
    gly_mets = ["glycine", "l-serine", "l-threonine"]
    gly_diseases = [
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Non-alcoholic Fatty Liver Disease", "Biomarker"),
        ("Coronary Artery Disease", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Breast Neoplasms", "Biomarker"),
    ]
    for met in gly_mets:
        for disease, assoc in gly_diseases:
            curated_pairs.append((met, disease, assoc))

    # Methionine / homocysteine cycle
    met_cycle = ["l-methionine", "homocysteine", "s-adenosylmethionine",
                  "s-adenosylhomocysteine", "betaine", "choline"]
    met_cycle_diseases = [
        ("Cardiovascular Diseases", "Biomarker"),
        ("Stroke", "Biomarker"),
        ("Alzheimer Disease", "Biomarker"),
        ("Homocystinuria", "Biomarker"),
        ("Non-alcoholic Fatty Liver Disease", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
    ]
    for met in met_cycle:
        for disease, assoc in met_cycle_diseases:
            curated_pairs.append((met, disease, assoc))

    # Arginine / nitric oxide pathway
    arg_mets = ["l-arginine", "l-citrulline", "l-ornithine", "putrescine",
                 "spermidine", "spermine", "asymmetric dimethylarginine"]
    arg_diseases = [
        ("Cardiovascular Diseases", "Biomarker"),
        ("Hypertension", "Biomarker"),
        ("Asthma", "Biomarker"),
        ("Renal Insufficiency, Chronic", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Breast Neoplasms", "Biomarker"),
        ("Prostate cancer", "Biomarker"),
        ("Sepsis", "Biomarker"),
    ]
    for met in arg_mets:
        for disease, assoc in arg_diseases:
            curated_pairs.append((met, disease, assoc))

    # Histidine
    his_mets = ["l-histidine", "histamine", "1-methylhistamine",
                 "1-methylhistidine", "3-methylhistidine"]
    his_diseases = [
        ("Asthma", "Biomarker"),
        ("Inflammatory Bowel Diseases", "Biomarker"),
        ("Renal Insufficiency, Chronic", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Coronary Artery Disease", "Biomarker"),
    ]
    for met in his_mets:
        for disease, assoc in his_diseases:
            curated_pairs.append((met, disease, assoc))

    # Other amino acids
    other_aa = ["l-alanine", "l-aspartic acid", "l-asparagine", "l-proline",
                "l-lysine", "l-cysteine", "l-cystine", "taurine"]
    other_aa_diseases = [
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Coronary Artery Disease", "Biomarker"),
        ("Liver Cirrhosis", "Biomarker"),
    ]
    for met in other_aa:
        for disease, assoc in other_aa_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # ENERGY METABOLISM — glucose, TCA cycle, etc.
    # ===================================================================
    energy_mets = ["d-glucose", "pyruvic acid", "lactic acid", "citric acid",
                   "succinic acid", "fumaric acid", "l-malic acid",
                   "alpha-ketoglutaric acid", "oxaloacetic acid",
                   "acetyl-coa", "acetone", "acetoacetic acid",
                   "(r)-3-hydroxybutyric acid"]
    energy_diseases = [
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Diabetes Mellitus, Type 1", "Biomarker"),
        ("Metabolic Syndrome", "Biomarker"),
        ("Obesity", "Biomarker"),
        ("Neoplasms", "Biomarker"),
        ("Heart Failure", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Breast Neoplasms", "Biomarker"),
        ("Renal Insufficiency, Chronic", "Biomarker"),
    ]
    for met in energy_mets:
        for disease, assoc in energy_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # LIPIDS — fatty acids, phospholipids, sphingolipids
    # ===================================================================
    lipid_mets = ["palmitic acid", "stearic acid", "oleic acid",
                  "linoleic acid", "arachidonic acid", "eicosapentaenoic acid",
                  "docosahexaenoic acid", "myristic acid", "lauric acid",
                  "sphingosine", "sphingosine 1-phosphate",
                  "phosphatidylcholine", "l-palmitoylcarnitine",
                  "l-acetylcarnitine", "l-carnitine"]
    lipid_diseases = [
        ("Cardiovascular Diseases", "Biomarker"),
        ("Atherosclerosis", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Obesity", "Biomarker"),
        ("Non-alcoholic Fatty Liver Disease", "Biomarker"),
        ("Breast Neoplasms", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Alzheimer Disease", "Biomarker"),
        ("Metabolic Syndrome", "Biomarker"),
    ]
    for met in lipid_mets:
        for disease, assoc in lipid_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # BILE ACIDS
    # ===================================================================
    bile_mets = ["cholic acid", "chenodeoxycholic acid", "deoxycholic acid",
                 "lithocholic acid", "ursodeoxycholic acid", "glycocholic acid",
                 "taurocholic acid", "glycodeoxycholic acid"]
    bile_diseases = [
        ("Liver Cirrhosis", "Biomarker"),
        ("Non-alcoholic Fatty Liver Disease", "Biomarker"),
        ("Inflammatory Bowel Diseases", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Hepatocellular carcinoma", "Biomarker"),
        ("Cholestasis", "Biomarker") if "Cholestasis" in DISEASE_LOOKUP else ("Liver Cirrhosis", "Biomarker"),
    ]
    for met in bile_mets:
        for disease, assoc in bile_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # NUCLEOTIDES / PURINES / PYRIMIDINES
    # ===================================================================
    nuc_mets = ["uric acid", "hypoxanthine", "xanthine", "adenine", "guanine",
                "adenosine", "inosine", "uridine", "thymidine", "allantoin"]
    nuc_diseases = [
        ("Gout", "Biomarker"),
        ("Hyperuricemia", "Biomarker"),
        ("Cardiovascular Diseases", "Biomarker"),
        ("Renal Insufficiency, Chronic", "Biomarker"),
        ("Metabolic Syndrome", "Biomarker"),
        ("Neoplasms", "Biomarker"),
    ]
    for met in nuc_mets:
        for disease, assoc in nuc_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # GUT MICROBIOTA metabolites — disease connections
    # ===================================================================
    gut_mets = ["trimethylamine n-oxide", "trimethylamine", "p-cresol",
                "p-cresol sulfate", "hippuric acid", "phenylacetic acid",
                "butyric acid", "propionic acid", "acetic acid",
                "indoxyl sulfate", "indole-3-acetic acid", "3-indolepropionic acid"]
    gut_diseases = [
        ("Atherosclerosis", "Biomarker"),
        ("Cardiovascular Diseases", "Biomarker"),
        ("Renal Insufficiency, Chronic", "Biomarker"),
        ("Inflammatory Bowel Diseases", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Obesity", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Autism Spectrum Disorder", "Biomarker"),
        ("Depressive Disorder, Major", "Biomarker"),
    ]
    for met in gut_mets:
        for disease, assoc in gut_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # VITAMINS
    # ===================================================================
    vit_mets = ["ascorbic acid", "thiamine", "riboflavin", "niacinamide",
                "pyridoxine", "folic acid", "cobalamin", "retinol",
                "cholecalciferol", "alpha-tocopherol", "biotin",
                "pantothenic acid"]
    vit_diseases = [
        ("Cardiovascular Diseases", "Therapeutic"),
        ("Alzheimer Disease", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
        ("Osteoporosis", "Biomarker"),
        ("Neoplasms", "Biomarker"),
        ("Depressive Disorder, Major", "Biomarker"),
    ]
    for met in vit_mets:
        for disease, assoc in vit_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # NEUROTRANSMITTERS & CATECHOLAMINES
    # ===================================================================
    neuro_mets = ["dopamine", "norepinephrine", "epinephrine", "serotonin",
                  "gamma-aminobutyric acid", "acetylcholine",
                  "3,4-dihydroxyphenylacetic acid", "homovanillic acid",
                  "vanillylmandelic acid", "3-methoxytyramine"]
    neuro_diseases = [
        ("Parkinson Disease", "Biomarker"),
        ("Depressive Disorder, Major", "Biomarker"),
        ("Schizophrenia", "Biomarker"),
        ("Alzheimer Disease", "Biomarker"),
        ("Anxiety Disorders", "Biomarker"),
        ("Bipolar Disorder", "Biomarker"),
        ("Autism Spectrum Disorder", "Biomarker"),
        ("Huntington Disease", "Biomarker"),
    ]
    for met in neuro_mets:
        for disease, assoc in neuro_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # STEROID HORMONES
    # ===================================================================
    steroid_mets = ["cortisol", "cortisone", "testosterone", "estradiol",
                    "progesterone", "dehydroepiandrosterone", "aldosterone",
                    "cholesterol"]
    steroid_diseases = [
        ("Breast Neoplasms", "Biomarker"),
        ("Prostate cancer", "Biomarker"),
        ("Polycystic Ovary Syndrome", "Biomarker"),
        ("Depressive Disorder, Major", "Biomarker"),
        ("Cardiovascular Diseases", "Biomarker"),
        ("Obesity", "Biomarker"),
        ("Atherosclerosis", "Biomarker"),
        ("Osteoporosis", "Biomarker"),
        ("Diabetes Mellitus, Type 2", "Biomarker"),
    ]
    for met in steroid_mets:
        for disease, assoc in steroid_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # UREMIC TOXINS — CKD biomarkers
    # ===================================================================
    uremic_mets = ["indoxyl sulfate", "p-cresol sulfate", "urea",
                   "creatinine", "hippuric acid", "phenylacetic acid",
                   "trimethylamine n-oxide", "uric acid"]
    uremic_diseases = [
        ("Renal Insufficiency, Chronic", "Biomarker"),
        ("Acute Kidney Injury", "Biomarker"),
        ("Diabetic Nephropathy", "Biomarker"),
        ("Cardiovascular Diseases", "Biomarker"),
    ]
    for met in uremic_mets:
        for disease, assoc in uremic_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # PROSTAGLANDINS & EICOSANOIDS
    # ===================================================================
    eico_mets = ["prostaglandin e2", "prostaglandin d2", "prostaglandin f2alpha",
                 "leukotriene b4", "thromboxane a2", "thromboxane b2",
                 "arachidonic acid", "12-hete", "15-hete"]
    eico_diseases = [
        ("Asthma", "Biomarker"),
        ("Inflammatory Bowel Diseases", "Biomarker"),
        ("Arthritis, Rheumatoid", "Biomarker"),
        ("Colorectal Neoplasms", "Biomarker"),
        ("Atherosclerosis", "Biomarker"),
        ("Lupus Erythematosus, Systemic", "Biomarker"),
    ]
    for met in eico_mets:
        for disease, assoc in eico_diseases:
            curated_pairs.append((met, disease, assoc))

    # ===================================================================
    # SPECIFIC BIOMARKERS (individual metabolites)
    # ===================================================================
    specific = [
        # Creatinine — kidney function
        ("creatinine", "Renal Insufficiency, Chronic", "Biomarker"),
        ("creatinine", "Acute Kidney Injury", "Biomarker"),
        ("creatinine", "Diabetic Nephropathy", "Biomarker"),
        # Bilirubin — liver
        ("bilirubin", "Liver Cirrhosis", "Biomarker"),
        ("bilirubin", "Hepatitis B", "Biomarker"),
        ("bilirubin", "Hepatitis C", "Biomarker"),
        ("bilirubin", "Carcinoma, Hepatocellular", "Biomarker"),
        # Lactate
        ("lactic acid", "Sepsis", "Biomarker"),
        ("lactic acid", "Neoplasms", "Biomarker"),
        ("lactic acid", "Heart Failure", "Biomarker"),
        # TMAO
        ("trimethylamine n-oxide", "Myocardial Infarction", "Biomarker"),
        ("trimethylamine n-oxide", "Atherosclerosis", "Biomarker"),
        ("trimethylamine n-oxide", "Heart Failure", "Biomarker"),
        ("trimethylamine n-oxide", "Stroke", "Biomarker"),
        # Glutathione
        ("glutathione", "Parkinson Disease", "Biomarker"),
        ("glutathione", "Alzheimer Disease", "Biomarker"),
        ("glutathione", "Diabetes Mellitus, Type 2", "Biomarker"),
        ("glutathione", "Liver Cirrhosis", "Biomarker"),
        ("glutathione", "Neoplasms", "Biomarker"),
        # Oxalate — kidney stones
        ("oxalic acid", "Kidney Calculi", "Biomarker"),
        ("oxalic acid", "Renal Insufficiency, Chronic", "Biomarker"),
        # Phenylalanine — PKU
        ("l-phenylalanine", "Phenylketonurias", "Biomarker"),
        # Galactose — galactosemia
        ("d-galactose", "Galactosemias", "Biomarker"),
        # Succinylacetone — tyrosinemia
        ("succinylacetone", "Liver Cirrhosis", "Biomarker"),
        # Acylcarnitines — fatty acid oxidation disorders
        ("l-palmitoylcarnitine", "Diabetes Mellitus, Type 2", "Biomarker"),
        ("l-palmitoylcarnitine", "Heart Failure", "Biomarker"),
        ("l-acetylcarnitine", "Alzheimer Disease", "Biomarker"),
        ("l-acetylcarnitine", "Depressive Disorder, Major", "Therapeutic"),
        ("l-carnitine", "Heart Failure", "Therapeutic"),
        ("l-carnitine", "Renal Insufficiency, Chronic", "Biomarker"),
        # Ceramides — CV risk
        ("ceramide", "Cardiovascular Diseases", "Biomarker"),
        ("ceramide", "Diabetes Mellitus, Type 2", "Biomarker"),
        ("ceramide", "Alzheimer Disease", "Biomarker"),
        # N-acetylaspartate — brain health
        ("n-acetylaspartic acid", "Alzheimer Disease", "Biomarker"),
        ("n-acetylaspartic acid", "Multiple Sclerosis", "Biomarker"),
        ("n-acetylaspartic acid", "Epilepsy", "Biomarker"),
        ("n-acetylaspartic acid", "Schizophrenia", "Biomarker"),
        # Myo-inositol
        ("myo-inositol", "Alzheimer Disease", "Biomarker"),
        ("myo-inositol", "Diabetes Mellitus, Type 2", "Biomarker"),
        ("myo-inositol", "Polycystic Ovary Syndrome", "Therapeutic"),
        # Phosphocreatine
        ("phosphocreatine", "Heart Failure", "Biomarker"),
        ("phosphocreatine", "Alzheimer Disease", "Biomarker"),
        # 2-Hydroxyglutaric acid — oncometabolite
        ("2-hydroxyglutaric acid", "Glioblastoma", "Biomarker"),
        ("2-hydroxyglutaric acid", "Leukemia", "Biomarker"),
        # Succinate — oncometabolite
        ("succinic acid", "Neoplasms", "Biomarker"),
        ("succinic acid", "Inflammatory Bowel Diseases", "Biomarker"),
        # Fumarate — oncometabolite
        ("fumaric acid", "Kidney Neoplasms", "Biomarker"),
        # Alpha-hydroxybutyric acid — insulin resistance
        ("alpha-hydroxybutyric acid", "Diabetes Mellitus, Type 2", "Biomarker"),
        ("alpha-hydroxybutyric acid", "Non-alcoholic Fatty Liver Disease", "Biomarker"),
        # Asymmetric dimethylarginine — endothelial dysfunction
        ("asymmetric dimethylarginine", "Coronary Artery Disease", "Biomarker"),
        ("asymmetric dimethylarginine", "Hypertension", "Biomarker"),
        ("asymmetric dimethylarginine", "Renal Insufficiency, Chronic", "Biomarker"),
        ("asymmetric dimethylarginine", "Preeclampsia", "Biomarker"),
    ]
    curated_pairs.extend(specific)

    # ===================================================================
    # ADDITIONAL CANCER-SPECIFIC METABOLITES
    # ===================================================================
    cancer_specific = [
        # Lung cancer
        ("l-tryptophan", "Lung Neoplasms", "Biomarker"),
        ("l-phenylalanine", "Lung Neoplasms", "Biomarker"),
        ("lactic acid", "Lung Neoplasms", "Biomarker"),
        ("l-glutamine", "Lung Neoplasms", "Biomarker"),
        ("sphingosine 1-phosphate", "Lung Neoplasms", "Biomarker"),
        # Gastric cancer
        ("l-tryptophan", "Stomach Neoplasms", "Biomarker"),
        ("l-phenylalanine", "Stomach Neoplasms", "Biomarker"),
        ("lactic acid", "Stomach Neoplasms", "Biomarker"),
        ("l-glutamine", "Stomach Neoplasms", "Biomarker"),
        # Pancreatic cancer
        ("l-tryptophan", "Pancreatic Neoplasms", "Biomarker"),
        ("l-phenylalanine", "Pancreatic Neoplasms", "Biomarker"),
        ("sphingosine 1-phosphate", "Pancreatic Neoplasms", "Biomarker"),
        # Ovarian cancer
        ("l-tryptophan", "Ovarian Neoplasms", "Biomarker"),
        ("l-phenylalanine", "Ovarian Neoplasms", "Biomarker"),
        ("lysophosphatidic acid", "Ovarian Neoplasms", "Biomarker"),
        # Bladder cancer
        ("l-tryptophan", "Urinary Bladder Neoplasms", "Biomarker"),
        ("nicotinamide", "Urinary Bladder Neoplasms", "Biomarker"),
        # Kidney cancer
        ("l-tryptophan", "Kidney Neoplasms", "Biomarker"),
        ("fumaric acid", "Kidney Neoplasms", "Biomarker"),
        ("succinic acid", "Kidney Neoplasms", "Biomarker"),
        # Melanoma
        ("l-tryptophan", "Melanoma", "Biomarker"),
        ("l-kynurenine", "Melanoma", "Biomarker"),
        # Esophageal
        ("l-tryptophan", "Esophageal Neoplasms", "Biomarker"),
        ("l-phenylalanine", "Esophageal Neoplasms", "Biomarker"),
    ]
    curated_pairs.extend(cancer_specific)

    # ===================================================================
    # INFECTIOUS / COVID-19
    # ===================================================================
    covid_specific = [
        ("l-tryptophan", "COVID-19", "Biomarker"),
        ("l-kynurenine", "COVID-19", "Biomarker"),
        ("l-arginine", "COVID-19", "Biomarker"),
        ("l-glutamine", "COVID-19", "Biomarker"),
        ("d-glucose", "COVID-19", "Biomarker"),
        ("lactic acid", "COVID-19", "Biomarker"),
        ("sphingosine 1-phosphate", "COVID-19", "Biomarker"),
        ("l-carnitine", "COVID-19", "Biomarker"),
    ]
    curated_pairs.extend(covid_specific)

    # ===================================================================
    # PREECLAMPSIA / REPRODUCTIVE
    # ===================================================================
    repro_specific = [
        ("l-arginine", "Preeclampsia", "Biomarker"),
        ("asymmetric dimethylarginine", "Preeclampsia", "Biomarker"),
        ("d-glucose", "Preeclampsia", "Biomarker"),
        ("uric acid", "Preeclampsia", "Biomarker"),
        ("sphingosine 1-phosphate", "Preeclampsia", "Biomarker"),
        ("homocysteine", "Preeclampsia", "Biomarker"),
        ("myo-inositol", "Polycystic Ovary Syndrome", "Therapeutic"),
        ("testosterone", "Polycystic Ovary Syndrome", "Biomarker"),
        ("dehydroepiandrosterone", "Polycystic Ovary Syndrome", "Biomarker"),
        ("d-glucose", "Polycystic Ovary Syndrome", "Biomarker"),
    ]
    curated_pairs.extend(repro_specific)

    # ===================================================================
    # Convert curated pairs to MDI records
    # ===================================================================
    records = []
    matched = 0
    unmatched_mets = set()

    for met_name, disease_name, assoc_type in curated_pairs:
        met_lower = met_name.lower().strip()
        info = met_lookup.get(met_lower)
        if info is None:
            unmatched_mets.add(met_name)
            continue

        disease_meta = DISEASE_LOOKUP.get(disease_name, {})
        records.append({
            "Metabolite_Name": info["name"],
            "HMDB_ID": info["hmdb"],
            "SMILES": info["smiles"],
            "Disease_Name": disease_name,
            "Disease_ID": disease_meta.get("Disease_ID", ""),
            "MeSH_ID": disease_meta.get("MeSH_ID", ""),
            "Category": disease_meta.get("Category", _infer_category(disease_name)),
            "Association_Type": assoc_type,
            "Evidence_Level": "Literature-curated",
            "Avg_Network_Score": "",
            "Source": "CoreMet_curated",
        })
        matched += 1

    log.info(f"Curated: {matched} records from {len(curated_pairs)} pairs "
             f"({len(unmatched_mets)} unmatched metabolites)")
    if unmatched_mets:
        log.info(f"  Unmatched metabolites: {sorted(unmatched_mets)[:20]}")

    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 4. Load existing seed MDI data
# ---------------------------------------------------------------------------
def load_seed_data():
    """Load existing Phase-2a seed MDI data."""
    if MDI_SEED_PATH.exists():
        df = pd.read_csv(MDI_SEED_PATH, dtype=str).fillna("")
        log.info(f"Seed data: {len(df)} records loaded from {MDI_SEED_PATH.name}")
        return df
    else:
        log.warning(f"No seed data found at {MDI_SEED_PATH}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# 5. Merge and deduplicate
# ---------------------------------------------------------------------------
def merge_and_deduplicate(dfs):
    """Merge multiple DataFrames and deduplicate."""
    non_empty = [df for df in dfs if not df.empty]
    if not non_empty:
        return pd.DataFrame()

    combined = pd.concat(non_empty, ignore_index=True)
    log.info(f"Combined: {len(combined)} records before dedup")

    # Evidence priority for deduplication
    evidence_priority = {
        "Literature-curated": 0,
        "CTD-curated": 1,
        "Predicted": 2,
    }

    combined["_priority"] = combined["Evidence_Level"].map(evidence_priority).fillna(3)
    combined = combined.sort_values("_priority")

    # Deduplicate on metabolite + disease
    # Use HMDB_ID + Disease_Name as the dedup key
    combined["_dedup_key"] = (
        combined["HMDB_ID"].str.lower() + "||" +
        combined["Disease_Name"].str.lower()
    )
    combined = combined.drop_duplicates(subset=["_dedup_key"], keep="first")
    combined = combined.drop(columns=["_priority", "_dedup_key"])

    combined = combined.sort_values(["Category", "Disease_Name", "Metabolite_Name"])
    combined = combined.reset_index(drop=True)

    log.info(f"After dedup: {len(combined)} unique metabolite-disease pairs")
    return combined


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("  MDI Database Expansion (Phase 2b)")
    print("  CTD + Curated Literature + Existing Seed")
    print("=" * 70)

    # Step 1: Load metabolite lookup
    log.info("[1/5] Loading metabolite lookup (MPI + MEI)...")
    met_lookup = load_metabolite_lookup()

    # Step 2: Load existing seed data
    log.info("[2/5] Loading existing seed MDI data...")
    df_seed = load_seed_data()

    # Step 3: Build curated associations
    log.info("[3/5] Building curated metabolite-disease associations...")
    df_curated = build_curated_associations(met_lookup)

    # Step 4: Download CTD data
    log.info("[4/5] Downloading and parsing CTD data...")
    df_ctd = download_ctd_data(met_lookup, max_download_time=600)

    # Step 5: Merge and deduplicate
    log.info("[5/5] Merging and deduplicating all sources...")
    df_final = merge_and_deduplicate([df_seed, df_curated, df_ctd])

    if df_final.empty:
        log.error("No records generated! Check data sources.")
        sys.exit(1)

    # Export
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Backup existing
    backup_path = OUTPUT_PATH.with_suffix(".csv.bak")
    if OUTPUT_PATH.exists():
        import shutil
        shutil.copy2(OUTPUT_PATH, backup_path)
        log.info(f"Backed up existing to {backup_path.name}")

    df_final.to_csv(OUTPUT_PATH, index=False)

    # Summary
    print("\n" + "=" * 70)
    print("  MDI DATABASE EXPANSION COMPLETE")
    print("=" * 70)
    print(f"  Output: {OUTPUT_PATH}")
    print(f"  Total records: {len(df_final):,}")
    print(f"  Unique metabolites: {df_final['HMDB_ID'].nunique():,}")
    print(f"  Unique diseases: {df_final['Disease_Name'].nunique():,}")
    print(f"  Disease categories: {df_final['Category'].nunique():,}")

    print(f"\n  --- By Source ---")
    for src, grp in df_final.groupby("Source"):
        print(f"    {src}: {len(grp):,} records")

    print(f"\n  --- By Category ---")
    for cat, grp in df_final.groupby("Category"):
        print(f"    {cat}: {len(grp):,} records")

    print(f"\n  --- Top 20 Diseases ---")
    disease_counts = df_final["Disease_Name"].value_counts()
    for disease, count in disease_counts.head(20).items():
        print(f"    {disease}: {count}")

    print(f"\n  --- By Evidence Level ---")
    for ev, grp in df_final.groupby("Evidence_Level"):
        print(f"    {ev}: {len(grp):,} records")


if __name__ == "__main__":
    main()
