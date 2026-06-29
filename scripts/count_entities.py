#!/usr/bin/env python3
"""Count unique entities across all CoreMet databases."""
import sys
sys.path.insert(0, "/Users/cheng.wang/Documents/mpi-web/MPI_web_server/mpi-vgae-web")

from app.config import Config
import pandas as pd

cfg = Config()
mpi = pd.read_csv(cfg.MPI_DB_PATH, dtype=str).fillna("")

from app.services.mdi_service import get_mdi_db
from app.services.mmi_service import get_mmi_db
from app.services.mdri_service import get_mdri_db
from app.services.mgi_service import get_mgi_db
from app.services.mgwas_service import get_mgwas_db

mdi = get_mdi_db()
mmi = get_mmi_db()
mdri = get_mdri_db()
mgi = get_mgi_db()
mgwas = get_mgwas_db()

# Metabolites - collect HMDB_ID + name pairs
met_map = {}  # hmdb_id -> name
for col_id, col_name, df in [
    ("HMDB ID", "Metabolite Name", mpi),
    ("HMDB_ID", "Metabolite_Name", mdi),
    ("HMDB_ID", "Metabolite_Name", mmi),
    ("HMDB_ID", "Metabolite_Name", mdri),
    ("HMDB_ID", "Metabolite_Name", mgi),
]:
    if col_id in df.columns:
        for _, row in df[[col_id, col_name]].drop_duplicates().iterrows():
            hid = str(row.get(col_id, "")).strip()
            name = str(row.get(col_name, "")).strip()
            if hid and hid.startswith("HMDB"):
                if hid not in met_map or (name and not met_map[hid]):
                    met_map[hid] = name

# Also metabolites from mGWAS
if "Metabolite" in mgwas.columns:
    for name in mgwas["Metabolite"].dropna().unique():
        name = str(name).strip()
        if name:
            key = f"MGWAS_{name}"
            if key not in met_map:
                met_map[key] = name

print(f"Metabolites: {len(met_map)}")

# Genes
gene_set = set()
for col, df in [("Gene Name", mpi), ("Gene_Symbol", mgi), ("Mapped_Gene", mgwas)]:
    if col in df.columns:
        gene_set.update(v.strip() for v in df[col].dropna().unique() if v and str(v).strip())
print(f"Genes: {len(gene_set)}")

# Proteins
prot_map = {}  # uniprot_id -> name
if "Uniprot ID" in mpi.columns and "Protein Name" in mpi.columns:
    for _, row in mpi[["Uniprot ID", "Protein Name"]].drop_duplicates().iterrows():
        uid = str(row.get("Uniprot ID", "")).strip()
        name = str(row.get("Protein Name", "")).strip()
        if uid:
            prot_map[uid] = name
print(f"Proteins: {len(prot_map)}")

# Diseases
disease_set = set()
if "Disease_Name" in mdi.columns:
    disease_set.update(v.strip() for v in mdi["Disease_Name"].dropna().unique() if v and str(v).strip())
print(f"Diseases: {len(disease_set)}")

# Microbes
microbe_set = set()
if "Microbe_Name" in mmi.columns:
    microbe_set.update(v.strip() for v in mmi["Microbe_Name"].dropna().unique() if v and str(v).strip())
print(f"Microbes: {len(microbe_set)}")

# Drugs
drug_map = {}
if "Drug_Name" in mdri.columns:
    for _, row in mdri[["Drug_Name"]].drop_duplicates().iterrows():
        name = str(row.get("Drug_Name", "")).strip()
        if name:
            drug_map[name] = ""
    if "DrugBank_ID" in mdri.columns:
        for _, row in mdri[["Drug_Name", "DrugBank_ID"]].drop_duplicates().iterrows():
            name = str(row.get("Drug_Name", "")).strip()
            dbid = str(row.get("DrugBank_ID", "")).strip()
            if name and dbid:
                drug_map[name] = dbid
print(f"Drugs: {len(drug_map)}")

# SNPs
snp_set = set()
if "rsID" in mgwas.columns:
    snp_set.update(v.strip() for v in mgwas["rsID"].dropna().unique() if v and str(v).strip())
print(f"SNPs: {len(snp_set)}")

total = len(met_map) + len(gene_set) + len(prot_map) + len(disease_set) + len(microbe_set) + len(drug_map) + len(snp_set)
print(f"TOTAL: {total}")
