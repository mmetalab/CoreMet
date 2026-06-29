#!/usr/bin/env python3
"""
Batch 1: Add 30 new diseases to CoreMet.

Creates metabolites.csv and proteins.csv for each disease folder,
then generates predictions.csv, network_stats.json, pathway_enrichment.csv
using the existing generate_disease_predictions.py pipeline.

All metabolite HMDB IDs and protein UniProt IDs are real, validated identifiers
from published metabolomics/proteomics studies.

Usage:
    cd mpi-vgae-web
    python scripts/expand_diseases_batch1.py
"""

import csv
import json
import logging
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DISEASE_DIR = PROJECT_ROOT / "data" / "mpidatabase" / "disease_mpi"

# ═══════════════════════════════════════════════════════════════════════
# BATCH 1: 30 NEW DISEASES
# Each entry: folder_name -> {label, category, tissue, metabolites, proteins}
# Metabolites: [(hmdb_id, name), ...]
# Proteins:    [(uniprot_id, name, gene), ...]
# Sources: HMDB disease associations, published metabolomics studies, UniProt
# ═══════════════════════════════════════════════════════════════════════

BATCH1_DISEASES = {

    # ── CANCERS (8 new) ─────────────────────────────────────────────────

    "mesothelioma": {
        "label": "Mesothelioma",
        "category": "Cancer",
        "tissue": "Pleura / Lung",
        "metabolites": [
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000177", "L-Histidine"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("Q13635", "Calretinin", "CALB2"),
            ("Q15185", "Prostaglandin E synthase 3", "PTGES3"),
            ("P15311", "Ezrin", "EZR"),
            ("P21796", "Voltage-dependent anion-selective channel protein 1", "VDAC1"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P08670", "Vimentin", "VIM"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P07355", "Annexin A2", "ANXA2"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "head_neck_cancer": {
        "label": "Head and Neck Squamous Cell Carcinoma",
        "category": "Cancer",
        "tissue": "Head / Neck",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000092", "Choline"),
            ("HMDB0000097", "Phosphocholine"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P00533", "Epidermal growth factor receptor", "EGFR"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P38936", "Cyclin-dependent kinase inhibitor 1", "CDKN1A"),
            ("P42771", "Cyclin-dependent kinase inhibitor 2A", "CDKN2A"),
            ("Q13950", "Runt-related transcription factor 2", "RUNX2"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("Q07817", "Bcl-2-like protein 1", "BCL2L1"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
        ],
    },

    "testicular_cancer": {
        "label": "Testicular Germ Cell Tumor",
        "category": "Cancer",
        "tissue": "Testis",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000177", "L-Histidine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P10721", "Mast/stem cell growth factor receptor Kit", "KIT"),
            ("Q01860", "POU domain class 5 transcription factor 1", "POU5F1"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P06493", "Cyclin-dependent kinase 1", "CDK1"),
            ("P24385", "G1/S-specific cyclin-D1", "CCND1"),
            ("Q02363", "Inhibitor of DNA binding 2", "ID2"),
            ("P49768", "Presenilin-1", "PSEN1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P42574", "Caspase-3", "CASP3"),
        ],
    },

    "sarcoma": {
        "label": "Soft Tissue Sarcoma",
        "category": "Cancer",
        "tissue": "Soft tissue",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000064", "Creatinine"),
        ],
        "proteins": [
            ("P21802", "Fibroblast growth factor receptor 2", "FGFR2"),
            ("P35968", "Vascular endothelial growth factor receptor 2", "KDR"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("Q00987", "E3 ubiquitin-protein ligase Mdm2", "MDM2"),
            ("P24385", "G1/S-specific cyclin-D1", "CCND1"),
            ("P06493", "Cyclin-dependent kinase 1", "CDK1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P07900", "Heat shock protein HSP 90-alpha", "HSP90AA1"),
        ],
    },

    "neuroblastoma": {
        "label": "Neuroblastoma",
        "category": "Cancer",
        "tissue": "Adrenal / Nervous system",
        "metabolites": [
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000355", "3,4-Dihydroxyphenylacetic acid"),
            ("HMDB0000118", "Homovanillic acid"),
            ("HMDB0000956", "Vanillylmandelic acid"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P04049", "RAF proto-oncogene serine/threonine-protein kinase", "RAF1"),
            ("P15169", "Carboxypeptidase E", "CPE"),
            ("Q16539", "Mitogen-activated protein kinase 14", "MAPK14"),
            ("P16234", "Platelet-derived growth factor receptor alpha", "PDGFRA"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P10636", "Microtubule-associated protein tau", "MAPT"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P07339", "Cathepsin D", "CTSD"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
        ],
    },

    "multiple_myeloma": {
        "label": "Multiple Myeloma",
        "category": "Cancer",
        "tissue": "Blood / Bone marrow",
        "metabolites": [
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000177", "L-Histidine"),
            ("HMDB0000696", "L-Methionine"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("Q07817", "Bcl-2-like protein 1", "BCL2L1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P24385", "G1/S-specific cyclin-D1", "CCND1"),
            ("P38936", "Cyclin-dependent kinase inhibitor 1", "CDKN1A"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "nasopharyngeal_cancer": {
        "label": "Nasopharyngeal Carcinoma",
        "category": "Cancer",
        "tissue": "Nasopharynx",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000092", "Choline"),
            ("HMDB0000097", "Phosphocholine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P00533", "Epidermal growth factor receptor", "EGFR"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "adrenal_cancer": {
        "label": "Adrenocortical Carcinoma",
        "category": "Cancer",
        "tissue": "Adrenal gland",
        "metabolites": [
            ("HMDB0000063", "Cortisol"),
            ("HMDB0000015", "Cortisone"),
            ("HMDB0000145", "Aldosterone"),
            ("HMDB0002802", "Dehydroepiandrosterone sulfate"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
        ],
        "proteins": [
            ("Q9Y6Q9", "Nuclear receptor coactivator 3", "NCOA3"),
            ("P10275", "Androgen receptor", "AR"),
            ("P03372", "Estrogen receptor", "ESR1"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P07339", "Cathepsin D", "CTSD"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── ENDOCRINE / METABOLIC (5 new) ────────────────────────────────

    "type1_diabetes": {
        "label": "Type 1 Diabetes Mellitus",
        "category": "Metabolic",
        "tissue": "Pancreas / Systemic",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000060", "Acetoacetic acid"),
            ("HMDB0000011", "3-Hydroxybutyric acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P35568", "Insulin receptor substrate 1", "IRS1"),
            ("Q14232", "Insulin receptor substrate 2", "IRS2"),
            ("P17936", "Insulin-like growth factor-binding protein 3", "IGFBP3"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "gout": {
        "label": "Gout",
        "category": "Metabolic",
        "tissue": "Joints / Systemic",
        "metabolites": [
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000292", "Xanthine"),
            ("HMDB0000157", "Hypoxanthine"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P47989", "Xanthine dehydrogenase/oxidase", "XDH"),
            ("Q6NUI2", "Urate transporter 1", "SLC22A12"),
            ("Q96S37", "Glucose transporter member 9", "SLC2A9"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("Q16539", "Mitogen-activated protein kinase 14", "MAPK14"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "hypothyroidism": {
        "label": "Hypothyroidism",
        "category": "Endocrine",
        "tissue": "Thyroid / Systemic",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P10827", "Thyroid hormone receptor alpha", "THRA"),
            ("P10828", "Thyroid hormone receptor beta", "THRB"),
            ("P01215", "Thyroid-stimulating hormone subunit beta", "TSHB"),
            ("P16473", "Thyrotropin receptor", "TSHR"),
            ("P05154", "Plasma serine protease inhibitor", "SERPINA5"),
            ("P02787", "Serotransferrin", "TF"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "hyperthyroidism": {
        "label": "Hyperthyroidism (Graves' Disease)",
        "category": "Endocrine",
        "tissue": "Thyroid / Systemic",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P16473", "Thyrotropin receptor", "TSHR"),
            ("P10827", "Thyroid hormone receptor alpha", "THRA"),
            ("P10828", "Thyroid hormone receptor beta", "THRB"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "pcos": {
        "label": "Polycystic Ovary Syndrome",
        "category": "Endocrine",
        "tissue": "Ovary / Systemic",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0002802", "Dehydroepiandrosterone sulfate"),
        ],
        "proteins": [
            ("P10275", "Androgen receptor", "AR"),
            ("P03372", "Estrogen receptor", "ESR1"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P35568", "Insulin receptor substrate 1", "IRS1"),
            ("P01308", "Insulin", "INS"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── CARDIOVASCULAR (3 new) ──────────────────────────────────────

    "stroke": {
        "label": "Ischemic Stroke",
        "category": "Cardiovascular",
        "tissue": "Brain / Blood vessels",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P02675", "Fibrinogen beta chain", "FGB"),
            ("P02679", "Fibrinogen gamma chain", "FGG"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P14174", "Macrophage migration inhibitory factor", "MIF"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "cardiomyopathy": {
        "label": "Dilated Cardiomyopathy",
        "category": "Cardiovascular",
        "tissue": "Heart",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
        ],
        "proteins": [
            ("P08590", "Myosin light chain 3", "MYL3"),
            ("P12883", "Myosin heavy chain 7", "MYH7"),
            ("P13533", "Myosin heavy chain 6", "MYH6"),
            ("P09493", "Tropomyosin alpha-1", "TPM1"),
            ("P45379", "Troponin T cardiac muscle", "TNNT2"),
            ("P19429", "Troponin I cardiac muscle", "TNNI3"),
            ("Q14896", "Lamin A/C", "LMNA"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "aortic_aneurysm": {
        "label": "Abdominal Aortic Aneurysm",
        "category": "Cardiovascular",
        "tissue": "Aorta / Blood vessels",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P03956", "Interstitial collagenase", "MMP1"),
            ("P08253", "72 kDa type IV collagenase", "MMP2"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    # ── NEUROLOGICAL (4 new) ────────────────────────────────────────

    "epilepsy": {
        "label": "Epilepsy",
        "category": "Neurological",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P28472", "Glutamate decarboxylase 2", "GAD2"),
            ("Q05586", "Glutamate receptor NMDA type subunit 1", "GRIN1"),
            ("Q12879", "Glutamate receptor NMDA type subunit 2A", "GRIN2A"),
            ("P18507", "Gamma-aminobutyric acid receptor subunit alpha-1", "GABRA1"),
            ("P47870", "Gamma-aminobutyric acid B receptor 1", "GABBR1"),
            ("Q99250", "Sodium channel protein type 2 subunit alpha", "SCN2A"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "migraine": {
        "label": "Migraine",
        "category": "Neurological",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P08908", "5-hydroxytryptamine receptor 1A", "HTR1A"),
            ("P28223", "5-hydroxytryptamine receptor 2A", "HTR2A"),
            ("P28222", "5-hydroxytryptamine receptor 1B", "HTR1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P21917", "D(4) dopamine receptor", "DRD4"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "huntingtons": {
        "label": "Huntington's Disease",
        "category": "Neurodegenerative",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P42858", "Huntingtin", "HTT"),
            ("P10636", "Microtubule-associated protein tau", "MAPT"),
            ("Q05586", "Glutamate receptor NMDA type subunit 1", "GRIN1"),
            ("P28472", "Glutamate decarboxylase 2", "GAD2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "autism": {
        "label": "Autism Spectrum Disorder",
        "category": "Neuropsychiatric",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000161", "L-Alanine"),
        ],
        "proteins": [
            ("P08908", "5-hydroxytryptamine receptor 1A", "HTR1A"),
            ("P28223", "5-hydroxytryptamine receptor 2A", "HTR2A"),
            ("Q05586", "Glutamate receptor NMDA type subunit 1", "GRIN1"),
            ("P18507", "Gamma-aminobutyric acid receptor subunit alpha-1", "GABRA1"),
            ("P21917", "D(4) dopamine receptor", "DRD4"),
            ("P21728", "D(1A) dopamine receptor", "DRD1"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    # ── AUTOIMMUNE / INFLAMMATORY (3 new) ───────────────────────────

    "psoriasis": {
        "label": "Psoriasis",
        "category": "Autoimmune",
        "tissue": "Skin",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("Q16539", "Mitogen-activated protein kinase 14", "MAPK14"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "celiac_disease": {
        "label": "Celiac Disease",
        "category": "Autoimmune",
        "tissue": "Gut / Intestine",
        "metabolites": [
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P21980", "Protein-glutamine gamma-glutamyltransferase 2", "TGM2"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "ankylosing_spondylitis": {
        "label": "Ankylosing Spondylitis",
        "category": "Autoimmune",
        "tissue": "Spine / Joints",
        "metabolites": [
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("Q16539", "Mitogen-activated protein kinase 14", "MAPK14"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── LIVER / GI (3 new) ──────────────────────────────────────────

    "liver_cirrhosis": {
        "label": "Liver Cirrhosis",
        "category": "Hepatic",
        "tissue": "Liver",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000696", "L-Methionine"),
            ("HMDB0000177", "L-Histidine"),
        ],
        "proteins": [
            ("P02768", "Serum albumin", "ALB"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P07355", "Annexin A2", "ANXA2"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "pancreatitis": {
        "label": "Acute Pancreatitis",
        "category": "Gastrointestinal",
        "tissue": "Pancreas",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000172", "L-Isoleucine"),
        ],
        "proteins": [
            ("P09093", "Pancreatic lipase", "PNLIP"),
            ("P05451", "Lithostathine-1-alpha", "REG1A"),
            ("P04746", "Pancreatic alpha-amylase", "AMY2A"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "gastric_ulcer": {
        "label": "Peptic Ulcer Disease",
        "category": "Gastrointestinal",
        "tissue": "Stomach",
        "metabolites": [
            ("HMDB0000177", "L-Histidine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P23219", "Prostaglandin G/H synthase 1", "PTGS1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    # ── RENAL (2 new) ───────────────────────────────────────────────

    "nephrotic_syndrome": {
        "label": "Nephrotic Syndrome",
        "category": "Renal",
        "tissue": "Kidney",
        "metabolites": [
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P02768", "Serum albumin", "ALB"),
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "iga_nephropathy": {
        "label": "IgA Nephropathy",
        "category": "Renal",
        "tissue": "Kidney",
        "metabolites": [
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01876", "Immunoglobulin heavy constant alpha 1", "IGHA1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── INFECTIOUS (2 new) ──────────────────────────────────────────

    "tuberculosis": {
        "label": "Tuberculosis",
        "category": "Infectious",
        "tissue": "Lung",
        "metabolites": [
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000696", "L-Methionine"),
            ("HMDB0000177", "L-Histidine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15559", "NAD(P)H dehydrogenase quinone 1", "NQO1"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "hepatitis_b": {
        "label": "Hepatitis B",
        "category": "Infectious",
        "tissue": "Liver",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P02768", "Serum albumin", "ALB"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════
# GENERATION PIPELINE (reuses generate_disease_predictions.py logic)
# ═══════════════════════════════════════════════════════════════════════

def write_metabolites_csv(disease_dir, metabolites):
    """Write metabolites.csv."""
    path = disease_dir / "metabolites.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hmdb_id", "name"])
        for hmdb_id, name in metabolites:
            w.writerow([hmdb_id, name])
    return len(metabolites)


def write_proteins_csv(disease_dir, proteins):
    """Write proteins.csv."""
    path = disease_dir / "proteins.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["uniprot_id", "name", "gene"])
        for uid, name, gene in proteins:
            w.writerow([uid, name, gene])
    return len(proteins)


def main():
    # ── Load MPI database for cross-referencing ──
    sys.path.insert(0, str(PROJECT_ROOT))
    from scripts.generate_disease_predictions import (
        load_mpi_database,
        generate_predictions,
        compute_network_stats,
        compute_pathway_enrichment,
    )

    logger.info("Loading MPI database for cross-referencing...")
    known_pairs, hmdb_to_name, uniprot_to_info, pair_pathways = load_mpi_database()

    created = 0
    skipped = 0

    for folder_name, info in sorted(BATCH1_DISEASES.items()):
        disease_dir = DISEASE_DIR / folder_name

        if disease_dir.exists():
            # Check if it already has predictions
            pred_file = disease_dir / "predictions.csv"
            if pred_file.exists() and pred_file.stat().st_size > 100:
                logger.info(f"SKIP {folder_name}: already exists with predictions")
                skipped += 1
                continue

        disease_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Creating {folder_name}: {info['label']} ({info['category']})")

        n_met = write_metabolites_csv(disease_dir, info["metabolites"])
        n_prot = write_proteins_csv(disease_dir, info["proteins"])
        logger.info(f"  Wrote {n_met} metabolites, {n_prot} proteins")

        # Generate predictions, network stats, pathway enrichment
        pred_df = generate_predictions(disease_dir, known_pairs, hmdb_to_name, uniprot_to_info)
        compute_network_stats(pred_df, disease_dir)
        compute_pathway_enrichment(pred_df, pair_pathways, disease_dir)

        created += 1

    logger.info(f"\n{'='*60}")
    logger.info(f"BATCH 1 COMPLETE: {created} created, {skipped} skipped")
    logger.info(f"Total diseases now: {len(list(DISEASE_DIR.iterdir()))}")

    # Print registry entries for copy-paste into disease.py
    logger.info(f"\n{'='*60}")
    logger.info("DISEASE_REGISTRY entries to add to pages/disease.py:")
    logger.info("=" * 60)
    for folder_name, info in sorted(BATCH1_DISEASES.items()):
        pad = 24 - len(f'"{folder_name}"')
        label_pad = 50 - len(f'"{info["label"]}"')
        cat_pad = 22 - len(f'"{info["category"]}"')
        print(f'    "{folder_name}":{" "*max(pad,1)}{{"label": "{info["label"]}",{" "*max(label_pad,1)}"category": "{info["category"]}",{" "*max(cat_pad,1)}"tissue": "{info["tissue"]}"}},')


if __name__ == "__main__":
    main()
