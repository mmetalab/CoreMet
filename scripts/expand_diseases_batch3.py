#!/usr/bin/env python3
"""
Batch 3: Add 30 more diseases to CoreMet (diseases 101–130).

Usage:
    cd mpi-vgae-web
    python scripts/expand_diseases_batch3.py
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
# BATCH 3: 30 NEW DISEASES
# ═══════════════════════════════════════════════════════════════════════

BATCH3_DISEASES = {

    # ── CANCER (6 new) ──────────────────────────────────────────────

    "retinoblastoma": {
        "label": "Retinoblastoma",
        "category": "Cancer",
        "tissue": "Eye",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P06400", "Retinoblastoma-associated protein", "RB1"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P24385", "G1/S-specific cyclin-D1", "CCND1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "gallbladder_cancer": {
        "label": "Gallbladder Carcinoma",
        "category": "Cancer",
        "tissue": "Gallbladder",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P00533", "Epidermal growth factor receptor", "EGFR"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P24385", "G1/S-specific cyclin-D1", "CCND1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "laryngeal_cancer": {
        "label": "Laryngeal Squamous Cell Carcinoma",
        "category": "Cancer",
        "tissue": "Larynx",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000092", "Choline"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P00533", "Epidermal growth factor receptor", "EGFR"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P42771", "Cyclin-dependent kinase inhibitor 2A", "CDKN2A"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P08670", "Vimentin", "VIM"),
        ],
    },

    "wilms_tumor": {
        "label": "Wilms Tumor (Nephroblastoma)",
        "category": "Cancer",
        "tissue": "Kidney",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
        ],
        "proteins": [
            ("P19544", "Wilms tumor protein", "WT1"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P24385", "G1/S-specific cyclin-D1", "CCND1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
        ],
    },

    "chronic_lymphocytic_leukemia": {
        "label": "Chronic Lymphocytic Leukemia",
        "category": "Cancer",
        "tissue": "Blood / Bone marrow",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("Q07817", "Bcl-2-like protein 1", "BCL2L1"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P38936", "Cyclin-dependent kinase inhibitor 1", "CDKN1A"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "thymic_cancer": {
        "label": "Thymoma",
        "category": "Cancer",
        "tissue": "Thymus",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
        ],
        "proteins": [
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P00533", "Epidermal growth factor receptor", "EGFR"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── OPHTHALMOLOGICAL (3 new) ────────────────────────────────────

    "glaucoma": {
        "label": "Glaucoma",
        "category": "Ophthalmological",
        "tissue": "Eye",
        "metabolites": [
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("Q99720", "Myocilin", "MYOC"),
            ("Q9NZR2", "Optineurin", "OPTN"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "macular_degeneration": {
        "label": "Age-Related Macular Degeneration",
        "category": "Ophthalmological",
        "tissue": "Eye / Retina",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P35968", "Vascular endothelial growth factor receptor 2", "KDR"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "diabetic_retinopathy": {
        "label": "Diabetic Retinopathy",
        "category": "Ophthalmological",
        "tissue": "Eye / Retina",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000243", "Pyruvic acid"),
        ],
        "proteins": [
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P35968", "Vascular endothelial growth factor receptor 2", "KDR"),
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── AUTOIMMUNE / DERMATOLOGICAL (4 new) ─────────────────────────

    "vitiligo": {
        "label": "Vitiligo",
        "category": "Autoimmune",
        "tissue": "Skin",
        "metabolites": [
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P14679", "Tyrosinase", "TYR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P07900", "Heat shock protein HSP 90-alpha", "HSP90AA1"),
        ],
    },

    "dermatomyositis": {
        "label": "Dermatomyositis",
        "category": "Autoimmune",
        "tissue": "Skin / Muscle",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "vasculitis": {
        "label": "Systemic Vasculitis",
        "category": "Autoimmune",
        "tissue": "Blood vessels",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "alopecia_areata": {
        "label": "Alopecia Areata",
        "category": "Autoimmune",
        "tissue": "Skin / Hair follicle",
        "metabolites": [
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P07900", "Heat shock protein HSP 90-alpha", "HSP90AA1"),
            ("P08670", "Vimentin", "VIM"),
        ],
    },

    # ── INFECTIOUS (2 new) ──────────────────────────────────────────

    "dengue": {
        "label": "Dengue Fever",
        "category": "Infectious",
        "tissue": "Blood / Systemic",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P02768", "Serum albumin", "ALB"),
        ],
    },

    "meningitis": {
        "label": "Bacterial Meningitis",
        "category": "Infectious",
        "tissue": "Brain / Meninges",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000177", "L-Histidine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P14174", "Macrophage migration inhibitory factor", "MIF"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P02768", "Serum albumin", "ALB"),
        ],
    },

    # ── GASTROINTESTINAL (3 new) ────────────────────────────────────

    "gallstones": {
        "label": "Cholelithiasis (Gallstones)",
        "category": "Gastrointestinal",
        "tissue": "Gallbladder / Liver",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P02768", "Serum albumin", "ALB"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
        ],
    },

    "diverticulitis": {
        "label": "Diverticulitis",
        "category": "Gastrointestinal",
        "tissue": "Colon",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "kidney_stones": {
        "label": "Nephrolithiasis (Kidney Stones)",
        "category": "Renal",
        "tissue": "Kidney",
        "metabolites": [
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000292", "Xanthine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P47989", "Xanthine dehydrogenase/oxidase", "XDH"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
        ],
    },

    # ── NEUROLOGICAL (3 new) ────────────────────────────────────────

    "neuropathy": {
        "label": "Diabetic Peripheral Neuropathy",
        "category": "Neurological",
        "tissue": "Peripheral nerves",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "dementia_lewy_body": {
        "label": "Dementia with Lewy Bodies",
        "category": "Neurodegenerative",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P37840", "Alpha-synuclein", "SNCA"),
            ("P10636", "Microtubule-associated protein tau", "MAPT"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "myotonic_dystrophy": {
        "label": "Myotonic Dystrophy",
        "category": "Neuromuscular",
        "tissue": "Muscle / Brain",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
        ],
    },

    # ── METABOLIC / SYSTEMIC (3 new) ────────────────────────────────

    "amyloidosis": {
        "label": "Systemic Amyloidosis",
        "category": "Metabolic",
        "tissue": "Systemic",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P02766", "Transthyretin", "TTR"),
            ("P02768", "Serum albumin", "ALB"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "fibromyalgia": {
        "label": "Fibromyalgia",
        "category": "Musculoskeletal",
        "tissue": "Muscle / Systemic",
        "metabolites": [
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P08908", "5-hydroxytryptamine receptor 1A", "HTR1A"),
            ("P28223", "5-hydroxytryptamine receptor 2A", "HTR2A"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "chronic_fatigue": {
        "label": "Chronic Fatigue Syndrome",
        "category": "Systemic",
        "tissue": "Systemic",
        "metabolites": [
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
            ("P07900", "Heat shock protein HSP 90-alpha", "HSP90AA1"),
        ],
    },

    # ── UROLOGICAL (2 new) ──────────────────────────────────────────

    "bph": {
        "label": "Benign Prostatic Hyperplasia",
        "category": "Urological",
        "tissue": "Prostate",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P10275", "Androgen receptor", "AR"),
            ("P07202", "Prostate-specific antigen", "KLK3"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P42574", "Caspase-3", "CASP3"),
        ],
    },

    "interstitial_cystitis": {
        "label": "Interstitial Cystitis",
        "category": "Urological",
        "tissue": "Bladder",
        "metabolites": [
            ("HMDB0000177", "L-Histidine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
        ],
    },

    # ── ENDOCRINE (1 new) ───────────────────────────────────────────

    "hyperaldosteronism": {
        "label": "Primary Hyperaldosteronism",
        "category": "Endocrine",
        "tissue": "Adrenal gland",
        "metabolites": [
            ("HMDB0000145", "Aldosterone"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P08235", "Mineralocorticoid receptor", "NR3C2"),
            ("P04150", "Glucocorticoid receptor", "NR3C1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01308", "Insulin", "INS"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },
}


# ═══════════════════════════════════════════════════════════════════════
# GENERATION HELPERS
# ═══════════════════════════════════════════════════════════════════════

def write_metabolites_csv(disease_dir, metabolites):
    path = disease_dir / "metabolites.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["hmdb_id", "name"])
        for hmdb_id, name in metabolites:
            w.writerow([hmdb_id, name])
    return len(metabolites)


def write_proteins_csv(disease_dir, proteins):
    path = disease_dir / "proteins.csv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["uniprot_id", "name", "gene"])
        for uid, name, gene in proteins:
            w.writerow([uid, name, gene])
    return len(proteins)


def main():
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

    for folder_name, info in sorted(BATCH3_DISEASES.items()):
        disease_dir = DISEASE_DIR / folder_name

        if disease_dir.exists():
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

        pred_df = generate_predictions(disease_dir, known_pairs, hmdb_to_name, uniprot_to_info)
        compute_network_stats(pred_df, disease_dir)
        compute_pathway_enrichment(pred_df, pair_pathways, disease_dir)

        created += 1

    logger.info(f"\n{'='*60}")
    logger.info(f"BATCH 3 COMPLETE: {created} created, {skipped} skipped")
    logger.info(f"Total diseases now: {len(list(DISEASE_DIR.iterdir()))}")


if __name__ == "__main__":
    main()
