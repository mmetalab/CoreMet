#!/usr/bin/env python3
"""
Batch 2: Add 30 more diseases to CoreMet (diseases 71–100).

Usage:
    cd mpi-vgae-web
    python scripts/expand_diseases_batch2.py
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
# BATCH 2: 30 NEW DISEASES
# ═══════════════════════════════════════════════════════════════════════

BATCH2_DISEASES = {

    # ── INFECTIOUS (4 new) ──────────────────────────────────────────

    "hepatitis_c": {
        "label": "Hepatitis C",
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
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000696", "L-Methionine"),
            ("HMDB0000243", "Pyruvic acid"),
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
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "hiv_aids": {
        "label": "HIV/AIDS",
        "category": "Infectious",
        "tissue": "Immune system / Systemic",
        "metabolites": [
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000684", "L-Kynurenine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000696", "L-Methionine"),
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P55211", "Caspase-9", "CASP9"),
            ("P10415", "Apoptosis regulator Bcl-2", "BCL2"),
            ("Q07817", "Bcl-2-like protein 1", "BCL2L1"),
            ("P04637", "Cellular tumor antigen p53", "TP53"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "sepsis": {
        "label": "Sepsis",
        "category": "Infectious",
        "tissue": "Systemic / Blood",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000684", "L-Kynurenine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P14174", "Macrophage migration inhibitory factor", "MIF"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P02768", "Serum albumin", "ALB"),
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "malaria": {
        "label": "Malaria",
        "category": "Infectious",
        "tissue": "Blood / Liver",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000177", "L-Histidine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P69905", "Hemoglobin subunit alpha", "HBA1"),
            ("P68871", "Hemoglobin subunit beta", "HBB"),
            ("P02768", "Serum albumin", "ALB"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── NEUROPSYCHIATRIC (3 new) ────────────────────────────────────

    "bipolar_disorder": {
        "label": "Bipolar Disorder",
        "category": "Neuropsychiatric",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P28223", "5-hydroxytryptamine receptor 2A", "HTR2A"),
            ("P08908", "5-hydroxytryptamine receptor 1A", "HTR1A"),
            ("P21917", "D(4) dopamine receptor", "DRD4"),
            ("P21728", "D(1A) dopamine receptor", "DRD1"),
            ("P31644", "5-hydroxytryptamine receptor 5A", "HTR5A"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "adhd": {
        "label": "Attention Deficit Hyperactivity Disorder",
        "category": "Neuropsychiatric",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P23975", "Sodium-dependent noradrenaline transporter", "SLC6A2"),
            ("Q01959", "Sodium-dependent dopamine transporter", "SLC6A3"),
            ("P21917", "D(4) dopamine receptor", "DRD4"),
            ("P21728", "D(1A) dopamine receptor", "DRD1"),
            ("P28223", "5-hydroxytryptamine receptor 2A", "HTR2A"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "anxiety": {
        "label": "Generalized Anxiety Disorder",
        "category": "Neuropsychiatric",
        "tissue": "Brain",
        "metabolites": [
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000112", "gamma-Aminobutyric acid"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P08908", "5-hydroxytryptamine receptor 1A", "HTR1A"),
            ("P28223", "5-hydroxytryptamine receptor 2A", "HTR2A"),
            ("P18507", "Gamma-aminobutyric acid receptor subunit alpha-1", "GABRA1"),
            ("P47870", "Gamma-aminobutyric acid B receptor 1", "GABBR1"),
            ("P21917", "D(4) dopamine receptor", "DRD4"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    # ── MUSCULOSKELETAL / METABOLIC (3 new) ─────────────────────────

    "osteoporosis": {
        "label": "Osteoporosis",
        "category": "Musculoskeletal",
        "tissue": "Bone",
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
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P03372", "Estrogen receptor", "ESR1"),
            ("O00300", "Tumor necrosis factor receptor superfamily member 11B", "TNFRSF11B"),
            ("O14788", "Tumor necrosis factor ligand superfamily member 11", "TNFSF11"),
            ("Q9Y6Q6", "Tumor necrosis factor ligand superfamily member 11", "TNFRSF11A"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "osteoarthritis": {
        "label": "Osteoarthritis",
        "category": "Musculoskeletal",
        "tissue": "Joints / Cartilage",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000158", "L-Tyrosine"),
        ],
        "proteins": [
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P03956", "Interstitial collagenase", "MMP1"),
            ("P08253", "72 kDa type IV collagenase", "MMP2"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "phenylketonuria": {
        "label": "Phenylketonuria",
        "category": "Metabolic",
        "tissue": "Liver / Brain",
        "metabolites": [
            ("HMDB0000159", "L-Phenylalanine"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000259", "Serotonin"),
            ("HMDB0000073", "Dopamine"),
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
        ],
        "proteins": [
            ("P00439", "Phenylalanine-4-hydroxylase", "PAH"),
            ("P17735", "Tyrosine aminotransferase", "TAT"),
            ("P07101", "Tyrosine 3-monooxygenase", "TH"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P42574", "Caspase-3", "CASP3"),
        ],
    },

    # ── OBSTETRIC / GYNECOLOGICAL (2 new) ───────────────────────────

    "preeclampsia": {
        "label": "Preeclampsia",
        "category": "Obstetric",
        "tissue": "Placenta / Systemic",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P17948", "Vascular endothelial growth factor receptor 1", "FLT1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P02768", "Serum albumin", "ALB"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "endometriosis": {
        "label": "Endometriosis",
        "category": "Gynecological",
        "tissue": "Uterus / Pelvis",
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
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P03372", "Estrogen receptor", "ESR1"),
            ("P04150", "Glucocorticoid receptor", "NR3C1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    # ── GENETIC / HEMATOLOGICAL (3 new) ─────────────────────────────

    "cystic_fibrosis": {
        "label": "Cystic Fibrosis",
        "category": "Genetic",
        "tissue": "Lung / Pancreas",
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
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P13569", "Cystic fibrosis transmembrane conductance regulator", "CFTR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "sickle_cell": {
        "label": "Sickle Cell Disease",
        "category": "Hematological",
        "tissue": "Blood",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("P69905", "Hemoglobin subunit alpha", "HBA1"),
            ("P68871", "Hemoglobin subunit beta", "HBB"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "hemophilia": {
        "label": "Hemophilia A",
        "category": "Hematological",
        "tissue": "Blood",
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
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P00451", "Coagulation factor VIII", "F8"),
            ("P00740", "Coagulation factor IX", "F9"),
            ("P00734", "Prothrombin", "F2"),
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P02675", "Fibrinogen beta chain", "FGB"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── METABOLIC / HEPATIC (3 new) ─────────────────────────────────

    "wilsons_disease": {
        "label": "Wilson's Disease",
        "category": "Metabolic",
        "tissue": "Liver / Brain",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000696", "L-Methionine"),
        ],
        "proteins": [
            ("P35670", "Copper-transporting ATPase 2", "ATP7B"),
            ("P00450", "Ceruloplasmin", "CP"),
            ("P02768", "Serum albumin", "ALB"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    "hemochromatosis": {
        "label": "Hereditary Hemochromatosis",
        "category": "Metabolic",
        "tissue": "Liver",
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
            ("HMDB0000243", "Pyruvic acid"),
        ],
        "proteins": [
            ("Q30201", "Hereditary hemochromatosis protein", "HFE"),
            ("P02787", "Serotransferrin", "TF"),
            ("P02792", "Ferritin light chain", "FTL"),
            ("P02794", "Ferritin heavy chain", "FTH1"),
            ("Q9NQ26", "Hepcidin", "HAMP"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "cushings_syndrome": {
        "label": "Cushing's Syndrome",
        "category": "Endocrine",
        "tissue": "Adrenal / Pituitary",
        "metabolites": [
            ("HMDB0000063", "Cortisol"),
            ("HMDB0000015", "Cortisone"),
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
        ],
        "proteins": [
            ("P04150", "Glucocorticoid receptor", "NR3C1"),
            ("P01189", "Pro-opiomelanocortin", "POMC"),
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    # ── AUTOIMMUNE (4 new) ──────────────────────────────────────────

    "myasthenia_gravis": {
        "label": "Myasthenia Gravis",
        "category": "Autoimmune",
        "tissue": "Neuromuscular junction",
        "metabolites": [
            ("HMDB0000895", "Acetylcholine"),
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
            ("P02708", "Acetylcholine receptor subunit alpha", "CHRNA1"),
            ("P11230", "Acetylcholine receptor subunit beta", "CHRNB1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P60174", "Triosephosphate isomerase", "TPI1"),
        ],
    },

    "sjogrens_syndrome": {
        "label": "Sjögren's Syndrome",
        "category": "Autoimmune",
        "tissue": "Exocrine glands",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000158", "L-Tyrosine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "scleroderma": {
        "label": "Systemic Sclerosis (Scleroderma)",
        "category": "Autoimmune",
        "tissue": "Skin / Systemic",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000159", "L-Phenylalanine"),
        ],
        "proteins": [
            ("P01137", "Transforming growth factor beta-1", "TGFB1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    "crohns_disease": {
        "label": "Crohn's Disease",
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
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000172", "L-Isoleucine"),
            ("HMDB0000177", "L-Histidine"),
        ],
        "proteins": [
            ("Q9HC29", "Nucleotide-binding oligomerization domain-containing protein 2", "NOD2"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── CARDIOVASCULAR (3 new) ──────────────────────────────────────

    "pulmonary_hypertension": {
        "label": "Pulmonary Arterial Hypertension",
        "category": "Cardiovascular",
        "tissue": "Lung / Heart",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000243", "Pyruvic acid"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000687", "L-Leucine"),
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P35968", "Vascular endothelial growth factor receptor 2", "KDR"),
            ("P29474", "Nitric oxide synthase brain", "NOS1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01137", "Transforming growth factor beta-1", "TGFB1"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "atrial_fibrillation": {
        "label": "Atrial Fibrillation",
        "category": "Cardiovascular",
        "tissue": "Heart",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000562", "Creatine"),
            ("HMDB0000254", "Succinic acid"),
            ("HMDB0000148", "L-Glutamic acid"),
            ("HMDB0000641", "L-Glutamine"),
            ("HMDB0000161", "L-Alanine"),
            ("HMDB0000883", "L-Valine"),
            ("HMDB0000929", "L-Tryptophan"),
            ("HMDB0000687", "L-Leucine"),
        ],
        "proteins": [
            ("P45379", "Troponin T cardiac muscle", "TNNT2"),
            ("P19429", "Troponin I cardiac muscle", "TNNI3"),
            ("Q14524", "Sodium channel protein type 5 subunit alpha", "SCN5A"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P02675", "Fibrinogen beta chain", "FGB"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "peripheral_artery": {
        "label": "Peripheral Artery Disease",
        "category": "Cardiovascular",
        "tissue": "Blood vessels / Limbs",
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
            ("HMDB0000562", "Creatine"),
        ],
        "proteins": [
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P05362", "Intercellular adhesion molecule 1", "ICAM1"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P02671", "Fibrinogen alpha chain", "FGA"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P35354", "Prostaglandin G/H synthase 2", "PTGS2"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── RESPIRATORY (1 new) ─────────────────────────────────────────

    "pulmonary_fibrosis": {
        "label": "Idiopathic Pulmonary Fibrosis",
        "category": "Respiratory",
        "tissue": "Lung",
        "metabolites": [
            ("HMDB0000190", "L-Lactic acid"),
            ("HMDB0000122", "D-Glucose"),
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
            ("P01137", "Transforming growth factor beta-1", "TGFB1"),
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P14780", "Matrix metalloproteinase-9", "MMP9"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    # ── RENAL (1 new) ───────────────────────────────────────────────

    "diabetic_nephropathy": {
        "label": "Diabetic Nephropathy",
        "category": "Renal",
        "tissue": "Kidney",
        "metabolites": [
            ("HMDB0000122", "D-Glucose"),
            ("HMDB0000064", "Creatinine"),
            ("HMDB0000289", "Uric acid"),
            ("HMDB0000190", "L-Lactic acid"),
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
            ("P15692", "Vascular endothelial growth factor A", "VEGFA"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01137", "Transforming growth factor beta-1", "TGFB1"),
            ("P02768", "Serum albumin", "ALB"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
        ],
    },

    # ── ENDOCRINE (2 new) ───────────────────────────────────────────

    "addisons_disease": {
        "label": "Addison's Disease",
        "category": "Endocrine",
        "tissue": "Adrenal gland",
        "metabolites": [
            ("HMDB0000063", "Cortisol"),
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
            ("P04150", "Glucocorticoid receptor", "NR3C1"),
            ("P08235", "Mineralocorticoid receptor", "NR3C2"),
            ("P01189", "Pro-opiomelanocortin", "POMC"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P01584", "Interleukin-1 beta", "IL1B"),
            ("P42574", "Caspase-3", "CASP3"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
        ],
    },

    "acromegaly": {
        "label": "Acromegaly",
        "category": "Endocrine",
        "tissue": "Pituitary",
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
            ("HMDB0000254", "Succinic acid"),
        ],
        "proteins": [
            ("P01241", "Somatotropin", "GH1"),
            ("P05019", "Insulin-like growth factor I", "IGF1"),
            ("P01308", "Insulin", "INS"),
            ("P06213", "Insulin receptor", "INSR"),
            ("P01375", "Tumor necrosis factor", "TNF"),
            ("P05231", "Interleukin-6", "IL6"),
            ("P14618", "Pyruvate kinase PKM", "PKM"),
            ("P04406", "Glyceraldehyde-3-phosphate dehydrogenase", "GAPDH"),
            ("P06733", "Alpha-enolase", "ENO1"),
            ("P40926", "Malate dehydrogenase, mitochondrial", "MDH2"),
        ],
    },

    # ── GI (1 new) ──────────────────────────────────────────────────

    "ulcerative_colitis": {
        "label": "Ulcerative Colitis",
        "category": "Autoimmune",
        "tissue": "Colon",
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
            ("HMDB0000177", "L-Histidine"),
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
}


# ═══════════════════════════════════════════════════════════════════════
# GENERATION HELPERS (same as batch 1)
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

    for folder_name, info in sorted(BATCH2_DISEASES.items()):
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
    logger.info(f"BATCH 2 COMPLETE: {created} created, {skipped} skipped")
    logger.info(f"Total diseases now: {len(list(DISEASE_DIR.iterdir()))}")


if __name__ == "__main__":
    main()
