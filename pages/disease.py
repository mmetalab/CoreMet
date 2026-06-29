"""
Disease Analysis page, interactive network explorer for pre-computed
disease-specific metabolite-protein interaction predictions.

Features:
  - Sidebar: disease selector, organism, confidence slider, network stats, export
  - Main: Cytoscape network graph (CoSE), stats row, hub tables, enrichment chart,
    predictions DataTable, node-click detail panel
"""

import io
import json

from dash import dcc, html, Input, Output, State, callback, dash_table, no_update, ctx
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import pandas as pd
import plotly.graph_objects as go

from components.page_header import make_page_header

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Full disease registry with metadata for filtering
DISEASE_REGISTRY = {
    # value: {label, category, tissue}
    # ── Cancers (20) ──
    "aml_leukemia":       {"label": "Acute Myeloid Leukemia",          "category": "Cancer",       "tissue": "Blood / Bone marrow"},
    "bladder_cancer":     {"label": "Bladder Cancer",                  "category": "Cancer",       "tissue": "Bladder"},
    "breast_cancer":      {"label": "Breast Cancer",                   "category": "Cancer",       "tissue": "Breast"},
    "cervical_cancer":    {"label": "Cervical Cancer",                 "category": "Cancer",       "tissue": "Cervix"},
    "cholangiocarcinoma": {"label": "Cholangiocarcinoma",              "category": "Cancer",       "tissue": "Bile duct / Liver"},
    "colorectal_cancer":  {"label": "Colorectal Cancer",               "category": "Cancer",       "tissue": "Colon / Rectum"},
    "endometrial_cancer": {"label": "Endometrial Cancer",              "category": "Cancer",       "tissue": "Uterus"},
    "esophageal_cancer":  {"label": "Esophageal Cancer",               "category": "Cancer",       "tissue": "Esophagus"},
    "gastric_cancer":     {"label": "Gastric Cancer",                  "category": "Cancer",       "tissue": "Stomach"},
    "glioblastoma":       {"label": "Glioblastoma",                    "category": "Cancer",       "tissue": "Brain"},
    "hcc":                {"label": "Hepatocellular Carcinoma (HCC)",  "category": "Cancer",       "tissue": "Liver"},
    "lung_cancer":        {"label": "Lung Cancer",                     "category": "Cancer",       "tissue": "Lung"},
    "lymphoma":           {"label": "Lymphoma",                        "category": "Cancer",       "tissue": "Lymph nodes"},
    "melanoma":           {"label": "Melanoma",                        "category": "Cancer",       "tissue": "Skin"},
    "oral_cancer":        {"label": "Oral Squamous Cell Carcinoma",    "category": "Cancer",       "tissue": "Oral cavity"},
    "ovarian_cancer":     {"label": "Ovarian Cancer",                  "category": "Cancer",       "tissue": "Ovary"},
    "pancreatic_cancer":  {"label": "Pancreatic Cancer",               "category": "Cancer",       "tissue": "Pancreas"},
    "thyroid_cancer":     {"label": "Papillary Thyroid Cancer",        "category": "Cancer",       "tissue": "Thyroid"},
    "prostate_cancer":    {"label": "Prostate Cancer",                 "category": "Cancer",       "tissue": "Prostate"},
    "renal_cell_carcinoma": {"label": "Renal Cell Carcinoma",          "category": "Cancer",       "tissue": "Kidney"},
    # ── Neurodegenerative / Neuropsychiatric (6) ──
    "alzheimers":         {"label": "Alzheimer's Disease",             "category": "Neurodegenerative", "tissue": "Brain"},
    "als":                {"label": "Amyotrophic Lateral Sclerosis",   "category": "Neurodegenerative", "tissue": "Brain / Spinal cord"},
    "depression":         {"label": "Major Depressive Disorder",       "category": "Neuropsychiatric",  "tissue": "Brain"},
    "multiple_sclerosis": {"label": "Multiple Sclerosis",              "category": "Autoimmune",   "tissue": "Brain / Spinal cord"},
    "parkinsons":         {"label": "Parkinson's Disease",             "category": "Neurodegenerative", "tissue": "Brain"},
    "schizophrenia":      {"label": "Schizophrenia",                   "category": "Neuropsychiatric",  "tissue": "Brain"},
    # ── Metabolic (4) ──
    "metabolic_syndrome": {"label": "Metabolic Syndrome",              "category": "Metabolic",    "tissue": "Systemic"},
    "nafld":              {"label": "Non-Alcoholic Fatty Liver Disease","category": "Metabolic",    "tissue": "Liver"},
    "obesity":            {"label": "Obesity",                         "category": "Metabolic",    "tissue": "Adipose / Systemic"},
    "t2_diabetes":        {"label": "Type 2 Diabetes",                 "category": "Metabolic",    "tissue": "Pancreas / Systemic"},
    # ── Cardiovascular (3) ──
    "atherosclerosis":    {"label": "Atherosclerosis",                 "category": "Cardiovascular","tissue": "Blood vessels"},
    "heart_failure":      {"label": "Heart Failure",                   "category": "Cardiovascular","tissue": "Heart"},
    "hypertension":       {"label": "Hypertension",                    "category": "Cardiovascular","tissue": "Blood vessels"},
    # ── Autoimmune / Inflammatory (3) ──
    "ibd":                {"label": "Inflammatory Bowel Disease",      "category": "Autoimmune",   "tissue": "Gut / Intestine"},
    "rheumatoid_arthritis": {"label": "Rheumatoid Arthritis",          "category": "Autoimmune",   "tissue": "Joints / Synovium"},
    "lupus":              {"label": "Systemic Lupus Erythematosus",    "category": "Autoimmune",   "tissue": "Systemic"},
    # ── Respiratory (2) ──
    "asthma":             {"label": "Asthma",                          "category": "Respiratory",  "tissue": "Lung"},
    "copd":               {"label": "Chronic Obstructive Pulmonary Disease", "category": "Respiratory", "tissue": "Lung"},
    # ── Other (2) ──
    "covid19":            {"label": "COVID-19",                        "category": "Infectious",   "tissue": "Lung / Systemic"},
    "chronic_kidney_disease": {"label": "Chronic Kidney Disease",      "category": "Renal",        "tissue": "Kidney"},

    # ═══════════════════════════════════════════════════════════════════
    # BATCH 1, 30 NEW DISEASES (auto-generated)
    # ═══════════════════════════════════════════════════════════════════
    # ── Cancers (8 new) ──
    "adrenal_cancer":        {"label": "Adrenocortical Carcinoma",               "category": "Cancer",           "tissue": "Adrenal gland"},
    "head_neck_cancer":      {"label": "Head and Neck Squamous Cell Carcinoma",  "category": "Cancer",           "tissue": "Head / Neck"},
    "mesothelioma":          {"label": "Mesothelioma",                           "category": "Cancer",           "tissue": "Pleura / Lung"},
    "multiple_myeloma":      {"label": "Multiple Myeloma",                       "category": "Cancer",           "tissue": "Blood / Bone marrow"},
    "nasopharyngeal_cancer": {"label": "Nasopharyngeal Carcinoma",               "category": "Cancer",           "tissue": "Nasopharynx"},
    "neuroblastoma":         {"label": "Neuroblastoma",                          "category": "Cancer",           "tissue": "Adrenal / Nervous system"},
    "sarcoma":               {"label": "Soft Tissue Sarcoma",                    "category": "Cancer",           "tissue": "Soft tissue"},
    "testicular_cancer":     {"label": "Testicular Germ Cell Tumor",             "category": "Cancer",           "tissue": "Testis"},
    # ── Endocrine (3 new) ──
    "hyperthyroidism":       {"label": "Hyperthyroidism (Graves' Disease)",      "category": "Endocrine",        "tissue": "Thyroid / Systemic"},
    "hypothyroidism":        {"label": "Hypothyroidism",                         "category": "Endocrine",        "tissue": "Thyroid / Systemic"},
    "pcos":                  {"label": "Polycystic Ovary Syndrome",              "category": "Endocrine",        "tissue": "Ovary / Systemic"},
    # ── Metabolic (2 new) ──
    "gout":                  {"label": "Gout",                                   "category": "Metabolic",        "tissue": "Joints / Systemic"},
    "type1_diabetes":        {"label": "Type 1 Diabetes Mellitus",               "category": "Metabolic",        "tissue": "Pancreas / Systemic"},
    # ── Cardiovascular (3 new) ──
    "aortic_aneurysm":       {"label": "Abdominal Aortic Aneurysm",             "category": "Cardiovascular",   "tissue": "Aorta / Blood vessels"},
    "cardiomyopathy":        {"label": "Dilated Cardiomyopathy",                "category": "Cardiovascular",   "tissue": "Heart"},
    "stroke":                {"label": "Ischemic Stroke",                        "category": "Cardiovascular",   "tissue": "Brain / Blood vessels"},
    # ── Neurological / Neuropsychiatric (4 new) ──
    "autism":                {"label": "Autism Spectrum Disorder",                "category": "Neuropsychiatric", "tissue": "Brain"},
    "epilepsy":              {"label": "Epilepsy",                               "category": "Neurological",     "tissue": "Brain"},
    "huntingtons":           {"label": "Huntington's Disease",                   "category": "Neurodegenerative","tissue": "Brain"},
    "migraine":              {"label": "Migraine",                               "category": "Neurological",     "tissue": "Brain"},
    # ── Autoimmune (3 new) ──
    "ankylosing_spondylitis":{"label": "Ankylosing Spondylitis",                "category": "Autoimmune",       "tissue": "Spine / Joints"},
    "celiac_disease":        {"label": "Celiac Disease",                         "category": "Autoimmune",       "tissue": "Gut / Intestine"},
    "psoriasis":             {"label": "Psoriasis",                              "category": "Autoimmune",       "tissue": "Skin"},
    # ── Hepatic / GI (3 new) ──
    "gastric_ulcer":         {"label": "Peptic Ulcer Disease",                   "category": "Gastrointestinal", "tissue": "Stomach"},
    "liver_cirrhosis":       {"label": "Liver Cirrhosis",                        "category": "Hepatic",          "tissue": "Liver"},
    "pancreatitis":          {"label": "Acute Pancreatitis",                     "category": "Gastrointestinal", "tissue": "Pancreas"},
    # ── Renal (2 new) ──
    "iga_nephropathy":       {"label": "IgA Nephropathy",                        "category": "Renal",            "tissue": "Kidney"},
    "nephrotic_syndrome":    {"label": "Nephrotic Syndrome",                     "category": "Renal",            "tissue": "Kidney"},
    # ── Infectious (2 new) ──
    "hepatitis_b":           {"label": "Hepatitis B",                            "category": "Infectious",       "tissue": "Liver"},
    "tuberculosis":          {"label": "Tuberculosis",                           "category": "Infectious",       "tissue": "Lung"},

    # ═══════════════════════════════════════════════════════════════════
    # BATCH 2, 30 NEW DISEASES (auto-generated)
    # ═══════════════════════════════════════════════════════════════════
    # ── Infectious (4 new) ──
    "hepatitis_c":           {"label": "Hepatitis C",                            "category": "Infectious",       "tissue": "Liver"},
    "hiv_aids":              {"label": "HIV/AIDS",                               "category": "Infectious",       "tissue": "Immune system / Systemic"},
    "malaria":               {"label": "Malaria",                                "category": "Infectious",       "tissue": "Blood / Liver"},
    "sepsis":                {"label": "Sepsis",                                 "category": "Infectious",       "tissue": "Systemic / Blood"},
    # ── Neuropsychiatric (3 new) ──
    "adhd":                  {"label": "Attention Deficit Hyperactivity Disorder","category": "Neuropsychiatric", "tissue": "Brain"},
    "anxiety":               {"label": "Generalized Anxiety Disorder",           "category": "Neuropsychiatric", "tissue": "Brain"},
    "bipolar_disorder":      {"label": "Bipolar Disorder",                       "category": "Neuropsychiatric", "tissue": "Brain"},
    # ── Musculoskeletal / Metabolic (3 new) ──
    "osteoarthritis":        {"label": "Osteoarthritis",                         "category": "Musculoskeletal",  "tissue": "Joints / Cartilage"},
    "osteoporosis":          {"label": "Osteoporosis",                           "category": "Musculoskeletal",  "tissue": "Bone"},
    "phenylketonuria":       {"label": "Phenylketonuria",                        "category": "Metabolic",        "tissue": "Liver / Brain"},
    # ── Obstetric / Gynecological (2 new) ──
    "endometriosis":         {"label": "Endometriosis",                          "category": "Gynecological",    "tissue": "Uterus / Pelvis"},
    "preeclampsia":          {"label": "Preeclampsia",                           "category": "Obstetric",        "tissue": "Placenta / Systemic"},
    # ── Genetic / Hematological (3 new) ──
    "cystic_fibrosis":       {"label": "Cystic Fibrosis",                        "category": "Genetic",          "tissue": "Lung / Pancreas"},
    "hemophilia":            {"label": "Hemophilia A",                            "category": "Hematological",    "tissue": "Blood"},
    "sickle_cell":           {"label": "Sickle Cell Disease",                    "category": "Hematological",    "tissue": "Blood"},
    # ── Metabolic / Hepatic (3 new) ──
    "hemochromatosis":       {"label": "Hereditary Hemochromatosis",             "category": "Metabolic",        "tissue": "Liver"},
    "wilsons_disease":       {"label": "Wilson's Disease",                       "category": "Metabolic",        "tissue": "Liver / Brain"},
    "cushings_syndrome":     {"label": "Cushing's Syndrome",                     "category": "Endocrine",        "tissue": "Adrenal / Pituitary"},
    # ── Autoimmune (5 new) ──
    "crohns_disease":        {"label": "Crohn's Disease",                        "category": "Autoimmune",       "tissue": "Gut / Intestine"},
    "myasthenia_gravis":     {"label": "Myasthenia Gravis",                      "category": "Autoimmune",       "tissue": "Neuromuscular junction"},
    "scleroderma":           {"label": "Systemic Sclerosis (Scleroderma)",       "category": "Autoimmune",       "tissue": "Skin / Systemic"},
    "sjogrens_syndrome":     {"label": "Sjögren's Syndrome",                     "category": "Autoimmune",       "tissue": "Exocrine glands"},
    "ulcerative_colitis":    {"label": "Ulcerative Colitis",                     "category": "Autoimmune",       "tissue": "Colon"},
    # ── Cardiovascular (3 new) ──
    "atrial_fibrillation":   {"label": "Atrial Fibrillation",                    "category": "Cardiovascular",   "tissue": "Heart"},
    "peripheral_artery":     {"label": "Peripheral Artery Disease",              "category": "Cardiovascular",   "tissue": "Blood vessels / Limbs"},
    "pulmonary_hypertension":{"label": "Pulmonary Arterial Hypertension",        "category": "Cardiovascular",   "tissue": "Lung / Heart"},
    # ── Respiratory (1 new) ──
    "pulmonary_fibrosis":    {"label": "Idiopathic Pulmonary Fibrosis",          "category": "Respiratory",      "tissue": "Lung"},
    # ── Renal (1 new) ──
    "diabetic_nephropathy":  {"label": "Diabetic Nephropathy",                   "category": "Renal",            "tissue": "Kidney"},
    # ── Endocrine (2 new) ──
    "acromegaly":            {"label": "Acromegaly",                             "category": "Endocrine",        "tissue": "Pituitary"},
    "addisons_disease":      {"label": "Addison's Disease",                      "category": "Endocrine",        "tissue": "Adrenal gland"},

    # ═══════════════════════ BATCH 3 (30 diseases) ═══════════════════════
    # ── Cancer (6 new) ──
    "retinoblastoma":                {"label": "Retinoblastoma",                          "category": "Cancer",           "tissue": "Eye"},
    "gallbladder_cancer":            {"label": "Gallbladder Carcinoma",                   "category": "Cancer",           "tissue": "Gallbladder"},
    "laryngeal_cancer":              {"label": "Laryngeal Squamous Cell Carcinoma",       "category": "Cancer",           "tissue": "Larynx"},
    "wilms_tumor":                   {"label": "Wilms Tumor (Nephroblastoma)",            "category": "Cancer",           "tissue": "Kidney"},
    "chronic_lymphocytic_leukemia":  {"label": "Chronic Lymphocytic Leukemia",            "category": "Cancer",           "tissue": "Blood / Bone marrow"},
    "thymic_cancer":                 {"label": "Thymoma",                                 "category": "Cancer",           "tissue": "Thymus"},
    # ── Ophthalmological (3 new) ──
    "glaucoma":                      {"label": "Glaucoma",                                "category": "Ophthalmological", "tissue": "Eye"},
    "macular_degeneration":          {"label": "Age-Related Macular Degeneration",        "category": "Ophthalmological", "tissue": "Eye / Retina"},
    "diabetic_retinopathy":          {"label": "Diabetic Retinopathy",                    "category": "Ophthalmological", "tissue": "Eye / Retina"},
    # ── Autoimmune / Dermatological (4 new) ──
    "vitiligo":                      {"label": "Vitiligo",                                "category": "Autoimmune",       "tissue": "Skin"},
    "dermatomyositis":               {"label": "Dermatomyositis",                         "category": "Autoimmune",       "tissue": "Skin / Muscle"},
    "vasculitis":                    {"label": "Systemic Vasculitis",                     "category": "Autoimmune",       "tissue": "Blood vessels"},
    "alopecia_areata":               {"label": "Alopecia Areata",                         "category": "Autoimmune",       "tissue": "Skin / Hair follicle"},
    # ── Infectious (2 new) ──
    "dengue":                        {"label": "Dengue Fever",                            "category": "Infectious",       "tissue": "Blood / Systemic"},
    "meningitis":                    {"label": "Bacterial Meningitis",                    "category": "Infectious",       "tissue": "Brain / Meninges"},
    # ── Gastrointestinal (3 new) ──
    "gallstones":                    {"label": "Cholelithiasis (Gallstones)",             "category": "Gastrointestinal", "tissue": "Gallbladder / Liver"},
    "diverticulitis":                {"label": "Diverticulitis",                          "category": "Gastrointestinal", "tissue": "Colon"},
    "gastroparesis":                 {"label": "Gastroparesis",                           "category": "Gastrointestinal", "tissue": "Stomach"},
    # ── Renal (1 new) ──
    "kidney_stones":                 {"label": "Nephrolithiasis (Kidney Stones)",         "category": "Renal",            "tissue": "Kidney"},
    # ── Neurological / Neurodegenerative (3 new) ──
    "neuropathy":                    {"label": "Diabetic Peripheral Neuropathy",          "category": "Neurological",     "tissue": "Peripheral nerves"},
    "dementia_lewy_body":            {"label": "Dementia with Lewy Bodies",              "category": "Neurodegenerative", "tissue": "Brain"},
    "myotonic_dystrophy":            {"label": "Myotonic Dystrophy",                     "category": "Neuromuscular",    "tissue": "Muscle / Brain"},
    # ── Metabolic / Systemic (3 new) ──
    "amyloidosis":                   {"label": "Systemic Amyloidosis",                   "category": "Metabolic",        "tissue": "Systemic"},
    "fibromyalgia":                  {"label": "Fibromyalgia",                            "category": "Musculoskeletal",  "tissue": "Muscle / Systemic"},
    "chronic_fatigue":               {"label": "Chronic Fatigue Syndrome",               "category": "Systemic",         "tissue": "Systemic"},
    # ── Urological (2 new) ──
    "bph":                           {"label": "Benign Prostatic Hyperplasia",            "category": "Urological",       "tissue": "Prostate"},
    "interstitial_cystitis":         {"label": "Interstitial Cystitis",                  "category": "Urological",       "tissue": "Bladder"},
    # ── Endocrine (1 new) ──
    "hyperaldosteronism":            {"label": "Primary Hyperaldosteronism",              "category": "Endocrine",        "tissue": "Adrenal gland"},
    # ── Hematological (2 new) ──
    "polycythemia_vera":             {"label": "Polycythemia Vera",                      "category": "Hematological",    "tissue": "Blood / Bone marrow"},
    "thalassemia":                   {"label": "Beta-Thalassemia",                       "category": "Hematological",    "tissue": "Blood"},
}

# Flat list for dropdown (built from registry)
DISEASE_OPTIONS = [
    {"label": info["label"], "value": key}
    for key, info in sorted(DISEASE_REGISTRY.items(), key=lambda x: x[1]["label"])
]

# Unique categories and tissues (for filter dropdowns)
CATEGORY_OPTIONS = sorted({info["category"] for info in DISEASE_REGISTRY.values()})
TISSUE_OPTIONS = sorted({info["tissue"] for info in DISEASE_REGISTRY.values()})

CYTOSCAPE_STYLESHEET = [
    # Metabolite nodes, blue circles
    {
        "selector": "node.metabolite",
        "style": {
            "background-color": "#3182ce",
            "label": "data(label)",
            "shape": "ellipse",
            "width": "mapData(degree, 1, 30, 28, 60)",
            "height": "mapData(degree, 1, 30, 28, 60)",
            "font-size": "9px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#2d3748",
            "text-margin-y": "4px",
            "border-width": "2px",
            "border-color": "#2b6cb0",
        },
    },
    # Protein nodes, red rounded-rectangles
    {
        "selector": "node.protein",
        "style": {
            "background-color": "#e53e3e",
            "label": "data(label)",
            "shape": "round-rectangle",
            "width": "mapData(degree, 1, 30, 32, 64)",
            "height": "mapData(degree, 1, 30, 24, 48)",
            "font-size": "9px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#2d3748",
            "text-margin-y": "4px",
            "border-width": "2px",
            "border-color": "#c53030",
        },
    },
    # Edges
    {
        "selector": "edge",
        "style": {
            "width": "mapData(score, 0.3, 1, 1, 4)",
            "line-color": "#a0aec0",
            "opacity": 0.6,
            "curve-style": "bezier",
        },
    },
    # Selected node highlight
    {
        "selector": ":selected",
        "style": {
            "border-width": "4px",
            "border-color": "#d69e2e",
            "background-color": "#ecc94b",
        },
    },
]

# ---------------------------------------------------------------------------
# Helper builders
# ---------------------------------------------------------------------------


def _make_stat_badge(icon_cls, value, label, color="#3182ce"):
    """Small stats badge for the overview row."""
    return html.Div(
        [
            html.Div(
                [
                    html.I(className=f"{icon_cls} me-1", style={"color": color}),
                    html.Span(str(value), style={"fontWeight": "700", "fontSize": "1.1rem"}),
                ],
                className="d-flex align-items-center",
            ),
            html.Div(label, style={"fontSize": "0.75rem", "color": "#718096"}),
        ],
        className="text-center",
        style={
            "padding": "12px 16px",
            "background": "#f7fafc",
            "borderRadius": "8px",
            "minWidth": "110px",
        },
    )


# ---------------------------------------------------------------------------
# LAYOUT
# ---------------------------------------------------------------------------

layout = html.Div(
    [
        html.Div(
            [
                make_page_header(
                    "Disease Network Explorer",
                    "Browse curated metabolite interaction networks for common diseases, "
                    "drawn from the CoreMet knowledge graph.",
                    [("Home", "/home"), ("Disease", None)],
                ),

                dbc.Row(
                    [
                        # ========== SIDEBAR (280 px) ==========
                        dbc.Col(
                            html.Div(
                                [
                                    # Category filter
                                    html.H6(
                                        [html.I(className="fas fa-tags me-2"), "Category"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dcc.Dropdown(
                                        id="disease-category-filter",
                                        options=[{"label": c, "value": c} for c in CATEGORY_OPTIONS],
                                        placeholder="All categories",
                                        clearable=True,
                                        className="mb-3",
                                        style={"fontSize": "0.85rem"},
                                    ),

                                    # Tissue filter
                                    html.H6(
                                        [html.I(className="fas fa-lungs me-2"), "Tissue / Organ"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dcc.Dropdown(
                                        id="disease-tissue-filter",
                                        options=[{"label": t, "value": t} for t in TISSUE_OPTIONS],
                                        placeholder="All tissues",
                                        clearable=True,
                                        className="mb-3",
                                        style={"fontSize": "0.85rem"},
                                    ),

                                    # Disease selector
                                    html.H6(
                                        [html.I(className="fas fa-heartbeat me-2"), "Disease"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dbc.Select(
                                        id="disease-select",
                                        options=[{"label": "-- Select a disease --", "value": ""}]
                                        + DISEASE_OPTIONS,
                                        value="hcc",
                                        className="mb-3",
                                    ),

                                    # Organism note
                                    html.Div(
                                        [
                                            html.I(className="fas fa-user me-1",
                                                   style={"color": "var(--cm-metabolite)"}),
                                            html.Em("Homo sapiens"),
                                            html.Span(", human disease data",
                                                       style={"color": "var(--cm-text-secondary)",
                                                              "fontSize": "0.8rem", "marginLeft": "4px"}),
                                        ],
                                        className="mb-3",
                                        style={"fontSize": "0.85rem"},
                                    ),

                                    # Confidence slider
                                    html.H6(
                                        [html.I(className="fas fa-sliders-h me-2"), "Confidence Threshold"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dcc.Slider(
                                        id="disease-confidence",
                                        min=0.0,
                                        max=1.0,
                                        step=0.05,
                                        value=0.30,
                                        marks={
                                            0: "0",
                                            0.25: "0.25",
                                            0.5: "0.5",
                                            0.75: "0.75",
                                            1: "1",
                                        },
                                        tooltip={"placement": "bottom", "always_visible": True},
                                    ),

                                    html.Hr(className="my-3"),

                                    # Network stats panel
                                    html.H6(
                                        [html.I(className="fas fa-chart-bar me-2"), "Network Stats"],
                                        className="fw-semibold mb-2",
                                    ),
                                    html.Div(id="disease-sidebar-stats", children=[
                                        html.P("Select a disease to view statistics.",
                                               style={"fontSize": "0.85rem", "color": "#a0aec0"}),
                                    ]),

                                    html.Hr(className="my-3"),

                                    # Export buttons
                                    html.H6(
                                        [html.I(className="fas fa-download me-2"), "Export"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dbc.ButtonGroup(
                                        [
                                            dbc.Button(
                                                [html.I(className="fas fa-file-csv me-1"), "CSV"],
                                                id="disease-export-csv",
                                                className="cm-btn-secondary",
                                                size="sm",
                                            ),
                                            dbc.Button(
                                                [html.I(className="fas fa-file-code me-1"), "JSON"],
                                                id="disease-export-json",
                                                className="cm-btn-secondary",
                                                size="sm",
                                            ),
                                            dbc.Button(
                                                [html.I(className="fas fa-project-diagram me-1"), "GraphML"],
                                                id="disease-export-graphml",
                                                className="cm-btn-secondary",
                                                size="sm",
                                            ),
                                            dbc.Button(
                                                [html.I(className="fas fa-image me-1"), "PNG"],
                                                id="disease-export-png",
                                                className="cm-btn-secondary",
                                                size="sm",
                                            ),
                                        ],
                                        className="d-flex flex-wrap",
                                    ),
                                    dcc.Download(id="disease-download"),
                                ],
                                className="cm-card",
                                style={"position": "sticky", "top": "80px"},
                            ),
                            width=3,
                            style={"minWidth": "280px", "maxWidth": "300px"},
                        ),

                        # ========== MAIN CONTENT ==========
                        dbc.Col(
                            [
                                # Hidden stores
                                dcc.Store(id="disease-data-store"),
                                dcc.Store(id="disease-elements-store"),

                                # Network graph
                                html.Div(
                                    [
                                        html.Div(
                                            [
                                                html.H5(
                                                    [html.I(className="fas fa-project-diagram me-2"),
                                                     "Interaction Network"],
                                                    className="cm-card-title d-inline",
                                                ),
                                                html.Span(
                                                    id="disease-graph-badge",
                                                    className="ms-2",
                                                    style={"fontSize": "0.8rem"},
                                                ),
                                            ],
                                            className="d-flex align-items-center mb-2",
                                        ),
                                        html.Div(
                                            [
                                                html.Span(
                                                    [
                                                        html.Span(
                                                            "",
                                                            style={
                                                                "display": "inline-block",
                                                                "width": "12px",
                                                                "height": "12px",
                                                                "borderRadius": "50%",
                                                                "background": "#3182ce",
                                                                "marginRight": "4px",
                                                            },
                                                        ),
                                                        "Metabolite",
                                                    ],
                                                    className="me-3",
                                                    style={"fontSize": "0.8rem"},
                                                ),
                                                html.Span(
                                                    [
                                                        html.Span(
                                                            "",
                                                            style={
                                                                "display": "inline-block",
                                                                "width": "12px",
                                                                "height": "12px",
                                                                "borderRadius": "4px",
                                                                "background": "#e53e3e",
                                                                "marginRight": "4px",
                                                            },
                                                        ),
                                                        "Protein",
                                                    ],
                                                    style={"fontSize": "0.8rem"},
                                                ),
                                            ],
                                            className="mb-2",
                                        ),
                                        cyto.Cytoscape(
                                            id="disease-cytoscape",
                                            layout={"name": "cose", "animate": True, "padding": 30},
                                            stylesheet=CYTOSCAPE_STYLESHEET,
                                            elements=[],
                                            style={"width": "100%", "height": "520px", "border": "1px solid #e2e8f0",
                                                    "borderRadius": "8px", "background": "#ffffff"},
                                            responsive=True,
                                        ),
                                    ],
                                    className="cm-card mb-3",
                                ),

                                # Stats badges row
                                dcc.Loading(
                                    id="disease-stats-loading",
                                    type="circle",
                                    color="#1a365d",
                                    children=html.Div(id="disease-stats-row", className="mb-3"),
                                ),

                                # Node detail panel (hidden by default)
                                html.Div(
                                    id="disease-node-detail",
                                    style={"display": "none"},
                                    className="cm-card mb-3",
                                ),

                                # Hub tables, two columns
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            html.Div(
                                                [
                                                    html.H5(
                                                        [html.I(className="fas fa-flask me-2"),
                                                         "Metabolite Hubs"],
                                                        className="cm-card-title",
                                                    ),
                                                    html.Div(id="disease-met-hubs"),
                                                ],
                                                className="cm-card",
                                            ),
                                            md=6,
                                        ),
                                        dbc.Col(
                                            html.Div(
                                                [
                                                    html.H5(
                                                        [html.I(className="fas fa-microscope me-2"),
                                                         "Protein Hubs"],
                                                        className="cm-card-title",
                                                    ),
                                                    html.Div(id="disease-prot-hubs"),
                                                ],
                                                className="cm-card",
                                            ),
                                            md=6,
                                        ),
                                    ],
                                    className="g-3 mb-3",
                                ),

                                # Pathway enrichment bar chart
                                html.Div(
                                    [
                                        html.H5(
                                            [html.I(className="fas fa-route me-2"),
                                             "Pathway Enrichment"],
                                            className="cm-card-title",
                                        ),
                                        dcc.Graph(
                                            id="disease-enrichment-chart",
                                            config={"displaylogo": False},
                                            style={"height": "400px"},
                                        ),
                                    ],
                                    className="cm-card mb-3",
                                ),

                                # Predictions DataTable
                                html.Div(
                                    [
                                        html.H5(
                                            [html.I(className="fas fa-table me-2"),
                                             "Predictions"],
                                            className="cm-card-title",
                                        ),
                                        html.Div(id="disease-predictions-table"),
                                    ],
                                    className="cm-card mb-3",
                                ),
                            ],
                        ),
                    ],
                    className="g-3",
                ),
            ],
            className="cm-page-container",
        ),

    ]
)


# ---------------------------------------------------------------------------
# CALLBACKS
# ---------------------------------------------------------------------------


@callback(
    Output("disease-select", "options"),
    Output("disease-tissue-filter", "options"),
    Input("disease-category-filter", "value"),
    Input("disease-tissue-filter", "value"),
)
def filter_disease_dropdown(category, tissue):
    """Filter the disease dropdown and tissue options based on category/tissue selection."""
    filtered = DISEASE_REGISTRY.items()
    if category:
        filtered = [(k, v) for k, v in filtered if v["category"] == category]
    # Build tissue options from current category-filtered set
    available_tissues = sorted({v["tissue"] for _, v in filtered})
    tissue_opts = [{"label": t, "value": t} for t in available_tissues]
    if tissue:
        filtered = [(k, v) for k, v in filtered if v["tissue"] == tissue]
    opts = [{"label": "-- Select a disease --", "value": ""}] + [
        {"label": v["label"], "value": k}
        for k, v in sorted(filtered, key=lambda x: x[1]["label"])
    ]
    return opts, tissue_opts


@callback(
    Output("disease-data-store", "data"),
    Output("disease-elements-store", "data"),
    Input("disease-select", "value"),
    Input("disease-confidence", "value"),
)
def load_disease_data(disease, confidence):
    """Load disease data and build Cytoscape elements whenever disease or threshold changes."""
    if not disease:
        return no_update, no_update

    from app.services.disease_service import DiseaseService

    svc = DiseaseService()
    data = svc.get_disease_data(disease)
    if not data:
        return no_update, no_update

    elements = svc.get_cytoscape_elements(disease, confidence_threshold=confidence or 0.3)

    # Serialise DataFrames
    store = {
        "disease": disease,
        "confidence": confidence,
        "network_stats": data.get("network_stats", {}),
    }
    # Predictions and enrichment stored as JSON strings
    preds = data.get("predictions", pd.DataFrame())
    if not preds.empty:
        store["predictions"] = preds.to_json(orient="split")
    enrich = data.get("enrichment", pd.DataFrame())
    if not enrich.empty:
        store["enrichment"] = enrich.to_json(orient="split")

    # Hub tables
    hubs = svc.get_hub_tables(disease)
    store["metabolite_hubs"] = hubs.get("metabolite_hubs", [])
    store["protein_hubs"] = hubs.get("protein_hubs", [])

    return store, elements


@callback(
    Output("disease-cytoscape", "elements"),
    Output("disease-graph-badge", "children"),
    Input("disease-elements-store", "data"),
)
def update_cytoscape(elements):
    """Push elements into the Cytoscape component."""
    if not elements:
        return [], ""
    n_nodes = sum(1 for e in elements if "source" not in e.get("data", {}))
    n_edges = sum(1 for e in elements if "source" in e.get("data", {}))
    badge = html.Span(
        f"{n_nodes} nodes, {n_edges} edges",
        className="cm-badge cm-badge-source",
    )
    return elements, badge


@callback(
    Output("disease-stats-row", "children"),
    Output("disease-sidebar-stats", "children"),
    Input("disease-data-store", "data"),
)
def update_stats(store):
    """Render stats badges and sidebar stats."""
    if not store:
        return "", html.P("No data.", style={"fontSize": "0.85rem", "color": "#a0aec0"})

    stats = store.get("network_stats", {})

    # Main stats row
    badges = dbc.Row(
        [
            dbc.Col(_make_stat_badge("fas fa-flask", stats.get("n_metabolites", 0), "Metabolites", "#3182ce"), width="auto"),
            dbc.Col(_make_stat_badge("fas fa-microscope", stats.get("n_proteins", 0), "Proteins", "#e53e3e"), width="auto"),
            dbc.Col(_make_stat_badge("fas fa-link", stats.get("n_edges", 0), "Edges", "#38a169"), width="auto"),
            dbc.Col(_make_stat_badge("fas fa-braille", stats.get("avg_degree", 0), "Avg Degree", "#d69e2e"), width="auto"),
            dbc.Col(_make_stat_badge("fas fa-compress-arrows-alt", stats.get("density", 0), "Density", "#805ad5"), width="auto"),
        ],
        className="g-2 justify-content-start",
    )

    # Sidebar stats
    sidebar = html.Div(
        [
            html.Div([html.Strong("Nodes: "), html.Span(str(stats.get("n_nodes", 0)))], style={"fontSize": "0.85rem"}),
            html.Div([html.Strong("Edges: "), html.Span(str(stats.get("n_edges", 0)))], style={"fontSize": "0.85rem"}),
            html.Div([html.Strong("Density: "), html.Span(f"{stats.get('density', 0):.4f}")], style={"fontSize": "0.85rem"}),
            html.Div([html.Strong("Components: "), html.Span(str(stats.get("n_components", 0)))], style={"fontSize": "0.85rem"}),
            html.Div([html.Strong("Avg Degree: "), html.Span(str(stats.get("avg_degree", 0)))], style={"fontSize": "0.85rem"}),
            html.Div([html.Strong("Threshold: "), html.Span(str(stats.get("score_threshold", 0.3)))], style={"fontSize": "0.85rem"}),
        ]
    )

    return badges, sidebar


@callback(
    Output("disease-met-hubs", "children"),
    Output("disease-prot-hubs", "children"),
    Input("disease-data-store", "data"),
)
def update_hub_tables(store):
    """Render hub tables."""
    if not store:
        empty = html.P("No data.", className="text-muted", style={"fontSize": "0.85rem"})
        return empty, empty

    def build_table(hubs, type_label):
        if not hubs:
            return html.P(f"No {type_label} hubs available.", className="text-muted", style={"fontSize": "0.85rem"})
        return dash_table.DataTable(
            data=hubs,
            columns=[
                {"name": type_label, "id": "node"},
                {"name": "Degree", "id": "degree"},
            ],
            sort_action="native",
            page_size=10,
            style_cell={"textAlign": "left", "fontSize": "0.85rem", "padding": "6px 10px", "fontFamily": "Arial, Helvetica, sans-serif"},
            style_header={"fontWeight": "600", "backgroundColor": "#1a365d", "color": "white", "fontFamily": "Arial, Helvetica, sans-serif"},
            style_data_conditional=[{"if": {"row_index": "odd"}, "backgroundColor": "#f7fafc"}],
        )

    met_table = build_table(store.get("metabolite_hubs", []), "Metabolite")
    prot_table = build_table(store.get("protein_hubs", []), "Protein")
    return met_table, prot_table


@callback(
    Output("disease-enrichment-chart", "figure"),
    Input("disease-data-store", "data"),
)
def update_enrichment_chart(store):
    """Render pathway enrichment horizontal bar chart."""
    fig = go.Figure()
    fig.update_layout(
        margin=dict(l=250, r=40, t=20, b=40),
        font=dict(family="Arial", size=11),
        plot_bgcolor="white",
        xaxis_title="Fold Enrichment",
    )

    if not store or "enrichment" not in store:
        fig.update_layout(
            annotations=[
                dict(text="No enrichment data available", x=0.5, y=0.5,
                     xref="paper", yref="paper", showarrow=False, font=dict(size=14, color="#a0aec0"))
            ]
        )
        return fig

    enrich_df = pd.read_json(io.StringIO(store["enrichment"]), orient="split")

    # Handle both column naming conventions
    col_map = {}
    for col in enrich_df.columns:
        lc = col.lower()
        if "fold" in lc and "enrich" in lc:
            col_map["fold"] = col
        if lc == "fdr":
            col_map["fdr"] = col
        if "pathway" in lc and "name" in lc:
            col_map["name"] = col

    fold_col = col_map.get("fold", "Fold_Enrichment")
    fdr_col = col_map.get("fdr", "FDR")
    name_col = col_map.get("name", "Pathway_Name")

    if fold_col not in enrich_df.columns or enrich_df.empty:
        fig.update_layout(
            annotations=[
                dict(text="No enriched pathways found", x=0.5, y=0.5,
                     xref="paper", yref="paper", showarrow=False, font=dict(size=14, color="#a0aec0"))
            ]
        )
        return fig

    # Filter significant and take top 20
    if fdr_col in enrich_df.columns:
        sig = enrich_df[enrich_df[fdr_col] <= 0.25].head(20)
    else:
        sig = enrich_df.head(20)

    if sig.empty:
        fig.update_layout(
            annotations=[
                dict(text="No significantly enriched pathways (FDR ≤ 0.25)", x=0.5, y=0.5,
                     xref="paper", yref="paper", showarrow=False, font=dict(size=14, color="#a0aec0"))
            ]
        )
        return fig

    colors = [
        "#38a169" if fdr < 0.01 else "#3182ce" if fdr < 0.05 else "#d69e2e" if fdr < 0.1 else "#a0aec0"
        for fdr in sig[fdr_col]
    ]

    fig.add_trace(
        go.Bar(
            y=sig[name_col].str[:50],
            x=sig[fold_col],
            orientation="h",
            marker_color=colors,
            text=[f"FDR={f:.2e}" for f in sig[fdr_col]],
            textposition="outside",
        )
    )
    fig.update_layout(height=max(300, len(sig) * 30))
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

    return fig


@callback(
    Output("disease-predictions-table", "children"),
    Input("disease-data-store", "data"),
)
def update_predictions_table(store):
    """Render full predictions DataTable."""
    if not store or "predictions" not in store:
        return html.P("Select a disease to view predictions.",
                       className="text-muted", style={"fontSize": "0.85rem"})

    preds_df = pd.read_json(io.StringIO(store["predictions"]), orient="split")
    if preds_df.empty:
        return html.P("No predictions available.", className="text-muted", style={"fontSize": "0.85rem"})

    # Round scores for display
    if "Prediction Score" in preds_df.columns:
        preds_df["Prediction Score"] = pd.to_numeric(
            preds_df["Prediction Score"], errors="coerce"
        ).round(5)

    return dash_table.DataTable(
        data=preds_df.to_dict("records"),
        columns=[{"name": c, "id": c} for c in preds_df.columns],
        sort_action="native",
        filter_action="native",
        page_size=20,
        style_cell={
            "textAlign": "left",
            "fontSize": "0.85rem",
            "padding": "6px 10px",
            "maxWidth": "200px",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
        },
        style_header={
            "fontWeight": "600",
            "backgroundColor": "#1a365d",
            "color": "white",
        },
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#f7fafc"},
        ],
    )


@callback(
    Output("disease-node-detail", "children"),
    Output("disease-node-detail", "style"),
    Input("disease-cytoscape", "tapNodeData"),
    State("disease-select", "value"),
    prevent_initial_call=True,
)
def show_node_detail(tap_data, disease):
    """Show detail panel when a node is clicked in the Cytoscape graph."""
    if not tap_data or not disease:
        return "", {"display": "none"}

    from app.services.disease_service import DiseaseService

    svc = DiseaseService()
    node_id = tap_data.get("id", "")
    detail = svc.get_node_detail(disease, node_id)

    node_type = detail.get("node_type", "unknown")
    type_color = "#3182ce" if node_type == "metabolite" else "#e53e3e"
    type_icon = "fas fa-flask" if node_type == "metabolite" else "fas fa-microscope"

    # Header
    header = html.Div(
        [
            html.I(className=f"{type_icon} me-2", style={"color": type_color}),
            html.Span(detail.get("name", node_id), style={"fontWeight": "700", "fontSize": "1.1rem"}),
            html.Span(
                f"  ({node_type})",
                style={"color": "#718096", "fontSize": "0.85rem", "marginLeft": "8px"},
            ),
        ],
        className="mb-2",
    )

    # Identifiers and links
    id_items = []
    if "hmdb_id" in detail:
        id_items.append(html.Div([html.Strong("HMDB: "), html.Span(detail["hmdb_id"])]))
    if "uniprot_id" in detail:
        id_items.append(html.Div([html.Strong("UniProt: "), html.Span(detail["uniprot_id"])]))
    if "gene" in detail:
        id_items.append(html.Div([html.Strong("Gene: "), html.Span(detail["gene"])]))
    if "external_url" in detail:
        id_items.append(
            html.Div(
                html.A(
                    [html.I(className="fas fa-external-link-alt me-1"), "View in external database"],
                    href=detail["external_url"],
                    target="_blank",
                    style={"fontSize": "0.85rem"},
                )
            )
        )
    id_items.append(html.Div([html.Strong("Degree: "), html.Span(str(detail.get("degree", 0)))]))

    # Top interactions mini-table
    interactions = detail.get("top_interactions", [])
    interactions_section = ""
    if interactions:
        interactions_section = html.Div(
            [
                html.H6("Top Interactions", className="fw-semibold mt-2 mb-1"),
                dash_table.DataTable(
                    data=interactions,
                    columns=[
                        {"name": "Partner", "id": "partner"},
                        {"name": "Score", "id": "score"},
                    ],
                    page_size=5,
                    style_cell={"textAlign": "left", "fontSize": "0.8rem", "padding": "4px 8px"},
                    style_header={"fontWeight": "600", "backgroundColor": "#edf2f7"},
                ),
            ]
        )

    content = html.Div([header, html.Div(id_items, style={"fontSize": "0.85rem"}), interactions_section])
    return content, {"display": "block"}


@callback(
    Output("disease-download", "data"),
    Input("disease-export-csv", "n_clicks"),
    Input("disease-export-json", "n_clicks"),
    Input("disease-export-graphml", "n_clicks"),
    Input("disease-export-png", "n_clicks"),
    State("disease-data-store", "data"),
    prevent_initial_call=True,
)
def handle_export(csv_clicks, json_clicks, graphml_clicks, png_clicks, store):
    """Handle export button clicks."""
    if not store or "predictions" not in store:
        return no_update

    trigger = ctx.triggered_id
    preds_df = pd.read_json(io.StringIO(store["predictions"]), orient="split")
    disease = store.get("disease", "disease")

    if trigger == "disease-export-csv":
        return dict(content=preds_df.to_csv(index=False), filename=f"{disease}_predictions.csv")

    elif trigger == "disease-export-json":
        return dict(content=preds_df.to_json(orient="records", indent=2), filename=f"{disease}_predictions.json")

    elif trigger == "disease-export-graphml":
        try:
            from app.services.export_service import export_graphml
            content = export_graphml(preds_df)
            return dict(content=content, filename=f"{disease}_network.graphml")
        except Exception:
            return dict(content=preds_df.to_csv(index=False), filename=f"{disease}_predictions.csv")

    return no_update


# -- Navbar shared callbacks -------------------------------------------------





# Export for main.py routing (navbar/footer handled globally)
page_content = layout
