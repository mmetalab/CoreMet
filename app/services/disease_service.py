"""
Disease analysis service — loads pre-computed disease MPI data and builds
Cytoscape-compatible network elements, hub tables, enrichment data,
and summary statistics.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any

import pandas as pd

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent.parent
DISEASE_DATA_DIR = PROJECT_ROOT / "data" / "mpidatabase" / "disease_mpi"

# Legacy disease catalogue — original 40 diseases (see pages/disease.py DISEASE_REGISTRY for full 130)
DISEASE_CATALOGUE = {
    # Cancers (20)
    "aml_leukemia": "Acute Myeloid Leukemia",
    "bladder_cancer": "Bladder Cancer",
    "breast_cancer": "Breast Cancer",
    "cervical_cancer": "Cervical Cancer",
    "cholangiocarcinoma": "Cholangiocarcinoma",
    "colorectal_cancer": "Colorectal Cancer",
    "endometrial_cancer": "Endometrial Cancer",
    "esophageal_cancer": "Esophageal Cancer",
    "gastric_cancer": "Gastric Cancer",
    "glioblastoma": "Glioblastoma",
    "hcc": "Hepatocellular Carcinoma (HCC)",
    "lung_cancer": "Lung Cancer",
    "lymphoma": "Lymphoma",
    "melanoma": "Melanoma",
    "oral_cancer": "Oral Squamous Cell Carcinoma",
    "ovarian_cancer": "Ovarian Cancer",
    "pancreatic_cancer": "Pancreatic Cancer",
    "thyroid_cancer": "Papillary Thyroid Cancer",
    "prostate_cancer": "Prostate Cancer",
    "renal_cell_carcinoma": "Renal Cell Carcinoma",
    # Neurodegenerative / Neuropsychiatric (6)
    "alzheimers": "Alzheimer's Disease",
    "als": "Amyotrophic Lateral Sclerosis (ALS)",
    "depression": "Major Depressive Disorder",
    "multiple_sclerosis": "Multiple Sclerosis",
    "parkinsons": "Parkinson's Disease",
    "schizophrenia": "Schizophrenia",
    # Metabolic (4)
    "metabolic_syndrome": "Metabolic Syndrome",
    "nafld": "Non-Alcoholic Fatty Liver Disease",
    "obesity": "Obesity",
    "t2_diabetes": "Type 2 Diabetes",
    # Cardiovascular (3)
    "atherosclerosis": "Atherosclerosis",
    "heart_failure": "Heart Failure",
    "hypertension": "Hypertension",
    # Autoimmune / Inflammatory (3)
    "ibd": "Inflammatory Bowel Disease",
    "rheumatoid_arthritis": "Rheumatoid Arthritis",
    "lupus": "Systemic Lupus Erythematosus",
    # Respiratory (2)
    "asthma": "Asthma",
    "copd": "Chronic Obstructive Pulmonary Disease",
    # Other (2)
    "covid19": "COVID-19",
    "chronic_kidney_disease": "Chronic Kidney Disease",
}


class DiseaseService:
    """Service for loading and querying pre-computed disease MPI datasets."""

    def __init__(self, data_dir: Optional[Path] = None):
        self._data_dir = data_dir or DISEASE_DATA_DIR
        self._cache: Dict[str, Dict[str, Any]] = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def get_diseases(self) -> List[Dict[str, str]]:
        """Return list of available diseases with label and key.

        Only diseases whose output directory actually exists on disk
        are included; the rest are listed but marked as unavailable.
        """
        results = []
        for key, label in DISEASE_CATALOGUE.items():
            available = (self._data_dir / key).is_dir()
            results.append({"key": key, "label": label, "available": available})
        return results

    def get_disease_data(self, disease_name: str) -> Dict[str, Any]:
        """Load all pre-computed artefacts for *disease_name*.

        Returns a dict with keys:
            metabolites   – DataFrame
            proteins      – DataFrame
            predictions   – DataFrame
            network_stats – dict
            enrichment    – DataFrame
        """
        if disease_name in self._cache:
            return self._cache[disease_name]

        disease_dir = self._data_dir / disease_name
        if not disease_dir.is_dir():
            logger.warning(f"Disease directory not found: {disease_dir}")
            return {}

        data: Dict[str, Any] = {}

        # Metabolites
        met_path = disease_dir / "metabolites.csv"
        data["metabolites"] = pd.read_csv(met_path) if met_path.exists() else pd.DataFrame()

        # Proteins
        prot_path = disease_dir / "proteins.csv"
        data["proteins"] = pd.read_csv(prot_path) if prot_path.exists() else pd.DataFrame()

        # Predictions
        pred_path = disease_dir / "predictions.csv"
        data["predictions"] = pd.read_csv(pred_path) if pred_path.exists() else pd.DataFrame()

        # Network stats
        stats_path = disease_dir / "network_stats.json"
        if stats_path.exists():
            with open(stats_path) as f:
                data["network_stats"] = json.load(f)
        else:
            data["network_stats"] = {}

        # Pathway enrichment
        enrich_path = disease_dir / "pathway_enrichment.csv"
        data["enrichment"] = pd.read_csv(enrich_path) if enrich_path.exists() else pd.DataFrame()

        self._cache[disease_name] = data
        return data

    def get_cytoscape_elements(
        self,
        disease_name: str,
        confidence_threshold: float = 0.3,
    ) -> List[Dict]:
        """Build Cytoscape.js-compatible elements list.

        Metabolites  -> blue circles         (classes: "metabolite")
        Proteins     -> red rounded-squares  (classes: "protein")
        Edges        -> score >= *confidence_threshold*
        """
        data = self.get_disease_data(disease_name)
        if not data:
            return []

        predictions = data.get("predictions", pd.DataFrame())
        if predictions.empty:
            # Fallback: show metabolite and protein nodes without edges
            logger.warning(f"No predictions for {disease_name}, building node-only elements")
            elements: List[Dict] = []
            metabolites = data.get("metabolites", pd.DataFrame())
            proteins = data.get("proteins", pd.DataFrame())
            for _, row in metabolites.iterrows():
                node_data = {
                    "id": str(row.get("name", row.get("hmdb_id", ""))),
                    "label": str(row.get("name", row.get("hmdb_id", ""))),
                    "node_type": "metabolite",
                    "degree": 0,
                }
                if pd.notna(row.get("hmdb_id")):
                    node_data["hmdb_id"] = str(row["hmdb_id"])
                elements.append({"data": node_data, "classes": "metabolite"})
            for _, row in proteins.iterrows():
                node_data = {
                    "id": str(row.get("name", row.get("uniprot_id", ""))),
                    "label": str(row.get("name", row.get("uniprot_id", ""))),
                    "node_type": "protein",
                    "degree": 0,
                }
                if pd.notna(row.get("uniprot_id")):
                    node_data["uniprot_id"] = str(row["uniprot_id"])
                if pd.notna(row.get("gene")):
                    node_data["gene"] = str(row["gene"])
                elements.append({"data": node_data, "classes": "protein"})
            return elements

        # Determine column names (handle both service-generated and seed-generated schemas)
        met_col = "Metabolite"
        prot_col = "Protein"
        score_col = "Prediction Score"
        hmdb_col = "HMDB_ID" if "HMDB_ID" in predictions.columns else None
        uniprot_col = "Uniprot_ID" if "Uniprot_ID" in predictions.columns else None
        gene_col = "Gene" if "Gene" in predictions.columns else None

        if score_col not in predictions.columns:
            return []

        scores = pd.to_numeric(predictions[score_col], errors="coerce")
        filtered = predictions[scores >= confidence_threshold].copy()

        if filtered.empty:
            return []

        elements: List[Dict] = []
        met_ids_seen = set()
        prot_ids_seen = set()

        # Degree accumulators
        degree: Dict[str, int] = {}

        # First pass: count degrees
        for _, row in filtered.iterrows():
            met = str(row[met_col])
            prot = str(row[prot_col])
            degree[met] = degree.get(met, 0) + 1
            degree[prot] = degree.get(prot, 0) + 1

        # Build metabolite nodes
        for _, row in filtered.iterrows():
            met = str(row[met_col])
            if met in met_ids_seen:
                continue
            met_ids_seen.add(met)
            node_data = {
                "id": met,
                "label": met,
                "node_type": "metabolite",
                "degree": degree.get(met, 0),
            }
            if hmdb_col and pd.notna(row.get(hmdb_col)):
                node_data["hmdb_id"] = str(row[hmdb_col])
            elements.append({"data": node_data, "classes": "metabolite"})

        # Build protein nodes
        for _, row in filtered.iterrows():
            prot = str(row[prot_col])
            if prot in prot_ids_seen:
                continue
            prot_ids_seen.add(prot)
            node_data = {
                "id": prot,
                "label": prot,
                "node_type": "protein",
                "degree": degree.get(prot, 0),
            }
            if uniprot_col and pd.notna(row.get(uniprot_col)):
                node_data["uniprot_id"] = str(row[uniprot_col])
            if gene_col and pd.notna(row.get(gene_col)):
                node_data["gene"] = str(row[gene_col])
            elements.append({"data": node_data, "classes": "protein"})

        # Build edges
        for _, row in filtered.iterrows():
            met = str(row[met_col])
            prot = str(row[prot_col])
            score = float(row[score_col])
            edge_data = {
                "id": f"{met}--{prot}",
                "source": met,
                "target": prot,
                "score": round(score, 5),
            }
            elements.append({"data": edge_data})

        return elements

    def get_hub_tables(self, disease_name: str, top_n: int = 10) -> Dict[str, List[Dict]]:
        """Return top metabolite and protein hubs from network_stats.

        Falls back to computing from predictions if stats file is missing.
        """
        data = self.get_disease_data(disease_name)
        stats = data.get("network_stats", {})

        met_hubs = stats.get("metabolite_hubs", [])[:top_n]
        prot_hubs = stats.get("protein_hubs", [])[:top_n]

        # Fallback: derive from predictions
        if not met_hubs or not prot_hubs:
            predictions = data.get("predictions", pd.DataFrame())
            if not predictions.empty and "Prediction Score" in predictions.columns:
                scores = pd.to_numeric(predictions["Prediction Score"], errors="coerce")
                above = predictions[scores >= 0.3]
                if not above.empty:
                    if not met_hubs:
                        met_counts = above["Metabolite"].value_counts().head(top_n)
                        met_hubs = [{"node": n, "degree": int(d)} for n, d in met_counts.items()]
                    if not prot_hubs:
                        prot_counts = above["Protein"].value_counts().head(top_n)
                        prot_hubs = [{"node": n, "degree": int(d)} for n, d in prot_counts.items()]

        return {"metabolite_hubs": met_hubs, "protein_hubs": prot_hubs}

    def get_enrichment_data(self, disease_name: str) -> pd.DataFrame:
        """Return pathway enrichment results for *disease_name*."""
        data = self.get_disease_data(disease_name)
        return data.get("enrichment", pd.DataFrame())

    def get_network_stats(self, disease_name: str) -> Dict[str, Any]:
        """Return network statistics dict."""
        data = self.get_disease_data(disease_name)
        return data.get("network_stats", {})

    def get_node_detail(self, disease_name: str, node_id: str) -> Dict[str, Any]:
        """Return detail for a single node: name, IDs, degree, top interactions."""
        data = self.get_disease_data(disease_name)
        predictions = data.get("predictions", pd.DataFrame())

        if predictions.empty:
            return {"id": node_id}

        # Find all edges involving this node
        is_met = predictions["Metabolite"].astype(str) == node_id
        is_prot = predictions["Protein"].astype(str) == node_id

        detail: Dict[str, Any] = {"id": node_id}

        if is_met.any():
            detail["node_type"] = "metabolite"
            detail["name"] = node_id
            edges = predictions[is_met].copy()
            # Find HMDB ID if available
            if "HMDB_ID" in edges.columns:
                hmdb = edges["HMDB_ID"].dropna().unique()
                if len(hmdb) > 0:
                    detail["hmdb_id"] = str(hmdb[0])
                    detail["external_url"] = f"https://hmdb.ca/metabolites/{hmdb[0]}"
            partner_col = "Protein"
        elif is_prot.any():
            detail["node_type"] = "protein"
            detail["name"] = node_id
            edges = predictions[is_prot].copy()
            if "Uniprot_ID" in edges.columns:
                uid = edges["Uniprot_ID"].dropna().unique()
                if len(uid) > 0:
                    detail["uniprot_id"] = str(uid[0])
                    detail["external_url"] = f"https://www.uniprot.org/uniprot/{uid[0]}"
            if "Gene" in edges.columns:
                gene = edges["Gene"].dropna().unique()
                if len(gene) > 0:
                    detail["gene"] = str(gene[0])
            partner_col = "Metabolite"
        else:
            return detail

        # Degree and top interactions
        scores = pd.to_numeric(edges["Prediction Score"], errors="coerce")
        edges = edges.assign(_score=scores).sort_values("_score", ascending=False)
        detail["degree"] = len(edges)
        detail["top_interactions"] = [
            {"partner": str(row[partner_col]), "score": round(float(row["_score"]), 5)}
            for _, row in edges.head(10).iterrows()
        ]

        return detail
