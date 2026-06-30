"""
Disease analysis service — loads pre-computed disease MPI data and builds
Cytoscape-compatible network elements, hub tables, enrichment data,
and summary statistics.
"""

import json
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Any

import pandas as pd

logger = logging.getLogger(__name__)

PROJECT_ROOT = Path(__file__).parent.parent.parent
DISEASE_DATA_DIR = PROJECT_ROOT / "data" / "mpidatabase" / "disease_mpi"
MAX_RELEASE_DISEASE_EDGES = 600

# Legacy disease catalogue — original 40 diseases (see pages/disease.py DISEASE_REGISTRY for the UI registry)
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

_DISEASE_ALIASES = {
    "alzheimers": ["Alzheimer Disease", "Alzheimer's Disease"],
    "breast_cancer": ["Breast Neoplasms", "Breast Cancer"],
    "colorectal_cancer": ["Colorectal cancer", "Colorectal Neoplasms"],
    "hcc": ["Hepatocellular carcinoma", "Carcinoma, Hepatocellular", "Liver Neoplasms"],
    "heart_failure": ["Heart Failure"],
    "hypertension": ["Hypertension"],
    "lung_cancer": ["Lung Neoplasms", "Lung Cancer"],
    "obesity": ["Obesity"],
    "parkinsons": ["Parkinson Disease", "Parkinsonian Disorders"],
    "schizophrenia": ["Schizophrenia"],
    "t2_diabetes": ["Type 2 diabetes mellitus", "Diabetes Mellitus, Type 2", "Type 2 Diabetes"],
    "type1_diabetes": ["Diabetes Mellitus, Type 1", "Type 1 Diabetes Mellitus"],
    "ulcerative_colitis": ["Ulcerative colitis", "Colitis, Ulcerative"],
    "crohns_disease": ["Crohn's disease", "Crohn Disease"],
}


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none", "null", "<na>", "n/a"} else text


def _norm(value: Any) -> str:
    text = _clean_text(value).lower()
    keep = [ch if ch.isalnum() else " " for ch in text]
    return " ".join("".join(keep).split())


def _valid_hmdb(value: Any) -> bool:
    text = _clean_text(value).upper()
    return text.startswith("HMDB") and text[4:].isdigit()


def _score_col(df: pd.DataFrame) -> str:
    return "Network Score" if "Network Score" in df.columns else "Prediction Score"


class DiseaseService:
    """Service for loading disease MPI datasets or release-backed fallback networks."""

    def __init__(self, data_dir: Optional[Path] = None):
        self._data_dir = data_dir or DISEASE_DATA_DIR
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._availability_cache: Dict[str, bool] = {}

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

    def has_release_data(self, disease_name: str, disease_label: str = "") -> bool:
        """Return whether release MDI/MPI tables can build a disease network."""
        cache_key = f"{disease_name}|{disease_label}"
        if cache_key in self._availability_cache:
            return self._availability_cache[cache_key]
        try:
            from app.services.mdi_service import get_mdi_db
            mdi = get_mdi_db()
            available = not self._match_mdi_for_disease(mdi, disease_name, disease_label, limit=1).empty
        except Exception as exc:
            logger.debug("Disease availability check failed for %s: %s", disease_name, exc)
            available = False
        self._availability_cache[cache_key] = available
        return available

    def get_disease_data(self, disease_name: str, disease_label: str = "") -> Dict[str, Any]:
        """Load all pre-computed artefacts for *disease_name*.

        Returns a dict with keys:
            metabolites   – DataFrame
            proteins      – DataFrame
            predictions   – DataFrame
            network_stats – dict
            enrichment    – DataFrame
        """
        cache_key = f"{disease_name}|{disease_label}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        disease_dir = self._data_dir / disease_name
        if not disease_dir.is_dir():
            data = self._build_release_disease_data(disease_name, disease_label)
            self._cache[cache_key] = data
            return data

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

        self._cache[cache_key] = data
        return data

    # ------------------------------------------------------------------
    # Release-table fallback
    # ------------------------------------------------------------------

    def _disease_terms(self, disease_name: str, disease_label: str = "") -> list[str]:
        label = disease_label or DISEASE_CATALOGUE.get(disease_name, disease_name.replace("_", " "))
        terms = [label, disease_name.replace("_", " ")]
        if "(" in label:
            terms.append(label.split("(")[0].strip())
        terms.extend(_DISEASE_ALIASES.get(disease_name, []))
        if label.lower().endswith(" cancer"):
            base = label[:-7].strip()
            terms.extend([f"{base} neoplasms", f"{base} carcinoma"])
        seen = set()
        clean_terms = []
        for term in terms:
            cleaned = _clean_text(term)
            normalized = _norm(cleaned)
            if len(normalized) >= 3 and normalized not in seen:
                clean_terms.append(cleaned)
                seen.add(normalized)
        return clean_terms

    def _match_mdi_for_disease(
        self,
        mdi: pd.DataFrame,
        disease_name: str,
        disease_label: str = "",
        limit: Optional[int] = None,
    ) -> pd.DataFrame:
        if mdi.empty or "Disease_Name" not in mdi.columns:
            return pd.DataFrame()

        disease_norm = mdi["Disease_Name"].astype("string").fillna("").map(_norm)
        source_cols = [c for c in ("Source", "source") if c in mdi.columns]
        source_norm = pd.Series("", index=mdi.index, dtype="string")
        for col in source_cols:
            source_norm = source_norm.str.cat(mdi[col].astype("string").fillna("").str.lower(), sep=" ")

        mask = pd.Series(False, index=mdi.index)
        marker = disease_name.lower()
        if marker:
            mask |= source_norm.str.contains(f"coremet_disease_mpi/{marker}", na=False, regex=False)
            mask |= source_norm.str.contains(f"coremet_case_study/{marker}", na=False, regex=False)

        for term in self._disease_terms(disease_name, disease_label):
            t_norm = _norm(term)
            if not t_norm:
                continue
            mask |= disease_norm.eq(t_norm)
            mask |= disease_norm.str.contains(t_norm, na=False, regex=False)

        result = mdi.loc[mask].copy()
        if limit:
            result = result.head(limit)
        return result

    def _build_release_disease_data(self, disease_name: str, disease_label: str = "") -> Dict[str, Any]:
        """Build a disease network from release MDI associations joined to MPI edges."""
        try:
            from app.services.mdi_service import get_mdi_db
            from app.services.mpi_service import get_mpi_db

            mdi = get_mdi_db()
            mdi_hits = self._match_mdi_for_disease(mdi, disease_name, disease_label)
            if mdi_hits.empty:
                logger.warning("No release MDI records found for disease: %s", disease_name)
                return {}

            mpi = get_mpi_db()
            if mpi.empty:
                logger.warning("MPI database unavailable for disease fallback: %s", disease_name)
                return {}

            hmdb_ids = {
                _clean_text(v) for v in mdi_hits.get("HMDB_ID", pd.Series(dtype=str))
                if _valid_hmdb(v)
            }
            met_names = {
                _norm(v) for v in mdi_hits.get("Metabolite_Name", pd.Series(dtype=str))
                if _clean_text(v)
            }
            if not hmdb_ids and not met_names:
                return {}

            mask = pd.Series(False, index=mpi.index)
            if hmdb_ids and "HMDB_ID" in mpi.columns:
                mask |= mpi["HMDB_ID"].astype("string").fillna("").isin(hmdb_ids)
            if met_names and "Metabolite_Name" in mpi.columns:
                mask |= mpi["Metabolite_Name"].astype("string").fillna("").map(_norm).isin(met_names)

            mpi_hits = mpi.loc[mask].copy()
            if mpi_hits.empty:
                logger.warning("No MPI edges found for disease metabolites: %s", disease_name)
                return {}

            score_map = self._build_disease_score_map(mdi_hits)
            mpi_hits["_disease_score"] = mpi_hits.apply(
                lambda row: self._disease_score_for_row(row, score_map), axis=1
            )
            if "confidence" in mpi_hits.columns:
                mpi_conf = pd.to_numeric(mpi_hits["confidence"], errors="coerce").fillna(0.85)
            else:
                mpi_conf = pd.Series(0.85, index=mpi_hits.index)
            network_score = ((mpi_conf + mpi_hits["_disease_score"]) / 2).clip(0, 1)
            mpi_hits["_network_score"] = network_score
            mpi_hits = mpi_hits.sort_values("_network_score", ascending=False).head(MAX_RELEASE_DISEASE_EDGES)

            predictions = pd.DataFrame({
                "Metabolite": mpi_hits.get("Metabolite_Name", "").astype(str),
                "HMDB_ID": mpi_hits.get("HMDB_ID", "").astype(str),
                "Protein": mpi_hits.get("Protein_Name", "").astype(str),
                "Uniprot_ID": mpi_hits.get("Uniprot_ID", "").astype(str),
                "Gene": mpi_hits.get("Gene_Name", "").astype(str),
                "Species": mpi_hits.get("Species", "").astype(str),
                "Pathway_Name": mpi_hits.get("Pathway_Name", "").astype(str),
                "Disease Evidence": mpi_hits["_disease_score"].round(3),
                "Network Score": mpi_hits["_network_score"].round(5),
                "MPI Source": mpi_hits.get("Evidence_Source", "").astype(str),
            }).drop_duplicates()

            metabolites = predictions[["Metabolite", "HMDB_ID"]].drop_duplicates()
            proteins = predictions[["Protein", "Uniprot_ID", "Gene"]].drop_duplicates()
            stats = self._compute_network_stats(predictions)
            enrichment = self._compute_pathway_enrichment(mpi_hits)

            return {
                "metabolites": metabolites,
                "proteins": proteins,
                "predictions": predictions,
                "network_stats": stats,
                "enrichment": enrichment,
                "source": "CoreMet release MDI+MPI",
                "mdi_records": mdi_hits.head(5000),
            }
        except Exception as exc:
            logger.exception("Failed to build release disease network for %s", disease_name)
            return {}

    def _build_disease_score_map(self, mdi_hits: pd.DataFrame) -> Dict[str, float]:
        scores = pd.to_numeric(mdi_hits.get("confidence", 0.75), errors="coerce").fillna(0.75)
        score_map: Dict[str, float] = {}
        for idx, row in mdi_hits.iterrows():
            score = float(scores.loc[idx])
            hmdb = _clean_text(row.get("HMDB_ID", ""))
            name = _norm(row.get("Metabolite_Name", ""))
            if hmdb:
                score_map[hmdb] = max(score_map.get(hmdb, 0), score)
            if name:
                score_map[name] = max(score_map.get(name, 0), score)
        return score_map

    def _disease_score_for_row(self, row: pd.Series, score_map: Dict[str, float]) -> float:
        hmdb = _clean_text(row.get("HMDB_ID", ""))
        name = _norm(row.get("Metabolite_Name", ""))
        return max(score_map.get(hmdb, 0), score_map.get(name, 0), 0.75)

    def _compute_network_stats(self, predictions: pd.DataFrame) -> Dict[str, Any]:
        score_col = _score_col(predictions)
        scores = pd.to_numeric(predictions[score_col], errors="coerce").fillna(0)
        active = predictions[scores >= 0.3]
        n_met = active["Metabolite"].nunique() if "Metabolite" in active.columns else 0
        n_prot = active["Protein"].nunique() if "Protein" in active.columns else 0
        n_edges = len(active)
        n_nodes = n_met + n_prot
        density = n_edges / (n_met * n_prot) if n_met and n_prot else 0
        avg_degree = round((2 * n_edges / n_nodes), 2) if n_nodes else 0
        met_hubs = [
            {"node": str(node), "degree": int(deg)}
            for node, deg in active["Metabolite"].value_counts().head(10).items()
        ] if "Metabolite" in active.columns else []
        prot_hubs = [
            {"node": str(node), "degree": int(deg)}
            for node, deg in active["Protein"].value_counts().head(10).items()
        ] if "Protein" in active.columns else []
        return {
            "n_metabolites": int(n_met),
            "n_proteins": int(n_prot),
            "n_edges": int(n_edges),
            "n_nodes": int(n_nodes),
            "density": round(float(density), 5),
            "avg_degree": avg_degree,
            "n_components": 1 if n_edges else 0,
            "score_threshold": 0.3,
            "metabolite_hubs": met_hubs,
            "protein_hubs": prot_hubs,
        }

    def _compute_pathway_enrichment(self, selected_mpi: pd.DataFrame) -> pd.DataFrame:
        if selected_mpi.empty or "Pathway_Name" not in selected_mpi.columns:
            return pd.DataFrame()
        selected_paths = selected_mpi["Pathway_Name"].astype("string").fillna("").str.strip()
        selected_paths = selected_paths[selected_paths.ne("")]
        if selected_paths.empty:
            return pd.DataFrame()
        try:
            from app.services.mpi_service import get_mpi_db
            universe = get_mpi_db()
            background = universe["Pathway_Name"].astype("string").fillna("").str.strip()
            background = background[background.ne("")]
        except Exception:
            background = selected_paths
        bg_counts = background.value_counts()
        sel_counts = selected_paths.value_counts()
        n = int(sel_counts.sum())
        N = int(bg_counts.sum())
        rows = []
        for pathway, k in sel_counts.head(50).items():
            K = int(bg_counts.get(pathway, k))
            fold = (k / n) / (K / N) if N and K and n else 1
            p_value = self._hypergeom_sf(int(k), N, K, n)
            rows.append({
                "Pathway_Name": pathway,
                "Count": int(k),
                "Background_Count": K,
                "Fold_Enrichment": round(float(fold), 3),
                "P_Value": p_value,
            })
        out = pd.DataFrame(rows)
        if out.empty:
            return out
        out = out.sort_values("P_Value")
        m = len(out)
        out["FDR"] = [
            min(float(p) * m / rank, 1.0)
            for rank, p in enumerate(out["P_Value"], start=1)
        ]
        return out.sort_values(["FDR", "P_Value"]).head(20)

    def _hypergeom_sf(self, k: int, N: int, K: int, n: int) -> float:
        try:
            from scipy.stats import hypergeom
            return float(hypergeom.sf(k - 1, N, K, n))
        except Exception:
            # Conservative fallback: no significance claim if scipy is unavailable.
            return 1.0

    def get_cytoscape_elements(
        self,
        disease_name: str,
        confidence_threshold: float = 0.3,
        disease_label: str = "",
    ) -> List[Dict]:
        """Build Cytoscape.js-compatible elements list.

        Metabolites  -> blue circles         (classes: "metabolite")
        Proteins     -> red rounded-squares  (classes: "protein")
        Edges        -> score >= *confidence_threshold*
        """
        data = self.get_disease_data(disease_name, disease_label)
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
        score_col = _score_col(predictions)
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

    def get_hub_tables(self, disease_name: str, top_n: int = 10, disease_label: str = "") -> Dict[str, List[Dict]]:
        """Return top metabolite and protein hubs from network_stats.

        Falls back to computing from predictions if stats file is missing.
        """
        data = self.get_disease_data(disease_name, disease_label)
        stats = data.get("network_stats", {})

        met_hubs = stats.get("metabolite_hubs", [])[:top_n]
        prot_hubs = stats.get("protein_hubs", [])[:top_n]

        # Fallback: derive from predictions
        if not met_hubs or not prot_hubs:
            predictions = data.get("predictions", pd.DataFrame())
            if not predictions.empty and _score_col(predictions) in predictions.columns:
                score_col = _score_col(predictions)
                scores = pd.to_numeric(predictions[score_col], errors="coerce")
                above = predictions[scores >= 0.3]
                if not above.empty:
                    if not met_hubs:
                        met_counts = above["Metabolite"].value_counts().head(top_n)
                        met_hubs = [{"node": n, "degree": int(d)} for n, d in met_counts.items()]
                    if not prot_hubs:
                        prot_counts = above["Protein"].value_counts().head(top_n)
                        prot_hubs = [{"node": n, "degree": int(d)} for n, d in prot_counts.items()]

        return {"metabolite_hubs": met_hubs, "protein_hubs": prot_hubs}

    def get_enrichment_data(self, disease_name: str, disease_label: str = "") -> pd.DataFrame:
        """Return pathway enrichment results for *disease_name*."""
        data = self.get_disease_data(disease_name, disease_label)
        return data.get("enrichment", pd.DataFrame())

    def get_network_stats(self, disease_name: str, disease_label: str = "") -> Dict[str, Any]:
        """Return network statistics dict."""
        data = self.get_disease_data(disease_name, disease_label)
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
        score_col = _score_col(edges)
        scores = pd.to_numeric(edges[score_col], errors="coerce")
        edges = edges.assign(_score=scores).sort_values("_score", ascending=False)
        detail["degree"] = len(edges)
        detail["top_interactions"] = [
            {"partner": str(row[partner_col]), "score": round(float(row["_score"]), 5)}
            for _, row in edges.head(10).iterrows()
        ]

        return detail
