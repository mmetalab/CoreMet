"""
CoreMet-AI — Graph Traversal & Path Finder

Executes structured graph queries over in-memory pandas DataFrames.
Supports multi-hop traversal across all 7 interaction layers.
Returns subgraphs (nodes + edges) and ranked paths.
"""

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# ── Layer definitions ─────────────────────────────────────────────────
LAYER_CONFIG = {
    "MPI": {
        "loader": "app.services.mpi_service.get_mpi_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "Protein_Name",
        "target_id": "Uniprot_ID",
        "target_type": "protein",
        "label": "Metabolite–Protein",
        "evidence_col": "Evidence_Source",
        "confidence_col": None,
        "pmid_col": None,
        "subtype_col": "interaction_subtype",
    },
    "MEI": {
        "loader": "app.services.mei_service.get_mei_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "Gene_Name",
        "target_id": "Uniprot_ID",
        "target_type": "protein",
        "label": "Metabolite–Enzyme",
        "evidence_col": "Evidence_Source",
        "confidence_col": None,
        "pmid_col": None,
        "subtype_col": None,
    },
    "MDI": {
        "loader": "app.services.mdi_service.get_mdi_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "Disease_Name",
        "target_id": "Disease_Name",
        "target_type": "disease",
        "label": "Metabolite–Disease",
        "evidence_col": "evidence_type",
        "confidence_col": "confidence",
        "pmid_col": "pmid",
        "subtype_col": "association_subtype",
    },
    "MMI": {
        "loader": "app.services.mmi_service.get_mmi_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "Microbe_Name",
        "target_id": "Taxonomy_ID",
        "target_type": "microbe",
        "label": "Metabolite–Microbe",
        "evidence_col": "Evidence_Level",
        "confidence_col": None,
        "pmid_col": "PMID",
        "subtype_col": "Relationship_Type",
    },
    "MDrI": {
        "loader": "app.services.mdri_service.get_mdri_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "Drug_Name",
        "target_id": "DrugBank_ID",
        "target_type": "drug",
        "label": "Metabolite–Drug",
        "evidence_col": "Evidence_Level",
        "confidence_col": None,
        "pmid_col": "PMID",
        "subtype_col": "Interaction_Type",
    },
    "MGI": {
        "loader": "app.services.mgi_service.get_mgi_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "Gene_Symbol",
        "target_id": "Gene_ID",
        "target_type": "gene",
        "label": "Metabolite–Gene",
        "evidence_col": "Source",
        "confidence_col": None,
        "pmid_col": "PMID",
        "subtype_col": "Interaction_Type",
    },
    "mGWAS": {
        "loader": "app.services.mgwas_service.get_mgwas_db",
        "met_col": "HMDB_ID",
        "met_name": "Metabolite_Name",
        "target_col": "rsID",
        "target_id": "rsID",
        "target_type": "snp",
        "label": "Metabolite–SNP",
        "evidence_col": "Source",
        "confidence_col": "P_Value",
        "pmid_col": "PMID",
        "subtype_col": "Trait",
    },
}

# Cache loaded DataFrames
_layer_cache: Dict[str, pd.DataFrame] = {}


def _get_layer_df(layer: str) -> pd.DataFrame:
    """Load a layer DataFrame (cached)."""
    if layer in _layer_cache:
        return _layer_cache[layer]
    cfg = LAYER_CONFIG[layer]
    module_path, func_name = cfg["loader"].rsplit(".", 1)
    import importlib
    mod = importlib.import_module(module_path)
    df = getattr(mod, func_name)()
    _layer_cache[layer] = df
    return df


def _parse_confidence(val) -> float:
    """Parse a confidence value to float, handling various formats."""
    if pd.isna(val) or val == "":
        return 0.5  # default moderate
    try:
        f = float(val)
        # P-values from mGWAS are small numbers; invert
        if f < 0.05:
            return min(1.0, -1 * (f - 1))  # small p-value → high confidence
        return min(1.0, f)
    except (ValueError, TypeError):
        # String evidence levels
        level_map = {
            "curated": 0.9, "experimental": 0.95, "predicted": 0.5,
            "inferred": 0.4, "text-mined": 0.3, "high": 0.9,
            "medium": 0.6, "low": 0.3,
        }
        return level_map.get(str(val).lower().strip(), 0.5)


# ── Core graph queries ────────────────────────────────────────────────

def get_metabolite_neighbors(
    hmdb_id: str,
    layers: Optional[List[str]] = None,
    max_per_layer: int = 50,
    metabolite_name: Optional[str] = None,
) -> List[dict]:
    """
    Get all neighbors of a metabolite across specified layers.
    Matches by HMDB_ID first, then by metabolite name as fallback.
    Returns a list of edge dicts.
    """
    if layers is None:
        layers = list(LAYER_CONFIG.keys())

    # Build a set of name variants for fallback matching
    name_variants = set()
    if metabolite_name:
        name_variants.add(metabolite_name.lower().strip())
        # Common acid/ate suffix conversions
        ml = metabolite_name.lower().strip()
        if ml.endswith("ic acid"):
            name_variants.add(ml[:-7] + "ate")  # butyric acid → butyrate
        elif ml.endswith("ate"):
            name_variants.add(ml[:-3] + "ic acid")  # butyrate → butyric acid
        if ml.endswith("ic acid"):
            base = ml[:-7]
            name_variants.add(base + "anoic acid")  # butyric acid → butanoic acid
        if ml.endswith("anoic acid"):
            base = ml[:-10]
            name_variants.add(base + "ic acid")  # butanoic acid → butyric acid
            name_variants.add(base + "ate")  # butanoic acid → butanoate

    edges = []
    for layer in layers:
        if layer not in LAYER_CONFIG:
            continue
        cfg = LAYER_CONFIG[layer]
        df = _get_layer_df(layer)
        if df.empty:
            continue

        # Primary match: by HMDB_ID
        mask = df[cfg["met_col"]].str.upper() == hmdb_id.upper()

        # Fallback match: by metabolite name variants
        if mask.sum() == 0 and name_variants:
            name_col = cfg["met_name"]
            if name_col in df.columns:
                mask = df[name_col].str.lower().str.strip().isin(name_variants)

        subset = df[mask].head(max_per_layer)
        for _, row in subset.iterrows():
            target_name = str(row.get(cfg["target_col"], ""))
            target_id = str(row.get(cfg["target_id"], target_name))
            met_name = str(row.get(cfg["met_name"], hmdb_id))
            conf_raw = row.get(cfg["confidence_col"], "") if cfg["confidence_col"] else ""
            evidence = str(row.get(cfg["evidence_col"], "")) if cfg["evidence_col"] else ""
            pmid = str(row.get(cfg["pmid_col"], "")) if cfg["pmid_col"] else ""
            subtype = str(row.get(cfg["subtype_col"], "")) if cfg["subtype_col"] else ""
            edges.append({
                "source": hmdb_id,
                "source_name": met_name,
                "source_type": "metabolite",
                "target": target_id,
                "target_name": target_name,
                "target_type": cfg["target_type"],
                "layer": layer,
                "layer_label": cfg["label"],
                "subtype": subtype,
                "confidence": _parse_confidence(conf_raw),
                "evidence": evidence,
                "pmid": pmid,
                "source_db": layer,
            })
    return edges


def get_entity_metabolites(
    entity_name: str,
    entity_type: str,
    layers: Optional[List[str]] = None,
    max_results: int = 100,
) -> List[dict]:
    """
    Find metabolites connected to a non-metabolite entity.
    Reverse lookup: given a disease/gene/drug/microbe/snp, find its metabolites.
    """
    if layers is None:
        layers = list(LAYER_CONFIG.keys())

    # Map entity type to which layers can find it
    type_to_layers = {
        "disease": ["MDI"],
        "gene": ["MGI"],
        "protein": ["MPI", "MEI"],
        "drug": ["MDrI"],
        "microbe": ["MMI"],
        "snp": ["mGWAS"],
    }
    relevant = [l for l in type_to_layers.get(entity_type, []) if l in layers]

    edges = []
    q = entity_name.lower().strip()

    for layer in relevant:
        cfg = LAYER_CONFIG[layer]
        df = _get_layer_df(layer)
        if df.empty:
            continue

        # Find rows matching the entity
        mask = df[cfg["target_col"]].str.lower().str.strip() == q
        if mask.sum() == 0:
            # Try target_id column
            mask = df[cfg["target_id"]].str.lower().str.strip() == q
        if mask.sum() == 0:
            # Fuzzy: prefix/genus match (e.g., "lactobacillus" matches "lactobacillus acidophilus")
            mask = df[cfg["target_col"]].str.lower().str.strip().str.startswith(q + " ")
        if mask.sum() == 0 and len(q) >= 4:
            # Fuzzy: substring match
            mask = df[cfg["target_col"]].str.lower().str.contains(q, na=False, regex=False)

        subset = df[mask].head(max_results)
        for _, row in subset.iterrows():
            hmdb = str(row.get(cfg["met_col"], ""))
            met_name = str(row.get(cfg["met_name"], hmdb))
            target_name = str(row.get(cfg["target_col"], ""))
            conf_raw = row.get(cfg["confidence_col"], "") if cfg["confidence_col"] else ""
            evidence = str(row.get(cfg["evidence_col"], "")) if cfg["evidence_col"] else ""
            pmid = str(row.get(cfg["pmid_col"], "")) if cfg["pmid_col"] else ""
            subtype = str(row.get(cfg["subtype_col"], "")) if cfg["subtype_col"] else ""
            edges.append({
                "source": hmdb,
                "source_name": met_name,
                "source_type": "metabolite",
                "target": entity_name,
                "target_name": target_name,
                "target_type": entity_type,
                "layer": layer,
                "layer_label": cfg["label"],
                "subtype": subtype,
                "confidence": _parse_confidence(conf_raw),
                "evidence": evidence,
                "pmid": pmid,
                "source_db": layer,
            })
    return edges


def _fuzzy_entity_match(name_a: str, name_b: str) -> bool:
    """Check if two entity names match, with genus-level and substring support."""
    a = name_a.lower().strip()
    b = name_b.lower().strip()
    if not a or not b:
        return False
    # Exact match
    if a == b:
        return True
    # Prefix/genus match (e.g., "lactobacillus" matches "lactobacillus acidophilus")
    if a.startswith(b + " ") or b.startswith(a + " "):
        return True
    # Substring match (for longer queries, min 4 chars)
    if len(a) >= 4 and a in b:
        return True
    if len(b) >= 4 and b in a:
        return True
    return False


def find_paths(
    source_name: str,
    source_type: str,
    source_id: str,
    target_name: str,
    target_type: str,
    target_id: str,
    layers: Optional[List[str]] = None,
    max_hops: int = 3,
    max_paths: int = 10,
    min_confidence: float = 0.0,
) -> List[List[dict]]:
    """
    Find paths between two entities through the metabolite-centered graph.

    Strategy:
    - If source is metabolite: 1-hop to targets, or 2-hop via intermediates
    - If target is metabolite: reverse lookup
    - If neither is metabolite: bridge through metabolites (2-hop or 3-hop)
    """
    if layers is None:
        layers = list(LAYER_CONFIG.keys())

    paths = []

    # ── Case 1: source is metabolite ──────────────────────────────────
    if source_type == "metabolite":
        hmdb = source_id
        neighbors = get_metabolite_neighbors(hmdb, layers=layers, max_per_layer=500,
                                             metabolite_name=source_name)

        # 1-hop: direct connection to target
        for edge in neighbors:
            if _fuzzy_entity_match(edge["target_name"], target_name):
                if edge["confidence"] >= min_confidence:
                    paths.append([edge])
            elif target_id and _fuzzy_entity_match(edge["target"], target_id):
                if edge["confidence"] >= min_confidence:
                    paths.append([edge])

        # 2-hop: metabolite → intermediate_entity → [shared_met] → target
        # Strategy: get target's metabolites, find intermediates from source's
        # neighbors that share metabolites with the target
        if len(paths) < max_paths and max_hops >= 2:
            # Get metabolites connected to target
            target_edges = get_entity_metabolites(
                target_name, target_type, layers=layers, max_results=200
            )
            target_hmdb_set = {e["source"] for e in target_edges if e["source"]}
            target_met_map = {e["source"]: e for e in target_edges}

            # Deduplicate intermediates and limit to top 20
            seen_int = set()
            unique_neighbors = []
            for edge1 in neighbors:
                int_key = edge1["target_name"].lower()
                if int_key not in seen_int:
                    seen_int.add(int_key)
                    unique_neighbors.append(edge1)
            unique_neighbors = unique_neighbors[:20]

            # For each intermediate, check if its metabolites overlap with target's
            seen_pairs = set()
            for edge1 in unique_neighbors:
                int_name = edge1["target_name"]
                int_type = edge1["target_type"]
                int_mets = get_entity_metabolites(
                    int_name, int_type, layers=layers, max_results=50
                )
                for edge2 in int_mets:
                    met_hmdb = edge2["source"]
                    if met_hmdb in target_hmdb_set:
                        pair_key = (int_name, met_hmdb)
                        if pair_key not in seen_pairs:
                            seen_pairs.add(pair_key)
                            tgt_edge = target_met_map[met_hmdb]
                            paths.append([edge1, edge2, tgt_edge])
                if len(paths) >= max_paths * 2:
                    break

    # ── Case 2: target is metabolite ──────────────────────────────────
    elif target_type == "metabolite":
        hmdb = target_id
        # Forward: get metabolite's neighbors, check if source is among them
        neighbors = get_metabolite_neighbors(hmdb, layers=layers, max_per_layer=500,
                                             metabolite_name=target_name)
        for edge in neighbors:
            if _fuzzy_entity_match(edge["target_name"], source_name):
                if edge["confidence"] >= min_confidence:
                    rev = {**edge,
                           "source": edge["target"], "source_name": edge["target_name"],
                           "source_type": edge["target_type"],
                           "target": edge["source"], "target_name": edge["source_name"],
                           "target_type": edge["source_type"]}
                    paths.append([rev])
            elif source_id and _fuzzy_entity_match(edge["target"], source_id):
                if edge["confidence"] >= min_confidence:
                    rev = {**edge,
                           "source": edge["target"], "source_name": edge["target_name"],
                           "source_type": edge["target_type"],
                           "target": edge["source"], "target_name": edge["source_name"],
                           "target_type": edge["source_type"]}
                    paths.append([rev])

        # Reverse: get source entity's metabolite connections, check if target met is among them
        if not paths:
            source_met_edges = get_entity_metabolites(
                source_name, source_type, layers=layers, max_results=500
            )
            tgt_name_lower = target_name.lower().strip()
            tgt_hmdb = hmdb.upper().strip()
            # Build name variants for target metabolite
            tgt_variants = {tgt_name_lower}
            if tgt_name_lower.endswith("ic acid"):
                tgt_variants.add(tgt_name_lower[:-7] + "ate")
            elif tgt_name_lower.endswith("ate"):
                tgt_variants.add(tgt_name_lower[:-3] + "ic acid")

            for edge in source_met_edges:
                met_name = edge["source_name"].lower().strip()
                met_hmdb = edge["source"].upper().strip()
                matched = (
                    met_hmdb == tgt_hmdb
                    or met_name in tgt_variants
                    or any(v in met_name for v in tgt_variants if len(v) >= 4)
                )
                if matched and edge["confidence"] >= min_confidence:
                    # Edge is: met → source_entity; reverse for display: source → met
                    paths.append([edge])

    # ── Case 3: neither is metabolite → bridge through metabolites ────
    else:
        # Source → metabolites
        source_edges = get_entity_metabolites(
            source_name, source_type, layers=layers, max_results=200
        )
        source_mets = {}
        for e in source_edges:
            hmdb = e["source"]
            met_name = e["source_name"]
            if hmdb and hmdb not in source_mets:
                source_mets[hmdb] = {"edge": e, "name": met_name}

        # Target → metabolites
        target_edges = get_entity_metabolites(
            target_name, target_type, layers=layers, max_results=200
        )
        target_mets = {}
        for e in target_edges:
            hmdb = e["source"]
            if hmdb and hmdb not in target_mets:
                target_mets[hmdb] = e

        # 2-hop: source → shared_metabolite → target
        shared = set(source_mets.keys()) & set(target_mets.keys())
        for hmdb in list(shared)[:max_paths * 3]:
            e1 = source_mets[hmdb]["edge"]
            e2 = target_mets[hmdb]
            avg_conf = (e1["confidence"] + e2["confidence"]) / 2
            if avg_conf >= min_confidence:
                paths.append([e1, e2])

        # 3-hop: source → met_A → intermediate → met_B → target
        if len(paths) < max_paths and max_hops >= 3:
            target_hmdb_set = set(target_mets.keys())
            for hmdb, met_info in list(source_mets.items())[:30]:
                met_name = met_info["name"]
                src_edge = met_info["edge"]
                # Get non-metabolite neighbors of this metabolite
                met_nbrs = get_metabolite_neighbors(
                    hmdb, layers=layers, max_per_layer=30,
                    metabolite_name=met_name,
                )
                for hop2 in met_nbrs:
                    # Direct match to target
                    if _fuzzy_entity_match(hop2["target_name"], target_name):
                        paths.append([src_edge, hop2])
                    # Check if this intermediate's metabolites overlap with target's
                    elif hop2["target_type"] != "metabolite":
                        int_mets = get_entity_metabolites(
                            hop2["target_name"], hop2["target_type"],
                            layers=layers, max_results=20,
                        )
                        for hop3_edge in int_mets:
                            if hop3_edge["source"] in target_hmdb_set:
                                tgt_edge = target_mets[hop3_edge["source"]]
                                paths.append([src_edge, hop2, tgt_edge])
                                break
                    if len(paths) >= max_paths * 3:
                        break
                if len(paths) >= max_paths * 3:
                    break

    # Sort by average confidence
    def path_score(p):
        if not p:
            return 0
        return sum(e["confidence"] for e in p) / len(p)

    paths.sort(key=path_score, reverse=True)
    return paths[:max_paths]


def build_subgraph(
    paths: List[List[dict]],
    extra_edges: Optional[List[dict]] = None,
) -> dict:
    """
    Convert a list of paths into a subgraph structure.
    Returns: {"nodes": [...], "edges": [...]}
    """
    nodes: Dict[str, dict] = {}
    edges: List[dict] = []
    edge_keys: Set[str] = set()

    def _add_node(name: str, ntype: str, nid: str):
        key = f"{ntype}:{name}"
        if key not in nodes:
            nodes[key] = {"id": key, "name": name, "type": ntype, "internal_id": nid}

    def _add_edge(edge: dict):
        src_key = f"{edge['source_type']}:{edge['source_name']}"
        tgt_key = f"{edge['target_type']}:{edge['target_name']}"
        ek = f"{src_key}→{tgt_key}:{edge['layer']}"
        if ek not in edge_keys:
            edge_keys.add(ek)
            edges.append({
                "source": src_key,
                "target": tgt_key,
                "layer": edge["layer"],
                "layer_label": edge["layer_label"],
                "subtype": edge.get("subtype", ""),
                "confidence": edge.get("confidence", 0.5),
                "evidence": edge.get("evidence", ""),
                "pmid": edge.get("pmid", ""),
                "source_db": edge.get("source_db", ""),
            })

    for path in paths:
        for edge in path:
            _add_node(edge["source_name"], edge["source_type"], edge["source"])
            _add_node(edge["target_name"], edge["target_type"], edge["target"])
            _add_edge(edge)

    if extra_edges:
        for edge in extra_edges:
            _add_node(edge["source_name"], edge["source_type"], edge["source"])
            _add_node(edge["target_name"], edge["target_type"], edge["target"])
            _add_edge(edge)

    return {"nodes": list(nodes.values()), "edges": edges}
