"""
CoreMet-AI — Path Ranking Engine

Scores and ranks candidate paths before presentation.
Each path is scored using:
  - average edge confidence
  - evidence type hierarchy (experimental > curated > inferred)
  - number of supporting databases
  - path length penalty
  - PMID support bonus
"""

import logging
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)

# ── Evidence hierarchy weights ────────────────────────────────────────
EVIDENCE_WEIGHTS = {
    "experimental": 1.0,
    "curated": 0.9,
    "literature": 0.8,
    "text-mined": 0.5,
    "predicted": 0.4,
    "inferred": 0.3,
    "computational": 0.3,
    "": 0.5,  # unknown
}

# ── Confidence labels ─────────────────────────────────────────────────
CONFIDENCE_LABELS = {
    (0.8, 1.01): "supported",
    (0.5, 0.8): "moderate evidence",
    (0.0, 0.5): "hypothesis-generating",
}


def _evidence_score(evidence_str: str) -> float:
    """Map an evidence type string to a numeric score."""
    if not evidence_str:
        return 0.5
    e = evidence_str.lower().strip()
    for key, val in EVIDENCE_WEIGHTS.items():
        if key in e:
            return val
    return 0.5


def score_path(path: List[dict]) -> dict:
    """
    Score a single path (list of edge dicts).

    Returns a dict with:
      - composite_score: float 0-1
      - avg_confidence: float
      - evidence_score: float
      - db_diversity: int (number of distinct source databases)
      - length_penalty: float
      - pmid_count: int
      - confidence_label: str
    """
    if not path:
        return {
            "composite_score": 0.0, "avg_confidence": 0.0,
            "evidence_score": 0.0, "db_diversity": 0,
            "length_penalty": 1.0, "pmid_count": 0,
            "confidence_label": "hypothesis-generating",
        }

    n = len(path)

    # Average edge confidence
    confidences = [e.get("confidence", 0.5) for e in path]
    avg_conf = sum(confidences) / n

    # Evidence type score (average across edges)
    ev_scores = [_evidence_score(e.get("evidence", "")) for e in path]
    avg_ev = sum(ev_scores) / n

    # Database diversity
    dbs = set(e.get("source_db", "") for e in path if e.get("source_db"))
    db_div = len(dbs)
    db_bonus = min(0.15, db_div * 0.05)

    # Path length penalty (shorter is better)
    length_penalty = 1.0 / (1.0 + 0.2 * (n - 1))

    # PMID support
    pmids = set()
    for e in path:
        pmid = e.get("pmid", "")
        if pmid and pmid != "nan":
            for p in str(pmid).split("|"):
                p = p.strip()
                if p and p != "nan":
                    pmids.add(p)
    pmid_bonus = min(0.1, len(pmids) * 0.02)

    # Composite score
    composite = (
        0.40 * avg_conf +
        0.25 * avg_ev +
        0.15 * length_penalty +
        0.10 * db_bonus / 0.15 +  # normalize to 0-1
        0.10 * min(1.0, pmid_bonus / 0.1)
    )
    composite = round(min(1.0, composite), 4)

    # Label
    label = "hypothesis-generating"
    for (lo, hi), lbl in CONFIDENCE_LABELS.items():
        if lo <= composite < hi:
            label = lbl
            break

    return {
        "composite_score": composite,
        "avg_confidence": round(avg_conf, 4),
        "evidence_score": round(avg_ev, 4),
        "db_diversity": db_div,
        "length_penalty": round(length_penalty, 4),
        "pmid_count": len(pmids),
        "confidence_label": label,
    }


def rank_paths(
    paths: List[List[dict]],
    top_n: int = 10,
    min_score: float = 0.0,
) -> List[dict]:
    """
    Score and rank a list of paths.

    Returns a list of dicts, each containing:
      - path: List[dict] (the edge list)
      - score: dict (from score_path)
      - rank: int
    """
    scored = []
    for path in paths:
        s = score_path(path)
        if s["composite_score"] >= min_score:
            scored.append({"path": path, "score": s})

    scored.sort(key=lambda x: x["score"]["composite_score"], reverse=True)

    for i, item in enumerate(scored[:top_n]):
        item["rank"] = i + 1

    return scored[:top_n]


def aggregate_evidence(ranked_paths: List[dict]) -> dict:
    """
    Compile evidence summary from ranked paths.

    Returns:
      - total_edges: int
      - total_pmids: int
      - evidence_distribution: dict
      - source_distribution: dict
      - layer_distribution: dict
      - confidence_distribution: dict
      - evidence_table: list of edge-level records
    """
    evidence_dist = {}
    source_dist = {}
    layer_dist = {}
    conf_buckets = {"high (≥0.8)": 0, "moderate (0.5–0.8)": 0, "low (<0.5)": 0}
    all_pmids = set()
    evidence_table = []

    for item in ranked_paths:
        for edge in item["path"]:
            # Evidence distribution
            ev = edge.get("evidence", "unknown") or "unknown"
            evidence_dist[ev] = evidence_dist.get(ev, 0) + 1

            # Source distribution
            src = edge.get("source_db", "unknown") or "unknown"
            source_dist[src] = source_dist.get(src, 0) + 1

            # Layer distribution
            layer = edge.get("layer", "unknown")
            layer_dist[layer] = layer_dist.get(layer, 0) + 1

            # Confidence buckets
            conf = edge.get("confidence", 0.5)
            if conf >= 0.8:
                conf_buckets["high (≥0.8)"] += 1
            elif conf >= 0.5:
                conf_buckets["moderate (0.5–0.8)"] += 1
            else:
                conf_buckets["low (<0.5)"] += 1

            # PMIDs
            pmid = edge.get("pmid", "")
            if pmid and pmid != "nan":
                for p in str(pmid).split("|"):
                    p = p.strip()
                    if p and p != "nan":
                        all_pmids.add(p)

            # Evidence table row
            evidence_table.append({
                "rank": item["rank"],
                "source_entity": edge.get("source_name", ""),
                "source_type": edge.get("source_type", ""),
                "target_entity": edge.get("target_name", ""),
                "target_type": edge.get("target_type", ""),
                "layer": layer,
                "subtype": edge.get("subtype", ""),
                "confidence": round(conf, 3),
                "evidence": ev,
                "pmid": edge.get("pmid", ""),
                "source_db": src,
            })

    return {
        "total_edges": len(evidence_table),
        "total_pmids": len(all_pmids),
        "evidence_distribution": evidence_dist,
        "source_distribution": source_dist,
        "layer_distribution": layer_dist,
        "confidence_distribution": conf_buckets,
        "evidence_table": evidence_table,
    }
