"""
CoreMet-AI — Query Orchestrator

Single entry point that coordinates the full AI pipeline:
  1. Parse query → structured plan
  2. Resolve entities → internal IDs
  3. Execute graph traversal → paths
  4. Rank paths → top N
  5. Aggregate evidence
  6. Summarize (grounded)
  7. Build subgraph for visualization

Returns a complete result object for UI rendering.
"""

import logging
import time
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def execute_ai_query(query: str, provider: str = "auto",
                     api_key: str = "", model: str = "") -> dict:
    """
    Execute the full CoreMet-AI pipeline.

    Args:
        query: Natural language question from the user.
        provider: LLM provider — "auto", "openai", "google", "template".
        api_key: Optional API key (overrides env).
        model: Optional model name (overrides env/default).

    Returns:
        {
            "status": "success" | "no_results" | "error",
            "query_plan": {...},
            "resolved_source": {...},
            "resolved_target": {...},
            "ranked_paths": [...],
            "evidence": {...},
            "summary": {...},
            "subgraph": {"nodes": [...], "edges": [...]},
            "timing_ms": int,
        }
    """
    t0 = time.time()

    try:
        # ── Stage 1: Parse query ──────────────────────────────────────
        from app.services.query_parser import parse_query
        plan = parse_query(query, provider=provider,
                           api_key=api_key, model=model)
        logger.info(f"AI query parsed: source='{plan.get('source_entity')}', "
                     f"target='{plan.get('target_entity')}', intent='{plan.get('intent')}'")

        # ── Stage 2: Resolve entities ─────────────────────────────────
        from app.services.entity_resolver import resolve_entity
        resolved_src = None
        resolved_tgt = None

        if plan.get("source_entity"):
            candidates = resolve_entity(
                plan["source_entity"],
                expected_type=plan.get("source_type") or None,
            )
            if candidates:
                resolved_src = candidates[0]
                plan["source_type"] = resolved_src["type"]
                plan["source_entity"] = resolved_src["name"]
                plan["source_id"] = resolved_src["id"]
                logger.info(f"Resolved source: {resolved_src['name']} "
                            f"({resolved_src['type']}, conf={resolved_src['confidence']})")

        if plan.get("target_entity"):
            candidates = resolve_entity(
                plan["target_entity"],
                expected_type=plan.get("target_type") or None,
            )
            if candidates:
                resolved_tgt = candidates[0]
                plan["target_type"] = resolved_tgt["type"]
                plan["target_entity"] = resolved_tgt["name"]
                plan["target_id"] = resolved_tgt["id"]
                logger.info(f"Resolved target: {resolved_tgt['name']} "
                            f"({resolved_tgt['type']}, conf={resolved_tgt['confidence']})")

        # ── Stage 3: Graph traversal ──────────────────────────────────
        from app.services.graph_traversal import (
            find_paths, get_metabolite_neighbors, get_entity_metabolites,
            build_subgraph,
        )

        paths = []
        extra_edges = []
        intent = plan.get("intent", "entity")

        if intent == "relation" and resolved_src and resolved_tgt:
            # Two-entity query: find connecting paths
            paths = find_paths(
                source_name=plan["source_entity"],
                source_type=plan["source_type"],
                source_id=plan.get("source_id", plan["source_entity"]),
                target_name=plan["target_entity"],
                target_type=plan["target_type"],
                target_id=plan.get("target_id", plan["target_entity"]),
                layers=plan.get("layers"),
                max_hops=plan.get("max_hops", 2),
                max_paths=plan.get("top_paths", 10),
                min_confidence=plan.get("confidence_threshold", 0.0),
            )
        elif resolved_src:
            # Single-entity query: explore neighbors
            if plan["source_type"] == "metabolite":
                edges = get_metabolite_neighbors(
                    plan.get("source_id", plan["source_entity"]),
                    layers=plan.get("layers"),
                    max_per_layer=30,
                    metabolite_name=plan.get("source_entity"),
                )
            else:
                edges = get_entity_metabolites(
                    plan["source_entity"],
                    plan["source_type"],
                    layers=plan.get("layers"),
                    max_results=50,
                )
            # Wrap each edge as a single-edge path
            paths = [[e] for e in edges[:plan.get("top_paths", 10)]]
            extra_edges = edges[plan.get("top_paths", 10):]

        if not paths:
            elapsed = int((time.time() - t0) * 1000)
            return {
                "status": "no_results",
                "query_plan": plan,
                "resolved_source": resolved_src,
                "resolved_target": resolved_tgt,
                "ranked_paths": [],
                "evidence": {"total_edges": 0, "total_pmids": 0,
                             "evidence_distribution": {}, "source_distribution": {},
                             "layer_distribution": {}, "confidence_distribution": {},
                             "evidence_table": []},
                "summary": {
                    "summary": (
                        f"No paths found connecting "
                        f"{plan.get('source_entity', 'the query entity')} "
                        f"{'to ' + plan.get('target_entity') if plan.get('target_entity') else ''} "
                        f"in CoreMet. Try broadening the search layers or "
                        f"reducing the confidence threshold."
                    ),
                    "method": "template",
                    "confidence_label": "no data",
                },
                "subgraph": {"nodes": [], "edges": []},
                "timing_ms": elapsed,
            }

        # ── Stage 4: Rank paths ───────────────────────────────────────
        from app.services.path_ranking import rank_paths, aggregate_evidence

        ranked = rank_paths(
            paths,
            top_n=plan.get("top_paths", 10),
            min_score=plan.get("confidence_threshold", 0.0),
        )

        # ── Stage 5: Aggregate evidence ───────────────────────────────
        evidence = aggregate_evidence(ranked)

        # ── Stage 6: Summarize ────────────────────────────────────────
        from app.services.ai_summarizer import summarize
        summary = summarize(plan, ranked, evidence,
                            provider=provider, api_key=api_key, model=model)

        # ── Stage 7: Build subgraph ───────────────────────────────────
        subgraph = build_subgraph(
            [item["path"] for item in ranked],
            extra_edges=extra_edges[:50],
        )

        elapsed = int((time.time() - t0) * 1000)

        return {
            "status": "success",
            "query_plan": plan,
            "resolved_source": resolved_src,
            "resolved_target": resolved_tgt,
            "ranked_paths": [
                {
                    "rank": item["rank"],
                    "score": item["score"],
                    "path": item["path"],
                    "path_description": _format_path_display(item["path"]),
                }
                for item in ranked
            ],
            "evidence": evidence,
            "summary": summary,
            "subgraph": subgraph,
            "timing_ms": elapsed,
        }

    except Exception as e:
        logger.error(f"AI query failed: {e}", exc_info=True)
        elapsed = int((time.time() - t0) * 1000)
        return {
            "status": "error",
            "error": str(e),
            "query_plan": {},
            "resolved_source": None,
            "resolved_target": None,
            "ranked_paths": [],
            "evidence": {"total_edges": 0, "total_pmids": 0,
                         "evidence_distribution": {}, "source_distribution": {},
                         "layer_distribution": {}, "confidence_distribution": {},
                         "evidence_table": []},
            "summary": {"summary": f"An error occurred: {e}", "method": "error",
                        "confidence_label": "error"},
            "subgraph": {"nodes": [], "edges": []},
            "timing_ms": elapsed,
        }


def _format_path_display(path: List[dict]) -> str:
    """Format a path as a concise display string."""
    if not path:
        return ""
    parts = []
    for edge in path:
        src = edge.get("source_name", "?")
        tgt = edge.get("target_name", "?")
        layer = edge.get("layer", "")
        parts.append(f"{src} →[{layer}]→ {tgt}")
    return " | ".join(parts)
