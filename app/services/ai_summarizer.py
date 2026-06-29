"""
CoreMet-AI — Grounded Summarizer

Generates short, structured explanations strictly grounded in retrieved paths.
Two modes:
  1. Template-based (always available, deterministic)
  2. LLM-enhanced (when OPENAI_API_KEY is set)

Every statement must be traceable to database records.
Outputs are labeled as: supported, moderate evidence, or hypothesis-generating.
"""

import logging
import os
from collections import Counter
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


def _format_path_description(path: List[dict]) -> str:
    """Format a single path as a readable sentence."""
    if not path:
        return ""
    parts = []
    for i, edge in enumerate(path):
        src = edge.get("source_name", "unknown")
        tgt = edge.get("target_name", "unknown")
        layer = edge.get("layer_label", edge.get("layer", ""))
        subtype = edge.get("subtype", "")
        if i == 0:
            if subtype:
                parts.append(f"{src} → {tgt} ({layer}, {subtype})")
            else:
                parts.append(f"{src} → {tgt} ({layer})")
        else:
            if subtype:
                parts.append(f"→ {tgt} ({layer}, {subtype})")
            else:
                parts.append(f"→ {tgt} ({layer})")
    return " ".join(parts)


def _dominant_layers(ranked_paths: List[dict]) -> List[str]:
    """Find the most common layers across top paths."""
    counter = Counter()
    for item in ranked_paths:
        for edge in item["path"]:
            counter[edge.get("layer_label", edge.get("layer", ""))] += 1
    return [layer for layer, _ in counter.most_common(3)]


def _dominant_intermediates(ranked_paths: List[dict]) -> List[str]:
    """Find the most common intermediate entities."""
    counter = Counter()
    for item in ranked_paths:
        for edge in item["path"]:
            # Count non-source, non-target entities
            for key in ["source_name", "target_name"]:
                name = edge.get(key, "")
                etype = edge.get(key.replace("_name", "_type"), "")
                if etype == "metabolite":
                    counter[f"{name} (metabolite)"] += 1
    return [name for name, _ in counter.most_common(5)]


def summarize_template(
    query_plan: dict,
    ranked_paths: List[dict],
    evidence: dict,
) -> str:
    """
    Template-based summarization. Always available.
    Produces a 3-6 sentence structured summary.
    """
    source = query_plan.get("source_entity", "the source entity")
    target = query_plan.get("target_entity", "")
    intent = query_plan.get("intent", "relation")

    n_paths = len(ranked_paths)
    n_edges = evidence.get("total_edges", 0)
    n_pmids = evidence.get("total_pmids", 0)

    if n_paths == 0:
        if target:
            return (
                f"No paths were found connecting {source} to {target} in CoreMet "
                f"within the specified traversal constraints. This may indicate that "
                f"the relationship is not documented in the current database, or that "
                f"a higher hop count or broader layer selection may be needed."
            )
        return (
            f"No significant interactions were found for {source} in CoreMet "
            f"with the current filters. Try broadening the search layers or "
            f"reducing the confidence threshold."
        )

    # Get scoring info
    top_score = ranked_paths[0]["score"] if ranked_paths else {}
    label = top_score.get("confidence_label", "hypothesis-generating")
    avg_conf = top_score.get("avg_confidence", 0)

    # Dominant layers and intermediates
    layers = _dominant_layers(ranked_paths)
    intermediates = _dominant_intermediates(ranked_paths)

    sentences = []

    # Sentence 1: Overall relationship
    if target and intent == "relation":
        sentences.append(
            f"CoreMet identifies {n_paths} path(s) connecting "
            f"**{source}** to **{target}**, supported by {n_edges} edges "
            f"across {len(evidence.get('layer_distribution', {}))} interaction layers."
        )
    else:
        sentences.append(
            f"CoreMet identifies {n_edges} interactions involving **{source}** "
            f"across {len(evidence.get('layer_distribution', {}))} biological layers."
        )

    # Sentence 2-3: Main mechanisms
    if layers:
        layer_str = ", ".join(layers[:2])
        sentences.append(
            f"The dominant interaction channels are {layer_str}."
        )

    if intermediates and target:
        int_str = ", ".join(intermediates[:3])
        sentences.append(
            f"Key metabolite intermediates bridging the connection include {int_str}."
        )

    # Sentence 4: Top path detail
    if ranked_paths:
        top_path = ranked_paths[0]["path"]
        path_desc = _format_path_description(top_path)
        sentences.append(f"The highest-ranked path is: {path_desc}.")

    # Sentence 5: Evidence and confidence
    if n_pmids > 0:
        sentences.append(
            f"These results are backed by {n_pmids} PubMed reference(s). "
            f"Overall confidence: **{label}** (average score: {avg_conf:.2f})."
        )
    else:
        sentences.append(
            f"Overall confidence: **{label}** (average score: {avg_conf:.2f})."
        )

    return " ".join(sentences)


# ── LLM-enhanced summarizer ──────────────────────────────────────────

_SUMMARIZER_SYSTEM_PROMPT = """You are a scientific summarizer for CoreMet, a metabolite-centered biological knowledge graph.

Given a user's query, the retrieved paths, and supporting evidence, write a concise explanation (3-6 sentences).

Structure:
- Sentence 1: Overall relationship between the entities
- Sentences 2-3: Main biological mechanisms found in the paths
- Last sentence: Confidence statement

CRITICAL RULES:
- Use ONLY information from the provided paths and evidence. Do NOT add external knowledge.
- Every claim must be traceable to a specific path or edge.
- Label the confidence as: "supported" (score ≥0.8), "moderate evidence" (0.5-0.8), or "hypothesis-generating" (<0.5).
- Be concise. No unnecessary words.
- Use scientific language appropriate for a biology researcher.
- Format entity names in **bold**."""


def _build_user_message(query_plan, ranked_paths, evidence):
    """Build the context message for LLM summarization."""
    path_descriptions = []
    for item in ranked_paths[:10]:
        desc = _format_path_description(item["path"])
        score = item["score"]
        path_descriptions.append(
            f"  Rank {item['rank']}: {desc} "
            f"[confidence: {score['avg_confidence']:.2f}, "
            f"evidence: {score['confidence_label']}, "
            f"PMIDs: {score['pmid_count']}]"
        )

    return f"""Query: {query_plan.get('raw_query', '')}
Source entity: {query_plan.get('source_entity', '')} ({query_plan.get('source_type', '')})
Target entity: {query_plan.get('target_entity', '')} ({query_plan.get('target_type', '')})

Retrieved paths ({len(ranked_paths)} total):
{chr(10).join(path_descriptions)}

Evidence summary:
- Total edges: {evidence.get('total_edges', 0)}
- Total PubMed references: {evidence.get('total_pmids', 0)}
- Layer distribution: {evidence.get('layer_distribution', {})}
- Evidence types: {evidence.get('evidence_distribution', {})}
- Confidence distribution: {evidence.get('confidence_distribution', {})}

Write a 3-6 sentence scientific summary grounded ONLY in the above data."""


def _summarize_openai(query_plan, ranked_paths, evidence,
                      api_key=None, model=None) -> Optional[str]:
    """Summarize using OpenAI (GPT-4o, GPT-4o-mini, etc.)."""
    api_key = api_key or os.environ.get("OPENAI_API_KEY", "")
    if not api_key:
        return None
    try:
        import openai
        client = openai.OpenAI(api_key=api_key)
        model = model or os.environ.get("OPENAI_MODEL", "gpt-4o-mini")
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": _SUMMARIZER_SYSTEM_PROMPT},
                {"role": "user", "content": _build_user_message(query_plan, ranked_paths, evidence)},
            ],
            temperature=0.1,
            max_tokens=400,
        )
        return response.choices[0].message.content.strip()
    except ImportError:
        logger.warning("openai package not installed")
        return None
    except Exception as e:
        logger.warning(f"OpenAI summarizer failed: {e}")
        return None


def _summarize_google(query_plan, ranked_paths, evidence,
                      api_key=None, model=None) -> Optional[str]:
    """Summarize using Google Gemini (free tier: gemini-2.0-flash)."""
    api_key = api_key or os.environ.get("GOOGLE_API_KEY", "")
    if not api_key:
        return None
    try:
        import google.generativeai as genai
        genai.configure(api_key=api_key)
        model_name = model or os.environ.get("GOOGLE_MODEL", "gemini-2.0-flash")
        gmodel = genai.GenerativeModel(
            model_name=model_name,
            system_instruction=_SUMMARIZER_SYSTEM_PROMPT,
        )
        response = gmodel.generate_content(
            _build_user_message(query_plan, ranked_paths, evidence),
            generation_config=genai.types.GenerationConfig(
                temperature=0.1,
                max_output_tokens=400,
            ),
        )
        return response.text.strip()
    except ImportError:
        logger.warning("google-generativeai package not installed")
        return None
    except Exception as e:
        logger.warning(f"Google summarizer failed: {e}")
        return None


def summarize_llm(
    query_plan: dict,
    ranked_paths: List[dict],
    evidence: dict,
    provider: str = "",
    api_key: str = "",
    model: str = "",
) -> Optional[str]:
    """
    LLM-enhanced summarization. Supports multiple providers.
    The LLM receives ONLY retrieved data — no external knowledge.
    """
    provider = provider or os.environ.get("AI_PROVIDER", "openai")

    if provider in ("google", "free", "gemini"):
        result = _summarize_google(query_plan, ranked_paths, evidence, api_key, model)
        if result:
            return result
        # Fall back to OpenAI if Google fails
        return _summarize_openai(query_plan, ranked_paths, evidence, api_key, model)
    elif provider == "openai":
        result = _summarize_openai(query_plan, ranked_paths, evidence, api_key, model)
        if result:
            return result
        return _summarize_google(query_plan, ranked_paths, evidence, api_key, model)
    else:
        # Try both
        result = _summarize_openai(query_plan, ranked_paths, evidence, api_key, model)
        if result:
            return result
        return _summarize_google(query_plan, ranked_paths, evidence, api_key, model)


def summarize(
    query_plan: dict,
    ranked_paths: List[dict],
    evidence: dict,
    provider: str = "",
    api_key: str = "",
    model: str = "",
) -> dict:
    """
    Generate a grounded summary. Uses LLM if available, else template.

    Returns:
        {"summary": str, "method": str, "confidence_label": str, "model": str}
    """
    provider = provider or os.environ.get("AI_PROVIDER", "openai")

    # Template mode — skip LLM entirely
    if provider == "template":
        template_summary = summarize_template(query_plan, ranked_paths, evidence)
        label = "hypothesis-generating"
        if ranked_paths:
            label = ranked_paths[0]["score"].get("confidence_label", label)
        return {"summary": template_summary, "method": "template", "confidence_label": label,
                "model": "template", "provider": "none"}

    # Try LLM first
    llm_summary = summarize_llm(query_plan, ranked_paths, evidence,
                                provider=provider, api_key=api_key, model=model)
    if llm_summary:
        label = "hypothesis-generating"
        if ranked_paths:
            label = ranked_paths[0]["score"].get("confidence_label", label)
        used_model = model or os.environ.get(
            "GOOGLE_MODEL" if provider in ("google", "free", "gemini") else "OPENAI_MODEL",
            "gpt-4o-mini"
        )
        return {"summary": llm_summary, "method": "llm", "confidence_label": label,
                "model": used_model, "provider": provider}

    # Fall back to template
    template_summary = summarize_template(query_plan, ranked_paths, evidence)
    label = "hypothesis-generating"
    if ranked_paths:
        label = ranked_paths[0]["score"].get("confidence_label", label)
    return {"summary": template_summary, "method": "template", "confidence_label": label,
            "model": "template", "provider": "none"}
