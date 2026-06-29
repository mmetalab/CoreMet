"""
CoreMet-AI — Query Parser

Converts natural language questions into structured query plans.
Two modes:
  1. Rule-based parser (always available, deterministic)
  2. LLM-enhanced parser (when OPENAI_API_KEY is set)

The output is always a JSON query plan — the contract between AI and the database.
"""

import json
import logging
import os
import re
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── Default query plan schema ─────────────────────────────────────────
DEFAULT_PLAN = {
    "source_entity": "",
    "source_type": "",
    "target_entity": "",
    "target_type": "",
    "max_hops": 2,
    "layers": ["MPI", "MEI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"],
    "evidence_filters": ["curated", "experimental"],
    "top_paths": 10,
    "confidence_threshold": 0.0,
}

# ── Intent patterns ───────────────────────────────────────────────────
_RELATION_PATTERNS = [
    # "How does X affect Y" / "How might X influence Y"
    (r"how\s+(?:does|might|could|do|can|would)\s+(.+?)\s+(?:affect|influence|impact|relate\s+to|interact\s+with|connect\s+to|modulate|regulate)\s+(.+?)[\?\.]?$", "relation"),
    # "What is the relationship between X and Y"
    (r"(?:what\s+is\s+the\s+)?relationship\s+between\s+(.+?)\s+and\s+(.+?)[\?\.]?$", "relation"),
    # "X and Y connection/link/association"
    (r"(.+?)\s+and\s+(.+?)\s+(?:connection|link|association|interaction|pathway)[\?\.]?$", "relation"),
    # "Connect X to Y"
    (r"connect\s+(.+?)\s+(?:to|with)\s+(.+?)[\?\.]?$", "relation"),
    # "Path from X to Y"
    (r"path\s+(?:from|between)\s+(.+?)\s+(?:to|and)\s+(.+?)[\?\.]?$", "relation"),
    # "Link between X and Y"
    (r"link\s+between\s+(.+?)\s+and\s+(.+?)[\?\.]?$", "relation"),
    # "Role of X in Y"
    (r"role\s+of\s+(.+?)\s+in\s+(.+?)[\?\.]?$", "relation"),
    # "X in Y"
    (r"(.+?)\s+in\s+(.+?)[\?\.]?$", "relation"),
]

_EXPLORE_PATTERNS = [
    # "What genes are connected to X" / "What diseases involve X"
    (r"what\s+(\w+)\s+(?:are|is)\s+(?:connected|linked|related|associated)\s+(?:to|with)\s+(.+?)[\?\.]?$", "explore"),
    # "What drugs interact with X" / "What genes interact with X"
    (r"what\s+(\w+)\s+(?:interact|interacts)\s+with\s+(.+?)[\?\.]?$", "explore"),
    # "Show me X neighbors/connections"
    (r"(?:show|list|find|get)\s+(?:me\s+)?(?:all\s+)?(?:the\s+)?(.+?)(?:'s)?\s+(?:neighbors|connections|interactions|links|partners)[\?\.]?$", "entity"),
    # "What interacts with X"
    (r"what\s+(?:interacts|connects)\s+with\s+(.+?)[\?\.]?$", "entity"),
    # "Tell me about X"
    (r"(?:tell\s+me\s+about|describe|summarize|explain)\s+(.+?)[\?\.]?$", "entity"),
]

# ── Layer inference keywords ──────────────────────────────────────────
_LAYER_KEYWORDS = {
    "gene": ["MGI"],
    "genetic": ["MGI", "mGWAS"],
    "protein": ["MPI", "MEI"],
    "enzyme": ["MEI"],
    "disease": ["MDI"],
    "cancer": ["MDI"],
    "drug": ["MDrI"],
    "pharmacol": ["MDrI"],
    "microb": ["MMI"],
    "gut": ["MMI"],
    "snp": ["mGWAS"],
    "gwas": ["mGWAS"],
    "variant": ["mGWAS"],
}


def _infer_layers(query: str) -> List[str]:
    """Infer which interaction layers are relevant from the query text."""
    q = query.lower()
    layers = set()
    for keyword, layer_list in _LAYER_KEYWORDS.items():
        if keyword in q:
            layers.update(layer_list)
    # If no specific layers inferred, use all
    return list(layers) if layers else list(DEFAULT_PLAN["layers"])


def _infer_hops(query: str) -> int:
    """Infer traversal depth from query."""
    q = query.lower()
    if "direct" in q or "1-hop" in q or "one hop" in q:
        return 1
    if "indirect" in q or "3-hop" in q or "three hop" in q:
        return 3
    return 2


def _clean_entity(text: str) -> str:
    """Clean extracted entity text."""
    text = text.strip().strip("\"'?.")
    # Remove common prefixes
    for prefix in ["the ", "a ", "an "]:
        if text.lower().startswith(prefix):
            text = text[len(prefix):]
    return text.strip()


def parse_query_rules(query: str) -> dict:
    """
    Rule-based query parser. Always available.
    Returns a structured query plan.
    """
    plan = dict(DEFAULT_PLAN)
    plan["raw_query"] = query
    plan["parser"] = "rule-based"

    q = query.strip()
    if not q:
        return plan

    # Try relation patterns (source → target)
    for pattern, intent in _RELATION_PATTERNS:
        m = re.search(pattern, q, re.IGNORECASE)
        if m:
            plan["source_entity"] = _clean_entity(m.group(1))
            plan["target_entity"] = _clean_entity(m.group(2))
            plan["intent"] = intent
            break

    # Try explore/entity patterns
    if not plan["source_entity"]:
        for pattern, intent in _EXPLORE_PATTERNS:
            m = re.search(pattern, q, re.IGNORECASE)
            if m:
                if intent == "explore" and m.lastindex >= 2:
                    plan["target_entity"] = _clean_entity(m.group(1))
                    plan["source_entity"] = _clean_entity(m.group(2))
                else:
                    plan["source_entity"] = _clean_entity(m.group(m.lastindex))
                plan["intent"] = intent
                break

    # Fallback: treat entire query as entity search
    if not plan["source_entity"]:
        plan["source_entity"] = _clean_entity(q)
        plan["intent"] = "entity"

    # Infer layers and hops
    plan["layers"] = _infer_layers(q)
    plan["max_hops"] = _infer_hops(q)

    return plan


# ── LLM-enhanced parser ──────────────────────────────────────────────

_LLM_SYSTEM_PROMPT = """You are a query parser for CoreMet, a metabolite-centered biological knowledge graph database.

Given a user's natural language question about metabolites, diseases, genes, proteins, drugs, microbes, or SNPs, extract a structured query plan as JSON.

The database has 8 interaction layers:
- MPI: Metabolite–Protein interactions (direct binding, transport, catalysis)
- MEI: Metabolite–Enzyme interactions (enzymatic reactions)
- MDI: Metabolite–Disease associations (biomarkers, causal links)
- MMI: Metabolite–Microbe interactions (production, consumption, biotransformation)
- MDrI: Metabolite–Drug interactions (drug targets, drug metabolism)
- MGI: Metabolite–Gene interactions (gene regulation, expression effects)
- mGWAS: Metabolite–SNP associations (genome-wide associations)

Entity types: metabolite, disease, gene, protein, drug, microbe, snp

Output ONLY valid JSON with these fields:
{
  "source_entity": "name of the source entity",
  "source_type": "one of: metabolite, disease, gene, protein, drug, microbe, snp, or empty",
  "target_entity": "name of the target entity, or empty if single-entity query",
  "target_type": "type of target, or empty",
  "max_hops": 2,
  "layers": ["list of relevant layers"],
  "evidence_filters": ["curated", "experimental"],
  "top_paths": 10,
  "confidence_threshold": 0.0,
  "intent": "relation or entity or explore"
}

Rules:
- For relation queries connecting two entities, ALWAYS include ALL layers to enable multi-hop bridging through metabolites: ["MPI", "MEI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"]. Metabolites are the central hub; paths often bridge through them.
- Only narrow layers if the user explicitly asks about a specific interaction type.
- If the query is about a single entity, set target_entity to empty string.
- Default max_hops is 2. Use 3 for complex indirect questions.
- Return ONLY the JSON object, no explanation."""


def _parse_openai(query: str, api_key: str = "", model: str = "") -> Optional[dict]:
    """Parse query using OpenAI."""
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
                {"role": "system", "content": _LLM_SYSTEM_PROMPT},
                {"role": "user", "content": query},
            ],
            temperature=0.0,
            max_tokens=500,
            response_format={"type": "json_object"},
        )
        content = response.choices[0].message.content
        plan = json.loads(content)
        plan["raw_query"] = query
        plan["parser"] = "llm"
        plan["model"] = model
        for key, default in DEFAULT_PLAN.items():
            if key not in plan:
                plan[key] = default
        return plan
    except ImportError:
        return None
    except Exception as e:
        logger.warning(f"OpenAI parser failed: {e}")
        return None


def _parse_google(query: str, api_key: str = "", model: str = "") -> Optional[dict]:
    """Parse query using Google Gemini."""
    api_key = api_key or os.environ.get("GOOGLE_API_KEY", "")
    if not api_key:
        return None
    try:
        import google.generativeai as genai
        genai.configure(api_key=api_key)
        model_name = model or os.environ.get("GOOGLE_MODEL", "gemini-2.0-flash")
        gmodel = genai.GenerativeModel(
            model_name=model_name,
            system_instruction=_LLM_SYSTEM_PROMPT,
        )
        response = gmodel.generate_content(
            query,
            generation_config=genai.types.GenerationConfig(
                temperature=0.0,
                max_output_tokens=500,
                response_mime_type="application/json",
            ),
        )
        content = response.text.strip()
        plan = json.loads(content)
        plan["raw_query"] = query
        plan["parser"] = "llm"
        plan["model"] = model_name
        for key, default in DEFAULT_PLAN.items():
            if key not in plan:
                plan[key] = default
        return plan
    except ImportError:
        return None
    except Exception as e:
        logger.warning(f"Google parser failed: {e}")
        return None


def parse_query_llm(query: str, provider: str = "",
                    api_key: str = "", model: str = "") -> Optional[dict]:
    """
    LLM-enhanced query parser. Supports OpenAI and Google Gemini.
    Falls back to None if unavailable.
    """
    provider = provider or os.environ.get("AI_PROVIDER", "openai")

    if provider in ("google", "free", "gemini"):
        result = _parse_google(query, api_key, model)
        if result:
            return result
        return _parse_openai(query, api_key, model)
    elif provider == "openai":
        result = _parse_openai(query, api_key, model)
        if result:
            return result
        return _parse_google(query, api_key, model)
    else:
        result = _parse_openai(query, api_key, model)
        if result:
            return result
        return _parse_google(query, api_key, model)


def parse_query(query: str, provider: str = "",
                api_key: str = "", model: str = "") -> dict:
    """
    Parse a natural language query into a structured plan.
    Uses LLM if available, otherwise falls back to rule-based.
    """
    if provider == "template":
        return parse_query_rules(query)
    llm_plan = parse_query_llm(query, provider=provider, api_key=api_key, model=model)
    if llm_plan:
        return llm_plan
    return parse_query_rules(query)
