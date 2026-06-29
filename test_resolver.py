"""Test entity resolution fixes."""
import app.services.entity_resolver as er
er._synonym_index = None  # Reset cache

from app.services.entity_resolver import resolve_entity

tests = [
    ("butyrate", None),
    ("butyric acid", None),
    ("colorectal cancer", "disease"),
    ("short-chain fatty acids", None),
    ("Lactobacillus", None),
    ("metformin", None),
    ("tryptophan", None),
    ("rs1260326", None),
    ("CYP1A2", None),
    ("depression", "disease"),
]

for query, etype in tests:
    results = resolve_entity(query, expected_type=etype)
    if results:
        r = results[0]
        print(f"{query:30s} -> {r['name']:30s} ({r['type']:12s}) id={r['id'][:20]} conf={r['confidence']}")
    else:
        print(f"{query:30s} -> NO MATCH")

# Clear caches for traversal tests
from app.services.graph_traversal import _layer_cache
_layer_cache.clear()
er._synonym_index = None

print("\n--- Full AI pipeline tests ---")
from app.services.ai_orchestrator import execute_ai_query

queries = [
    "How does butyrate influence colorectal cancer?",
    "Tell me about metformin",
    "What is the relationship between CYP1A2 and obesity?",
    "Connect rs1260326 to metabolic syndrome",
    "Role of Lactobacillus in short-chain fatty acids",
    "What drugs interact with glutathione?",
]

for q in queries:
    er._synonym_index = None
    _layer_cache.clear()
    r = execute_ai_query(q)
    src = r['query_plan'].get('source_entity', '?')
    srctype = r['query_plan'].get('source_type', '?')
    tgt = r['query_plan'].get('target_entity', '')
    tgttype = r['query_plan'].get('target_type', '')
    npaths = len(r.get('ranked_paths', []))
    nnodes = len(r.get('subgraph', {}).get('nodes', []))
    nedges = r.get('evidence', {}).get('total_edges', 0)
    conf = r.get('summary', {}).get('confidence_label', '?')
    ms = r.get('timing_ms', 0)
    tgt_str = f" -> {tgt} ({tgttype})" if tgt else ""
    print(f"Q: {q[:50]:50s} | {src} ({srctype}){tgt_str} | paths={npaths} nodes={nnodes} edges={nedges} [{conf}] {ms}ms")

print("\nDONE")
