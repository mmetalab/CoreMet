"""Quick test of the CoreMet-AI pipeline."""
import sys
import json

from app.services.ai_orchestrator import execute_ai_query

# Test 1: Relation query
print("=== Test 1: Relation query ===")
r = execute_ai_query("How does butyrate influence colorectal cancer?")
print(f"Status: {r['status']}")
print(f"Source: {r['query_plan'].get('source_entity')} ({r['query_plan'].get('source_type')})")
print(f"Target: {r['query_plan'].get('target_entity')} ({r['query_plan'].get('target_type')})")
print(f"Intent: {r['query_plan'].get('intent')}")
print(f"Resolved source: {r.get('resolved_source')}")
print(f"Resolved target: {r.get('resolved_target')}")
print(f"Ranked paths: {len(r.get('ranked_paths', []))}")
print(f"Subgraph nodes: {len(r.get('subgraph', {}).get('nodes', []))}")
print(f"Subgraph edges: {len(r.get('subgraph', {}).get('edges', []))}")
print(f"Evidence edges: {r.get('evidence', {}).get('total_edges', 0)}")
print(f"PMIDs: {r.get('evidence', {}).get('total_pmids', 0)}")
print(f"Summary: {r.get('summary', {}).get('summary', '')[:300]}")
print(f"Confidence: {r.get('summary', {}).get('confidence_label')}")
print(f"Timing: {r.get('timing_ms')} ms")
print()

# Test 2: Entity query
print("=== Test 2: Entity query ===")
r2 = execute_ai_query("Tell me about metformin")
print(f"Status: {r2['status']}")
print(f"Source: {r2['query_plan'].get('source_entity')} ({r2['query_plan'].get('source_type')})")
print(f"Ranked paths: {len(r2.get('ranked_paths', []))}")
print(f"Subgraph nodes: {len(r2.get('subgraph', {}).get('nodes', []))}")
print(f"Summary: {r2.get('summary', {}).get('summary', '')[:200]}")
print(f"Timing: {r2.get('timing_ms')} ms")
print()

# Test 3: Gene query
print("=== Test 3: Gene-disease query ===")
r3 = execute_ai_query("What is the relationship between CYP1A2 and obesity?")
print(f"Status: {r3['status']}")
print(f"Source: {r3['query_plan'].get('source_entity')} ({r3['query_plan'].get('source_type')})")
print(f"Target: {r3['query_plan'].get('target_entity')} ({r3['query_plan'].get('target_type')})")
print(f"Ranked paths: {len(r3.get('ranked_paths', []))}")
print(f"Subgraph nodes: {len(r3.get('subgraph', {}).get('nodes', []))}")
print(f"Summary: {r3.get('summary', {}).get('summary', '')[:200]}")
print(f"Timing: {r3.get('timing_ms')} ms")
print()

# Test 4: SNP query
print("=== Test 4: SNP query ===")
r4 = execute_ai_query("Connect rs1260326 to metabolic syndrome")
print(f"Status: {r4['status']}")
print(f"Source: {r4['query_plan'].get('source_entity')} ({r4['query_plan'].get('source_type')})")
print(f"Ranked paths: {len(r4.get('ranked_paths', []))}")
print(f"Subgraph nodes: {len(r4.get('subgraph', {}).get('nodes', []))}")
print(f"Summary: {r4.get('summary', {}).get('summary', '')[:200]}")
print(f"Timing: {r4.get('timing_ms')} ms")
print()

# Test 5: Microbe query
print("=== Test 5: Microbe query ===")
r5 = execute_ai_query("Role of Lactobacillus in short-chain fatty acids")
print(f"Status: {r5['status']}")
print(f"Ranked paths: {len(r5.get('ranked_paths', []))}")
print(f"Subgraph nodes: {len(r5.get('subgraph', {}).get('nodes', []))}")
print(f"Summary: {r5.get('summary', {}).get('summary', '')[:200]}")
print(f"Timing: {r5.get('timing_ms')} ms")

print("\n=== ALL TESTS COMPLETE ===")
