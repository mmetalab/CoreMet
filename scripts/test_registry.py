#!/usr/bin/env python3
"""Test the CoreMet ID registry."""
import sys, os, time
sys.path.insert(0, "/Users/cheng.wang/Documents/mpi-web/MPI_web_server/mpi-vgae-web")
os.chdir("/Users/cheng.wang/Documents/mpi-web/MPI_web_server/mpi-vgae-web")

# Delete old registry to force rebuild
reg_file = "data/coremetdb_entity_registry.json"
if os.path.exists(reg_file):
    os.remove(reg_file)

from app.services.entity_registry import build_registry, lookup_id, get_entity, get_type_stats

t0 = time.time()
build_registry()
print(f"Registry built in {time.time()-t0:.1f}s")
print(f"Stats: {get_type_stats()}")

# Test lookups
tests = [
    ("Butyric acid", "metabolite"),
    ("HMDB0000039", None),
    ("ALB", "gene"),
    ("Albumin", "protein"),
    ("Colorectal cancer", "disease"),
    ("Metformin", "drug"),
    ("Lactobacillus", "microbe"),
    ("rs1260326", "snp"),
    ("LDHA", "gene"),
    ("P02768", None),
]
for name, etype in tests:
    cid = lookup_id(name, etype)
    ent = get_entity(cid) if cid else None
    print(f"  {name:25s} ({etype or 'any':12s}) -> {cid or 'NOT FOUND':12s}  {ent['name'] if ent else ''}")

# Verify file was saved
import os
print(f"\nRegistry file: {reg_file} ({os.path.getsize(reg_file)/1024/1024:.1f} MB)")
