#!/usr/bin/env python3
"""Build HMDB → PubChem CID / ChEBI ID mapping via PubChem REST API.

Strategy:
  1. For each unique (HMDB_ID, SMILES) pair, resolve SMILES → PubChem CID
  2. Batch CID → ChEBI via PubChem xrefs/RegistryID endpoint
  3. Save mapping CSV: hmdb_id, metabolite_name, smiles, pubchem_cid, chebi_id

Rate limiting: PubChem allows ~5 requests/second for individual lookups.
We use asyncio + aiohttp for concurrent but throttled requests.
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

# Project paths
BASE_DIR = Path(__file__).resolve().parent.parent
DB_PATH = BASE_DIR / "data" / "mpidatabase" / "MPIDB_v2.csv"
OUTPUT_DIR = BASE_DIR / "data" / "mappings"
OUTPUT_CSV = OUTPUT_DIR / "hmdb_xref_mapping.csv"
CHECKPOINT_JSON = OUTPUT_DIR / "hmdb_xref_checkpoint.json"

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
REQUEST_DELAY = 0.22  # ~4.5 req/sec, safely under PubChem's 5/sec limit
XREF_BATCH_SIZE = 80  # CIDs per batch xref request


def load_unique_metabolites(db_path: str) -> List[dict]:
    """Load unique (HMDB ID, SMILES, Name) from MPIDB."""
    seen = set()
    results = []
    with open(db_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            hmdb = row.get("HMDB ID", "").strip()
            if not hmdb or hmdb in seen:
                continue
            seen.add(hmdb)
            results.append({
                "hmdb_id": hmdb,
                "name": row.get("Metabolite Name", "").strip(),
                "smiles": row.get("SMILES", "").strip(),
            })
    return results


def load_checkpoint() -> dict:
    """Load checkpoint of already-resolved mappings."""
    if CHECKPOINT_JSON.exists():
        with open(CHECKPOINT_JSON) as f:
            return json.load(f)
    return {}


def save_checkpoint(mapping: dict):
    """Save checkpoint."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(CHECKPOINT_JSON, "w") as f:
        json.dump(mapping, f, indent=1)


def smiles_to_cid(smiles: str) -> Optional[int]:
    """Resolve SMILES → PubChem CID via PubChem REST API."""
    if not smiles:
        return None
    encoded = urllib.parse.quote(smiles, safe="")
    url = f"{PUBCHEM_BASE}/smiles/{encoded}/property/InChIKey/JSON"
    req = urllib.request.Request(url)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
            return data["PropertyTable"]["Properties"][0]["CID"]
    except Exception:
        return None


def name_to_cid(name: str) -> Optional[int]:
    """Fallback: resolve metabolite name → PubChem CID."""
    if not name:
        return None
    encoded = urllib.parse.quote(name, safe="")
    url = f"{PUBCHEM_BASE}/name/{encoded}/property/InChIKey/JSON"
    req = urllib.request.Request(url)
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
            return data["PropertyTable"]["Properties"][0]["CID"]
    except Exception:
        return None


def batch_cid_to_chebi(cids: list[int]) -> Dict[int, str]:
    """Batch lookup CID → ChEBI ID via PubChem xrefs."""
    if not cids:
        return {}
    cid_str = ",".join(str(c) for c in cids)
    url = f"{PUBCHEM_BASE}/cid/{cid_str}/xrefs/RegistryID/JSON"
    req = urllib.request.Request(url)
    result = {}
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
            for info in data.get("InformationList", {}).get("Information", []):
                cid = info["CID"]
                regs = info.get("RegistryID", [])
                chebi_matches = [r for r in regs if r.upper().startswith("CHEBI:")]
                if chebi_matches:
                    result[cid] = chebi_matches[0]
    except Exception:
        pass
    return result


def main():
    print(f"Loading metabolites from {DB_PATH}...")
    metabolites = load_unique_metabolites(str(DB_PATH))
    print(f"Found {len(metabolites)} unique HMDB IDs")

    # Load checkpoint
    checkpoint = load_checkpoint()
    print(f"Checkpoint has {len(checkpoint)} resolved entries")

    # Phase 1: SMILES → PubChem CID
    print("\n=== Phase 1: SMILES → PubChem CID ===")
    resolved = 0
    failed = 0
    skipped = 0
    for i, met in enumerate(metabolites):
        hmdb = met["hmdb_id"]
        if hmdb in checkpoint and checkpoint[hmdb].get("pubchem_cid"):
            skipped += 1
            continue

        # Try SMILES first, then name fallback
        cid = smiles_to_cid(met["smiles"])
        if cid is None:
            cid = name_to_cid(met["name"])
        time.sleep(REQUEST_DELAY)

        if cid is not None:
            checkpoint[hmdb] = {
                "name": met["name"],
                "smiles": met["smiles"],
                "pubchem_cid": cid,
                "chebi_id": None,
            }
            resolved += 1
        else:
            checkpoint[hmdb] = {
                "name": met["name"],
                "smiles": met["smiles"],
                "pubchem_cid": None,
                "chebi_id": None,
            }
            failed += 1

        if (i + 1) % 50 == 0:
            save_checkpoint(checkpoint)
            total_with_cid = sum(1 for v in checkpoint.values() if v.get("pubchem_cid"))
            print(f"  [{i+1}/{len(metabolites)}] resolved={total_with_cid}, failed={failed}")

    save_checkpoint(checkpoint)
    total_with_cid = sum(1 for v in checkpoint.values() if v.get("pubchem_cid"))
    print(f"\nPhase 1 complete: {total_with_cid} CIDs resolved, {failed} failed, {skipped} skipped")

    # Phase 2: Batch CID → ChEBI
    print("\n=== Phase 2: Batch CID → ChEBI ===")
    cid_to_hmdb = {}
    for hmdb, info in checkpoint.items():
        cid = info.get("pubchem_cid")
        if cid and not info.get("chebi_id"):
            cid_to_hmdb[cid] = hmdb

    cids_needing_chebi = list(cid_to_hmdb.keys())
    print(f"  {len(cids_needing_chebi)} CIDs need ChEBI resolution")

    chebi_found = 0
    for batch_start in range(0, len(cids_needing_chebi), XREF_BATCH_SIZE):
        batch = cids_needing_chebi[batch_start:batch_start + XREF_BATCH_SIZE]
        chebi_map = batch_cid_to_chebi(batch)
        for cid, chebi in chebi_map.items():
            hmdb = cid_to_hmdb[cid]
            checkpoint[hmdb]["chebi_id"] = chebi
            chebi_found += 1
        time.sleep(REQUEST_DELAY * 2)

        if (batch_start // XREF_BATCH_SIZE + 1) % 5 == 0:
            save_checkpoint(checkpoint)
            print(f"  Batch {batch_start // XREF_BATCH_SIZE + 1}: {chebi_found} ChEBI IDs found so far")

    save_checkpoint(checkpoint)
    print(f"\nPhase 2 complete: {chebi_found} ChEBI IDs resolved")

    # Phase 3: Write final CSV
    print("\n=== Phase 3: Writing output CSV ===")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["hmdb_id", "metabolite_name", "smiles", "pubchem_cid", "chebi_id"])
        for hmdb, info in sorted(checkpoint.items()):
            writer.writerow([
                hmdb,
                info.get("name", ""),
                info.get("smiles", ""),
                info.get("pubchem_cid", ""),
                info.get("chebi_id", ""),
            ])

    total_cid = sum(1 for v in checkpoint.values() if v.get("pubchem_cid"))
    total_chebi = sum(1 for v in checkpoint.values() if v.get("chebi_id"))
    print(f"\nFinal mapping saved to {OUTPUT_CSV}")
    print(f"  Total: {len(checkpoint)} metabolites")
    print(f"  PubChem CID: {total_cid} ({100*total_cid/len(checkpoint):.1f}%)")
    print(f"  ChEBI ID: {total_chebi} ({100*total_chebi/len(checkpoint):.1f}%)")


if __name__ == "__main__":
    main()
