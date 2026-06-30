#!/usr/bin/env python3
"""
compute_db_stats.py — Single source of truth for CoreMet statistics.

Reads the v1.0 release database CSVs and writes data/coremetdb_stats.json.
Every downstream artifact (manuscript tables, publication figures, web home-page
stats) must read its numbers from this JSON so they can never drift again.

Usage:
    conda activate mpi-vgae
    python scripts/compute_db_stats.py
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
OUT_PATH = DATA_DIR / "coremetdb_stats.json"
CHUNK = 200_000

# Curated disease-network explorer count (distinct from raw MDI disease count).
# Kept explicit so the two are never conflated again.
CURATED_DISEASE_NETWORKS = 112

# Per-database configuration. Column names differ across files (MPI uses spaces).
# metabolite = HMDB id column; target = the partner-entity id column.
DB_CONFIG = {
    "MPI": {
        "label": "Metabolite–Protein",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mpi.csv",
        "metabolite": "HMDB ID",
        "target": "Uniprot ID",
        "target_label": "proteins",
        "organism": "Species",
        "pathway": "Pathway_Name",
        "evidence": "Evidence_Source",
        "pmid": None,
        "dedup": ["HMDB ID", "Uniprot ID", "Species"],
        # Curated-only: keep rows with a real Evidence_Source (KEGG/Rhea/original)
        # or Rhea 'reaction_participant'; drop the 112,743 unsourced 'pathway_participant'
        # co-membership edges (no source/species/protein name). Note: curated KEGG rows
        # are ALSO labelled pathway_participant, so Evidence_Source is the discriminator.
        "keep_cols": ["Evidence_Source", "interaction_subtype"],
        "keep_mask": lambda c: c["Evidence_Source"].notna() | (c["interaction_subtype"] == "reaction_participant"),
        "source": "KEGG, Rhea",
    },
    "MEI": {
        "label": "Metabolite–Enzyme",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mei.csv",
        "metabolite": "HMDB_ID",
        "target": "EC_Number",
        "target_label": "EC numbers",
        "organism": "Species",
        "pathway": "Pathway_Name",
        "evidence": "Evidence_Source",
        "pmid": None,
        "dedup": ["HMDB_ID", "EC_Number", "Species"],
        "source": "KEGG, Rhea",
    },
    "MDI": {
        "label": "Metabolite–Disease",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mdi.csv",
        "metabolite": "HMDB_ID",
        "target": "Disease_Name",
        "target_label": "diseases",
        "organism": None,
        "evidence": "Evidence_Level",
        "pmid": "pmid",
        "dedup": ["HMDB_ID", "Disease_Name"],
        "source": "CTD",
    },
    "MMI": {
        "label": "Metabolite–Microbe",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mmi.csv",
        "metabolite": "HMDB_ID",
        "target": "Microbe_Name",
        "target_label": "microbes",
        "organism": "Organism",
        "evidence": "Evidence_Level",
        "pmid": "PMID",
        "dedup": ["HMDB_ID", "Microbe_Name", "Relationship_Type"],
        "source": "gutMGene, AGORA2",
    },
    "MDrI": {
        "label": "Metabolite–Drug",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mdri.csv",
        "metabolite": "HMDB_ID",
        "target": "Drug_Name",
        "target_label": "drugs",
        "organism": None,
        "evidence": "Evidence_Level",
        "pmid": "PMID",
        "dedup": ["HMDB_ID", "Drug_Name", "Interaction_Type"],
        "source": "DrugBank",
    },
    "MGI": {
        "label": "Metabolite–Gene",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mgi.csv",
        "metabolite": "HMDB_ID",
        "target": "Gene_Symbol",
        "target_label": "genes",
        "organism": "Organism",
        "evidence": None,
        "pmid": "PMID",
        "dedup": ["HMDB_ID", "Gene_Symbol", "Organism", "Interaction_Type"],
        "source": "CTD",
    },
    "mGWAS": {
        "label": "Metabolite–SNP",
        "path": DATA_DIR / "databases" / "release" / "coremetdb_mgwas.csv",
        "metabolite": "HMDB_ID",
        "target": "rsID",
        "target_label": "SNPs",
        "organism": None,
        "evidence": None,
        "pmid": "PMID",
        "dedup": ["HMDB_ID", "rsID"],
        "source": "GWAS Catalog",
    },
}


def _clean(series: pd.Series) -> pd.Series:
    """Drop NaN / blank / common null sentinels and return string values."""
    s = series.dropna().astype(str).str.strip()
    return s[~s.isin({"", "nan", "NaN", "None", "-", "NA", "n/a"})]


def process_db(name: str, cfg: dict) -> tuple[dict, set]:
    path = cfg["path"]
    if not path.exists():
        raise FileNotFoundError(f"{name}: missing {path}")

    usecols = [cfg["metabolite"], cfg["target"]]
    for key in ("organism", "pathway", "evidence", "pmid"):
        col = cfg.get(key)
        if col:
            usecols.append(col)
    dedup_cols = cfg.get("dedup") or []
    usecols.extend(dedup_cols)
    keep_mask = cfg.get("keep_mask")
    usecols.extend(cfg.get("keep_cols") or [])
    usecols = list(dict.fromkeys(usecols))  # dedupe, preserve order

    n_rows = 0
    n_excluded = 0
    metabolites: set[str] = set()
    targets: set[str] = set()
    organisms: set[str] = set()
    pathways: set[str] = set()
    evidence = Counter()
    pmid_present = 0
    edge_keys: set[str] = set()

    reader = pd.read_csv(
        path, usecols=usecols, dtype=str, chunksize=CHUNK,
        keep_default_na=True, on_bad_lines="skip", low_memory=False,
    )
    for chunk in reader:
        if keep_mask is not None:
            before = len(chunk)
            chunk = chunk[keep_mask(chunk)]
            n_excluded += before - len(chunk)
        n_rows += len(chunk)
        metabolites.update(_clean(chunk[cfg["metabolite"]]))
        targets.update(_clean(chunk[cfg["target"]]))
        if cfg.get("organism"):
            organisms.update(_clean(chunk[cfg["organism"]]))
        if cfg.get("pathway"):
            for cell in _clean(chunk[cfg["pathway"]]):
                for p in str(cell).split(";"):
                    p = p.strip()
                    if p:
                        pathways.add(p)
        if cfg.get("evidence"):
            evidence.update(_clean(chunk[cfg["evidence"]]))
        if cfg.get("pmid"):
            pmid_present += len(_clean(chunk[cfg["pmid"]]))
        if dedup_cols:
            keys = chunk[dedup_cols].fillna("").astype(str)
            edge_keys.update(keys.agg("\x1f".join, axis=1).tolist())

    n_edges = len(edge_keys) if dedup_cols else n_rows
    stat = {
        "label": cfg["label"],
        "source": cfg["source"],
        "interactions": n_edges,
        "raw_rows": n_rows,
        "excluded_rows": n_excluded if keep_mask is not None else 0,
        "dedup_ratio": round(n_edges / n_rows, 4) if n_rows else None,
        "metabolites": len(metabolites),
        "target_label": cfg["target_label"],
        "targets": len(targets),
        "pmid_coverage": round(pmid_present / n_rows, 4) if (cfg.get("pmid") and n_rows) else None,
    }
    if cfg.get("organism"):
        stat["organisms"] = len(organisms)
    if cfg.get("pathway"):
        stat["pathways"] = len(pathways)
    if cfg.get("evidence"):
        stat["evidence_distribution"] = dict(evidence.most_common())

    print(f"  {name:6s} edges={n_edges:>9,}  (raw={n_rows:>9,}, keep={stat['dedup_ratio']})  "
          f"metabolites={stat['metabolites']:>6,}  {cfg['target_label']}={stat['targets']:>7,}")
    return stat, metabolites


def main() -> None:
    print("Computing CoreMet statistics from v1.0 release CSVs ...")
    per_db: dict[str, dict] = {}
    metabolite_sets: dict[str, set] = {}

    for name, cfg in DB_CONFIG.items():
        stat, mets = process_db(name, cfg)
        per_db[name] = stat
        metabolite_sets[name] = mets

    total_interactions = sum(d["interactions"] for d in per_db.values())

    # Cross-type metabolite overlap: how many DB-types each metabolite appears in.
    appearances = Counter()
    for mets in metabolite_sets.values():
        appearances.update(mets)
    overlap_hist = Counter(appearances.values())  # {n_types: n_metabolites}
    union_metabolites = len(appearances)
    in_ge_3 = sum(c for k, c in overlap_hist.items() if k >= 3)
    in_all = overlap_hist.get(len(DB_CONFIG), 0)

    # Weighted overall PMID coverage across DBs that carry PMIDs.
    pmid_rows = sum(d["interactions"] for n, d in per_db.items() if d["pmid_coverage"] is not None)
    pmid_hits = sum(int(d["pmid_coverage"] * d["interactions"]) for d in per_db.values()
                    if d["pmid_coverage"] is not None)

    stats = {
        "version": "v1.0",
        "totals": {
            "interactions": total_interactions,
            "database_types": len(DB_CONFIG),
            "unique_metabolites_union": union_metabolites,
            "curated_disease_networks": CURATED_DISEASE_NETWORKS,
            "pmid_coverage_overall": round(pmid_hits / pmid_rows, 4) if pmid_rows else None,
        },
        "databases": per_db,
        "cross_type_overlap": {
            "metabolites_by_db_count": {str(k): overlap_hist[k] for k in sorted(overlap_hist)},
            "metabolites_in_ge_3_types": in_ge_3,
            "metabolites_in_all_types": in_all,
        },
        "_note": (
            "Generated by scripts/compute_db_stats.py for the CoreMet v1.0 release. "
            "Interaction counts match the release CSVs; raw_rows/excluded_rows capture "
            "source-lineage filtering where available. Do not hand-edit. "
            "'curated_disease_networks' (112) is the disease-network explorer count "
            "and is DISTINCT from MDI 'diseases' (raw target count)."
        ),
    }

    OUT_PATH.write_text(json.dumps(stats, indent=2))
    print(f"\nTotal interactions: {total_interactions:,}")
    print(f"Union metabolites:  {union_metabolites:,}  (>=3 types: {in_ge_3:,}, all types: {in_all})")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
