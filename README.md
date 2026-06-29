# CoreMet

**A metabolite-centered cross-domain interaction knowledge graph.**

CoreMet integrates and harmonizes seven types of metabolite interaction into one
interoperable, evidence-aware web resource. Starting from any metabolite (or a disease,
drug, gene, or variant) you can traverse outward across proteins, enzymes, diseases,
microbes, drugs, genes, and genomic variants in a single query.

| | |
|---|---|
| **Interactions** | 1,952,688 (7 curated interaction types) |
| **Metabolites** | 30,674 (HMDB-anchored; 60 present in all 7 layers) |
| **Provenance** | 5-level evidence hierarchy, normalized confidence, PubMed refs |
| **Access** | Web UI, 15 REST endpoints, bulk CSV download, no login |
| **Stack** | Python, Dash/Flask, deployed on Render |

| Layer | Type | Interactions | Source |
|-------|------|-------------:|--------|
| MPI | Metabolite–Protein | 38,061 | KEGG, Rhea |
| MEI | Metabolite–Enzyme | 47,551 | KEGG, Rhea |
| MDI | Metabolite–Disease | 82,882 | CTD, HMDB |
| MMI | Metabolite–Microbe | 77,605 | gutMGene, AGORA2 |
| MDrI | Metabolite–Drug | 3,500 | DrugBank |
| MGI | Metabolite–Gene | 1,658,745 | CTD |
| mGWAS | Metabolite–SNP | 44,344 | GWAS Catalog |

## Features

- **Browse** each interaction layer with evidence-aware filters (evidence type, confidence, species) and CSV export.
- **Search** any metabolite, protein, gene, disease, microbe, drug, or SNP.
- **Cross-layer profile** page aggregating every interaction for a metabolite in one view.
- **Network** exploration with an interactive multi-type Cytoscape graph.
- **Downloads** of clean per-layer edge lists whose counts match the publication exactly.
- **REST API** for programmatic access.
- **Prediction (beta)**, an optional graph-neural-network tool for hypothesis generation. The curated database is the core resource.

## Repository layout

```
app/            Dash app factory, config, services (data access per layer)
pages/          Page layouts (home, database browser, profile, network, downloads, API docs)
pages/modules/  Per-layer browse modules (mpi, mei, mdi, mmi, mdri, mgi, mgwas)
components/     Navbar, footer, page header, stat cards
scripts/        Reproducible data pipeline, figure/stat generators, fetch_data.py
assets/         CSS theme and static assets
api/            REST API blueprint
tests/          pytest suite
data/           NOT in git (hosted on Zenodo); only coremetdb_stats.json is committed
render.yaml     Render deployment blueprint
```

## Run locally

```bash
conda create -n coremet python=3.10 -y && conda activate coremet
pip install -r requirements.txt

# Fetch the runtime data (Zenodo URL from DEPLOY_RENDER.md):
DATA_BUNDLE_URL="https://zenodo.org/records/21032647/files/coremet_runtime_core.zip?download=1" \
    python scripts/fetch_data.py

python run.py            # serves on http://localhost:8080
```

The full `data/` tree (~15 GB of raw sources) is **not** committed. Only the small
`data/coremetdb_stats.json` is tracked; the runtime data (release CSVs plus the entity
registry, about 380 MB) is downloaded by `scripts/fetch_data.py` from Zenodo. See
**[DATA_README.md](DATA_README.md)** for the data dictionary.

## Deploy

See **[DEPLOY_RENDER.md](DEPLOY_RENDER.md)** for the full GitHub, Zenodo, and Render
protocol with exact build/start commands, environment variables, plan, and parameters.

## Reproducing the data and figures

```bash
python scripts/compute_db_stats.py        # data/coremetdb_stats.json (single source of truth)
python scripts/build_release_csvs.py      # data/databases/release/*.csv (clean, paper-matching)
python scripts/generate_database_figures.py
```

## License

Code: MIT (see [LICENSE](LICENSE)). Data: CC BY 4.0.

## Citation

A formal citation will be available once the CoreMet manuscript is published
(in preparation, 2026). Please link to the resource URL in the meantime.
