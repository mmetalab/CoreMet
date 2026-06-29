# CoreMet — Feature Documentation

**Version**: 1.0  
**Last updated**: 2026-04-04  
**Server**: Python 3.9 + Dash 2.18 + dash-bootstrap-components + dash-cytoscape

---

## 1. Architecture Overview

CoreMet is a metabolite-centered knowledge graph web server integrating seven interaction databases spanning molecular, microbial, genetic, disease, and pharmacologic layers. The system is built as a single-page Dash application with server-side rendering and client-side interactivity.

### Technology Stack

| Component | Technology |
|-----------|-----------|
| Framework | Dash 2.18 (Plotly) |
| UI components | dash-bootstrap-components 1.x |
| Graph visualization | dash-cytoscape 1.0 |
| Data backend | Pandas DataFrames (in-memory) |
| AI pipeline | OpenAI / Google Gemini API |
| Styling | Custom CSS (Arial, Nature Journal style) |
| Deployment | Gunicorn, Render.com ready |

### Entity-Type Color System

All visualizations, badges, filters, and legends use a consistent color code:

| Entity Type | Color | Hex |
|------------|-------|-----|
| Metabolite | Orange | `#e27a3f` |
| Protein | Blue | `#3182ce` |
| Gene | Gold | `#d69e2e` |
| Disease | Red | `#e53e3e` |
| Microbe | Green | `#38a169` |
| Drug | Purple | `#805ad5` |
| SNP | Teal | `#319795` |

---

## 2. Navigation & Page Structure

**Navbar order**: CoreMet (logo) | CoreMet-AI | Search | Explore | Downloads | API | Help | Cite

### Page Inventory

| Route | Page | Purpose |
|-------|------|---------|
| `/home` | Homepage | Hero + search bar + stats + action cards |
| `/coremetai` | CoreMet-AI | Natural language graph query with AI |
| `/search` | Search | Global multi-entity search |
| `/explore` | Explore Hub | Browse by entity type or interaction layer |
| `/explore/<type>` | Entity Browse | Paginated list of entities of one type |
| `/mpi` | MPI Module | Metabolite–Protein Interactions table |
| `/mdi` | MDI Module | Metabolite–Disease Interactions table |
| `/mmi` | MMI Module | Metabolite–Microbe Interactions table |
| `/mdri` | MDrI Module | Metabolite–Drug Interactions table |
| `/mgi` | MGI Module | Metabolite–Gene Interactions table |
| `/mgwas` | mGWAS Module | Metabolite–SNP Associations table |
| `/metabolite?id=` | Metabolite Detail | Entity page with 5 tabs |
| `/disease-detail?name=` | Disease Detail | Entity page with 5 tabs |
| `/gene?name=` | Gene Detail | Entity page with 5 tabs |
| `/protein?name=` | Protein Detail | Entity page with 5 tabs |
| `/drug?name=` | Drug Detail | Entity page with 5 tabs |
| `/microbe?name=` | Microbe Detail | Entity page with 5 tabs |
| `/snp?name=` | SNP Detail | Entity page with 5 tabs |
| `/enrichment` | Enrichment Analysis | 5-type enrichment (pathway, disease, microbe, drug, gene) |
| `/network` | Network Viewer | Multi-type Cytoscape graph |
| `/downloads` | Downloads | Versioned dataset download cards |
| `/api-docs` | API Documentation | REST API reference with examples |
| `/help` | Help & About | Getting started, FAQ, entity definitions, citation |

---

## 3. Homepage

### Components (top to bottom)

1. **Hero block** — Title "CoreMet" with subtitle describing the knowledge graph
2. **Search bar** — Large autocomplete search across all 7 entity types with placeholder text
3. **Example chips** — Clickable examples: butyrate, metformin, rs1260326 (link to `/search?q=...`)
4. **Statistics cards** (8) — Total interactions, metabolites, genes, proteins, microbes, diseases, drugs, SNPs
5. **Action panels** (3) — "Search Entities", "Explore Networks", "Query with AI" with links
6. **Interaction layer section** — 7 colored chips linking to MPI/MDI/MMI/MDrI/MGI/mGWAS modules
7. **Citation card** — Copy-ready citation block

---

## 4. Search

### Features

- **Server-side pre-rendering**: When navigating to `/search?q=butyrate`, results are rendered server-side (no double-refresh)
- **Multi-entity search**: Searches across metabolites, diseases, genes, proteins, drugs, microbes, SNPs simultaneously
- **Grouped results**: Results grouped by entity type with color-coded badges and icons
- **Direct links**: Each result links to its entity detail page (`/metabolite?id=...`, `/gene?name=...`, etc.)
- **Example chips**: Quick-search examples for common queries

### Search Coverage

| Entity Type | Source Database | ID Column |
|------------|----------------|-----------|
| Metabolite | MPI | HMDB ID |
| Disease | MDI | Disease_Name |
| Gene | MGI | Gene_Symbol |
| Protein | MPI | Protein Name |
| Drug | MDrI | Drug_Name |
| Microbe | MMI | Microbe_Name |
| SNP | mGWAS | rsID |

---

## 5. Entity Detail Pages

All entity types share an identical layout structure for consistency.

### Page Header

- Primary name (e.g., "Butyric acid")
- Canonical identifier (e.g., HMDB0000039)
- Entity type badge (colored)
- External links (HMDB, KEGG, PubChem, ChEBI, UniProt, NCBI Gene, dbSNP, DrugBank — context-dependent)

### Summary Cards Row

Compact cards showing connected entity counts per layer with color-coded left borders:
- Linked proteins, genes, diseases, microbes, drugs, SNPs
- Total evidence count

### Five Tabs

#### Overview Tab
- Basic identifiers and cross-references
- Layer distribution bar chart (neighbors per interaction type)
- Top connected entities per layer

#### Network Tab
- **Interactive Cytoscape graph** centered on the entity
- **2-hop cross-layer traversal**: Center entity → Ring 1 (metabolites) → Ring 2 (genes, proteins, drugs, diseases, microbes, SNPs)
- `cose` layout triggered on tab activation (fixes hidden-tab rendering)
- Color-coded nodes by entity type
- Solid edges for hop 1, dashed for hop 2
- Clickable for node details

#### Interactions Tab
- Structured DataTable of all edges involving the entity
- Grouped by interaction layer
- Columns: target entity, target type, edge type, subtype, confidence, evidence type, source, PMID
- **Entity links**: Clickable metabolite/gene/disease/etc. names link to their detail pages
- **Source links**: Database names link to external source URLs
- **PMID links**: PubMed IDs link to PubMed
- Markdown rendering via `presentation: "markdown"` in DataTable

#### Evidence Tab
- Evidence type distribution chart (Plotly donut)
- Source database distribution chart
- Layer distribution chart
- Confidence histogram
- Supporting reference count
- Raw evidence table

#### Downloads Tab
- Export current entity interactions as CSV (one button per layer)
- API endpoint reference for programmatic access

---

## 6. Interaction Layer Modules

Each of the 7 interaction databases has a dedicated browse page (`/mpi`, `/mdi`, ..., `/mgwas`).

### Common Features (module_factory.py)

- **Statistics cards**: Total interactions, unique entities, organisms, etc.
- **Visualizations**: 3–5 Plotly charts (donut, bar, heatmap, histogram)
- **Filters sidebar**:
  - Free-text search
  - Evidence type dropdown
  - Confidence slider
  - **"Hide predicted/inferred" toggle** (filters out `inferred`, `computational`, `predicted`)
  - Organism dropdown (when applicable)
  - **Extra biological filters**: disease category, pathway, tissue, chromosome, experimental method, etc.
- **DataTable** with:
  - Entity links (entity names → entity detail pages)
  - Source links (database names → external URLs with version info)
  - PMID links (→ PubMed)
  - Pagination, sorting, filtering
  - CSV export
- **Reset filters button**

### Per-Module Data Sources (Verified)

| Module | Source Database | URL | Version |
|--------|---------------|-----|---------|
| MPI | KEGG | https://www.genome.jp/kegg/reaction/ | Release 2025 |
| MPI | Rhea | https://www.rhea-db.org/ | v139 (2024) |
| MDI | CTD | https://ctdbase.org/ | 2025 release |
| MDI | HMDB | https://hmdb.ca/ | v5.0 (2024) |
| MMI | gutMGene | http://bio-annotation.cn/gutmgene/ | v2.0 (2023) |
| MMI | AGORA2 | https://www.vmh.life/#microbe | VMH 2023 |
| MDrI | DrugBank | https://go.drugbank.com/releases/latest | v5.1 (2024) |
| MGI | CTD | https://ctdbase.org/ | 2025 release |
| mGWAS | GWAS Catalog | https://www.ebi.ac.uk/gwas/ | 2025-03-17, GRCh38.p14, dbSNP 156 |
| mGWAS | Shin et al. | https://pubmed.ncbi.nlm.nih.gov/24816252/ | Nat Genet 46:543–550 (2014) |
| mGWAS | Yin et al. | https://pubmed.ncbi.nlm.nih.gov/38120091/ | Nat Genet (2024) |

### Evidence Classification per Database

| Module | Rows | Experimental | Curated | Computational | Predicted |
|--------|------|-------------|---------|---------------|-----------|
| MPI | 37,596 | — | 37,596 | — | — |
| MDI | 82,882 | — | 82,879 | — | 3 |
| MMI | 83,149 | 2,323 | — | 80,826 | — |
| MDrI | 3,500 | — | 3,500 | — | — |
| MGI | 1,658,745 | 1,474,517 | 184,228 | — | — |
| mGWAS | 44,344 | 44,344 | — | — | — |

### Biological Filters per Module

| Module | Extra Filters |
|--------|---------------|
| MPI | Pathway, EC Number, Evidence Source, Organism (10 species) |
| MDI | Disease, Association Type, Evidence Level |
| MMI | Microbe Rank, Tissue, Relationship Type, Experimental Method, Organism |
| MDrI | Interaction Type, Tissue, Evidence Level |
| MGI | Interaction Type (132 terms), Interaction Actions, Organism (7 species) |
| mGWAS | Chromosome (23), Source, Trait, Trait Category, Significance Tier |

---

## 7. Explore Hub

### Three Browse Modes

1. **Browse by entity type** — 7 cards (Metabolites, Genes, Proteins, Diseases, Microbes, Drugs, SNPs) linking to `/explore/<type>` pages with paginated entity lists
2. **Browse by interaction layer** — 7 cards with descriptions linking to `/mpi`, `/mdi`, etc.
3. **Network explorer** — Link to `/network` for open-ended graph exploration

### Entity Browse Pages (`/explore/<type>`)

Each entity type has a paginated browse page showing:
- Entity name (linked to detail page)
- Canonical ID
- Number of total interactions
- Top interaction count

---

## 8. CoreMet-AI

### Overview

CoreMet-AI is a query compiler for biological reasoning, not a chatbot. Users write natural language questions, and the system translates them into structured graph traversals, retrieves evidence-grounded subgraphs, ranks results, and generates strictly grounded explanations.

### 7-Stage Pipeline

| Stage | Function | Implementation |
|-------|----------|---------------|
| 1. Intent parsing | NL → structured JSON query plan | Rule-based patterns + LLM enhancement |
| 2. Entity resolution | Names → internal database IDs | Synonym index across all 7 DBs, fuzzy matching, genus-level microbe matching |
| 3. Graph traversal | Execute subgraph retrieval | 1–3 hop paths across 8 layers (MPI + MEI + MDI + MMI + MDrI + MGI + mGWAS) |
| 4. Path ranking | Score and select top paths | Composite: confidence × evidence type × DB diversity × path length × PMID coverage |
| 5. Evidence aggregation | Compile provenance | Source database, confidence labels, supporting references |
| 6. LLM summarization | Generate grounded explanation | OpenAI / Google Gemini / template fallback |
| 7. UI rendering | Three synchronized panels | Answer + Query Plan + Subgraph & Evidence |

### AI Model Selection

Users can select from:

| Model | Provider | Tier | Requirement |
|-------|----------|------|-------------|
| Gemini 2.0 Flash | Google | Free | `GOOGLE_API_KEY` in `.env` or entered in UI |
| GPT-4o-mini | OpenAI | API key | `OPENAI_API_KEY` |
| GPT-4o | OpenAI | API key | `OPENAI_API_KEY` |
| GPT-4-turbo | OpenAI | API key | `OPENAI_API_KEY` |
| Gemini Pro | Google | API key | `GOOGLE_API_KEY` |
| Template (no LLM) | — | Free | No API key needed |

### UI Layout

- **Model selection bar** — Dropdown + API key input field
- **Query input** — Large text box with example prompts
- **3-panel results**:
  - Center: AI-generated summary (grounded, 3–6 sentences)
  - Right: Parsed query plan (source entity, target entity, layers, hops, filters)
  - Bottom: Interactive Cytoscape subgraph + sortable evidence table

### API Endpoint

```
POST /api/v1/ai/query
Content-Type: application/json

{
  "query": "How might butyrate influence colorectal cancer?",
  "provider": "openai",
  "model": "gpt-4o-mini"
}
```

---

## 9. Downloads & API

### Downloads Page

- Versioned download cards for each database layer
- Format: CSV
- Record counts shown per dataset
- CC BY 4.0 licensing
- ML resources section (embeddings, decoders)

### REST API

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/health` | GET | Health check |
| `/api/v1/species` | GET | Available organisms |
| `/api/v1/autocomplete?q=` | GET | Entity autocomplete |
| `/api/v1/mmi/stats` | GET | MMI statistics |
| `/api/v1/mgi/stats` | GET | MGI statistics |
| `/api/v1/mgwas/stats` | GET | mGWAS statistics |
| `/api/v1/ai/query` | POST | CoreMet-AI query |

---

## 10. Help & Documentation

Includes:
- Getting started guide
- Definition of each interaction type (MPI, MEI, MDI, MMI, MDrI, MGI, mGWAS)
- Evidence hierarchy explanation (experimental > curated > computational > predicted)
- Confidence score methodology
- Identifier mapping (HMDB, InChIKey, SMILES, UniProt, Gene_ID, rsID, DrugBank_ID)
- FAQ
- Contact and citation instructions

---

## 11. Design Principles

1. **Journal-like aesthetic** — Arial font, white background, minimal shadows, flat cards
2. **Entity-type color consistency** — Same colors in graphs, chips, filters, legends, statistics
3. **Every page answers**: "Given this entity, what are its strongest cross-layer biological connections?"
4. **Evidence transparency** — Every interaction has an evidence type, source database, and (where available) PMID
5. **Hide predicted toggle** — Users can filter out inferred/computational interactions for higher-confidence views
6. **Source link verification** — Every source database links to its official URL with version/year

---

## 12. Configuration

### Environment Variables (`.env`)

```
OPENAI_API_KEY=sk-...          # OpenAI API key for CoreMet-AI
OPENAI_MODEL=gpt-4o-mini       # Default OpenAI model
GOOGLE_API_KEY=AIza...          # Google API key for free tier
GOOGLE_MODEL=gemini-2.0-flash  # Default Google model
```

### Starting the Server

```bash
cd /path/to/mpi-vgae-web
conda activate mpi-vgae
python run.py
```

Server runs at http://localhost:8080

---

## 13. Next To-Do List (Ranked by Priority)

### CRITICAL (before Paper B submission)

| # | Task | Effort | Impact |
|---|------|--------|--------|
| 1 | Real web screenshots for Paper B figures (incl. CoreMet-AI demo) | 2 hrs | Very high |
| 2 | Figure regeneration — Paper B Fig 1 & 2 with final numbers | 2 hrs | Very high |
| 3 | Expand PubChem/ChEBI cross-refs via UniChem/PUG batch mapping | 1 day | High |
| 4 | Obtain Google Gemini free API key for CoreMet-AI free tier | 10 min | High |

### HIGH PRIORITY

| # | Task | Effort | Impact |
|---|------|--------|--------|
| 5 | Add MiMeDB for MMI (manual download from mimedb.org) | 2 hrs | High |
| 6 | Add MarkerDB clinical biomarkers for MDI | 2 hrs | High |
| 7 | mGWAS expansion — Chen et al. 2023 UK Biobank (Nature Genetics) | 1 day | High |
| 8 | mGWAS expansion — Schlosser et al. 2020 kidney metabolite GWAS | 4 hrs | Medium-high |
| 9 | DisGeNET for MDI — additional disease associations | 1 day | Medium-high |

### MEDIUM PRIORITY

| # | Task | Effort | Impact |
|---|------|--------|--------|
| 10 | SABIO-RK kinetics for MEI (fix SBML parsing) | 4 hrs | Medium |
| 11 | BRENDA enzyme kinetics for MEI | 2 days | Medium |
| 12 | STITCH chemical-gene for MGI validation | 1 day | Medium |
| 13 | SMPDB disease pathways for MDI context | 4 hrs | Medium |
| 14 | Foundation model retrain on expanded databases | 2 hrs | Medium |

### LOW PRIORITY (post-submission)

| # | Task | Effort | Impact |
|---|------|--------|--------|
| 15 | Docker deployment config | 1 day | Low |
| 16 | Metabolite chemical class annotations | 4 hrs | Low |
| 17 | ArangoDB migration for graph traversal performance | 3 days | Low |
| 18 | User authentication and session management | 2 days | Low |
| 19 | Batch query API endpoint | 1 day | Low |
| 20 | Curriculum learning for MDrI improvement | 1 day | Low |

---

## 14. Test Suite

231 tests covering:
- API endpoints
- Page rendering (all routes)
- Data services (MPI/MEI/MDI/MMI/MDrI/MGI/mGWAS)
- Enrichment analysis
- Molecular validation
- Foundation model inference
- Export functionality

Run tests:
```bash
conda activate mpi-vgae
cd mpi-vgae-web
pytest tests/ -v --tb=short
```
