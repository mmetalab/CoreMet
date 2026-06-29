# Changelog

All notable changes to CoreMet are documented in this file.

---

## [1.5.0] — 2026-03-04

### Added — Database Expansions
- **MDrI expansion** — 97 → **3,500** metabolite–drug interactions via DrugBank/HMDB cross-referencing (3,283 DrugBank cross-ref + 120 enzyme bridge)
- **MDI expansion (CTD)** — 7,045 → **55,921** metabolite–disease interactions via CTD Chemical–Disease (DirectEvidence only, 15,647 CAS-mapped metabolites)
- **MMI expansion (AGORA2)** — 2,323 → **83,149** metabolite–microbe interactions via 818 AGORA2 genome-scale metabolic models (exchange reactions with HMDB/KEGG/PubChem/ChEBI annotations)
- **Total interaction count**: 37,596 MPI + 196,998 MEI + 55,921 MDI + 83,149 MMI + 3,500 MDrI = **377,164**

### Added — MDrI Integration
- **Network page**: Drug node style (teal diamond), MDrI edge style, "Drug" query type radio, "Drug (MDrI)" toggle checkbox, drug legend, "Metformin" example
- **Profile page**: MDrI search in `_search_all_databases()`, drug type badges, drug pie colour (#319795), MDrI section card with Drug_Name/DrugBank_ID/Interaction_Type columns
- **Drug Enrichment tab**: Fisher's exact test + BH FDR against MDrI database, PK/PD bar chart, summary stats, results table, CSV download
- **MDrI REST API**: `/api/v1/mdri/stats` and `/api/v1/mdri/search` endpoints with query parameters (q, metabolite, drug, interaction_type, evidence, limit)

### Added — Reports
- **Unified PDF report**: `EnrichmentPDF` class in `pdf_report_service.py` — CoreMet-branded enrichment report combining pathway, disease, microbe, and drug analyses; "Download Unified PDF Report" button on enrichment page

### Added — Expansion Scripts
- `scripts/expand_mdri.py` — MDrI expansion from DrugBank XML + HMDB XML
- `scripts/expand_mdi_ctd.py` — MDI expansion from CTD Chemical–Disease TSV
- `scripts/expand_mmi_agora2.py` — MMI expansion from 818 AGORA2 SBML models

### Changed
- `requirements.txt` — Added `fpdf2 == 2.8.7`, `lxml == 6.0.2`

### Removed
- `mpivgae.py` — Deprecated legacy VGAE script
- `pages/interactions.py` — Deprecated interactions page
- `pages/docking.py` — Deprecated docking page

---

## [1.4.0] — 2026-03-03

### Added
- **Gallery page** — New `/gallery` route showcasing three curated case studies:
  - *Hepatocellular Carcinoma (HCC)* — VGAE prediction & network analysis (1,635 predictions, density 0.41)
  - *Papillary Thyroid Cancer* — Pathway enrichment focus (922 predictions, FDR-significant pathways)
  - *Alzheimer's Disease* — Neurodegenerative MPI network (119 edges, density 0.44)
  - Cross-disease comparison card highlighting unique metabolic fingerprints
  - Interactive hub-degree bar charts and enrichment visualizations loaded from live data
- **Gallery** added to navbar (between Profile and Docs)

### Changed
- **Disease panels expanded to 130 diseases** across 22 categories and 73 tissues:
  - Batch 1 (+30): mesothelioma, head/neck cancer, testicular cancer, sarcoma, neuroblastoma, multiple myeloma, nasopharyngeal cancer, adrenal cancer, type 1 diabetes, gout, hypothyroidism, hyperthyroidism, PCOS, stroke, cardiomyopathy, aortic aneurysm, epilepsy, migraine, Huntington's, autism, psoriasis, celiac disease, ankylosing spondylitis, liver cirrhosis, pancreatitis, gastric ulcer, nephrotic syndrome, IgA nephropathy, tuberculosis, hepatitis B
  - Batch 2 (+30): hepatitis C, HIV/AIDS, sepsis, malaria, bipolar disorder, ADHD, anxiety, osteoporosis, osteoarthritis, phenylketonuria, preeclampsia, endometriosis, cystic fibrosis, sickle cell, hemophilia, Wilson's, hemochromatosis, Cushing's, myasthenia gravis, Sjögren's, scleroderma, Crohn's, pulmonary hypertension, atrial fibrillation, peripheral artery disease, pulmonary fibrosis, diabetic nephropathy, Addison's, acromegaly, ulcerative colitis
  - Batch 3 (+30): retinoblastoma, gallbladder cancer, laryngeal cancer, Wilms tumor, CLL, thymoma, glaucoma, macular degeneration, diabetic retinopathy, vitiligo, dermatomyositis, vasculitis, alopecia areata, dengue, meningitis, gallstones, diverticulitis, gastroparesis, kidney stones, neuropathy, Lewy body dementia, myotonic dystrophy, amyloidosis, fibromyalgia, chronic fatigue syndrome, BPH, interstitial cystitis, hyperaldosteronism, polycythemia vera, thalassemia
- **New disease categories**: Ophthalmological, Urological, Neuromuscular, Systemic, Gynecological, Obstetric, Genetic, Hepatic, Hematological (expanded from original 10 → 22)
- **About page** — Updated platform statistics: 130 diseases, 22 categories, 73 tissues
- **Home page** — `disease_panels` stat updated from 40 → 130; added `disease_categories` = 22
- **Documentation page** — Disease Explorer tutorial updated to reference 130 panels
- **Help page** — MDI description updated to 7,045 associations across 130 disease panels
- **README.md** — Disease Networks feature line updated to 130 diseases across 22 categories
- **disease_service.py** — Comment updated to clarify legacy 40-disease catalogue vs full 130-entry registry

---

## [1.3.0] — 2026-03-03

### Added
- **MDrI (Metabolite–Drug Interaction) module** — 97 curated drug–metabolite interactions
- **Documentation page** (`/documentation`) — Step-by-step tutorials for all modules
- **About page** (`/about`) — Contact, citation, data sources, licensing

### Changed
- Entity chart replaced Barpolar with horizontal bar (log scale, 7 entity types)
- Network page supports multi-type interaction queries
- Disease page dynamic sidebar labels for MPI/MDI/MMI switching

---

## [1.2.0] — 2026-03-03

### Added
- **Disease expansion Batch 0** — Original 40 disease panels with VGAE predictions
- Interactive Plotly charts on home page (composition donut, category bar, treemap, word cloud)
- Highlight numbers strip, coverage cards, workflow steps, feature grid

### Changed
- Home page redesign with ~180 lines of custom CSS
- All Plotly charts standardized to Arial font

---

## [1.1.0] — 2026-03-03

### Added
- MEI (Metabolite–Enzyme Interaction) database: 196,998 interactions
- MDI (Metabolite–Disease Interaction) database: 7,045 interactions
- MMI (Metabolite–Microbe Interaction) database: 2,323 interactions
- Database page v2 with tabbed MPI/MDI/MMI/MDrI views

---

## [1.0.0] — 2026-03-02

### Added
- Initial release of CoreMet web platform
- MPI database: 37,596 metabolite–protein interactions across 10 organisms
- VGAE-based interaction prediction engine
- Pathway enrichment analysis (Fisher's exact test, BH FDR correction)
- Disease Explorer with network visualization (Cytoscape)
- Metabolite Profile page with external database links
