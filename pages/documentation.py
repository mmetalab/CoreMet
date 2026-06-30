"""
Documentation page, module-by-module tutorial for the CoreMet platform.
Covers: Database, Predict, Disease, Enrichment, Network, Profile.
"""

from dash import html, dcc
import dash_bootstrap_components as dbc

from components.page_header import make_page_header


# ---------------------------------------------------------------------------
# Helper: collapsible tutorial section
# ---------------------------------------------------------------------------

def _tutorial_section(section_id, icon, title, steps, tips=None):
    """Build a numbered tutorial section with optional tips."""
    step_items = []
    for i, (step_title, step_desc) in enumerate(steps, 1):
        step_items.append(
            html.Div([
                dbc.Row([
                    dbc.Col(
                        html.Div(
                            str(i),
                            style={
                                "width": "32px", "height": "32px", "borderRadius": "50%",
                                "background": "linear-gradient(135deg, #1a365d, #38b2ac)",
                                "color": "#fff", "display": "flex", "alignItems": "center",
                                "justifyContent": "center", "fontWeight": "700",
                                "fontSize": "0.9rem", "flexShrink": "0",
                            },
                        ),
                        width="auto",
                    ),
                    dbc.Col([
                        html.H6(step_title, style={"fontWeight": "600", "marginBottom": "2px"}),
                        html.P(step_desc, className="text-muted mb-0", style={"fontSize": "0.85rem"}),
                    ]),
                ], className="align-items-start"),
            ], style={"marginBottom": "12px"})
        )

    children = [
        html.Div([
            html.Div(html.I(className=f"fas {icon}"), className="coverage-icon icon-mpi",
                     style={"marginRight": "12px"}),
            html.H4(title, style={"fontWeight": "600", "marginBottom": "0"}),
        ], style={"display": "flex", "alignItems": "center", "marginBottom": "16px"}),
        *step_items,
    ]

    if tips:
        children.append(html.Div([
            html.Div([
                html.I(className="fas fa-lightbulb me-2", style={"color": "#d69e2e"}),
                html.Strong("Tips"),
            ], style={"marginBottom": "6px"}),
            html.Ul([html.Li(t, style={"fontSize": "0.85rem", "color": "#718096"}) for t in tips]),
        ], style={
            "background": "rgba(214,158,46,0.08)", "borderRadius": "8px",
            "padding": "12px 16px", "marginTop": "8px",
        }))

    return html.Div(children, className="cm-card", style={"marginBottom": "20px"})


# ---------------------------------------------------------------------------
# LAYOUT
# ---------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        make_page_header(
            title="Documentation",
            subtitle="Step-by-step tutorials for every module of CoreMet.",
            breadcrumb_items=[("Home", "/home"), ("Documentation", None)],
        ),

        # ── Table of Contents ──────────────────────────────────
        html.Div([
            html.H5([html.I(className="fas fa-list me-2"), "Table of Contents"],
                     style={"fontWeight": "600", "marginBottom": "12px"}),
            html.Div([
                html.A("1. Database Browser", href="#sec-database",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("2. Interaction Prediction", href="#sec-predict",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("3. Disease Explorer", href="#sec-disease",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("4. Pathway Enrichment", href="#sec-enrichment",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("5. Network Explorer", href="#sec-network",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("6. Metabolite Profile", href="#sec-profile",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("7. REST API", href="#sec-api",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("8. Batch Search", href="#sec-batch",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
                html.A("9. Data Formats & IDs", href="#sec-formats",
                       style={"display": "block", "marginBottom": "4px", "fontSize": "0.9rem"}),
            ]),
        ], className="cm-card", style={"marginBottom": "24px"}),

        # ── 1. Database ────────────────────────────────────────
        html.Div(id="sec-database"),
        _tutorial_section(
            "database", "fa-database", "1. Database Browser",
            [
                ("Select Interaction Type",
                 "Choose between Metabolite–Protein (MPI), Metabolite–Disease (MDI), "
                 "Metabolite–Microbe (MMI), or Metabolite–Drug (MDrI) using the radio buttons "
                 "at the top. The sidebar filters and table columns will adapt automatically."),
                ("Filter with the Sidebar",
                 "Use the left sidebar to narrow results. For MPI: filter by organism, pathway, "
                 "and evidence source. For MDI: filter by disease category, disease name, and "
                 "evidence level. For MMI: filter by organism, tissue, and relationship type. "
                 "For MDrI: filter by interaction type, tissue, and evidence level."),
                ("Search Globally",
                 "Type any keyword in the search bar, metabolite name, HMDB ID, UniProt ID, "
                 "gene name, disease, drug name, or pathway. The search is case-insensitive and "
                 "filters across all visible columns."),
                ("Sort & Paginate",
                 "Click any column header to sort ascending/descending. The table uses server-side "
                 "pagination, navigate pages with the pagination controls."),
                ("Download Results",
                 "Click 'Download CSV' to export the currently filtered dataset as a CSV file "
                 "for offline analysis."),
            ],
            tips=[
                "HMDB IDs, UniProt IDs, and DrugBank IDs are clickable hyperlinks to their respective databases.",
                "The mini donut chart on the sidebar shows the distribution of the current filter group.",
                "Switch interaction types frequently, the database covers 38K+ MPI, 48K+ MEI, 83K+ MDI, 78K+ MMI, 3.5K MDrI, 1.66M MGI, and 44K+ mGWAS records.",
            ],
        ),

        # ── 2. Predict ────────────────────────────────────────
        html.Div(id="sec-predict"),
        _tutorial_section(
            "predict", "fa-project-diagram", "2. Optional Prediction",
            [
                ("Prepare Your Input",
                 "Prepare a list of metabolites (HMDB IDs or names) and proteins (UniProt IDs, "
                 "gene names, or protein names). You can type them directly or upload a CSV/TXT file."),
                ("Upload Metabolites",
                 "Paste your metabolite list in the left text area. One identifier per line. "
                 "The system auto-detects HMDB IDs (e.g., HMDB0000001) and resolves names to IDs."),
                ("Upload Proteins",
                 "Paste your protein list in the right text area. Supports UniProt IDs "
                 "(e.g., P04217), gene names (e.g., ALB), or full protein names."),
                ("Run Prediction",
                 "Click 'Predict Interactions' to submit an optional beta job. Prediction is for "
                 "hypothesis generation and is separate from the curated database content."),
                ("View & Export Results",
                 "Results appear in an interactive table. Sort by score, filter by metabolite/protein, "
                 "and download as CSV. Predictions above 0.5 are considered high-confidence."),
            ],
            tips=[
                "Use the example data buttons to quickly test with pre-loaded metabolite and protein lists.",
                "Higher prediction scores indicate stronger model support, not curated evidence.",
                "For database-only workflows, use Browse, Search, Downloads, or the REST API without running prediction.",
            ],
        ),

        # ── 3. Disease ────────────────────────────────────────
        html.Div(id="sec-disease"),
        _tutorial_section(
            "disease", "fa-heartbeat", "3. Disease Explorer",
            [
                ("Select a Disease",
                 "Use the left sidebar to pick from 112 release-backed disease network panels. "
                 "Filter by category (Cancer, Neurodegenerative, Metabolic, Ophthalmological, etc.) or tissue type."),
                ("View Category & Tissue Filters",
                 "Use the 'Category' and 'Tissue' dropdowns above the disease list to narrow "
                 "the sidebar. For example, select 'Cancer' to only see oncology panels."),
                ("Explore the Network",
                 "The main panel shows an interactive Cytoscape network of metabolite–protein "
                 "interactions for the selected disease. Nodes are colour-coded by type."),
                ("Check Pathway Enrichment",
                 "Below the network, a bar chart shows pathway enrichment results (Fisher's exact test). "
                 "Bars are coloured by FDR significance, darker = more significant."),
                ("Review Disease Enrichment",
                 "The disease enrichment tab shows which other diseases share metabolite signatures "
                 "with the selected panel."),
            ],
            tips=[
                "Hover over network nodes to see metabolite/protein details.",
                "The enrichment charts support zooming and panning, drag to zoom, double-click to reset.",
                "Each disease panel contains curated metabolite biomarkers from published literature.",
            ],
        ),

        # ── 4. Enrichment ─────────────────────────────────────
        html.Div(id="sec-enrichment"),
        _tutorial_section(
            "enrichment", "fa-chart-bar", "4. Pathway Enrichment",
            [
                ("Input Your Gene/Protein Set",
                 "Enter a list of proteins or genes from your experiment. The enrichment "
                 "tool maps them to KEGG pathways and computes over-representation."),
                ("Select Parameters",
                 "Choose the organism, significance threshold (FDR cutoff), and minimum pathway "
                 "size. Default parameters work well for most analyses."),
                ("Run Enrichment",
                 "Click 'Run Enrichment' to execute Fisher's exact test against curated pathway "
                 "databases (KEGG). Results include p-values, FDR-adjusted q-values, and odds ratios."),
                ("Interpret Results",
                 "The bar chart displays the top enriched pathways ranked by -log10(FDR). "
                 "A DataTable below provides full details including pathway IDs and gene lists."),
            ],
            tips=[
                "Pathways with FDR ≤ 0.05 are generally considered statistically significant.",
                "Download the enrichment table as CSV for further analysis in R or Python.",
                "The enrichment results can be cross-referenced with the Disease Explorer findings.",
            ],
        ),

        # ── 5. Network ────────────────────────────────────────
        html.Div(id="sec-network"),
        _tutorial_section(
            "network", "fa-atom", "5. Network Explorer",
            [
                ("Search for an Entity",
                 "Type any metabolite, protein, disease, microbe, or drug name in the search bar. "
                 "The autocomplete dropdown suggests matches from all entity types."),
                ("View Ego-Graph",
                 "The query node appears at the center with all its direct interaction partners. "
                 "Edges are colour-coded: blue for MPI, red for MDI, purple for MMI, teal for MDrI."),
                ("Explore Interactively",
                 "Click any neighbor node to re-center the network on that entity. "
                 "Use the mouse wheel to zoom. Drag nodes to rearrange the layout."),
                ("Filter by Edge Type",
                 "Toggle MPI, MDI, MMI, and MDrI edges on/off using the checkboxes below the graph "
                 "to focus on specific interaction types."),
                ("Export",
                 "Download the current network as a PNG image or the edge list as CSV."),
            ],
            tips=[
                "Works for all entity types, try searching for a disease name to see its metabolite network.",
                "The layout algorithm (CoSE) automatically arranges nodes to minimise edge crossing.",
                "Multi-type queries let you see how metabolites bridge proteins, diseases, and microbes.",
            ],
        ),

        # ── 6. Metabolite Profile ─────────────────────────────
        html.Div(id="sec-profile"),
        _tutorial_section(
            "profile", "fa-flask", "6. Metabolite Profile",
            [
                ("Select a Metabolite",
                 "Search by HMDB ID or metabolite name. The profile page shows comprehensive "
                 "information about the selected metabolite."),
                ("View Interaction Summary",
                 "See all proteins, diseases, microbes, and drugs associated with this metabolite. "
                 "Counts and top partners are displayed in summary cards."),
                ("Explore Cross-References",
                 "External links to HMDB, KEGG, PubChem, and ChEBI are provided where available."),
                ("Check SMILES & Structure",
                 "If available, the SMILES string and a 2D structure rendering are shown."),
            ],
            tips=[
                "The profile page aggregates data from all four databases (MPI, MDI, MMI, MDrI).",
                "Use this as a starting point to investigate a metabolite's biological role.",
            ],
        ),

        # ── 7. REST API ───────────────────────────────────────
        html.Div(id="sec-api"),
        html.Div([
            html.Div([
                html.Div(html.I(className="fas fa-code"), className="coverage-icon icon-org",
                         style={"marginRight": "12px"}),
                html.H4("7. REST API", style={"fontWeight": "600", "marginBottom": "0"}),
            ], style={"display": "flex", "alignItems": "center", "marginBottom": "16px"}),

            html.P("CoreMet exposes a RESTful API for programmatic access; prediction is optional and beta.",
                   style={"fontSize": "0.9rem", "color": "#718096", "marginBottom": "12px"}),

            html.Div([
                html.Table([
                    html.Thead(html.Tr([
                        html.Th("Endpoint", style={"padding": "8px 12px", "fontWeight": "600"}),
                        html.Th("Method", style={"padding": "8px 12px", "fontWeight": "600"}),
                        html.Th("Description", style={"padding": "8px 12px", "fontWeight": "600"}),
                    ]), style={"backgroundColor": "#1a365d", "color": "#f7fafc"}),
                    html.Tbody([
                        html.Tr([html.Td("/api/v1/predict", style={"padding": "8px 12px"}),
                                 html.Td("POST", style={"padding": "8px 12px"}),
                                 html.Td("Optional beta metabolite–protein prediction job",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/v1/database/search", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Search the MPI database by keyword or structured filters",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/v1/species", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("List all available species",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/v1/results/<id>", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Retrieve prediction results by job ID",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/v1/health", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Health check endpoint",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/mmi/stats", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Get MMI database statistics",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/mmi/search", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Search MMI interactions by metabolite or microbe",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/mdri/stats", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Get MDrI database statistics",
                                         style={"padding": "8px 12px"})]),
                        html.Tr([html.Td("/api/mdri/search", style={"padding": "8px 12px"}),
                                 html.Td("GET", style={"padding": "8px 12px"}),
                                 html.Td("Search MDrI interactions by metabolite or drug",
                                         style={"padding": "8px 12px"})]),
                    ]),
                ], style={
                    "width": "100%", "borderCollapse": "collapse",
                    "border": "1px solid var(--cm-border)", "borderRadius": "8px",
                    "fontSize": "0.85rem",
                }),
            ]),

            html.Div([
                html.Div([
                    html.I(className="fas fa-terminal me-2", style={"color": "#38a169"}),
                    html.Strong("Example: cURL prediction request"),
                ], style={"marginBottom": "8px"}),
                html.Pre(
                    'curl -X POST https://www.coremet.org/api/v1/predict \\\n'
                    '  -H "Content-Type: application/json" \\\n'
                    '  -d \'{"metabolites":[{"name":"D-Glucose","hmdb_id":"HMDB0000122","smiles":"C(C1C(C(C(C(O1)O)O)O)O)O"}],"proteins":[{"uniprot_id":"P35557","name":"Glucokinase","gene":"GCK","organism":"Homo sapiens","sequence":"M"}],"organism":"Homo sapiens"}\'',
                    style={
                        "background": "#1a202c", "color": "#e2e8f0",
                        "padding": "12px 16px", "borderRadius": "8px",
                        "fontSize": "0.82rem", "overflow": "auto",
                    },
                ),
            ], style={"marginTop": "16px"}),
        ], className="cm-card", style={"marginBottom": "20px"}),

        # ── 8. Batch Search ──────────────────────────────────────
        html.Div(id="sec-batch"),
        _tutorial_section(
            "batch", "fa-list-ol", "8. Batch Search",
            [
                ("Navigate to Batch Search",
                 "Click 'Batch' in the navbar or go to /batch-search."),
                ("Enter Your Metabolites",
                 "Paste up to 200 HMDB IDs or metabolite names in the text area, one per line. "
                 "You can mix HMDB IDs (e.g., HMDB0000122) and common names (e.g., Glucose)."),
                ("Run the Search",
                 "Click 'Search' to look up each metabolite across all five databases (MPI, MEI, MDI, MMI, MDrI). "
                 "A summary table shows the number of interactions found per database."),
                ("Explore Results",
                 "Metabolite names and HMDB IDs in the results table are clickable, "
                 "click any to open its full detail page. Not-found metabolites are highlighted in red."),
                ("Download CSV",
                 "Click 'Download CSV' to export the batch results for offline analysis."),
            ],
            tips=[
                "Use 'Load Example' to try with 7 sample metabolites before pasting your own list.",
                "Maximum 200 metabolites per batch to keep response times reasonable.",
            ],
        ),

        # ── 9. Data Formats ─────────────────────────────────────
        html.Div(id="sec-formats"),
        html.Div([
            html.Div([
                html.Div(html.I(className="fas fa-file-alt"), className="coverage-icon icon-mdi",
                         style={"marginRight": "12px"}),
                html.H4("9. Data Formats & Identifiers", style={"fontWeight": "600", "marginBottom": "0"}),
            ], style={"display": "flex", "alignItems": "center", "marginBottom": "16px"}),

            html.P("CoreMet accepts multiple identifier formats:", style={"fontSize": "0.9rem", "color": "#718096"}),

            dbc.Row([
                dbc.Col([
                    html.H6("Metabolites", style={"fontWeight": "600"}),
                    html.Ul([
                        html.Li("HMDB ID, e.g., HMDB0000001"),
                        html.Li("Common name, e.g., Glucose, L-Alanine"),
                        html.Li("KEGG ID, e.g., C00031 (for MMI)"),
                    ], style={"fontSize": "0.85rem", "color": "#718096"}),
                ], md=3),
                dbc.Col([
                    html.H6("Proteins", style={"fontWeight": "600"}),
                    html.Ul([
                        html.Li("UniProt ID, e.g., P04217"),
                        html.Li("Gene name, e.g., ALB, INS"),
                        html.Li("Protein name, e.g., Serum albumin"),
                    ], style={"fontSize": "0.85rem", "color": "#718096"}),
                ], md=3),
                dbc.Col([
                    html.H6("Diseases & Microbes", style={"fontWeight": "600"}),
                    html.Ul([
                        html.Li("Disease name, e.g., Breast Cancer"),
                        html.Li("MeSH ID, e.g., D001943"),
                        html.Li("Microbe name, e.g., Lactobacillus"),
                    ], style={"fontSize": "0.85rem", "color": "#718096"}),
                ], md=3),
                dbc.Col([
                    html.H6("Drugs", style={"fontWeight": "600"}),
                    html.Ul([
                        html.Li("DrugBank ID, e.g., DB00945"),
                        html.Li("Drug name, e.g., Aspirin, Metformin"),
                        html.Li("Interaction type, Pharmacokinetic / Pharmacodynamic"),
                    ], style={"fontSize": "0.85rem", "color": "#718096"}),
                ], md=3),
            ]),
        ], className="cm-card", style={"marginBottom": "20px"}),

    ], className="cm-page-container"),
])

page_content = layout
