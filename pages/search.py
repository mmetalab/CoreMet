"""
Search page, main entry point for entity search.
Global search with mixed-entity autocomplete and search results display.

Route: /search  or  /search?q=butyrate
"""

import os
from pathlib import Path
from urllib.parse import parse_qs

from dash import html, dcc, Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

from components.page_header import make_page_header

# ── Batch example datasets ──────────────────────────────────────────────
_BATCH_EXAMPLES_DIR = Path(__file__).parent.parent / "data" / "batch_examples"

BATCH_EXAMPLES = {
    "gut-brain": {
        "label": "Gut–Brain Axis",
        "desc": "20 entities: SCFAs, neurotransmitters, gut microbes, neurological diseases, and related genes",
        "file": "gut_brain_axis.csv",
        "icon": "fas fa-brain",
    },
    "hcc": {
        "label": "HCC Biomarkers",
        "desc": "20 entities: amino acids, oncogenes, liver diseases, and targeted therapies",
        "file": "hcc_biomarkers.csv",
        "icon": "fas fa-lungs",
    },
    "t2d": {
        "label": "Type 2 Diabetes",
        "desc": "20 entities: BCAAs, glucose, antidiabetic drugs, diabetes SNPs, and key genes",
        "file": "t2d_metabolic_panel.csv",
        "icon": "fas fa-vial",
    },
    "scfa": {
        "label": "SCFA Network",
        "desc": "20 entities: short-chain fatty acids, producing microbes, GI diseases, and receptors",
        "file": "scfa_network.csv",
        "icon": "fas fa-bacterium",
    },
    "crc": {
        "label": "Colorectal Multi-omics",
        "desc": "20 entities: metabolites, microbes, tumor suppressors, and chemotherapy drugs",
        "file": "colorectal_multiomics.csv",
        "icon": "fas fa-dna",
    },
}


def _load_example(key):
    """Load an example dataset and return entity names as newline-separated string."""
    info = BATCH_EXAMPLES.get(key)
    if not info:
        return ""
    fpath = _BATCH_EXAMPLES_DIR / info["file"]
    if not fpath.exists():
        return ""
    lines = fpath.read_text().strip().split("\n")
    # Skip header if present
    if lines and lines[0].lower().startswith("entity"):
        lines = lines[1:]
    return "\n".join(l.strip() for l in lines if l.strip())

# Entity type → (color, icon)
ENTITY_STYLE = {
    "metabolite": ("#e27a3f", "fas fa-atom"),
    "protein":    ("#3182ce", "fas fa-cubes"),
    "gene":       ("#d69e2e", "fas fa-dna"),
    "disease":    ("#e53e3e", "fas fa-heartbeat"),
    "microbe":    ("#38a169", "fas fa-bacterium"),
    "drug":       ("#805ad5", "fas fa-pills"),
    "snp":        ("#319795", "fas fa-map-marker-alt"),
}


def _entity_badge(etype):
    """Small colored chip for entity type."""
    color = ENTITY_STYLE.get(etype, ("#718096", "fas fa-circle"))[0]
    return html.Span(etype.capitalize(), className="search-result-badge",
                     style={"color": color, "background": f"{color}14",
                            "border": f"1px solid {color}30"})


def build_search_page(search_qs=None):
    """Build the search page, optionally with a pre-filled query."""
    query = ""
    if search_qs:
        params = parse_qs(search_qs.lstrip("?"))
        query = params.get("q", [""])[0]

    # Pre-render results server-side if query is present (no double-refresh)
    initial_results = None
    if query:
        results = _search_all(query)
        initial_results = _render_results(query, results)

    return html.Div([
        html.Div([
            # ── Header ──
            html.Div([
                html.H1("Search CoreMet", style={
                    "fontSize": "1.8rem", "fontWeight": "700", "color": "#1a202c",
                    "marginBottom": "8px",
                }),
                html.P(
                    "Find any metabolite, disease, gene, protein, drug, microbe, or SNP.",
                    style={"fontSize": "0.95rem", "color": "#718096", "marginBottom": "24px"},
                ),
            ], style={"textAlign": "center", "paddingTop": "32px"}),

            # ── Tabs: Single Search | Batch Search ──
            dbc.Tabs([
                dbc.Tab(label="Single Search", tab_id="tab-single", children=[
                    html.Div([
                        # Search bar
                        html.Div([
                            dbc.InputGroup([
                                dbc.Input(
                                    id="search-input",
                                    placeholder="Search a metabolite, disease, gene, protein, drug, microbe, or SNP…",
                                    value=query,
                                    type="text",
                                    debounce=True,
                                    style={"fontSize": "1rem", "padding": "12px 16px",
                                           "borderRadius": "8px 0 0 8px"},
                                ),
                                dbc.Button(
                                    html.I(className="fas fa-search"),
                                    id="search-btn",
                                    color="primary",
                                    style={"borderRadius": "0 8px 8px 0", "padding": "12px 20px"},
                                ),
                            ], style={"maxWidth": "640px", "margin": "0 auto"}),
                        ], style={"marginBottom": "24px", "marginTop": "24px"}),

                        # Example chips
                        html.Div([
                            html.Span("Examples: ", style={"fontSize": "0.8rem", "color": "#a0aec0"}),
                            html.A("Butyrate", href="/search?q=butyrate", className="cm-example-chip"),
                            html.A("LDHA", href="/search?q=LDHA", className="cm-example-chip"),
                            html.A("Metformin", href="/search?q=metformin", className="cm-example-chip"),
                            html.A("Colorectal cancer", href="/search?q=colorectal+cancer",
                                   className="cm-example-chip"),
                            html.A("rs1260326", href="/search?q=rs1260326", className="cm-example-chip"),
                        ], style={
                            "display": "flex", "justifyContent": "center", "alignItems": "center",
                            "gap": "4px", "marginBottom": "32px", "flexWrap": "wrap",
                        }),

                        # Results
                        html.Div(initial_results, id="search-results-container"),
                    ]),
                ]),

                dbc.Tab(label="Batch Search (up to 200)", tab_id="tab-batch", children=[
                    html.Div([
                        html.P(
                            "Paste up to 200 entity names (one per line) to resolve, map, and visualize their connections.",
                            style={"fontSize": "0.9rem", "color": "#718096", "marginTop": "20px",
                                   "marginBottom": "12px", "textAlign": "center"},
                        ),

                        # ── Example datasets ──
                        html.Div([
                            html.Span("Load example: ", style={
                                "fontSize": "0.8rem", "color": "#718096", "fontWeight": "600",
                                "marginRight": "8px",
                            }),
                            *[
                                dbc.Button(
                                    [html.I(className=f"{info['icon']} me-1"), info["label"]],
                                    id=f"batch-example-{key}",
                                    color="outline-secondary", size="sm",
                                    title=info["desc"],
                                    style={"fontSize": "0.75rem", "borderRadius": "20px",
                                           "margin": "2px 4px", "padding": "4px 12px",
                                           "color": "#4a5568", "borderColor": "#cbd5e0",
                                           "backgroundColor": "#f7fafc"},
                                )
                                for key, info in BATCH_EXAMPLES.items()
                            ],
                        ], style={
                            "display": "flex", "flexWrap": "wrap", "justifyContent": "center",
                            "alignItems": "center", "marginBottom": "12px",
                        }),

                        # ── Template download ──
                        html.Div([
                            dbc.Button(
                                [html.I(className="fas fa-download me-1"), "Download CSV Template"],
                                id="batch-template-btn", color="link", size="sm",
                                style={"fontSize": "0.78rem", "color": "#3182ce"},
                            ),
                            dcc.Download(id="batch-template-download"),
                        ], style={"textAlign": "center", "marginBottom": "16px"}),

                        # ── Text area ──
                        dbc.Textarea(
                            id="batch-input",
                            placeholder="Butyrate\nLDHA\nMetformin\nColorectal cancer\nrs1260326\nLactobacillus\nAlbumin",
                            style={"height": "180px", "fontSize": "0.9rem", "fontFamily": "Arial, sans-serif",
                                   "borderRadius": "8px", "marginBottom": "12px"},
                        ),

                        # ── Action buttons ──
                        html.Div([
                            dbc.Button([html.I(className="fas fa-search me-2"), "Resolve & Map Entities"],
                                       id="batch-search-btn", color="primary", size="md",
                                       style={"borderRadius": "8px"}),
                            dbc.Button([html.I(className="fas fa-eraser me-1"), "Clear"],
                                       id="batch-clear-btn", color="secondary", size="md", outline=True,
                                       style={"borderRadius": "8px", "marginLeft": "8px"}),
                            html.Span(id="batch-count-label",
                                      style={"fontSize": "0.8rem", "color": "#a0aec0", "marginLeft": "12px"}),
                        ], style={"textAlign": "center", "marginBottom": "24px"}),

                        # Loading spinner for batch results
                        dcc.Loading(
                            id="batch-loading",
                            type="circle",
                            children=html.Div(id="batch-results-container"),
                        ),

                        # Hidden download trigger
                        dcc.Download(id="batch-results-download"),
                        dcc.Store(id="batch-results-csv-store"),
                    ], style={"maxWidth": "900px", "margin": "0 auto"}),
                ]),
            ], id="search-tabs", active_tab="tab-single" if not query else "tab-single",
               style={"marginBottom": "24px"}),

            # Hidden stores
            dcc.Store(id="search-initial-query", data=""),
            dcc.Interval(id="search-auto-trigger", interval=9999999, max_intervals=0, disabled=True),

        ], className="cm-page-container"),
    ])


# Static fallback for routing
page_content = build_search_page()


# --------------------------------------------------------------------------
# CALLBACKS
# --------------------------------------------------------------------------

@callback(
    Output("search-results-container", "children"),
    Input("search-btn", "n_clicks"),
    Input("search-input", "n_submit"),
    Input("search-auto-trigger", "n_intervals"),
    State("search-input", "value"),
    State("search-initial-query", "data"),
    prevent_initial_call=True,
)
def run_search(n_clicks, n_submit, n_intervals, query, initial_query):
    """Search across all entity types and display results."""
    from dash import ctx
    if ctx.triggered_id == "search-auto-trigger" and initial_query:
        query = initial_query
    if not query or len(query.strip()) < 2:
        return html.Div("Enter at least 2 characters to search.", style={
            "textAlign": "center", "color": "#a0aec0", "padding": "48px",
        })

    query = query.strip()
    results = _search_all(query)

    # Enrich results with CoreMet IDs
    try:
        from app.services.entity_registry import lookup_id
        for r in results:
            cid = lookup_id(r["name"], r["type"])
            if not cid and r.get("id"):
                cid = lookup_id(r["id"])
            r["coremetdb_id"] = cid or ""
    except Exception:
        pass

    return _render_results(query, results)


def _render_results(query, results):
    """Render search results as grouped cards (used by both server-side and callback)."""
    if not results:
        return html.Div([
            html.Div(html.I(className="fas fa-search", style={"fontSize": "2rem", "color": "#e2e8f0"})),
            html.P(f"No results found for \"{query}\".", style={
                "color": "#718096", "marginTop": "12px",
            }),
        ], style={"textAlign": "center", "padding": "48px"})

    # Group by entity type
    grouped = {}
    for r in results:
        grouped.setdefault(r["type"], []).append(r)

    cards = []
    for etype, items in grouped.items():
        color = ENTITY_STYLE.get(etype, ("#718096", "fas fa-circle"))[0]
        icon = ENTITY_STYLE.get(etype, ("#718096", "fas fa-circle"))[1]

        rows = []
        for item in items[:20]:
            coremetdb_id = item.get("coremetdb_id", "")
            ext_id = item.get("id", "")
            href = item.get("href", "#")

            # CoreMet ID cell, make it a clickable link
            id_children = []
            if coremetdb_id:
                id_children.append(html.A(
                    coremetdb_id, href=href,
                    className="result-id",
                    style={"textDecoration": "none"},
                ))
            if ext_id and ext_id != coremetdb_id:
                id_children.append(html.Span(
                    ext_id, className="result-ext-id",
                ))

            rows.append(
                html.A(
                    dbc.Row([
                        dbc.Col(
                            html.Div([
                                html.I(className=icon, style={"color": color, "marginRight": "10px"}),
                                html.Span(item["name"], className="result-name"),
                            ], style={"display": "flex", "alignItems": "center"}),
                            md=4,
                        ),
                        dbc.Col(html.Div(id_children), md=4),
                        dbc.Col(_entity_badge(etype), md=2),
                        dbc.Col(
                            html.I(className="fas fa-chevron-right", style={"color": "#e2e8f0"}),
                            md=2, className="text-end",
                        ),
                    ], className="align-items-center"),
                    href=href,
                    style={
                        "display": "block", "padding": "12px 16px", "textDecoration": "none",
                        "borderBottom": "1px solid #f0f0f0", "transition": "background 0.1s",
                    },
                    className="search-result-row",
                )
            )

        cards.append(html.Div([
            html.Div([
                html.Span(style={
                    "width": "8px", "height": "8px", "borderRadius": "50%",
                    "background": color, "display": "inline-block", "marginRight": "8px",
                }),
                html.Span(f"{etype.capitalize()}s", className="group-label"),
                html.Span(f" ({len(items)})", className="group-count"),
            ], className="search-group-header"),
            *rows,
        ], className="cm-card", style={"marginBottom": "16px", "padding": "0", "overflow": "hidden"}))

    return html.Div([
        html.Div(f"{len(results)} results for \"{query}\"", style={
            "fontSize": "0.85rem", "color": "#718096", "marginBottom": "16px",
        }),
        *cards,
    ])


def _search_all(query):
    """Search across all databases for matching entities."""
    results = []
    q = query.lower()

    # Common metabolite synonyms for improved matching
    _SYNONYMS = {
        "butyrate": "butyric acid", "acetate": "acetic acid", "pyruvate": "pyruvic acid",
        "lactate": "lactic acid", "succinate": "succinic acid", "fumarate": "fumaric acid",
        "citrate": "citric acid", "malate": "malic acid", "glutamate": "glutamic acid",
        "aspartate": "aspartic acid", "oxalate": "oxalic acid", "formate": "formic acid",
        "propionate": "propionic acid", "valerate": "valeric acid",
        "albumin": "ALB", "serum albumin": "ALB",
    }
    # Also try synonym form when searching
    alt_q = _SYNONYMS.get(q, "").lower()

    # Search metabolites (MPI)
    try:
        from app.config import Config
        import pandas as pd
        cfg = Config()
        mpi = pd.read_csv(cfg.MPI_DB_PATH, usecols=["HMDB ID", "Metabolite Name"], dtype=str)
        mpi = mpi.dropna(subset=["Metabolite Name"]).drop_duplicates(subset=["HMDB ID"])
        mask = mpi["Metabolite Name"].str.lower().str.contains(q, na=False)
        if alt_q:
            mask = mask | mpi["Metabolite Name"].str.lower().str.contains(alt_q, na=False)
        matches = mpi[mask].head(10)
        for _, row in matches.iterrows():
            results.append({
                "type": "metabolite", "name": row["Metabolite Name"],
                "id": row["HMDB ID"], "href": f"/metabolite?id={row['HMDB ID']}",
            })
    except Exception:
        pass

    from urllib.parse import quote_plus

    # Search genes (MGI), before diseases to avoid substring collisions (e.g. ALB vs Albuminuria)
    try:
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if mgi is not None and "Gene_Symbol" in mgi.columns:
            mask = mgi["Gene_Symbol"].str.lower().str.contains(q, na=False)
            if alt_q:
                mask = mask | mgi["Gene_Symbol"].str.lower().str.contains(alt_q, na=False)
            gene_matches = mgi[mask]
            for name in gene_matches["Gene_Symbol"].unique()[:10]:
                gid = ""
                row = gene_matches[gene_matches["Gene_Symbol"] == name].iloc[0]
                if "Gene_ID" in gene_matches.columns:
                    gid = str(row.get("Gene_ID", ""))
                results.append({
                    "type": "gene", "name": name, "id": gid,
                    "href": f"/gene?name={quote_plus(name)}",
                })
    except Exception:
        pass

    # Search proteins (MPI)
    try:
        from app.config import Config as Cfg2
        cfg2 = Cfg2()
        mpi2 = pd.read_csv(cfg2.MPI_DB_PATH, usecols=["Uniprot ID", "Protein Name"], dtype=str)
        mpi2 = mpi2.dropna(subset=["Protein Name"]).drop_duplicates(subset=["Protein Name"])
        prot_matches = mpi2[mpi2["Protein Name"].str.lower().str.contains(q, na=False)].head(10)
        for _, row in prot_matches.iterrows():
            results.append({
                "type": "protein", "name": row["Protein Name"],
                "id": str(row.get("Uniprot ID", "")),
                "href": f"/protein?name={quote_plus(row['Protein Name'])}",
            })
    except Exception:
        pass

    # Search diseases (MDI)
    try:
        from app.services.mdi_service import get_mdi_db
        mdi = get_mdi_db()
        if mdi is not None and "Disease_Name" in mdi.columns:
            disease_matches = mdi[mdi["Disease_Name"].str.lower().str.contains(q, na=False)]
            for name in disease_matches["Disease_Name"].unique()[:10]:
                did = ""
                row = disease_matches[disease_matches["Disease_Name"] == name].iloc[0]
                if "Disease_ID" in disease_matches.columns:
                    did = str(row.get("Disease_ID", ""))
                results.append({
                    "type": "disease", "name": name, "id": did,
                    "href": f"/disease-detail?name={quote_plus(name)}",
                })
    except Exception:
        pass

    # Search drugs (MDrI)
    try:
        from app.services.mdri_service import get_mdri_db
        mdri = get_mdri_db()
        if mdri is not None and "Drug_Name" in mdri.columns:
            drug_matches = mdri[mdri["Drug_Name"].str.lower().str.contains(q, na=False)]
            for name in drug_matches["Drug_Name"].unique()[:10]:
                did = ""
                row = drug_matches[drug_matches["Drug_Name"] == name].iloc[0]
                if "DrugBank_ID" in drug_matches.columns:
                    did = str(row.get("DrugBank_ID", ""))
                results.append({
                    "type": "drug", "name": name, "id": did,
                    "href": f"/drug?name={quote_plus(name)}",
                })
    except Exception:
        pass

    # Search microbes (MMI)
    try:
        from app.services.mmi_service import get_mmi_db
        mmi = get_mmi_db()
        if mmi is not None and "Microbe_Name" in mmi.columns:
            mic_matches = mmi[mmi["Microbe_Name"].str.lower().str.contains(q, na=False)]
            for name in mic_matches["Microbe_Name"].unique()[:10]:
                tid = ""
                row = mic_matches[mic_matches["Microbe_Name"] == name].iloc[0]
                if "Taxonomy_ID" in mic_matches.columns:
                    tid = str(row.get("Taxonomy_ID", ""))
                results.append({
                    "type": "microbe", "name": name, "id": tid,
                    "href": f"/microbe?name={quote_plus(name)}",
                })
    except Exception:
        pass

    # Search SNPs (mGWAS)
    try:
        from app.services.mgwas_service import get_mgwas_db
        mgwas = get_mgwas_db()
        if mgwas is not None and "rsID" in mgwas.columns:
            snp_matches = mgwas[mgwas["rsID"].str.lower().str.contains(q, na=False)]
            for name in snp_matches["rsID"].unique()[:10]:
                results.append({
                    "type": "snp", "name": name, "id": name,
                    "href": f"/snp?name={quote_plus(name)}",
                })
    except Exception:
        pass

    # Enrich all results with CoreMet IDs
    try:
        from app.services.entity_registry import lookup_id
        for r in results:
            cid = lookup_id(r["name"], r["type"])
            if not cid and r.get("id"):
                cid = lookup_id(r["id"], r["type"])
            r["coremetdb_id"] = cid or ""
    except Exception:
        pass

    return results


# --------------------------------------------------------------------------
# BATCH SEARCH
# --------------------------------------------------------------------------

# Entity type colors consistent with site theme
_TYPE_COLORS = {
    "metabolite": "#e27a3f", "protein": "#3182ce", "gene": "#d69e2e",
    "disease": "#e53e3e", "microbe": "#38a169", "drug": "#805ad5", "snp": "#319795",
}


def _resolve_single(name):
    """Resolve a single entity name to type and ID. Returns best match."""
    name = name.strip()
    if not name or len(name) < 2:
        return None

    results = _search_all(name)
    if not results:
        return None

    nl = name.lower()
    # 1. Exact match (highest priority)
    exact = [r for r in results if r["name"].lower() == nl]
    if exact:
        return exact[0]
    # 2. Metabolite with exact or startswith match
    met = [r for r in results if r["type"] == "metabolite"]
    met_starts = [r for r in met if r["name"].lower().startswith(nl)]
    if met_starts:
        return min(met_starts, key=lambda r: len(r["name"]))
    # 3. Shortest name among all results (closest match)
    return min(results, key=lambda r: len(r["name"]))


def _find_edges_between(resolved_entities):
    """Find all database edges connecting the resolved entities."""
    from urllib.parse import quote_plus
    edges = []
    names_lower = {r["name"].lower() for r in resolved_entities}
    ids_lower = {r.get("id", "").lower() for r in resolved_entities if r.get("id")}

    def _matches(val):
        if not val:
            return False
        return str(val).lower() in names_lower or str(val).lower() in ids_lower

    # MPI/MEI: Metabolite ↔ Protein/Gene
    try:
        from app.config import Config
        import pandas as pd
        cfg = Config()
        mpi = pd.read_csv(cfg.MPI_DB_PATH, dtype=str).fillna("")
        for _, row in mpi.iterrows():
            src_match = _matches(row.get("Metabolite Name")) or _matches(row.get("HMDB ID"))
            tgt_match = _matches(row.get("Protein Name")) or _matches(row.get("Gene Name"))
            if src_match and tgt_match:
                edges.append({
                    "source": row.get("Metabolite Name", ""),
                    "target": row.get("Protein Name") or row.get("Gene Name", ""),
                    "layer": "MPI", "source_type": "metabolite", "target_type": "protein",
                })
    except Exception:
        pass

    # MDI: Metabolite ↔ Disease
    try:
        from app.services.mdi_service import get_mdi_db
        mdi = get_mdi_db()
        if mdi is not None:
            for _, row in mdi.iterrows():
                src_match = _matches(row.get("Metabolite_Name")) or _matches(row.get("HMDB_ID"))
                tgt_match = _matches(row.get("Disease_Name"))
                if src_match and tgt_match:
                    edges.append({
                        "source": str(row.get("Metabolite_Name", "")),
                        "target": str(row.get("Disease_Name", "")),
                        "layer": "MDI", "source_type": "metabolite", "target_type": "disease",
                    })
                    if len(edges) > 500:
                        break
    except Exception:
        pass

    # MMI: Metabolite ↔ Microbe
    try:
        from app.services.mmi_service import get_mmi_db
        mmi = get_mmi_db()
        if mmi is not None:
            for _, row in mmi.iterrows():
                src_match = _matches(row.get("Metabolite_Name")) or _matches(row.get("HMDB_ID"))
                tgt_match = _matches(row.get("Microbe_Name"))
                if src_match and tgt_match:
                    edges.append({
                        "source": str(row.get("Metabolite_Name", "")),
                        "target": str(row.get("Microbe_Name", "")),
                        "layer": "MMI", "source_type": "metabolite", "target_type": "microbe",
                    })
                    if len(edges) > 500:
                        break
    except Exception:
        pass

    # MDrI: Metabolite ↔ Drug
    try:
        from app.services.mdri_service import get_mdri_db
        mdri = get_mdri_db()
        if mdri is not None:
            for _, row in mdri.iterrows():
                src_match = _matches(row.get("Metabolite_Name")) or _matches(row.get("HMDB_ID"))
                tgt_match = _matches(row.get("Drug_Name"))
                if src_match and tgt_match:
                    edges.append({
                        "source": str(row.get("Metabolite_Name", "")),
                        "target": str(row.get("Drug_Name", "")),
                        "layer": "MDrI", "source_type": "metabolite", "target_type": "drug",
                    })
                    if len(edges) > 500:
                        break
    except Exception:
        pass

    # MGI: Metabolite ↔ Gene
    try:
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if mgi is not None:
            for _, row in mgi.head(200000).iterrows():
                src_match = _matches(row.get("Metabolite_Name")) or _matches(row.get("HMDB_ID"))
                tgt_match = _matches(row.get("Gene_Symbol"))
                if src_match and tgt_match:
                    edges.append({
                        "source": str(row.get("Metabolite_Name", "")),
                        "target": str(row.get("Gene_Symbol", "")),
                        "layer": "MGI", "source_type": "metabolite", "target_type": "gene",
                    })
                    if len(edges) > 500:
                        break
    except Exception:
        pass

    # mGWAS: Metabolite ↔ SNP
    try:
        from app.services.mgwas_service import get_mgwas_db
        mgwas = get_mgwas_db()
        if mgwas is not None:
            for _, row in mgwas.iterrows():
                src_match = _matches(row.get("Metabolite_Name")) or _matches(row.get("HMDB_ID"))
                tgt_match = _matches(row.get("rsID"))
                if src_match and tgt_match:
                    edges.append({
                        "source": str(row.get("Metabolite_Name", "")),
                        "target": str(row.get("rsID", "")),
                        "layer": "mGWAS", "source_type": "metabolite", "target_type": "snp",
                    })
                    if len(edges) > 500:
                        break
    except Exception:
        pass

    # Deduplicate edges
    seen = set()
    unique = []
    for e in edges:
        key = (e["source"].lower(), e["target"].lower(), e["layer"])
        if key not in seen:
            seen.add(key)
            unique.append(e)
    return unique[:500]


def _build_batch_network(resolved, edges):
    """Build Cytoscape elements for batch search results."""
    nodes = {}
    for r in resolved:
        nid = r["name"]
        if nid not in nodes:
            nodes[nid] = {"type": r["type"], "href": r.get("href", "#")}

    # Nodes from edges (may include entities not in the input)
    for e in edges:
        for key in ["source", "target"]:
            nid = e[key]
            if nid not in nodes:
                nodes[nid] = {"type": e[f"{key}_type"], "href": "#"}

    elements = []
    for name, info in nodes.items():
        elements.append({
            "data": {
                "id": name, "label": name[:25] + ("…" if len(name) > 25 else ""),
                "type": info["type"], "href": info["href"],
            },
            "classes": info["type"],
        })

    for i, e in enumerate(edges):
        elements.append({
            "data": {
                "id": f"e{i}", "source": e["source"], "target": e["target"],
                "layer": e["layer"],
            },
            "classes": e["layer"],
        })

    return elements


# ── Example dataset loaders ──────────────────────────────────────────────

@callback(
    Output("batch-input", "value", allow_duplicate=True),
    [Input(f"batch-example-{key}", "n_clicks") for key in BATCH_EXAMPLES],
    prevent_initial_call=True,
)
def load_batch_example(*args):
    """Load example dataset into batch textarea."""
    from dash import ctx
    if not ctx.triggered_id:
        return no_update
    key = ctx.triggered_id.replace("batch-example-", "")
    return _load_example(key)


@callback(
    Output("batch-template-download", "data"),
    Input("batch-template-btn", "n_clicks"),
    prevent_initial_call=True,
)
def download_template(n_clicks):
    """Download batch search CSV template."""
    template = "entity_name\nButyric acid\nLDHA\nMetformin\nColorectal cancer\nrs1260326\nLactobacillus\nAlbumin\n"
    return dict(content=template, filename="coremetdb_batch_template.csv")


@callback(
    Output("batch-input", "value", allow_duplicate=True),
    Input("batch-clear-btn", "n_clicks"),
    prevent_initial_call=True,
)
def clear_batch(n_clicks):
    """Clear the batch text area."""
    return ""


@callback(
    Output("batch-results-download", "data"),
    Input("batch-download-csv-btn", "n_clicks"),
    State("batch-results-csv-store", "data"),
    prevent_initial_call=True,
)
def download_batch_csv(n_clicks, csv_data):
    """Download batch search results as CSV, only on explicit click."""
    if not n_clicks or not csv_data:
        return no_update
    return dict(content=csv_data, filename="coremetdb_batch_results.csv")


@callback(
    Output("batch-results-container", "children"),
    Output("batch-results-csv-store", "data"),
    Input("batch-search-btn", "n_clicks"),
    State("batch-input", "value"),
    prevent_initial_call=True,
)
def run_batch_search(n_clicks, text):
    """Resolve batch entities, find connections, build network + table."""
    if not text or not text.strip():
        return html.Div("Paste entity names (one per line) above.", style={
            "textAlign": "center", "color": "#a0aec0", "padding": "48px"}), no_update

    # Also handle comma-separated input
    raw = text.strip()
    if "\n" not in raw and "," in raw:
        lines = [l.strip() for l in raw.split(",") if l.strip()]
    else:
        lines = [l.strip() for l in raw.split("\n") if l.strip()]
    if len(lines) > 200:
        lines = lines[:200]

    # Resolve all entities
    resolved = []
    unresolved = []

    # Lazy import CoreMet ID registry
    try:
        from app.services.entity_registry import lookup_id
        _has_registry = True
    except Exception:
        _has_registry = False

    for name in lines:
        result = _resolve_single(name)
        if result:
            result["input"] = name  # track original input
            # Add CoreMet ID
            if _has_registry:
                cid = lookup_id(result["name"], result["type"])
                if not cid and result.get("id"):
                    cid = lookup_id(result["id"])
                result["coremetdb_id"] = cid or ""
            resolved.append(result)
        else:
            unresolved.append(name)

    if not resolved:
        return html.Div([
            html.I(className="fas fa-exclamation-triangle",
                   style={"fontSize": "1.5rem", "color": "#ed8936"}),
            html.P("Could not resolve any of the input entities.", style={
                "color": "#718096", "marginTop": "8px"}),
        ], style={"textAlign": "center", "padding": "48px"}), no_update

    # Find edges between resolved entities
    edges = _find_edges_between(resolved)

    # Summary stats
    type_counts = {}
    for r in resolved:
        type_counts[r["type"]] = type_counts.get(r["type"], 0) + 1

    stat_chips = []
    for t, c in sorted(type_counts.items(), key=lambda x: -x[1]):
        color = _TYPE_COLORS.get(t, "#718096")
        stat_chips.append(html.Span(f"{c} {t}s", style={
            "display": "inline-block", "padding": "4px 12px", "borderRadius": "20px",
            "fontSize": "0.75rem", "fontWeight": "600", "color": color,
            "background": f"{color}14", "border": f"1px solid {color}30", "margin": "2px",
        }))

    # Build result sections
    sections = []

    # 1. Summary bar
    sections.append(html.Div([
        html.Div([
            html.Span(f"{len(resolved)} resolved", style={
                "fontWeight": "700", "color": "#2d3748", "marginRight": "8px"}),
            html.Span(f"/ {len(lines)} input", style={"color": "#a0aec0", "marginRight": "16px"}),
            html.Span(f"{len(edges)} connections found", style={
                "fontWeight": "600", "color": "#3182ce"}),
        ], style={"marginBottom": "8px"}),
        html.Div(stat_chips),
        html.Div(f"{len(unresolved)} unresolved: {', '.join(unresolved[:10])}" +
                 ("…" if len(unresolved) > 10 else ""),
                 style={"fontSize": "0.75rem", "color": "#a0aec0", "marginTop": "4px"})
        if unresolved else None,
    ], className="cm-card", style={"padding": "16px", "marginBottom": "16px"}))

    # 2. Network visualization
    if edges:
        elements = _build_batch_network(resolved, edges)
        # Cytoscape stylesheet
        cyto_style = [
            {"selector": "node", "style": {
                "label": "data(label)", "text-valign": "center",
                "font-size": "9px", "font-family": "Arial, sans-serif",
                "width": 28, "height": 28, "text-outline-width": 2,
                "text-outline-color": "#fff", "color": "#2d3748",
            }},
            {"selector": "edge", "style": {
                "width": 1.5, "curve-style": "bezier",
                "target-arrow-shape": "triangle", "target-arrow-color": "#cbd5e0",
                "line-color": "#cbd5e0", "opacity": 0.7,
            }},
        ]
        for etype, color in _TYPE_COLORS.items():
            cyto_style.append({"selector": f".{etype}", "style": {
                "background-color": color, "border-color": color, "border-width": 2,
            }})
        # Layer-based edge colors
        layer_colors = {"MPI": "#3182ce", "MEI": "#3182ce", "MDI": "#e53e3e",
                        "MMI": "#38a169", "MDrI": "#805ad5", "MGI": "#d69e2e",
                        "mGWAS": "#319795"}
        for layer, color in layer_colors.items():
            cyto_style.append({"selector": f".{layer}", "style": {
                "line-color": color, "target-arrow-color": color,
            }})

        # Legend
        legend_items = []
        for etype, color in _TYPE_COLORS.items():
            if any(r["type"] == etype for r in resolved):
                legend_items.append(html.Span([
                    html.Span(style={"width": "10px", "height": "10px", "borderRadius": "50%",
                                     "background": color, "display": "inline-block",
                                     "marginRight": "4px"}),
                    etype.capitalize(),
                ], style={"marginRight": "12px", "fontSize": "0.75rem", "color": "#4a5568"}))

        sections.append(html.Div([
            html.H6([html.I(className="fas fa-project-diagram me-2"), "Interaction Network"],
                     style={"fontWeight": "600", "color": "#1a202c", "marginBottom": "8px"}),
            html.Div(legend_items, style={"marginBottom": "8px"}),
            cyto.Cytoscape(
                id="batch-cytoscape",
                elements=elements,
                layout={"name": "cose", "animate": True, "animationDuration": 500,
                         "nodeRepulsion": 8000, "idealEdgeLength": 80,
                         "gravity": 0.3, "padding": 30},
                stylesheet=cyto_style,
                style={"width": "100%", "height": "480px", "background": "#fafbfc",
                       "border": "1px solid #e2e8f0", "borderRadius": "8px"},
                responsive=False,
            ),
        ], className="cm-card", style={"padding": "16px", "marginBottom": "16px"}))
    else:
        sections.append(html.Div([
            html.P("No direct connections found between the input entities in CoreMet.",
                   style={"color": "#a0aec0", "textAlign": "center", "padding": "32px",
                          "fontSize": "0.9rem"}),
        ], className="cm-card", style={"marginBottom": "16px"}))

    # 3. Entity table
    table_rows = []
    for r in resolved:
        color = _TYPE_COLORS.get(r["type"], "#718096")
        # Count connections for this entity
        conn = sum(1 for e in edges if e["source"].lower() == r["name"].lower()
                   or e["target"].lower() == r["name"].lower())
        cid = r.get("coremetdb_id", "")
        table_rows.append(html.Tr([
            html.Td(html.A(r["name"], href=r.get("href", "#"), style={
                "fontWeight": "600", "color": "#3182ce", "textDecoration": "underline"})),
            html.Td(html.Span(cid, style={"fontSize": "0.78rem", "color": "#3182ce",
                                           "fontWeight": "600"}) if cid else ""),
            html.Td(r.get("id", ""), style={"fontSize": "0.8rem", "color": "#718096"}),
            html.Td(html.Span(r["type"].capitalize(), style={
                "padding": "2px 10px", "borderRadius": "20px", "fontSize": "0.7rem",
                "fontWeight": "600", "color": color, "background": f"{color}14",
                "border": f"1px solid {color}30",
            })),
            html.Td(str(conn), style={"fontWeight": "600", "color": "#3182ce",
                                       "textAlign": "center"}),
        ]))

    sections.append(html.Div([
        html.H6([html.I(className="fas fa-table me-2"), f"Resolved Entities ({len(resolved)})"],
                 style={"fontWeight": "600", "color": "#1a202c", "marginBottom": "8px"}),
        html.Div([
            html.Table([
                html.Thead(html.Tr([
                    html.Th("Entity", style={"padding": "8px 12px", "fontWeight": "600",
                                             "fontSize": "0.8rem", "color": "#4a5568",
                                             "borderBottom": "2px solid #e2e8f0"}),
                    html.Th("CoreMet ID", style={"padding": "8px 12px", "fontWeight": "600",
                                                  "fontSize": "0.8rem", "color": "#4a5568",
                                                  "borderBottom": "2px solid #e2e8f0"}),
                    html.Th("External ID", style={"padding": "8px 12px", "fontWeight": "600",
                                          "fontSize": "0.8rem", "color": "#4a5568",
                                          "borderBottom": "2px solid #e2e8f0"}),
                    html.Th("Type", style={"padding": "8px 12px", "fontWeight": "600",
                                            "fontSize": "0.8rem", "color": "#4a5568",
                                            "borderBottom": "2px solid #e2e8f0"}),
                    html.Th("Connections", style={"padding": "8px 12px", "fontWeight": "600",
                                                   "fontSize": "0.8rem", "color": "#4a5568",
                                                   "borderBottom": "2px solid #e2e8f0",
                                                   "textAlign": "center"}),
                ])),
                html.Tbody(table_rows),
            ], style={"width": "100%", "borderCollapse": "collapse"}),
        ], style={"overflowX": "auto"}),
    ], className="cm-card", style={"padding": "16px", "marginBottom": "16px"}))

    # 4. Download button
    sections.append(html.Div([
        dbc.Button(
            [html.I(className="fas fa-file-csv me-2"), "Download Results as CSV"],
            id="batch-download-csv-btn", color="success", size="sm", outline=True,
            style={"borderRadius": "8px"},
        ),
    ], style={"textAlign": "center", "marginBottom": "24px"}))

    # Build CSV string for download
    csv_lines = ["input_name,resolved_name,coremetdb_id,entity_type,external_id,connections"]
    for r in resolved:
        conn = sum(1 for e in edges if e["source"].lower() == r["name"].lower()
                   or e["target"].lower() == r["name"].lower())
        csv_lines.append(f"{r.get('input', r['name'])},{r['name']},{r.get('coremetdb_id', '')},"
                         f"{r['type']},{r.get('id', '')},{conn}")
    if unresolved:
        for u in unresolved:
            csv_lines.append(f"{u},,unresolved,,0")
    csv_data = "\n".join(csv_lines)

    return html.Div(sections), csv_data
