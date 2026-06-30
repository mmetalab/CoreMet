"""
Network Visualization page, interactive multi-type Cytoscape ego-graph.

Users enter a query for any node type (metabolite, protein, disease, microbe,
drug, gene, or SNP); the page searches six databases (MPI, MDI, MMI, MDrI,
MGI, mGWAS) and renders a colour-coded force-directed network:
  - Blue edges/nodes   : Metabolite–Protein   (MPI)
  - Red edges/nodes    : Metabolite–Disease   (MDI)
  - Purple edges/nodes : Metabolite–Microbe   (MMI)
  - Teal edges/nodes   : Metabolite–Drug      (MDrI)
  - Gold edges/nodes   : Metabolite–Gene      (MGI)
  - Violet edges/nodes : Metabolite–SNP       (mGWAS)

Interactive controls: type toggles, organism filter, node-click detail.
"""

from dash import dcc, html, Input, Output, State, callback, no_update, ctx
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

from components.page_header import make_page_header

# ---------------------------------------------------------------------------
# Cytoscape stylesheet, multi-type colour-coded
# ---------------------------------------------------------------------------

NETWORK_STYLESHEET = [
    # ── Centre metabolite node ──
    {
        "selector": "node.metabolite",
        "style": {
            "background-color": "#ed8936",
            "label": "data(label)",
            "shape": "ellipse",
            "width": "mapData(degree, 0, 50, 52, 96)",
            "height": "mapData(degree, 0, 50, 52, 96)",
            "font-size": "9px",
            "font-weight": "700",
            "text-valign": "center",
            "text-halign": "center",
            "color": "#1a202c",
            "text-max-width": "80px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "3px",
            "border-color": "#c05621",
            "z-index": 10,
        },
    },
    # ── Protein nodes (blue) ──
    {
        "selector": "node.protein",
        "style": {
            "background-color": "#3182ce",
            "label": "data(label)",
            "shape": "round-rectangle",
            "width": "mapData(degree, 1, 10, 30, 56)",
            "height": "mapData(degree, 1, 10, 22, 42)",
            "font-size": "6.5px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#2d3748",
            "text-margin-y": "3px",
            "text-max-width": "60px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "2px",
            "border-color": "#2b6cb0",
        },
    },
    # ── Disease nodes (red) ──
    {
        "selector": "node.disease",
        "style": {
            "background-color": "#e53e3e",
            "label": "data(label)",
            "shape": "triangle",
            "width": "mapData(degree, 1, 5, 34, 56)",
            "height": "mapData(degree, 1, 5, 34, 56)",
            "font-size": "6.5px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#742a2a",
            "text-margin-y": "3px",
            "text-max-width": "65px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "2px",
            "border-color": "#c53030",
        },
    },
    # ── Microbe nodes (purple) ──
    {
        "selector": "node.microbe",
        "style": {
            "background-color": "#805ad5",
            "label": "data(label)",
            "shape": "hexagon",
            "width": "mapData(degree, 1, 10, 32, 56)",
            "height": "mapData(degree, 1, 10, 32, 56)",
            "font-size": "6.5px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#44337a",
            "text-margin-y": "3px",
            "text-max-width": "60px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "2px",
            "border-color": "#6b46c1",
        },
    },
    # ── Drug nodes (teal) ──
    {
        "selector": "node.drug",
        "style": {
            "background-color": "#319795",
            "label": "data(label)",
            "shape": "diamond",
            "width": "mapData(degree, 1, 10, 32, 56)",
            "height": "mapData(degree, 1, 10, 32, 56)",
            "font-size": "6.5px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#234e52",
            "text-margin-y": "3px",
            "text-max-width": "60px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "2px",
            "border-color": "#2c7a7b",
        },
    },
    # ── Gene nodes (gold) ──
    {
        "selector": "node.gene",
        "style": {
            "background-color": "#d69e2e",
            "label": "data(label)",
            "shape": "round-pentagon",
            "width": "mapData(degree, 1, 10, 32, 56)",
            "height": "mapData(degree, 1, 10, 32, 56)",
            "font-size": "6.5px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#744210",
            "text-margin-y": "3px",
            "text-max-width": "60px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "2px",
            "border-color": "#b7791f",
        },
    },
    # ── SNP nodes (violet) ──
    {
        "selector": "node.snp",
        "style": {
            "background-color": "#9f7aea",
            "label": "data(label)",
            "shape": "star",
            "width": "mapData(degree, 1, 10, 28, 48)",
            "height": "mapData(degree, 1, 10, 28, 48)",
            "font-size": "6.5px",
            "text-valign": "bottom",
            "text-halign": "center",
            "color": "#553c9a",
            "text-margin-y": "3px",
            "text-max-width": "60px",
            "text-wrap": "ellipsis",
            "text-overflow-wrap": "anywhere",
            "text-outline-color": "#fff",
            "text-outline-width": "1px",
            "border-width": "2px",
            "border-color": "#805ad5",
        },
    },
    # ── Edges per type ──
    {
        "selector": "edge.mpi-edge",
        "style": {
            "line-color": "#3182ce",
            "width": 2,
            "opacity": 0.65,
            "curve-style": "bezier",
            "target-arrow-shape": "none",
        },
    },
    {
        "selector": "edge.mdi-edge",
        "style": {
            "line-color": "#e53e3e",
            "width": 2.5,
            "opacity": 0.7,
            "curve-style": "bezier",
            "target-arrow-shape": "none",
        },
    },
    {
        "selector": "edge.mmi-edge",
        "style": {
            "line-color": "#805ad5",
            "width": 2,
            "opacity": 0.65,
            "curve-style": "bezier",
            "target-arrow-shape": "none",
        },
    },
    {
        "selector": "edge.mdri-edge",
        "style": {
            "line-color": "#319795",
            "width": 2,
            "opacity": 0.65,
            "curve-style": "bezier",
            "target-arrow-shape": "none",
        },
    },
    {
        "selector": "edge.mgi-edge",
        "style": {
            "line-color": "#d69e2e",
            "width": 2,
            "opacity": 0.65,
            "curve-style": "bezier",
            "target-arrow-shape": "none",
        },
    },
    {
        "selector": "edge.mgwas-edge",
        "style": {
            "line-color": "#9f7aea",
            "width": 2,
            "opacity": 0.65,
            "curve-style": "bezier",
            "target-arrow-shape": "none",
        },
    },
    # ── Selection highlight ──
    {
        "selector": ":selected",
        "style": {
            "border-width": "4px",
            "border-color": "#d69e2e",
            "background-color": "#ecc94b",
        },
    },
]

# ---------------------------------------------------------------------------
# LEGEND HELPER
# ---------------------------------------------------------------------------

def _legend_item(colour, shape_label, label):
    """Small colour + label legend entry."""
    shapes = {
        "circle": {"borderRadius": "50%"},
        "rectangle": {"borderRadius": "4px"},
        "diamond": {"transform": "rotate(45deg)", "borderRadius": "2px"},
        "triangle": {"clipPath": "polygon(50% 0%, 0% 100%, 100% 100%)"},
        "hexagon": {"borderRadius": "3px"},
    }
    swatch = {
        "display": "inline-block",
        "width": "14px",
        "height": "14px",
        "background": colour,
        "marginRight": "6px",
        "verticalAlign": "middle",
        **shapes.get(shape_label, {}),
    }
    return html.Span(
        [html.Span(style=swatch), label],
        style={"marginRight": "16px", "fontSize": "0.82rem", "whiteSpace": "nowrap"},
    )


def _edge_legend(colour, label):
    """Line swatch + label."""
    return html.Span(
        [
            html.Span(style={
                "display": "inline-block",
                "width": "22px",
                "height": "3px",
                "background": colour,
                "marginRight": "6px",
                "verticalAlign": "middle",
                "borderRadius": "2px",
            }),
            label,
        ],
        style={"marginRight": "16px", "fontSize": "0.82rem", "whiteSpace": "nowrap"},
    )


# ---------------------------------------------------------------------------
# EXAMPLE BUTTONS
# ---------------------------------------------------------------------------

EXAMPLE_QUERIES = [
    {"name": "Glucose", "type": "metabolite"},
    {"name": "Dopamine", "type": "metabolite"},
    {"name": "Tryptophan", "type": "metabolite"},
    {"name": "TP53", "type": "protein"},
    {"name": "CYP1A2", "type": "gene"},
    {"name": "Alzheimer", "type": "disease"},
    {"name": "Lactobacillus", "type": "microbe"},
    {"name": "Cholesterol", "type": "metabolite"},
    {"name": "Metformin", "type": "drug"},
]

DEFAULT_QUERY = EXAMPLE_QUERIES[0]


# ---------------------------------------------------------------------------
# LAYOUT
# ---------------------------------------------------------------------------

page_content = html.Div(
    [
        html.Div(
            [
                make_page_header(
                    "Network Visualization",
                    "Explore multi-type interaction networks, search by metabolite, protein, disease, or microbe.",
                    [("Home", "/home"), ("Network", None)],
                ),

                dbc.Row(
                    [
                        # ======== SIDEBAR (300 px) ========
                        dbc.Col(
                            html.Div(
                                [
                                    # Search bar
                                    html.H6(
                                        [html.I(className="fas fa-search me-2"), "Search Query"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dbc.RadioItems(
                                        id="net-query-type",
                                        options=[
                                            {"label": "Metabolite", "value": "metabolite"},
                                            {"label": "Protein", "value": "protein"},
                                            {"label": "Disease", "value": "disease"},
                                            {"label": "Microbe", "value": "microbe"},
                                            {"label": "Drug", "value": "drug"},
                                            {"label": "Gene", "value": "gene"},
                                            {"label": "SNP", "value": "snp"},
                                        ],
                                        value="metabolite",
                                        inline=True,
                                        className="mb-2",
                                        style={"fontSize": "0.85rem"},
                                    ),
                                    dbc.InputGroup(
                                        [
                                            dbc.Input(
                                                id="net-query-input",
                                                placeholder="e.g. Glucose, TP53, Alzheimer …",
                                                type="text",
                                                value=DEFAULT_QUERY["name"],
                                                debounce=True,
                                            ),
                                            dbc.Button(
                                                html.I(className="fas fa-project-diagram"),
                                                id="net-search-btn",
                                                color="primary",
                                                n_clicks=0,
                                            ),
                                        ],
                                        className="mb-3",
                                    ),

                                    # Examples
                                    html.Div(
                                        [
                                            html.Small("Try: ", className="text-muted"),
                                            *[
                                                html.A(
                                                    f"{ex['name']}",
                                                    id={"type": "net-example", "index": i},
                                                    href="#",
                                                    className="me-2",
                                                    style={"fontSize": "0.82rem"},
                                                    title=ex["type"],
                                                )
                                                for i, ex in enumerate(EXAMPLE_QUERIES[:5])
                                            ],
                                        ],
                                        className="mb-3",
                                    ),

                                    html.Hr(),

                                    # Interaction-type toggles
                                    html.H6(
                                        [html.I(className="fas fa-layer-group me-2"), "Interaction Types"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dbc.Checklist(
                                        id="net-type-filter",
                                        options=[
                                            {"label": html.Span([
                                                html.Span(style={"display": "inline-block", "width": "12px",
                                                                  "height": "12px", "borderRadius": "50%",
                                                                  "background": "#3182ce", "marginRight": "6px",
                                                                  "verticalAlign": "middle"}),
                                                "Protein (MPI)",
                                            ]), "value": "MPI"},
                                            {"label": html.Span([
                                                html.Span(style={"display": "inline-block", "width": "12px",
                                                                  "height": "12px", "borderRadius": "50%",
                                                                  "background": "#e53e3e", "marginRight": "6px",
                                                                  "verticalAlign": "middle"}),
                                                "Disease (MDI)",
                                            ]), "value": "MDI"},
                                            {"label": html.Span([
                                                html.Span(style={"display": "inline-block", "width": "12px",
                                                                  "height": "12px", "borderRadius": "50%",
                                                                  "background": "#805ad5", "marginRight": "6px",
                                                                  "verticalAlign": "middle"}),
                                                "Microbe (MMI)",
                                            ]), "value": "MMI"},
                                            {"label": html.Span([
                                                html.Span(style={"display": "inline-block", "width": "12px",
                                                                  "height": "12px", "borderRadius": "50%",
                                                                  "background": "#319795", "marginRight": "6px",
                                                                  "verticalAlign": "middle"}),
                                                "Drug (MDrI)",
                                            ]), "value": "MDrI"},
                                            {"label": html.Span([
                                                html.Span(style={"display": "inline-block", "width": "12px",
                                                                  "height": "12px", "borderRadius": "50%",
                                                                  "background": "#d69e2e", "marginRight": "6px",
                                                                  "verticalAlign": "middle"}),
                                                "Gene (MGI)",
                                            ]), "value": "MGI"},
                                            {"label": html.Span([
                                                html.Span(style={"display": "inline-block", "width": "12px",
                                                                  "height": "12px", "borderRadius": "50%",
                                                                  "background": "#9f7aea", "marginRight": "6px",
                                                                  "verticalAlign": "middle"}),
                                                "SNP (mGWAS)",
                                            ]), "value": "mGWAS"},
                                        ],
                                        value=["MPI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"],
                                        className="mb-3",
                                        switch=True,
                                    ),

                                    html.Hr(),

                                    # Organism filter
                                    html.H6(
                                        [html.I(className="fas fa-globe me-2"), "Organism Filter"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dcc.Dropdown(
                                        id="net-organism-filter",
                                        placeholder="All organisms",
                                        clearable=True,
                                        options=[],  # populated by callback
                                        className="mb-3",
                                        style={"fontSize": "0.85rem"},
                                    ),

                                    html.Hr(),

                                    # Layout selector
                                    html.H6(
                                        [html.I(className="fas fa-th me-2"), "Layout"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dbc.RadioItems(
                                        id="net-layout-select",
                                        options=[
                                            {"label": "Force-directed (CoSE)", "value": "cose"},
                                            {"label": "Concentric", "value": "concentric"},
                                            {"label": "Circle", "value": "circle"},
                                            {"label": "Breadth-first", "value": "breadthfirst"},
                                            {"label": "Grid", "value": "grid"},
                                        ],
                                        value="cose",
                                        className="mb-3",
                                        inline=False,
                                    ),

                                    html.Hr(),

                                    # Max nodes slider
                                    html.H6(
                                        [html.I(className="fas fa-sliders-h me-2"), "Max Nodes"],
                                        className="fw-semibold mb-2",
                                    ),
                                    dcc.Slider(
                                        id="net-max-nodes",
                                        min=20,
                                        max=500,
                                        step=20,
                                        value=200,
                                        marks={20: "20", 100: "100", 200: "200", 300: "300", 500: "500"},
                                        tooltip={"placement": "bottom"},
                                    ),

                                    html.Hr(),

                                    # Stats panel
                                    html.Div(id="net-stats-panel"),
                                ],
                                className="cm-card net-sidebar",
                                style={
                                    "position": "sticky", "top": "80px",
                                    "maxHeight": "calc(100vh - 100px)",
                                    "overflowY": "auto",
                                },
                            ),
                            xs=12, md=3,
                        ),

                        # ======== MAIN AREA ========
                        dbc.Col(
                            [
                                # Legend row
                                html.Div(
                                    [
                                        html.Span("Nodes: ", className="fw-semibold",
                                                   style={"fontSize": "0.82rem"}),
                                        _legend_item("#ed8936", "circle", "Metabolite"),
                                        _legend_item("#3182ce", "rectangle", "Protein"),
                                        _legend_item("#e53e3e", "triangle", "Disease"),
                                        _legend_item("#805ad5", "hexagon", "Microbe"),
                                        _legend_item("#319795", "diamond", "Drug"),
                                        _legend_item("#d69e2e", "hexagon", "Gene"),
                                        _legend_item("#9f7aea", "circle", "SNP"),
                                        html.Span("│ ", style={"margin": "0 6px", "color": "#cbd5e0"}),
                                        html.Span("Edges: ", className="fw-semibold",
                                                   style={"fontSize": "0.82rem"}),
                                        _edge_legend("#3182ce", "MPI"),
                                        _edge_legend("#e53e3e", "MDI"),
                                        _edge_legend("#805ad5", "MMI"),
                                        _edge_legend("#319795", "MDrI"),
                                        _edge_legend("#d69e2e", "MGI"),
                                        _edge_legend("#9f7aea", "mGWAS"),
                                    ],
                                    style={
                                        "display": "flex",
                                        "flexWrap": "wrap",
                                        "alignItems": "center",
                                        "padding": "10px 16px",
                                        "background": "#f7fafc",
                                        "borderRadius": "8px",
                                        "marginBottom": "12px",
                                        "border": "1px solid #e2e8f0",
                                    },
                                ),

                                # Cytoscape graph
                                dcc.Loading(
                                    id="net-graph-loading",
                                    type="circle",
                                    color="#1a365d",
                                    children=html.Div(
                                        [
                                            cyto.Cytoscape(
                                                id="net-cytoscape",
                                                layout={"name": "cose", "animate": True,
                                                        "padding": 50,
                                                        "nodeRepulsion": 12000,
                                                        "idealEdgeLength": 100,
                                                        "nodeDimensionsIncludeLabels": True},
                                                stylesheet=NETWORK_STYLESHEET,
                                                elements=[],
                                                style={
                                                    "width": "100%",
                                                    "height": "600px",
                                                    "border": "1px solid #e2e8f0",
                                                    "borderRadius": "8px",
                                                    "background": "#ffffff",
                                                },
                                                responsive=True,
                                            ),
                                        ],
                                        className="cm-card mb-3",
                                    ),
                                ),

                                # Node detail panel (appears on click)
                                html.Div(id="net-node-detail", style={"display": "none"},
                                         className="cm-card mb-3"),

                                # Type-breakdown stat badges
                                html.Div(id="net-type-badges", className="mb-3"),
                            ],
                            xs=12, md=9,
                        ),
                    ],
                    className="g-3",
                ),
            ],
            className="cm-page-container",
        ),

        # Hidden stores
        dcc.Store(id="net-elements-store"),
        dcc.Store(id="net-stats-store"),
    ]
)


# ---------------------------------------------------------------------------
# CALLBACKS
# ---------------------------------------------------------------------------

@callback(
    Output("net-organism-filter", "options"),
    Input("url", "pathname"),
)
def populate_organism_options(_pathname):
    """Load organism options once on page load."""
    try:
        from app.services.network_service import get_available_organisms
        orgs = get_available_organisms()
        return [{"label": o, "value": o} for o in orgs]
    except Exception:
        return []


@callback(
    Output("net-query-input", "value"),
    Output("net-query-type", "value"),
    Input({"type": "net-example", "index": 0}, "n_clicks"),
    Input({"type": "net-example", "index": 1}, "n_clicks"),
    Input({"type": "net-example", "index": 2}, "n_clicks"),
    Input({"type": "net-example", "index": 3}, "n_clicks"),
    Input({"type": "net-example", "index": 4}, "n_clicks"),
    prevent_initial_call=True,
)
def set_example_query(*args):
    """Fill input and set query type when an example link is clicked."""
    triggered = ctx.triggered_id
    if triggered and isinstance(triggered, dict):
        idx = triggered.get("index", 0)
        if 0 <= idx < len(EXAMPLE_QUERIES):
            ex = EXAMPLE_QUERIES[idx]
            return ex["name"], ex["type"]
    return no_update, no_update


@callback(
    Output("net-elements-store", "data"),
    Output("net-stats-store", "data"),
    Input("url", "pathname"),
    Input("net-search-btn", "n_clicks"),
    Input("net-query-input", "n_submit"),
    Input({"type": "net-example", "index": 0}, "n_clicks"),
    Input({"type": "net-example", "index": 1}, "n_clicks"),
    Input({"type": "net-example", "index": 2}, "n_clicks"),
    Input({"type": "net-example", "index": 3}, "n_clicks"),
    Input({"type": "net-example", "index": 4}, "n_clicks"),
    State("net-query-input", "value"),
    State("net-query-type", "value"),
    State("net-type-filter", "value"),
    State("net-organism-filter", "value"),
    State("net-max-nodes", "value"),
)
def run_network_search(
    pathname, _n_clicks, _n_submit,
    _ex0, _ex1, _ex2, _ex3, _ex4,
    query, query_type, types, organism, max_nodes,
):
    """Search all databases and build network elements."""
    triggered = ctx.triggered_id
    if triggered == "url" and pathname not in ("/network", "/Network"):
        return no_update, no_update

    if isinstance(triggered, dict) and triggered.get("type") == "net-example":
        idx = triggered.get("index", 0)
        if 0 <= idx < len(EXAMPLE_QUERIES):
            ex = EXAMPLE_QUERIES[idx]
            query = ex["name"]
            query_type = ex["type"]

    if not query:
        query = DEFAULT_QUERY["name"]
        query_type = DEFAULT_QUERY["type"]

    if not query or len(query.strip()) < 2:
        return no_update, no_update

    if triggered in (None, "url"):
        # Keep the automatic first render lightweight on small Render instances.
        # Users can still run all-layer searches via the search button/toggles.
        types = ["MPI"]
        max_nodes = min(max_nodes or 100, 100)

    from app.services.network_service import build_network_elements

    result = build_network_elements(
        query,
        query_type=query_type or "metabolite",
        include_types=types or ["MPI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"],
        organism_filter=organism,
        max_nodes=max_nodes or 200,
    )

    return result.get("elements", []), result.get("stats", {})


@callback(
    Output("net-cytoscape", "elements"),
    Output("net-cytoscape", "layout"),
    Input("net-elements-store", "data"),
    Input("net-layout-select", "value"),
)
def update_graph(elements, layout_name):
    """Push elements to Cytoscape and update layout."""
    layout_opts = {
        "name": layout_name or "cose",
        "animate": True,
        "padding": 50,
        "nodeDimensionsIncludeLabels": True,
    }
    if layout_name == "cose":
        layout_opts.update({"nodeRepulsion": 12000, "idealEdgeLength": 100})
    elif layout_name == "concentric":
        layout_opts.update({
            "concentric": "function(node){ return node.data('degree') || 0; }",
            "levelWidth": "function(nodes){ return 2; }",
        })
    elif layout_name == "breadthfirst":
        layout_opts.update({"directed": False, "spacingFactor": 1.2})

    return elements or [], layout_opts


@callback(
    Output("net-stats-panel", "children"),
    Output("net-type-badges", "children"),
    Input("net-stats-store", "data"),
)
def update_stats(stats):
    """Render stats sidebar and type badges."""
    if not stats:
        return (
            html.Div(
                [
                    html.I(className="fas fa-info-circle me-2 text-muted"),
                    html.Span("Enter a metabolite name to build the network.",
                              className="text-muted"),
                ],
                style={"padding": "12px"},
            ),
            [],
        )

    tc = stats.get("type_counts", {})
    total_n = stats.get("total_nodes", 0)
    total_e = stats.get("total_edges", 0)
    q_label = stats.get("query_label", stats.get("metabolite_label", ""))
    q_type = stats.get("query_type", "metabolite")
    icon_cls = {"metabolite": "fas fa-atom", "protein": "fas fa-dna",
                "disease": "fas fa-heartbeat", "microbe": "fas fa-bacterium"}.get(q_type, "fas fa-atom")

    # Sidebar stats
    sidebar = html.Div([
        html.Div(
            [
                html.I(className=f"{icon_cls} me-2", style={"color": "#ed8936"}),
                html.Strong(q_label),
            ],
            className="mb-2",
            style={"fontSize": "1rem"},
        ),
        html.Div(
            f"{total_n:,} nodes  ·  {total_e:,} edges",
            style={"fontSize": "0.85rem", "color": "#718096", "marginBottom": "8px"},
        ),
        *[
            html.Div(
                [
                    html.Span(
                        style={
                            "display": "inline-block",
                            "width": "10px",
                            "height": "10px",
                            "borderRadius": "50%",
                            "background": {"MPI": "#3182ce",
                                           "MDI": "#e53e3e", "MMI": "#805ad5",
                                           "MDrI": "#319795", "MGI": "#d69e2e",
                                           "mGWAS": "#9f7aea"}.get(t, "#999"),
                            "marginRight": "6px",
                        },
                    ),
                    html.Span(f"{t}: ", style={"fontWeight": "600"}),
                    html.Span(f"{c:,} interactions"),
                ],
                style={"fontSize": "0.85rem", "marginBottom": "4px"},
            )
            for t, c in tc.items() if c > 0
        ],
    ])

    # Badges below graph
    badge_colours = {"MPI": "primary", "MDI": "danger", "MMI": "info", "MDrI": "success",
                      "MGI": "warning", "mGWAS": "secondary"}
    badges = html.Div(
        [
            dbc.Badge(
                f"{t}  {c:,}",
                color=badge_colours.get(t, "secondary"),
                className="me-2",
                pill=True,
                style={"fontSize": "0.85rem"},
            )
            for t, c in tc.items() if c > 0
        ],
        className="d-flex flex-wrap",
    )

    return sidebar, badges


@callback(
    Output("net-node-detail", "children"),
    Output("net-node-detail", "style"),
    Input("net-cytoscape", "tapNodeData"),
)
def show_node_detail(data):
    """Show a detail card when a node is clicked."""
    if not data:
        return no_update, {"display": "none"}

    node_type = data.get("node_type", "unknown")
    label = data.get("label", "")
    degree = data.get("degree", 0)
    nid = data.get("id", "")

    icon_map = {
        "metabolite": "fas fa-atom",
        "protein": "fas fa-dna",
        "disease": "fas fa-heartbeat",
        "microbe": "fas fa-bacterium",
        "drug": "fas fa-pills",
        "gene": "fas fa-dna",
        "snp": "fas fa-map-marker-alt",
    }
    colour_map = {
        "metabolite": "#ed8936",
        "protein": "#3182ce",
        "disease": "#e53e3e",
        "microbe": "#805ad5",
        "drug": "#319795",
        "gene": "#d69e2e",
        "snp": "#9f7aea",
    }

    icon = icon_map.get(node_type, "fas fa-circle")
    colour = colour_map.get(node_type, "#718096")
    type_display = node_type.capitalize()

    # Build identifier info
    id_parts = nid.split(":", 1) if ":" in nid else ("", nid)
    id_value = id_parts[1] if len(id_parts) > 1 else nid

    children = [
        html.Div(
            [
                html.I(className=f"{icon} me-2", style={"color": colour, "fontSize": "1.2rem"}),
                html.Strong(label, style={"fontSize": "1rem"}),
                dbc.Badge(type_display, className="ms-2", style={"background": colour}),
                html.A(
                    html.I(className="fas fa-times"),
                    id="net-close-detail",
                    style={"marginLeft": "auto", "cursor": "pointer", "color": "#a0aec0"},
                ),
            ],
            className="d-flex align-items-center mb-2",
        ),
        html.Div(
            [
                html.Span("ID: ", style={"fontWeight": "600", "color": "#718096", "fontSize": "0.82rem"}),
                html.Code(id_value, style={"fontSize": "0.82rem"}),
                html.Span(f"  ·  {degree} connection{'s' if degree != 1 else ''}",
                          style={"color": "#718096", "fontSize": "0.82rem", "marginLeft": "8px"}),
            ],
        ),
    ]

    # External links for known types
    if node_type == "protein" and id_value.startswith(("P", "Q", "A", "O")):
        children.append(html.Div(
            html.A(f"UniProt: {id_value}", href=f"https://www.uniprot.org/uniprot/{id_value}",
                   target="_blank", style={"fontSize": "0.82rem"}),
            className="mt-1",
        ))
    elif node_type == "metabolite":
        children.append(html.Div(
            html.A("Search in Database →", href=f"/database?q={label}",
                   style={"fontSize": "0.82rem"}),
            className="mt-1",
        ))
    elif node_type == "disease":
        children.append(html.Div(
            html.A("View Disease Details →", href="/disease",
                   style={"fontSize": "0.82rem"}),
            className="mt-1",
        ))
    elif node_type == "drug" and id_value.startswith("DB"):
        children.append(html.Div(
            html.A(f"DrugBank: {id_value}", href=f"https://go.drugbank.com/drugs/{id_value}",
                   target="_blank", style={"fontSize": "0.82rem"}),
            className="mt-1",
        ))
    elif node_type == "gene":
        children.append(html.Div(
            html.A(f"NCBI Gene: {id_value}", href=f"https://www.ncbi.nlm.nih.gov/gene/?term={id_value}",
                   target="_blank", style={"fontSize": "0.82rem"}),
            className="mt-1",
        ))
    elif node_type == "snp" and id_value.startswith("rs"):
        children.append(html.Div(
            html.A(f"dbSNP: {id_value}", href=f"https://www.ncbi.nlm.nih.gov/snp/{id_value}",
                   target="_blank", style={"fontSize": "0.82rem"}),
            className="mt-1",
        ))

    return children, {"display": "block"}
