"""
CoreMet-AI, Interactive AI Query Page

A dedicated AI-assisted graph exploration interface.
Users submit natural language questions; the system returns:
  1. Grounded text summary (center)
  2. Parsed query plan (right, transparency layer)
  3. Interactive subgraph + evidence table (bottom)
"""

import json

import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
from dash import Input, Output, State, callback, dcc, html, no_update

# ── Entity type colors (consistent with site) ────────────────────────
TYPE_COLORS = {
    "metabolite": "#e27a3f",
    "protein": "#3182ce",
    "gene": "#d69e2e",
    "disease": "#e53e3e",
    "microbe": "#38a169",
    "drug": "#805ad5",
    "snp": "#319795",
}

# ── Example queries ───────────────────────────────────────────────────
EXAMPLE_QUERIES = [
    "How does butyrate influence colorectal cancer?",
    "What genes are connected to metformin?",
    "Relationship between tryptophan and depression",
    "Role of Lactobacillus in short-chain fatty acids",
    "Connect rs1260326 to metabolic syndrome",
    "What drugs interact with glutathione?",
]

# ── Cytoscape stylesheet ─────────────────────────────────────────────
_CYTO_STYLE = [
    {
        "selector": "node",
        "style": {
            "label": "data(label)",
            "font-family": "Arial, Helvetica, sans-serif",
            "font-size": "11px",
            "text-valign": "center",
            "text-halign": "center",
            "background-color": "data(color)",
            "color": "#1a202c",
            "width": "40px",
            "height": "40px",
            "border-width": "2px",
            "border-color": "#fff",
            "text-outline-width": 1,
            "text-outline-color": "#fff",
            "text-wrap": "ellipsis",
            "text-max-width": "80px",
        },
    },
    {
        "selector": "node[?is_source]",
        "style": {
            "width": "55px",
            "height": "55px",
            "font-size": "13px",
            "font-weight": "bold",
            "border-width": "3px",
            "border-color": "#2d3748",
        },
    },
    {
        "selector": "node[?is_target]",
        "style": {
            "width": "50px",
            "height": "50px",
            "font-size": "12px",
            "font-weight": "bold",
            "border-width": "3px",
            "border-color": "#e53e3e",
            "border-style": "dashed",
        },
    },
    {
        "selector": "edge",
        "style": {
            "width": 2,
            "line-color": "#a0aec0",
            "target-arrow-color": "#a0aec0",
            "target-arrow-shape": "triangle",
            "curve-style": "bezier",
            "opacity": 0.7,
            "font-family": "Arial, Helvetica, sans-serif",
            "font-size": "9px",
        },
    },
    {
        "selector": "edge[?is_ranked]",
        "style": {
            "width": 3,
            "line-color": "#2b6cb0",
            "target-arrow-color": "#2b6cb0",
            "opacity": 1.0,
            "label": "data(layer)",
        },
    },
    {
        "selector": ":selected",
        "style": {
            "border-color": "#2b6cb0",
            "border-width": "4px",
        },
    },
]


def _build_query_input():
    """Build the query input section with model selection."""
    return html.Div([
        html.Div([
            html.H2("CoreMet-AI", style={"fontWeight": "700", "color": "#1a202c",
                                          "fontFamily": "Arial, Helvetica, sans-serif"}),
            html.P(
                "Ask a metabolite-centered biological question. "
                "The system retrieves evidence-grounded paths from CoreMet, "
                "ranks them, and generates a traceable summary.",
                style={"color": "#4a5568", "fontSize": "0.95rem",
                       "fontFamily": "Arial, Helvetica, sans-serif"},
            ),
        ], style={"marginBottom": "1rem"}),

        # ── Model selection row ───────────────────────────────────
        html.Div([
            dbc.Row([
                dbc.Col([
                    html.Label("AI Model", style={
                        "fontWeight": "600", "fontSize": "0.85rem",
                        "color": "#4a5568", "fontFamily": "Arial, Helvetica, sans-serif",
                        "marginBottom": "4px",
                    }),
                    dbc.Select(
                        id="ai-model-select",
                        options=[
                            {"label": "Google Gemini Flash (free)", "value": "google-free"},
                            {"label": "OpenAI GPT-4o-mini", "value": "openai-gpt-4o-mini"},
                            {"label": "OpenAI GPT-4o", "value": "openai-gpt-4o"},
                            {"label": "OpenAI GPT-4-turbo", "value": "openai-gpt-4-turbo"},
                            {"label": "Google Gemini Pro", "value": "google-gemini-pro"},
                            {"label": "Template (no LLM)", "value": "template"},
                        ],
                        value="google-free",
                        style={
                            "fontFamily": "Arial, Helvetica, sans-serif",
                            "fontSize": "0.9rem",
                        },
                    ),
                ], md=4),
                dbc.Col([
                    html.Label("API Key (optional for free tier)", style={
                        "fontWeight": "600", "fontSize": "0.85rem",
                        "color": "#4a5568", "fontFamily": "Arial, Helvetica, sans-serif",
                        "marginBottom": "4px",
                    }),
                    dbc.Input(
                        id="ai-api-key-input",
                        type="password",
                        placeholder="sk-... or AIza... (stored locally, never sent to our server)",
                        style={
                            "fontFamily": "Arial, Helvetica, sans-serif",
                            "fontSize": "0.9rem",
                        },
                    ),
                ], md=8),
            ], className="mb-2"),
            html.Div([
                html.I(className="fas fa-info-circle me-1", style={"color": "#a0aec0"}),
                html.Span(
                    "Free tier uses Google Gemini Flash, no API key required. "
                    "For OpenAI models, paste your key above. "
                    "Keys are used client-side only and never stored on our server.",
                    style={"fontSize": "0.8rem", "color": "#a0aec0",
                           "fontFamily": "Arial, Helvetica, sans-serif"},
                ),
            ], style={"marginBottom": "0.75rem"}),
        ], style={
            "backgroundColor": "#f7fafc", "border": "1px solid #e2e8f0",
            "borderRadius": "8px", "padding": "1rem",
            "marginBottom": "1rem",
        }),

        # Query input
        dbc.InputGroup([
            dbc.Textarea(
                id="ai-query-input",
                placeholder="Ask a metabolite-centered biological question...",
                style={
                    "resize": "none", "height": "60px",
                    "fontFamily": "Arial, Helvetica, sans-serif",
                    "fontSize": "1rem", "borderRadius": "8px 0 0 8px",
                },
                debounce=True,
            ),
            dbc.Button(
                [html.I(className="fas fa-search me-2"), "Query"],
                id="ai-submit-btn",
                color="primary",
                style={
                    "borderRadius": "0 8px 8px 0", "minWidth": "100px",
                    "fontFamily": "Arial, Helvetica, sans-serif",
                    "fontWeight": "600",
                },
            ),
        ], style={"marginBottom": "0.75rem"}),

        # Example queries
        html.Div([
            html.Span("Try: ", style={"color": "#718096", "fontSize": "0.85rem",
                                       "fontFamily": "Arial, Helvetica, sans-serif"}),
            *[
                html.Span(
                    q,
                    id={"type": "ai-example", "index": i},
                    className="cm-example-chip",
                    n_clicks=0,
                    style={
                        "display": "inline-block", "padding": "2px 10px",
                        "margin": "2px 4px", "borderRadius": "12px",
                        "fontSize": "0.8rem", "cursor": "pointer",
                        "backgroundColor": "#edf2f7", "color": "#2d3748",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "transition": "background 0.2s",
                    },
                )
                for i, q in enumerate(EXAMPLE_QUERIES)
            ],
        ]),

        # Loading spinner
        dcc.Loading(
            id="ai-loading",
            type="dot",
            color="#3182ce",
            children=html.Div(id="ai-results-container"),
        ),

        # Hidden store for results
        dcc.Store(id="ai-result-store"),
    ], className="cm-page-container", style={"maxWidth": "1200px", "margin": "0 auto",
                                               "padding": "2rem 1.5rem"})


def _build_results_section():
    """Build the results display area (initially hidden)."""
    return html.Div([
        # ── Top row: Summary (center) + Query Plan (right) ────────
        dbc.Row([
            # Summary panel
            dbc.Col([
                html.Div([
                    html.H5("Answer", style={
                        "fontWeight": "700", "marginBottom": "0.5rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "borderBottom": "2px solid #3182ce", "paddingBottom": "0.5rem",
                    }),
                    html.Div(id="ai-summary-text", style={
                        "fontSize": "0.95rem", "lineHeight": "1.7",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "color": "#2d3748",
                    }),
                    html.Div(id="ai-confidence-badge", style={"marginTop": "0.5rem"}),
                    html.Div(id="ai-timing-info", style={
                        "fontSize": "0.75rem", "color": "#a0aec0", "marginTop": "0.5rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                    }),
                ], style={
                    "backgroundColor": "#fff", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1.25rem",
                    "minHeight": "150px",
                }),
            ], md=7),

            # Query plan panel (transparency layer)
            dbc.Col([
                html.Div([
                    html.H5("Query Plan", style={
                        "fontWeight": "700", "marginBottom": "0.5rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "borderBottom": "2px solid #38a169", "paddingBottom": "0.5rem",
                    }),
                    html.Div(id="ai-query-plan", style={
                        "fontSize": "0.85rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                    }),
                ], style={
                    "backgroundColor": "#f7fafc", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1.25rem",
                    "minHeight": "150px",
                }),
            ], md=5),
        ], className="mb-3"),

        # ── Bottom: Graph + Evidence ──────────────────────────────
        dbc.Row([
            # Graph panel
            dbc.Col([
                html.Div([
                    html.Div([
                        html.H5("Subgraph", style={
                            "fontWeight": "700", "fontFamily": "Arial, Helvetica, sans-serif",
                            "display": "inline-block",
                        }),
                        # Layer filter checkboxes
                        html.Div(id="ai-layer-filters", style={
                            "display": "inline-block", "marginLeft": "1rem",
                        }),
                    ], style={"marginBottom": "0.5rem",
                              "borderBottom": "2px solid #805ad5", "paddingBottom": "0.5rem"}),
                    cyto.Cytoscape(
                        id="ai-cytoscape",
                        layout={"name": "cose", "animate": True, "animationDuration": 500,
                                "nodeRepulsion": 8000, "idealEdgeLength": 80,
                                "gravity": 0.3, "padding": 30},
                        style={"width": "100%", "height": "420px",
                               "border": "1px solid #e2e8f0", "borderRadius": "8px",
                               "backgroundColor": "#fafafa"},
                        stylesheet=_CYTO_STYLE,
                        elements=[],
                        responsive=False,
                        autoRefreshLayout=False,
                    ),
                    # Node/edge detail drawer
                    html.Div(id="ai-selection-detail", style={
                        "marginTop": "0.5rem", "fontSize": "0.85rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "color": "#4a5568", "minHeight": "40px",
                    }),
                ], style={
                    "backgroundColor": "#fff", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1rem",
                }),
            ], md=7),

            # Evidence panel
            dbc.Col([
                html.Div([
                    html.H5("Evidence", style={
                        "fontWeight": "700", "fontFamily": "Arial, Helvetica, sans-serif",
                        "borderBottom": "2px solid #e53e3e", "paddingBottom": "0.5rem",
                        "marginBottom": "0.5rem",
                    }),
                    html.Div(id="ai-evidence-stats", style={
                        "marginBottom": "0.5rem",
                    }),
                    html.Div(
                        id="ai-evidence-table-container",
                        style={"maxHeight": "380px", "overflowY": "auto"},
                    ),
                    dbc.Button(
                        [html.I(className="fas fa-download me-2"), "Export CSV"],
                        id="ai-export-btn",
                        size="sm",
                        color="secondary",
                        outline=True,
                        style={"marginTop": "0.5rem",
                               "fontFamily": "Arial, Helvetica, sans-serif"},
                    ),
                    dcc.Download(id="ai-download-csv"),
                ], style={
                    "backgroundColor": "#fff", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1rem",
                }),
            ], md=5),
        ]),
    ])


# ── Page layout ───────────────────────────────────────────────────────

page_content = html.Div([
    _build_query_input(),
    html.Div(id="ai-results-wrapper", style={"padding": "0 1.5rem", "maxWidth": "1200px",
                                               "margin": "0 auto"}),
])


# ── Callbacks ─────────────────────────────────────────────────────────

@callback(
    Output("ai-query-input", "value"),
    Input({"type": "ai-example", "index": 0}, "n_clicks"),
    Input({"type": "ai-example", "index": 1}, "n_clicks"),
    Input({"type": "ai-example", "index": 2}, "n_clicks"),
    Input({"type": "ai-example", "index": 3}, "n_clicks"),
    Input({"type": "ai-example", "index": 4}, "n_clicks"),
    Input({"type": "ai-example", "index": 5}, "n_clicks"),
    prevent_initial_call=True,
)
def fill_example(*clicks):
    """Fill the query input with an example query when clicked."""
    from dash import ctx
    triggered = ctx.triggered_id
    if triggered and isinstance(triggered, dict):
        idx = triggered.get("index", 0)
        if 0 <= idx < len(EXAMPLE_QUERIES):
            return EXAMPLE_QUERIES[idx]
    return no_update


@callback(
    Output("ai-result-store", "data"),
    Output("ai-results-wrapper", "children"),
    Input("ai-submit-btn", "n_clicks"),
    State("ai-query-input", "value"),
    State("ai-model-select", "value"),
    State("ai-api-key-input", "value"),
    prevent_initial_call=True,
)
def run_ai_query(n_clicks, query, model_selection, user_api_key):
    """Execute the AI pipeline and render results."""
    if not query or not query.strip():
        return no_update, no_update

    # Parse model selection → provider + model
    provider = "template"
    model = ""
    api_key = user_api_key or ""

    if model_selection == "google-free":
        provider = "google"
        model = "gemini-2.0-flash"
        # Use env key if user didn't paste one
        if not api_key:
            import os
            api_key = os.environ.get("GOOGLE_API_KEY", "")
    elif model_selection and model_selection.startswith("openai-"):
        provider = "openai"
        model = model_selection.replace("openai-", "")
        if not api_key:
            import os
            api_key = os.environ.get("OPENAI_API_KEY", "")
    elif model_selection and model_selection.startswith("google-"):
        provider = "google"
        model = model_selection.replace("google-", "")
        if not api_key:
            import os
            api_key = os.environ.get("GOOGLE_API_KEY", "")
    elif model_selection == "template":
        provider = "template"

    from app.services.ai_orchestrator import execute_ai_query
    result = execute_ai_query(query.strip(), provider=provider,
                              api_key=api_key, model=model)

    # Build the results UI
    children = _render_results(result)
    return result, children


def _render_results(result: dict):
    """Render the full results section from an AI pipeline result."""
    status = result.get("status", "error")
    plan = result.get("query_plan", {})
    summary_data = result.get("summary", {})
    evidence = result.get("evidence", {})
    subgraph = result.get("subgraph", {"nodes": [], "edges": []})
    ranked_paths = result.get("ranked_paths", [])
    timing = result.get("timing_ms", 0)

    # ── Summary text ──────────────────────────────────────────────
    summary_text = summary_data.get("summary", "No results.")
    summary_method = summary_data.get("method", "template")
    confidence_label = summary_data.get("confidence_label", "")

    conf_color = {"supported": "#38a169", "moderate evidence": "#d69e2e",
                  "hypothesis-generating": "#e53e3e", "no data": "#a0aec0",
                  "error": "#e53e3e"}.get(confidence_label, "#a0aec0")

    summary_div = html.Div([
        dcc.Markdown(summary_text, style={
            "fontFamily": "Arial, Helvetica, sans-serif",
            "fontSize": "0.95rem", "lineHeight": "1.7",
        }),
    ])

    confidence_badge = html.Span([
        html.Span(confidence_label.upper(), style={
            "backgroundColor": conf_color, "color": "#fff",
            "padding": "2px 10px", "borderRadius": "12px",
            "fontSize": "0.75rem", "fontWeight": "600",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }),
        html.Span(f"  via {summary_method}", style={
            "fontSize": "0.75rem", "color": "#a0aec0", "marginLeft": "0.5rem",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }),
        html.Span(f"  ({summary_data.get('model', 'template')})", style={
            "fontSize": "0.75rem", "color": "#a0aec0",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }) if summary_data.get("model") and summary_data["model"] != "template" else None,
    ])

    timing_div = html.Span(f"Completed in {timing} ms")

    # ── Query plan ────────────────────────────────────────────────
    plan_items = []
    if plan.get("source_entity"):
        plan_items.append(_plan_row("Source", plan["source_entity"],
                                     plan.get("source_type", "")))
    if plan.get("target_entity"):
        plan_items.append(_plan_row("Target", plan["target_entity"],
                                     plan.get("target_type", "")))
    plan_items.append(_plan_row("Intent", plan.get("intent", "-"), ""))
    plan_items.append(_plan_row("Max hops", str(plan.get("max_hops", 2)), ""))
    layers = plan.get("layers", [])
    plan_items.append(_plan_row("Layers", ", ".join(layers), ""))
    plan_items.append(_plan_row("Parser", plan.get("parser", "rule-based"), ""))

    # Resolution info
    src_res = result.get("resolved_source")
    tgt_res = result.get("resolved_target")
    if src_res:
        plan_items.append(_plan_row(
            "Resolved source",
            f"{src_res['name']} ({src_res['type']})",
            f"conf: {src_res['confidence']}"
        ))
    if tgt_res:
        plan_items.append(_plan_row(
            "Resolved target",
            f"{tgt_res['name']} ({tgt_res['type']})",
            f"conf: {tgt_res['confidence']}"
        ))

    plan_div = html.Div(plan_items)

    # ── Cytoscape graph ───────────────────────────────────────────
    elements = _build_cytoscape_elements(subgraph, plan, ranked_paths)

    # Layer filter chips
    layer_counts = evidence.get("layer_distribution", {})
    layer_chips = [
        html.Span(
            f"{layer} ({count})",
            style={
                "display": "inline-block", "padding": "1px 8px",
                "margin": "0 3px", "borderRadius": "10px",
                "fontSize": "0.75rem", "backgroundColor": "#edf2f7",
                "fontFamily": "Arial, Helvetica, sans-serif",
            },
        )
        for layer, count in layer_counts.items()
    ]

    # ── Evidence stats ────────────────────────────────────────────
    ev_stats = html.Div([
        html.Span(f"{evidence.get('total_edges', 0)} edges", style={
            "fontWeight": "600", "marginRight": "1rem",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }),
        html.Span(f"{evidence.get('total_pmids', 0)} PMIDs", style={
            "fontWeight": "600", "marginRight": "1rem",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }),
        html.Span(f"{len(ranked_paths)} ranked paths", style={
            "fontWeight": "600",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }),
    ], style={"fontSize": "0.85rem", "color": "#4a5568"})

    # ── Evidence table ────────────────────────────────────────────
    ev_table_rows = []
    for row in evidence.get("evidence_table", [])[:50]:
        type_color = TYPE_COLORS.get(row.get("target_type", ""), "#a0aec0")
        ev_table_rows.append(
            html.Tr([
                html.Td(str(row.get("rank", "")), style={"width": "30px"}),
                html.Td(row.get("source_entity", ""), style={"fontWeight": "500"}),
                html.Td("→"),
                html.Td(row.get("target_entity", ""), style={"fontWeight": "500"}),
                html.Td(
                    html.Span(row.get("layer", ""), style={
                        "backgroundColor": "#edf2f7", "padding": "1px 6px",
                        "borderRadius": "8px", "fontSize": "0.75rem",
                    }),
                ),
                html.Td(row.get("subtype", "")[:30]),
                html.Td(f"{row.get('confidence', 0):.2f}"),
                html.Td(
                    html.A(
                        row.get("pmid", "")[:15],
                        href=f"https://pubmed.ncbi.nlm.nih.gov/{row.get('pmid', '')}",
                        target="_blank",
                        style={"color": "#3182ce", "textDecoration": "none"},
                    ) if row.get("pmid") and row["pmid"] != "nan" else "-"
                ),
            ], style={"fontSize": "0.8rem", "fontFamily": "Arial, Helvetica, sans-serif"})
        )

    ev_table = html.Table([
        html.Thead(html.Tr([
            html.Th("#"), html.Th("Source"), html.Th(""), html.Th("Target"),
            html.Th("Layer"), html.Th("Subtype"), html.Th("Conf."), html.Th("PMID"),
        ], style={"fontSize": "0.75rem", "color": "#718096", "borderBottom": "2px solid #e2e8f0",
                  "fontFamily": "Arial, Helvetica, sans-serif"})),
        html.Tbody(ev_table_rows),
    ], style={"width": "100%", "borderCollapse": "collapse"}) if ev_table_rows else html.P(
        "No evidence records.", style={"color": "#a0aec0", "fontStyle": "italic"}
    )

    # ── Assemble ──────────────────────────────────────────────────
    return html.Div([
        # Top row
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.H5("Answer", style={
                        "fontWeight": "700", "marginBottom": "0.5rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "borderBottom": "2px solid #3182ce", "paddingBottom": "0.5rem",
                    }),
                    summary_div,
                    html.Div([confidence_badge], style={"marginTop": "0.5rem"}),
                    html.Div([timing_div], style={
                        "fontSize": "0.75rem", "color": "#a0aec0", "marginTop": "0.5rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                    }),
                ], style={
                    "backgroundColor": "#fff", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1.25rem",
                }),
            ], md=7),
            dbc.Col([
                html.Div([
                    html.H5("Query Plan", style={
                        "fontWeight": "700", "marginBottom": "0.5rem",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "borderBottom": "2px solid #38a169", "paddingBottom": "0.5rem",
                    }),
                    plan_div,
                ], style={
                    "backgroundColor": "#f7fafc", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1.25rem",
                }),
            ], md=5),
        ], className="mb-3"),

        # Bottom row
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.Div([
                        html.H5("Subgraph", style={
                            "fontWeight": "700",
                            "fontFamily": "Arial, Helvetica, sans-serif",
                            "display": "inline-block",
                        }),
                        html.Span(layer_chips, style={"marginLeft": "1rem"}),
                    ], style={"marginBottom": "0.5rem",
                              "borderBottom": "2px solid #805ad5", "paddingBottom": "0.5rem"}),
                    cyto.Cytoscape(
                        id="ai-cytoscape",
                        layout={"name": "cose", "animate": True, "animationDuration": 500,
                                "nodeRepulsion": 8000, "idealEdgeLength": 80,
                                "gravity": 0.3, "padding": 30},
                        style={"width": "100%", "height": "420px",
                               "border": "1px solid #e2e8f0", "borderRadius": "8px",
                               "backgroundColor": "#fafafa"},
                        stylesheet=_CYTO_STYLE,
                        elements=elements,
                        responsive=False,
                        autoRefreshLayout=False,
                    ),
                    # Legend
                    html.Div([
                        html.Span([
                            html.Span("", style={
                                "display": "inline-block", "width": "10px", "height": "10px",
                                "borderRadius": "50%", "backgroundColor": color,
                                "marginRight": "4px", "verticalAlign": "middle",
                            }),
                            html.Span(etype.capitalize(), style={
                                "fontSize": "0.75rem", "marginRight": "12px",
                                "fontFamily": "Arial, Helvetica, sans-serif",
                                "verticalAlign": "middle",
                            }),
                        ])
                        for etype, color in TYPE_COLORS.items()
                    ], style={"marginTop": "0.5rem", "textAlign": "center"}),
                ], style={
                    "backgroundColor": "#fff", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1rem",
                }),
            ], md=7),
            dbc.Col([
                html.Div([
                    html.H5("Evidence", style={
                        "fontWeight": "700",
                        "fontFamily": "Arial, Helvetica, sans-serif",
                        "borderBottom": "2px solid #e53e3e", "paddingBottom": "0.5rem",
                        "marginBottom": "0.5rem",
                    }),
                    ev_stats,
                    html.Div(ev_table, style={"maxHeight": "380px", "overflowY": "auto",
                                               "marginTop": "0.5rem"}),
                ], style={
                    "backgroundColor": "#fff", "border": "1px solid #e2e8f0",
                    "borderRadius": "8px", "padding": "1rem",
                }),
            ], md=5),
        ]),
    ], style={"marginTop": "1.5rem"})


def _plan_row(label: str, value: str, detail: str):
    """Build a single row in the query plan display."""
    return html.Div([
        html.Span(f"{label}: ", style={
            "fontWeight": "600", "color": "#4a5568",
            "fontFamily": "Arial, Helvetica, sans-serif",
            "fontSize": "0.85rem",
        }),
        html.Span(value, style={
            "color": "#1a202c", "fontFamily": "Arial, Helvetica, sans-serif",
            "fontSize": "0.85rem",
        }),
        html.Span(f"  {detail}", style={
            "color": "#a0aec0", "fontSize": "0.75rem", "marginLeft": "0.5rem",
            "fontFamily": "Arial, Helvetica, sans-serif",
        }) if detail else None,
    ], style={"marginBottom": "4px"})


def _build_cytoscape_elements(subgraph: dict, plan: dict, ranked_paths: list) -> list:
    """Convert subgraph dict to Cytoscape elements."""
    elements = []
    source_name = plan.get("source_entity", "").lower()
    target_name = plan.get("target_entity", "").lower()

    # Collect ranked edge keys for highlighting
    ranked_edge_keys = set()
    for item in ranked_paths:
        for edge in item.get("path", []):
            src_key = f"{edge.get('source_type', '')}:{edge.get('source_name', '')}"
            tgt_key = f"{edge.get('target_type', '')}:{edge.get('target_name', '')}"
            ranked_edge_keys.add(f"{src_key}→{tgt_key}")

    for node in subgraph.get("nodes", []):
        ntype = node.get("type", "metabolite")
        nname = node.get("name", "")
        is_source = nname.lower() == source_name
        is_target = nname.lower() == target_name
        elements.append({
            "data": {
                "id": node["id"],
                "label": nname[:20],
                "color": TYPE_COLORS.get(ntype, "#a0aec0"),
                "type": ntype,
                "is_source": is_source,
                "is_target": is_target,
                "full_name": nname,
            },
        })

    for edge in subgraph.get("edges", []):
        ek = f"{edge['source']}→{edge['target']}"
        elements.append({
            "data": {
                "source": edge["source"],
                "target": edge["target"],
                "layer": edge.get("layer", ""),
                "subtype": edge.get("subtype", ""),
                "confidence": edge.get("confidence", 0),
                "pmid": edge.get("pmid", ""),
                "is_ranked": ek in ranked_edge_keys,
            },
        })

    return elements


# ── Export callback ───────────────────────────────────────────────────

@callback(
    Output("ai-download-csv", "data"),
    Input("ai-export-btn", "n_clicks"),
    State("ai-result-store", "data"),
    prevent_initial_call=True,
)
def export_evidence_csv(n_clicks, result):
    """Export the evidence table as CSV."""
    if not result:
        return no_update
    evidence = result.get("evidence", {})
    table = evidence.get("evidence_table", [])
    if not table:
        return no_update

    import pandas as pd
    df = pd.DataFrame(table)
    return dcc.send_data_frame(df.to_csv, "coremetai_evidence.csv", index=False)
