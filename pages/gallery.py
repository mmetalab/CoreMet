"""
Gallery page, curated case studies demonstrating CoreMet workflows.

Shows real results from three disease case studies:
  • Hepatocellular Carcinoma (HCC), GraphSAGE-based prediction
  • Papillary Thyroid Cancer, pathway enrichment focus
  • Alzheimer's Disease, neurodegenerative MPI network
"""

import json
from pathlib import Path

from dash import html, dcc
import dash_bootstrap_components as dbc
import plotly.graph_objects as go

from components.page_header import make_page_header

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DISEASE_DIR = PROJECT_ROOT / "data" / "mpidatabase" / "disease_mpi"
CASE_DIR = PROJECT_ROOT / "data" / "case_studies"

# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------
_PAL = {
    "navy":   "#1a365d",
    "teal":   "#38b2ac",
    "red":    "#e53e3e",
    "blue":   "#3182ce",
    "purple": "#805ad5",
    "orange": "#dd6b20",
    "green":  "#38a169",
    "grey":   "#718096",
}

_CARD_STYLE = {
    "background": "#fff",
    "borderRadius": "12px",
    "border": "1px solid rgba(26,54,93,0.08)",
    "boxShadow": "0 1px 3px rgba(26,54,93,0.06)",
    "padding": "24px",
    "marginBottom": "24px",
}


# ---------------------------------------------------------------------------
# Helper, load live disease panel data from disk
# ---------------------------------------------------------------------------

def _load_panel(disease_key):
    """Load network_stats.json from a disease panel folder."""
    path = DISEASE_DIR / disease_key / "network_stats.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


def _load_case(disease_key):
    """Load results.json from case_studies if available."""
    path = CASE_DIR / disease_key / "results.json"
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


# ---------------------------------------------------------------------------
# Stat badge helper
# ---------------------------------------------------------------------------

def _stat_badge(icon, value, label, color):
    return html.Div([
        html.I(className=f"fas {icon}", style={"color": color, "fontSize": "1.2rem"}),
        html.Div([
            html.Div(str(value), style={"fontWeight": "700", "fontSize": "1.1rem", "color": _PAL["navy"]}),
            html.Div(label, style={"fontSize": "0.75rem", "color": _PAL["grey"]}),
        ], style={"marginLeft": "10px"}),
    ], style={"display": "flex", "alignItems": "center", "minWidth": "120px"})


# ---------------------------------------------------------------------------
# Case study card builder
# ---------------------------------------------------------------------------

def _case_study_card(
    idx, title, disease_key, category, tissue,
    description, highlights, key_findings, color,
):
    panel = _load_panel(disease_key)
    case = _load_case(disease_key)

    # ---------- stat badges ----------
    if panel:
        n_met = panel.get("n_metabolites", "-")
        n_prot = panel.get("n_proteins", "-")
        n_edges = panel.get("n_edges", "-")
        density = panel.get("density", 0)
    else:
        n_met = n_prot = n_edges = "-"
        density = 0

    stats_row = html.Div([
        _stat_badge("fa-flask", n_met, "Metabolites", color),
        _stat_badge("fa-dna", n_prot, "Proteins", color),
        _stat_badge("fa-project-diagram", n_edges, "Edges", color),
        _stat_badge("fa-chart-line", f"{density:.2f}" if isinstance(density, float) else density, "Density", color),
    ], style={"display": "flex", "flexWrap": "wrap", "gap": "20px", "marginBottom": "16px"})

    # ---------- hub network mini-chart ----------
    fig_hub = None
    if panel:
        met_hubs = panel.get("metabolite_hubs", [])[:5]
        prot_hubs = panel.get("protein_hubs", [])[:5]
        all_hubs = (
            [(h["node"], h["degree"], "Metabolite") for h in met_hubs]
            + [(h["node"], h["degree"], "Protein") for h in prot_hubs]
        )
        all_hubs.sort(key=lambda x: x[1], reverse=True)
        top10 = all_hubs[:10]
        if top10:
            names = [h[0][:25] + ("…" if len(h[0]) > 25 else "") for h in top10]
            degrees = [h[1] for h in top10]
            types = [h[2] for h in top10]
            bar_colors = [_PAL["teal"] if t == "Metabolite" else _PAL["blue"] for t in types]
            fig_hub = go.Figure(go.Bar(
                y=names[::-1], x=degrees[::-1],
                orientation="h",
                marker_color=bar_colors[::-1],
                hovertemplate="<b>%{y}</b><br>Degree: %{x}<extra></extra>",
            ))
            fig_hub.update_layout(
                font=dict(family="Arial, Helvetica, sans-serif"),
                height=260, margin=dict(l=10, r=10, t=30, b=10),
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                xaxis=dict(title="Degree", gridcolor="rgba(0,0,0,0.06)"),
                yaxis=dict(automargin=True),
                title=dict(text="Top Hub Nodes", font=dict(size=12, color=_PAL["grey"])),
            )

    # ---------- enrichment mini-chart ----------
    fig_enrich = None
    if case:
        enrich = case.get("enrichment", [])
        # filter to named pathways only
        named = [e for e in enrich if not e["Pathway_Name"].startswith("ec") and not e["Pathway_Name"].startswith("map0")]
        named.sort(key=lambda e: e["Fold_Enrichment"], reverse=True)
        top_pw = named[:8]
        if top_pw:
            pw_names = [p["Pathway_Name"][:30] + ("…" if len(p["Pathway_Name"]) > 30 else "") for p in top_pw]
            pw_fe = [p["Fold_Enrichment"] for p in top_pw]
            pw_pvals = [p["P_value"] for p in top_pw]
            fig_enrich = go.Figure(go.Bar(
                y=pw_names[::-1], x=pw_fe[::-1],
                orientation="h",
                marker_color=[_PAL["purple"]] * len(pw_names),
                hovertemplate="<b>%{y}</b><br>Fold enrichment: %{x:.1f}<extra></extra>",
            ))
            fig_enrich.update_layout(
                font=dict(family="Arial, Helvetica, sans-serif"),
                height=240, margin=dict(l=10, r=10, t=30, b=10),
                paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                xaxis=dict(title="Fold Enrichment", gridcolor="rgba(0,0,0,0.06)"),
                yaxis=dict(automargin=True),
                title=dict(text="Top Enriched Pathways", font=dict(size=12, color=_PAL["grey"])),
            )

    # ---------- assemble card ----------
    header = html.Div([
        html.Div([
            html.Span(f"Case Study {idx}", style={
                "background": color, "color": "#fff", "padding": "3px 10px",
                "borderRadius": "20px", "fontSize": "0.75rem", "fontWeight": "600",
                "marginRight": "10px",
            }),
            html.Span(category, style={
                "background": "rgba(26,54,93,0.06)", "color": _PAL["navy"],
                "padding": "3px 10px", "borderRadius": "20px",
                "fontSize": "0.75rem", "fontWeight": "500", "marginRight": "8px",
            }),
            html.Span(tissue, style={
                "background": "rgba(56,178,172,0.08)", "color": _PAL["teal"],
                "padding": "3px 10px", "borderRadius": "20px",
                "fontSize": "0.75rem", "fontWeight": "500",
            }),
        ], style={"marginBottom": "10px"}),
        html.H4(title, style={"fontWeight": "700", "color": _PAL["navy"], "marginBottom": "6px"}),
        html.P(description, style={"fontSize": "0.9rem", "color": _PAL["grey"], "lineHeight": "1.65"}),
    ])

    # highlights list
    hl_list = html.Ul([
        html.Li(h, style={"fontSize": "0.85rem", "color": "#2d3748", "marginBottom": "4px"})
        for h in highlights
    ], style={"paddingLeft": "18px", "marginBottom": "16px"})

    # findings callout
    findings_box = html.Div([
        html.Div([
            html.I(className="fas fa-lightbulb me-2", style={"color": _PAL["orange"]}),
            html.Strong("Key Findings", style={"color": _PAL["navy"]}),
        ], style={"marginBottom": "8px", "display": "flex", "alignItems": "center"}),
        html.Ul([
            html.Li(f, style={"fontSize": "0.85rem", "color": "#2d3748", "marginBottom": "4px"})
            for f in key_findings
        ], style={"paddingLeft": "18px", "marginBottom": "0"}),
    ], style={
        "background": "rgba(221,107,32,0.05)", "borderLeft": f"3px solid {_PAL['orange']}",
        "padding": "14px 18px", "borderRadius": "0 8px 8px 0", "marginBottom": "16px",
    })

    # charts row
    charts = []
    if fig_hub:
        charts.append(dbc.Col(dcc.Graph(figure=fig_hub, config={"displayModeBar": False}), md=6))
    if fig_enrich:
        charts.append(dbc.Col(dcc.Graph(figure=fig_enrich, config={"displayModeBar": False}), md=6))
    elif fig_hub:
        # fill second col with "Explore" CTA
        charts.append(dbc.Col(html.Div([
            html.I(className="fas fa-arrow-right me-2"),
            html.A("Explore this disease panel →", href="/disease",
                   style={"fontWeight": "600", "color": color, "fontSize": "0.9rem", "textDecoration": "none"}),
        ], style={"display": "flex", "alignItems": "center", "justifyContent": "center", "height": "100%"}), md=6))

    chart_row = dbc.Row(charts, className="g-3") if charts else html.Div()

    # CTA button
    cta = html.Div([
        html.A([
            html.I(className="fas fa-external-link-alt me-2"),
            f"Open {title.split('-')[0].strip()} in Disease Explorer",
        ], href="/disease", style={
            "background": color, "color": "#fff", "padding": "8px 20px",
            "borderRadius": "8px", "fontSize": "0.85rem", "fontWeight": "600",
            "textDecoration": "none", "display": "inline-flex", "alignItems": "center",
        }),
    ], style={"marginTop": "16px"})

    return html.Div([
        header,
        stats_row,
        hl_list,
        findings_box,
        chart_row,
        cta,
    ], style=_CARD_STYLE)


# ---------------------------------------------------------------------------
# Workflow summary card
# ---------------------------------------------------------------------------

def _workflow_card():
    steps = [
        ("fa-database", "Database Query", "Search 377K+ curated interactions across MPI, MEI, MDI, MMI, and MDrI modules."),
        ("fa-heartbeat", "Disease Panel", "Select from 130 pre-computed disease panels spanning 22 categories."),
        ("fa-brain", "ML Prediction", "Run GraphSAGE-based prediction for novel metabolite–protein interactions."),
        ("fa-chart-bar", "Enrichment", "Identify enriched KEGG pathways with Fisher's exact test (BH-corrected)."),
        ("fa-project-diagram", "Network Viz", "Visualise interaction networks with hub detection and community analysis."),
    ]
    cols = [
        dbc.Col(html.Div([
            html.Div(html.I(className=f"fas {icon}", style={"color": _PAL["teal"], "fontSize": "1.3rem"}),
                     style={"marginBottom": "8px"}),
            html.Div(title, style={"fontWeight": "600", "fontSize": "0.85rem", "color": _PAL["navy"], "marginBottom": "4px"}),
            html.Div(desc, style={"fontSize": "0.78rem", "color": _PAL["grey"], "lineHeight": "1.5"}),
        ], style={"textAlign": "center", "padding": "12px"}), md=True)
        for icon, title, desc in steps
    ]
    return html.Div([
        html.H5([
            html.I(className="fas fa-route me-2", style={"color": _PAL["teal"]}),
            "CoreMet Discovery Workflow",
        ], style={"fontWeight": "700", "color": _PAL["navy"], "marginBottom": "16px"}),
        html.P(
            "Each case study below follows this general workflow, from curated database query "
            "through prediction, enrichment, and network analysis.",
            style={"fontSize": "0.9rem", "color": _PAL["grey"], "marginBottom": "16px", "lineHeight": "1.65"},
        ),
        dbc.Row(cols, className="g-2"),
    ], style=_CARD_STYLE)


# ---------------------------------------------------------------------------
# Cross-disease comparison card
# ---------------------------------------------------------------------------

def _cross_disease_card():
    path = CASE_DIR / "cross_disease_comparison.json"
    if not path.exists():
        return html.Div()
    with open(path) as f:
        cmp = json.load(f)

    diseases = cmp.get("diseases", [])
    shared_met = cmp.get("shared_hub_metabolites", [])
    shared_prot = cmp.get("shared_hub_proteins", [])

    return html.Div([
        html.Div([
            html.Span("Bonus", style={
                "background": _PAL["orange"], "color": "#fff", "padding": "3px 10px",
                "borderRadius": "20px", "fontSize": "0.75rem", "fontWeight": "600",
                "marginRight": "10px",
            }),
            html.H5("Cross-Disease Comparison", style={
                "fontWeight": "700", "color": _PAL["navy"], "marginBottom": "0", "display": "inline",
            }),
        ], style={"marginBottom": "14px"}),
        html.P([
            f"Comparing hub nodes across {len(diseases)} diseases ",
            html.Span(f"({', '.join(d.replace('_', ' ').title() for d in diseases)})",
                      style={"fontStyle": "italic"}),
            " reveals the extent of metabolic overlap between distinct pathologies.",
        ], style={"fontSize": "0.9rem", "color": _PAL["grey"], "lineHeight": "1.65"}),
        dbc.Row([
            dbc.Col([
                html.Div([
                    html.I(className="fas fa-flask me-2", style={"color": _PAL["teal"]}),
                    html.Strong(f"{len(shared_met)} shared hub metabolites"),
                ], style={"marginBottom": "8px"}),
                html.P(
                    ", ".join(shared_met) if shared_met else "No shared hub metabolites, each disease has a unique metabolic fingerprint.",
                    style={"fontSize": "0.85rem", "color": "#2d3748"},
                ),
            ], md=6),
            dbc.Col([
                html.Div([
                    html.I(className="fas fa-dna me-2", style={"color": _PAL["blue"]}),
                    html.Strong(f"{len(shared_prot)} shared hub proteins"),
                ], style={"marginBottom": "8px"}),
                html.P(
                    ", ".join(shared_prot) if shared_prot else "No shared hub proteins, distinct molecular mechanisms across conditions.",
                    style={"fontSize": "0.85rem", "color": "#2d3748"},
                ),
            ], md=6),
        ]),
        html.Div([
            html.I(className="fas fa-lightbulb me-2", style={"color": _PAL["orange"]}),
            html.Span(
                "This unique-hub pattern supports disease-specific metabolomic signatures "
                "and validates the discriminatory power of the CoreMet prediction engine.",
                style={"fontSize": "0.85rem", "color": "#2d3748"},
            ),
        ], style={
            "background": "rgba(221,107,32,0.05)", "borderLeft": f"3px solid {_PAL['orange']}",
            "padding": "12px 16px", "borderRadius": "0 8px 8px 0", "marginTop": "12px",
        }),
    ], style=_CARD_STYLE)


# ---------------------------------------------------------------------------
# LAYOUT
# ---------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        make_page_header(
            title="Gallery",
            subtitle="Curated case studies demonstrating metabolomics discovery with CoreMet.",
            breadcrumb_items=[("Home", "/home"), ("Gallery", None)],
        ),

        # Workflow overview
        _workflow_card(),

        # ── Case Study 1: HCC ─────────────────────────────────
        _case_study_card(
            idx=1,
            title="Hepatocellular Carcinoma, MPI Prediction & Network Analysis",
            disease_key="hcc",
            category="Cancer",
            tissue="Liver",
            description=(
                "Hepatocellular carcinoma (HCC) is the most common primary liver malignancy. "
                "Using CoreMet's GraphSAGE-based prediction engine, we predicted 1,635 metabolite–protein "
                "interactions from 20 HCC-associated metabolites and 15 key liver-cancer proteins. "
                "Network analysis reveals central metabolic hubs driving the Warburg effect."
            ),
            highlights=[
                "20 metabolites × 15 proteins → 300 candidate pairs, 1,635 total GraphSAGE-scored predictions",
                "Network: 35 nodes, 123 edges, density 0.41, a tightly connected module",
                "D-Glucose (degree 10) and RAC-alpha Ser/Thr kinase (degree 13) are top hubs",
                "Pathway enrichment implicates bile acid biosynthesis and taurine metabolism",
            ],
            key_findings=[
                "L-Lactic acid ↔ PKM axis recapitulates the Warburg effect in HCC metabolism",
                "Glucose-6-phosphatase appears as a high-degree hub, consistent with disrupted gluconeogenesis",
                "Top enriched pathways, primary bile acid biosynthesis (FE = 25.4) and taurine metabolism (FE = 23.3), are established HCC biomarker routes",
            ],
            color=_PAL["red"],
        ),

        # ── Case Study 2: Thyroid cancer ──────────────────────
        _case_study_card(
            idx=2,
            title="Papillary Thyroid Cancer, Pathway Enrichment Focus",
            disease_key="thyroid_cancer",
            category="Cancer",
            tissue="Thyroid",
            description=(
                "Papillary thyroid cancer (PTC) is the most prevalent endocrine malignancy. "
                "CoreMet predicted 922 interactions across a large network of 607 nodes and 855 edges. "
                "Fisher's exact pathway enrichment (BH-corrected FDR < 0.05) identified amino acid "
                "and cysteine metabolism as top dysregulated pathways."
            ),
            highlights=[
                "922 significant predictions from the GraphSAGE model (prediction threshold ≥ 0.3)",
                "Large network: 607 nodes, 855 edges across thyroid-specific MPI pairs",
                "L-Glutamic acid is the dominant metabolite hub (degree 202)",
                "Multiple pathways reach genome-wide significance (FDR < 0.05)",
            ],
            key_findings=[
                "Alanine/aspartate/glutamate metabolism is the #1 enriched pathway (FE = 5.49, FDR = 4.6 × 10⁻¹⁶)",
                "Cysteine and methionine metabolism ranks #2 (FE = 4.49, FDR = 4.6 × 10⁻¹⁶), consistent with oxidative stress in thyroid tumors",
                "RAS pathway proteins (KRAS, BRAF) connect to amino acid metabolites, a known PTC driver axis",
            ],
            color=_PAL["purple"],
        ),

        # ── Case Study 3: Alzheimer's ────────────────────────
        _case_study_card(
            idx=3,
            title="Alzheimer's Disease, Neurodegenerative MPI Network",
            disease_key="alzheimers",
            category="Neurodegenerative",
            tissue="Brain",
            description=(
                "Alzheimer's disease (AD) is the leading cause of dementia worldwide. "
                "CoreMet's disease panel includes 18 AD-associated metabolites and 15 key "
                "neurodegeneration proteins. Network analysis reveals a dense interaction core "
                "(density 0.44) centred on neurotransmitter and amyloid-processing pathways."
            ),
            highlights=[
                "33 nodes (18 metabolites + 15 proteins), 119 edges, density 0.44",
                "Serotonin (degree 11) and Dopamine (degree 8) are top metabolite hubs, both neurotransmitters",
                "BDNF (degree 12) and Presenilin-1 (degree 12) are top protein hubs",
                "APP (Amyloid Precursor Protein) connects to 10 metabolites",
            ],
            key_findings=[
                "The serotonin ↔ BDNF hub link matches the serotonin–neurotrophin hypothesis of AD",
                "Presenilin-1 and APP form a densely connected sub-module, recapitulating the amyloid cascade",
                "Uric acid (degree 10) emerges as a protective metabolite hub, consistent with its known antioxidant role in neurodegeneration",
            ],
            color=_PAL["blue"],
        ),

        # ── Explore CTA ───────────────────────────────────────
        html.Div([
            html.H5("Start Your Own Discovery", style={"fontWeight": "700", "color": _PAL["navy"], "marginBottom": "8px"}),
            html.P(
                "CoreMet provides 130 pre-computed disease panels, GraphSAGE-based interaction prediction, "
                "and pathway enrichment, ready for your metabolomics research.",
                style={"fontSize": "0.9rem", "color": _PAL["grey"], "marginBottom": "16px"},
            ),
            html.Div([
                html.A([html.I(className="fas fa-heartbeat me-2"), "Disease Explorer"],
                       href="/disease", className="btn", style={
                           "background": _PAL["red"], "color": "#fff", "borderRadius": "8px",
                           "padding": "10px 20px", "fontWeight": "600", "fontSize": "0.85rem",
                           "textDecoration": "none", "marginRight": "12px",
                       }),
                html.A([html.I(className="fas fa-brain me-2"), "Predict Interactions"],
                       href="/predict", className="btn", style={
                           "background": _PAL["teal"], "color": "#fff", "borderRadius": "8px",
                           "padding": "10px 20px", "fontWeight": "600", "fontSize": "0.85rem",
                           "textDecoration": "none", "marginRight": "12px",
                       }),
                html.A([html.I(className="fas fa-database me-2"), "Browse Database"],
                       href="/database", className="btn", style={
                           "background": _PAL["blue"], "color": "#fff", "borderRadius": "8px",
                           "padding": "10px 20px", "fontWeight": "600", "fontSize": "0.85rem",
                           "textDecoration": "none",
                       }),
            ], style={"display": "flex", "flexWrap": "wrap", "gap": "8px"}),
        ], style={**_CARD_STYLE, "textAlign": "center"}),

    ], className="cm-page-container"),
])

page_content = layout
