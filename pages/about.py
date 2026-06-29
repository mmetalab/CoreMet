"""
About page, contact information, team, citation, and licensing.
"""

from dash import html, dcc
import dash_bootstrap_components as dbc

from components.page_header import make_page_header


# ---------------------------------------------------------------------------
# Helper: contact / info card
# ---------------------------------------------------------------------------

def _info_card(icon, title, children):
    return html.Div([
        html.Div([
            html.Div(html.I(className=f"fas {icon}"), className="coverage-icon icon-mpi",
                     style={"marginRight": "12px"}),
            html.H5(title, style={"fontWeight": "600", "marginBottom": "0"}),
        ], style={"display": "flex", "alignItems": "center", "marginBottom": "16px"}),
        *children,
    ], className="cm-card", style={"marginBottom": "20px"})


# ---------------------------------------------------------------------------
# LAYOUT
# ---------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        make_page_header(
            title="About",
            subtitle="Contact information, team, citation, and data licensing.",
            breadcrumb_items=[("Home", "/home"), ("About", None)],
        ),

        # ── Overview ───────────────────────────────────────────
        _info_card("fa-info-circle", "About CoreMet", [
            html.P([
                "CoreMet is a comprehensive web platform for exploring metabolite–protein, "
                "metabolite–disease, metabolite–microbe, and metabolite–drug interactions. "
                "It integrates curated data from multiple public databases and provides "
                "GraphSAGE-based machine learning predictions for novel interaction discovery."
            ], style={"fontSize": "0.9rem", "color": "#718096", "lineHeight": "1.65"}),
            html.P([
                "The platform currently hosts ",
                html.Strong("377,164"),
                " curated interactions spanning ",
                html.Strong("10 organisms"),
                ", ",
                html.Strong("130 diseases"),
                " across ",
                html.Strong("22 categories"),
                ", ",
                html.Strong("3,500 drug interactions"),
                ", and ",
                html.Strong("100+ tissues"),
                "."
            ], style={"fontSize": "0.9rem", "color": "#718096", "lineHeight": "1.65"}),
        ]),

        # ── Contact ────────────────────────────────────────────
        _info_card("fa-envelope", "Contact", [
            dbc.Row([
                dbc.Col([
                    html.Div([
                        html.H6("Cheng Wang, Ph.D.", style={"fontWeight": "600", "marginBottom": "4px"}),
                        html.P("Developer & Lead Researcher", className="text-muted",
                               style={"fontSize": "0.85rem", "marginBottom": "4px"}),
                        html.P([
                            html.I(className="fas fa-envelope me-2",
                                   style={"color": "#38b2ac", "fontSize": "0.85rem"}),
                            html.A("cheng.wang@example.edu",
                                   href="mailto:cheng.wang@example.edu",
                                   style={"fontSize": "0.85rem", "color": "#1a365d"}),
                        ], className="mb-1"),
                    ], style={
                        "background": "rgba(26,54,93,0.03)", "borderRadius": "8px",
                        "padding": "16px", "border": "1px solid rgba(26,54,93,0.08)",
                    }),
                ], md=6),
                dbc.Col([
                    html.Div([
                        html.H6("Lab PI / Corresponding Author", style={"fontWeight": "600", "marginBottom": "4px"}),
                        html.P("Faculty of Health Sciences", className="text-muted",
                               style={"fontSize": "0.85rem", "marginBottom": "4px"}),
                        html.P([
                            html.I(className="fas fa-envelope me-2",
                                   style={"color": "#38b2ac", "fontSize": "0.85rem"}),
                            html.A("pi@example.edu",
                                   href="mailto:pi@example.edu",
                                   style={"fontSize": "0.85rem", "color": "#1a365d"}),
                        ], className="mb-1"),
                    ], style={
                        "background": "rgba(26,54,93,0.03)", "borderRadius": "8px",
                        "padding": "16px", "border": "1px solid rgba(26,54,93,0.08)",
                    }),
                ], md=6),
            ]),
        ]),

        # ── Data Sources ──────────────────────────────────────
        _info_card("fa-database", "Data Sources", [
            html.P("CoreMet integrates data curated from the following public databases:",
                   style={"fontSize": "0.9rem", "color": "#718096", "marginBottom": "12px"}),
            dbc.Row([
                dbc.Col([
                    html.H6("Metabolite–Protein (MPI & MEI)", style={"fontWeight": "600", "fontSize": "0.85rem"}),
                    html.Ul([
                        html.Li(html.A("HMDB", href="https://hmdb.ca", target="_blank")),
                        html.Li(html.A("UniProt", href="https://www.uniprot.org", target="_blank")),
                        html.Li(html.A("STITCH", href="http://stitch.embl.de", target="_blank")),
                        html.Li(html.A("KEGG", href="https://www.genome.jp/kegg/", target="_blank")),
                        html.Li(html.A("Rhea", href="https://www.rhea-db.org", target="_blank")),
                    ], style={"fontSize": "0.85rem"}),
                ], md=3),
                dbc.Col([
                    html.H6("Metabolite–Disease (MDI)", style={"fontWeight": "600", "fontSize": "0.85rem"}),
                    html.Ul([
                        html.Li(html.A("HMDB", href="https://hmdb.ca", target="_blank")),
                        html.Li(html.A("DisGeNET", href="https://www.disgenet.org", target="_blank")),
                        html.Li(html.A("CTD", href="http://ctdbase.org", target="_blank")),
                    ], style={"fontSize": "0.85rem"}),
                ], md=3),
                dbc.Col([
                    html.H6("Metabolite–Microbe (MMI)", style={"fontWeight": "600", "fontSize": "0.85rem"}),
                    html.Ul([
                        html.Li(html.A("gutMGene", href="http://bio-annotation.cn/gutmgene/", target="_blank")),
                        html.Li(html.A("KEGG", href="https://www.genome.jp/kegg/", target="_blank")),
                        html.Li(html.A("HMDB", href="https://hmdb.ca", target="_blank")),
                    ], style={"fontSize": "0.85rem"}),
                ], md=3),
                dbc.Col([
                    html.H6("Metabolite–Drug (MDrI)", style={"fontWeight": "600", "fontSize": "0.85rem"}),
                    html.Ul([
                        html.Li(html.A("DrugBank", href="https://go.drugbank.com", target="_blank")),
                        html.Li(html.A("HMDB", href="https://hmdb.ca", target="_blank")),
                        html.Li(html.A("CTD", href="http://ctdbase.org", target="_blank")),
                        html.Li(html.A("PubMed", href="https://pubmed.ncbi.nlm.nih.gov", target="_blank")),
                    ], style={"fontSize": "0.85rem"}),
                ], md=3),
            ]),
        ]),

        # ── Source Code / Repository ──────────────────────────
        _info_card("fa-code-branch", "Source Code & Repository", [
            html.P("CoreMet is open-source software. The code and curated datasets are available on GitHub:",
                   style={"fontSize": "0.9rem", "color": "#718096", "marginBottom": "12px"}),
            dbc.Row([
                dbc.Col([
                    html.A([
                        html.I(className="fab fa-github me-2", style={"fontSize": "1.2rem"}),
                        "CoreMet Web Platform",
                    ], href="https://github.com/cormet/cormet-web", target="_blank",
                       className="btn",
                       style={
                           "background": "#1a365d", "color": "#fff", "borderRadius": "8px",
                           "padding": "10px 20px", "fontSize": "0.9rem", "fontWeight": "500",
                           "textDecoration": "none", "display": "inline-flex", "alignItems": "center",
                       }),
                ], width="auto"),
                dbc.Col([
                    html.A([
                        html.I(className="fas fa-book me-2", style={"fontSize": "1.2rem"}),
                        "Documentation",
                    ], href="/documentation",
                       className="btn",
                       style={
                           "background": "#38b2ac", "color": "#fff", "borderRadius": "8px",
                           "padding": "10px 20px", "fontSize": "0.9rem", "fontWeight": "500",
                           "textDecoration": "none", "display": "inline-flex", "alignItems": "center",
                       }),
                ], width="auto"),
            ], className="g-2"),
        ]),

        # ── License / Disclaimer ──────────────────────────────
        _info_card("fa-balance-scale", "License & Disclaimer", [
            html.P([
                "CoreMet is released under the ",
                html.A("MIT License", href="https://opensource.org/licenses/MIT", target="_blank"),
                ". The curated interaction data is provided for academic research purposes. "
                "We make no warranties regarding the completeness or accuracy of the data. "
                "Users should independently verify critical findings."
            ], style={"fontSize": "0.85rem", "color": "#718096", "lineHeight": "1.65"}),
            html.P([
                "Medical disclaimer: Information provided on this platform is not intended "
                "to be a substitute for professional medical advice, diagnosis, or treatment."
            ], style={"fontSize": "0.85rem", "color": "#718096", "lineHeight": "1.65", "marginBottom": "0"}),
        ]),

    ], className="cm-page-container"),
])

page_content = layout
