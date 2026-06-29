"""
Downloads page, clean, practical database download interface.
Presents CoreMet as reusable infrastructure with versioned download cards.

Route: /downloads
"""

import json
from pathlib import Path

from dash import html, dcc
import dash_bootstrap_components as dbc

from components.page_header import make_page_header


def _load_counts():
    """Record counts from the canonical stats JSON (single source of truth)."""
    fallback = {"total": "1,952,688", "MPI": "38,061", "MEI": "47,551", "MDI": "82,882",
                "MMI": "77,605", "MDrI": "3,500", "MGI": "1,658,745", "mGWAS": "44,344"}
    try:
        s = json.loads((Path(__file__).parent.parent / "data" / "coremetdb_stats.json").read_text())
        c = {k: f"{v['interactions']:,}" for k, v in s["databases"].items()}
        c["total"] = f"{s['totals']['interactions']:,}"
        return c
    except Exception:
        return fallback


_C = _load_counts()


def _download_card(name, description, fmt, records, version, href, color="#1a365d", disabled=False):
    """Single download card with metadata."""
    return dbc.Col(
        html.Div([
            html.Div([
                html.Span(style={
                    "width": "8px", "height": "8px", "borderRadius": "50%",
                    "background": color, "display": "inline-block", "marginRight": "8px",
                }),
                html.Span(name, style={"fontWeight": "600", "fontSize": "0.95rem", "color": "#1a202c"}),
            ], style={"marginBottom": "8px"}),
            html.P(description, style={"fontSize": "0.8rem", "color": "#718096", "marginBottom": "12px", "lineHeight": "1.5"}),
            html.Div([
                html.Div([
                    html.Span("Format: ", style={"color": "#a0aec0", "fontSize": "0.75rem"}),
                    html.Span(fmt, style={"fontSize": "0.75rem", "fontWeight": "500"}),
                ]),
                html.Div([
                    html.Span("Records: ", style={"color": "#a0aec0", "fontSize": "0.75rem"}),
                    html.Span(records, style={"fontSize": "0.75rem", "fontWeight": "500"}),
                ]),
                html.Div([
                    html.Span("Version: ", style={"color": "#a0aec0", "fontSize": "0.75rem"}),
                    html.Span(version, style={"fontSize": "0.75rem", "fontWeight": "500"}),
                ]),
            ], style={"display": "flex", "gap": "16px", "flexWrap": "wrap", "marginBottom": "12px"}),
            dbc.Button([
                html.I(className="fas fa-download me-2"),
                "Download" if not disabled else "Optional",
            ], href=None if disabled else href, size="sm", outline=True, color="primary",
               disabled=disabled, style={"fontSize": "0.8rem"}),
        ], className="cm-card", style={"padding": "20px", "height": "100%"}),
        md=6, lg=4, className="mb-3",
    )


layout = html.Div([
    html.Div([
        make_page_header(
            title="Downloads",
            subtitle="Download CoreMet datasets for local analysis.",
            breadcrumb_items=[("Home", "/home"), ("Downloads", None)],
        ),

        # ── License & update info ──
        html.Div([
            html.P([
                "All CoreMet datasets are released under the ",
                html.A("CC BY 4.0", href="https://creativecommons.org/licenses/by/4.0/",
                       target="_blank", style={"fontWeight": "600"}),
                " license. Data is updated quarterly. Last update: ",
                html.Strong("March 2026"),
                ".",
            ], style={"fontSize": "0.9rem", "color": "#4a5568", "marginBottom": "4px"}),
            html.P(
                "When using these datasets, please cite the CoreMet publication.",
                style={"fontSize": "0.8rem", "color": "#a0aec0"},
            ),
        ], style={
            "background": "#f8f9fa", "borderRadius": "8px", "padding": "16px 20px",
            "border": "1px solid #e2e8f0", "marginBottom": "32px",
        }),

        # ── Full database downloads ──
        html.H3("Full database", style={
            "fontSize": "1.15rem", "fontWeight": "600", "color": "#1a202c",
            "marginBottom": "16px",
        }),

        dbc.Row([
            _download_card(
                "Complete Edge List", "All curated, deduplicated interactions across all seven layers.",
                "CSV", _C["total"], "v3.0",
                "/api/v1/download/full-edges", "#1a365d",
            ),
            _download_card(
                "Node Metadata", "All entities with identifiers, aliases, and annotations.",
                "JSON", "~50,000", "v3.0",
                "/api/v1/download/node-metadata", "#1a365d",
            ),
            _download_card(
                "Schema Documentation", "Data dictionary and column definitions.",
                "Markdown", "-", "v3.0",
                "/api/v1/download/schema", "#1a365d",
            ),
        ], className="g-3 mb-4"),

        # ── Per-layer downloads ──
        html.H3("By interaction layer", style={
            "fontSize": "1.15rem", "fontWeight": "600", "color": "#1a202c",
            "marginBottom": "16px",
        }),

        dbc.Row([
            _download_card(
                "MPI: Metabolite–Protein", "Curated metabolite–protein binding/catalysis interactions.",
                "CSV", _C["MPI"], "v3.0",
                "/api/v1/download/mpi", "#3182ce",
            ),
            _download_card(
                "MEI: Metabolite–Enzyme", "Metabolite–enzyme (EC) catalytic relationships.",
                "CSV", _C["MEI"], "v3.0",
                "/api/v1/download/mei", "#00a3c4",
            ),
            _download_card(
                "MDI: Metabolite–Disease", "Metabolite–disease phenotype associations.",
                "CSV", _C["MDI"], "v3.0",
                "/api/v1/download/mdi", "#e53e3e",
            ),
            _download_card(
                "MMI: Metabolite–Microbe", "Gut microbe metabolite production/consumption.",
                "CSV", _C["MMI"], "v3.0",
                "/api/v1/download/mmi", "#38a169",
            ),
            _download_card(
                "MDrI: Metabolite–Drug", "Drug–metabolite pharmacokinetic interactions.",
                "CSV", _C["MDrI"], "v3.0",
                "/api/v1/download/mdri", "#805ad5",
            ),
            _download_card(
                "MGI: Metabolite–Gene", "Chemical–gene interactions from CTD.",
                "CSV", _C["MGI"], "v3.0",
                "/api/v1/download/mgi", "#d69e2e",
            ),
            _download_card(
                "mGWAS: Metabolite–SNP", "GWAS-derived SNP–metabolite links.",
                "CSV", _C["mGWAS"], "v3.0",
                "/api/v1/download/mgwas", "#319795",
            ),
        ], className="g-3 mb-4"),

        # ── ML / Embeddings ──
        html.H3("Machine learning resources", style={
            "fontSize": "1.15rem", "fontWeight": "600", "color": "#1a202c",
            "marginBottom": "16px",
        }),

        dbc.Row([
            _download_card(
                "Node Embeddings", "128-dim GraphSAGE embeddings for all entities.",
                "NPZ", "~50,000", "v3.0",
                "/api/v1/download/embeddings", "#718096", disabled=True,
            ),
            _download_card(
                "Pre-trained Model", "GraphSAGE model weights for transfer learning.",
                "PT", "-", "v3.0",
                "/api/v1/download/model", "#718096", disabled=True,
            ),
        ], className="g-3 mb-5"),

    ], className="cm-page-container"),
])

page_content = layout
