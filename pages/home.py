"""
Landing page, clean, journal-like homepage for CoreMet.
Hero + search bar + stats + "What you can do" + architecture overview.
"""

from functools import lru_cache
from urllib.parse import quote_plus

from dash import html, Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc


# ---------------------------------------------------------------------------
# Compute live database stats at import time
# ---------------------------------------------------------------------------
def _compute_stats():
    """Load canonical database statistics from data/coremetdb_stats.json.

    This is the single source of truth produced by scripts/compute_db_stats.py;
    the manuscript, the publication figures, and this page all read the same
    numbers. Hard-coded values below are only a fallback if the JSON is missing.
    """
    import json
    from pathlib import Path

    _data_dir = Path(__file__).parent.parent / "data"
    stats = {
        "mpi_total": 38061, "mei_total": 47551, "mdi_total": 82882, "mmi_total": 77605,
        "mdri_total": 3500, "mgi_total": 1658745, "mgwas_total": 44344,
        "metabolites": 30674, "proteins": 22257, "ecs": 2406, "diseases": 3595,
        "microbes": 1262, "drugs": 2162, "snps": 33623, "genes": 50164,
        "organisms": 10, "pathways": 296, "total_all": 1952688,
    }
    try:
        s = json.loads((_data_dir / "coremetdb_stats.json").read_text())
        db = s["databases"]
        stats.update({
            "mpi_total": db["MPI"]["interactions"], "mei_total": db["MEI"]["interactions"],
            "mdi_total": db["MDI"]["interactions"], "mmi_total": db["MMI"]["interactions"],
            "mdri_total": db["MDrI"]["interactions"], "mgi_total": db["MGI"]["interactions"],
            "mgwas_total": db["mGWAS"]["interactions"],
            "metabolites": s["totals"]["unique_metabolites_union"],
            "proteins": db["MPI"]["targets"], "ecs": db["MEI"]["targets"],
            "diseases": db["MDI"]["targets"], "microbes": db["MMI"]["targets"],
            "drugs": db["MDrI"]["targets"], "snps": db["mGWAS"]["targets"],
            "genes": db["MGI"]["targets"],
            "organisms": db["MPI"].get("organisms", 10),
            "pathways": db["MPI"].get("pathways", 296),
            "total_all": s["totals"]["interactions"],
        })
    except Exception:
        pass
    return stats


_STATS = _compute_stats()


# ---------------------------------------------------------------------------
# Entity color mapping (consistent everywhere)
# ---------------------------------------------------------------------------
ENTITY_COLORS = {
    "metabolite": "#e27a3f",
    "protein":    "#3182ce",
    "gene":       "#d69e2e",
    "disease":    "#e53e3e",
    "microbe":    "#38a169",
    "drug":       "#805ad5",
    "snp":        "#319795",
}


# ---------------------------------------------------------------------------
# Helper: stat card
# ---------------------------------------------------------------------------
def _stat_card(value, label, color, icon):
    return html.Div([
        html.Div(html.I(className=icon), style={
            "fontSize": "1.3rem", "color": color, "marginBottom": "4px",
        }),
        html.Div(f"{value:,}", style={
            "fontSize": "1.5rem", "fontWeight": "700", "color": "#1a202c", "lineHeight": "1.2",
        }),
        html.Div(label, style={
            "fontSize": "0.75rem", "color": "#718096", "fontWeight": "500",
        }),
    ], style={
        "textAlign": "center", "padding": "20px 12px", "flex": "1 1 0",
        "minWidth": "120px",
    })


# ---------------------------------------------------------------------------
# Helper: action card
# ---------------------------------------------------------------------------
def _action_card(icon, title, desc, href, color, disabled=False):
    opacity = "0.5" if disabled else "1"
    cursor = "not-allowed" if disabled else "pointer"
    return dbc.Col(
        html.A(
            html.Div([
                html.Div(html.I(className=icon), style={
                    "width": "48px", "height": "48px", "borderRadius": "10px",
                    "display": "flex", "alignItems": "center", "justifyContent": "center",
                    "fontSize": "1.25rem", "color": color,
                    "background": f"{color}14", "marginBottom": "16px",
                }),
                html.H4(title, style={
                    "fontSize": "1rem", "fontWeight": "600", "color": "#1a202c",
                    "marginBottom": "8px",
                }),
                html.P(desc, style={
                    "fontSize": "0.85rem", "color": "#718096", "margin": "0",
                    "lineHeight": "1.55",
                }),
            ], className="cm-card", style={
                "height": "100%", "padding": "28px 24px", "textAlign": "left",
                "opacity": opacity,
            }),
            href=href if not disabled else None,
            style={"textDecoration": "none", "cursor": cursor},
        ),
        md=4, className="mb-3",
    )


# ---------------------------------------------------------------------------
# Helper: layer chip
# ---------------------------------------------------------------------------
def _layer_chip(label, count, color, href):
    return html.A(
        html.Div([
            html.Span(style={
                "width": "10px", "height": "10px", "borderRadius": "50%",
                "background": color, "display": "inline-block", "marginRight": "8px",
                "flexShrink": "0",
            }),
            html.Span(label, style={"fontWeight": "600", "fontSize": "0.85rem", "marginRight": "8px"}),
            html.Span(f"{count:,}", style={
                "fontSize": "0.8rem", "color": "#718096", "marginLeft": "auto",
            }),
        ], style={
            "display": "flex", "alignItems": "center", "padding": "10px 16px",
            "background": "#fff", "border": "1px solid #e2e8f0", "borderRadius": "8px",
            "transition": "border-color 0.15s",
        }),
        href=href, style={"textDecoration": "none"},
    )


@lru_cache(maxsize=1)
def _load_home_metabolite_index():
    """Small cached metabolite index for homepage suggestions."""
    try:
        from app.config import Config
        import pandas as pd

        cfg = Config()
        df = pd.read_csv(cfg.MPI_DB_PATH, usecols=["HMDB ID", "Metabolite Name"], dtype=str)
        df = df.fillna("").drop_duplicates(subset=["HMDB ID"])
        df = df[df["Metabolite Name"].str.strip().ne("")]
        return df[["HMDB ID", "Metabolite Name"]]
    except Exception:
        return None


def _home_search_href(query: str) -> str:
    return f"/search?q={quote_plus(query.strip())}"


# --------------------------------------------------------------------------
# LAYOUT
# --------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        # ══════════════════════════════════════════════════════════
        # HERO SECTION, title + subtitle + search bar + examples
        # ══════════════════════════════════════════════════════════
        html.Div([
            html.H1("CoreMet", style={
                "fontSize": "2.8rem", "fontWeight": "700", "color": "#1a202c",
                "marginBottom": "12px", "letterSpacing": "-0.02em",
            }),
            html.P(
                "An integrated metabolite-centered knowledge graph for biomedical and "
                "disease-focused research, connecting metabolites to genes, proteins, "
                "diseases, gut microbiota, drugs, and genetic variants.",
                style={
                    "fontSize": "1.1rem", "color": "#4a5568", "maxWidth": "700px",
                    "margin": "0 auto 32px", "lineHeight": "1.6",
                },
            ),

            # ── Central search bar ──
            html.Div([
                dbc.InputGroup([
                    dbc.Input(
                        id="home-search-input",
                        placeholder="Search a metabolite, disease, gene, protein, drug, microbe, or SNP…",
                        type="text",
                        debounce=False,
                        className="home-search-input",
                    ),
                    dbc.Button(
                        html.I(className="fas fa-search"),
                        id="home-search-btn",
                        color="primary",
                        className="home-search-button",
                        n_clicks=0,
                    ),
                ], className="home-search-control"),
                html.Div(id="home-search-suggestions", className="home-search-suggestions"),
            ], className="home-search-wrap"),

            # ── Example prompts ──
            html.Div([
                html.Span("Try: ", style={"fontSize": "0.8rem", "color": "#a0aec0", "marginRight": "8px"}),
                html.A("butyrate", href="/search?q=butyrate", className="cm-example-chip"),
                html.A("metformin", href="/search?q=metformin", className="cm-example-chip"),
                html.A("rs1260326", href="/search?q=rs1260326", className="cm-example-chip"),
            ], style={"display": "flex", "justifyContent": "center", "alignItems": "center", "gap": "4px"}),

        ], style={
            "textAlign": "center", "padding": "56px 24px 40px",
        }),

        # ══════════════════════════════════════════════════════════
        # STATS BAR
        # ══════════════════════════════════════════════════════════
        html.Div([
            html.Div([
                _stat_card(_STATS["total_all"], "Total Interactions", "#1a365d", "fas fa-link"),
                _stat_card(_STATS["metabolites"], "Metabolites", ENTITY_COLORS["metabolite"], "fas fa-atom"),
                _stat_card(_STATS["genes"], "Genes", ENTITY_COLORS["gene"], "fas fa-dna"),
                _stat_card(_STATS["proteins"], "Proteins", ENTITY_COLORS["protein"], "fas fa-cubes"),
                _stat_card(_STATS["diseases"], "Diseases", ENTITY_COLORS["disease"], "fas fa-heartbeat"),
                _stat_card(_STATS["microbes"], "Microbes", ENTITY_COLORS["microbe"], "fas fa-bacterium"),
                _stat_card(_STATS["drugs"], "Drugs", ENTITY_COLORS["drug"], "fas fa-pills"),
                _stat_card(_STATS["snps"], "SNPs", ENTITY_COLORS["snp"], "fas fa-map-marker-alt"),
            ], style={
                "display": "flex", "flexWrap": "wrap", "justifyContent": "center",
                "borderTop": "1px solid #e2e8f0", "borderBottom": "1px solid #e2e8f0",
                "paddingTop": "8px", "paddingBottom": "8px",
            }),
        ], style={"marginBottom": "48px"}),

        # ══════════════════════════════════════════════════════════
        # WHAT YOU CAN DO, 3 action panels
        # ══════════════════════════════════════════════════════════
        html.Div([
            html.H2("What you can do", style={
                "fontSize": "1.5rem", "fontWeight": "700", "color": "#1a202c",
                "marginBottom": "8px", "textAlign": "center",
            }),
            html.P("Search, explore, and download the metabolite knowledge graph.", style={
                "fontSize": "0.9rem", "color": "#718096", "textAlign": "center",
                "marginBottom": "32px",
            }),
        ]),

        dbc.Row([
            _action_card(
                "fas fa-search", "Search entities",
                "Find any metabolite and see all connected diseases, genes, proteins, "
                "drugs, microbes, and genetic variants with evidence sources.",
                "/search", ENTITY_COLORS["metabolite"],
            ),
            _action_card(
                "fas fa-project-diagram", "Explore networks",
                "Discover multi-omics connections between metabolites and disease "
                "mechanisms through interactive cross-layer graph visualization.",
                "/explore", ENTITY_COLORS["protein"],
            ),
            _action_card(
                "fas fa-download", "Download & reuse",
                "Export any interaction layer with full provenance, or pull the whole "
                "knowledge graph via bulk CSV and a documented REST API for your pipelines.",
                "/downloads", ENTITY_COLORS["gene"],
            ),
        ], className="g-3 mb-5"),

        # ══════════════════════════════════════════════════════════
        # KNOWLEDGE GRAPH ARCHITECTURE, 7 interaction layers
        # ══════════════════════════════════════════════════════════
        html.Div([
            html.H2("Seven interaction layers", style={
                "fontSize": "1.5rem", "fontWeight": "700", "color": "#1a202c",
                "marginBottom": "8px", "textAlign": "center",
            }),
            html.P(
                "CoreMet integrates seven curated interaction layers into a unified "
                "knowledge graph for biomedical hypothesis generation, disease biomarker "
                "discovery, and multi-omics pathway analysis.",
                style={"fontSize": "0.9rem", "color": "#718096", "textAlign": "center",
                       "marginBottom": "32px", "maxWidth": "620px", "margin": "0 auto 32px"},
            ),
        ]),

        html.Div([
            dbc.Row([
                dbc.Col(_layer_chip("MPI: Metabolite–Protein", _STATS["mpi_total"],
                                    ENTITY_COLORS["protein"], "/mpi"), md=6, lg=4, className="mb-2"),
                dbc.Col(_layer_chip("MEI: Metabolite–Enzyme", _STATS["mei_total"],
                                    "#00a3c4", "/mei"), md=6, lg=4, className="mb-2"),
                dbc.Col(_layer_chip("MDI: Metabolite–Disease", _STATS["mdi_total"],
                                    ENTITY_COLORS["disease"], "/mdi"), md=6, lg=4, className="mb-2"),
                dbc.Col(_layer_chip("MMI: Metabolite–Microbe", _STATS["mmi_total"],
                                    ENTITY_COLORS["microbe"], "/mmi"), md=6, lg=4, className="mb-2"),
                dbc.Col(_layer_chip("MDrI: Metabolite–Drug", _STATS["mdri_total"],
                                    ENTITY_COLORS["drug"], "/mdri"), md=6, lg=4, className="mb-2"),
                dbc.Col(_layer_chip("MGI: Metabolite–Gene", _STATS["mgi_total"],
                                    ENTITY_COLORS["gene"], "/mgi"), md=6, lg=4, className="mb-2"),
                dbc.Col(_layer_chip("mGWAS: Metabolite–SNP", _STATS["mgwas_total"],
                                    ENTITY_COLORS["snp"], "/mgwas"), md=6, lg=4, className="mb-2"),
            ], className="g-2"),
        ], style={"maxWidth": "900px", "margin": "0 auto", "marginBottom": "48px"}),

    ], className="cm-page-container"),
])


# --------------------------------------------------------------------------
# CALLBACKS
# --------------------------------------------------------------------------

@callback(
    Output("home-search-suggestions", "children"),
    Input("home-search-input", "value"),
    prevent_initial_call=True,
)
def home_search_autocomplete(search_value):
    """Lightweight homepage suggestions; search still works for any free text."""
    if not search_value or len(search_value.strip()) < 2:
        return None

    query = search_value.strip()
    q = query.lower()
    rows = [
        html.A([
            html.I(className="fas fa-search me-2"),
            html.Span(f'Search all CoreMet for "{query}"'),
        ], href=_home_search_href(query), className="home-search-suggestion-row primary"),
    ]

    index = _load_home_metabolite_index()
    if index is not None and not index.empty:
        name_lower = index["Metabolite Name"].str.lower()
        hmdb_lower = index["HMDB ID"].str.lower()
        starts = index[name_lower.str.startswith(q, na=False) | hmdb_lower.str.startswith(q, na=False)]
        contains = index[
            (name_lower.str.contains(q, na=False, regex=False) | hmdb_lower.str.contains(q, na=False, regex=False))
            & ~index.index.isin(starts.index)
        ]
        match_indexes = list(starts.head(4).index)
        if len(match_indexes) < 4:
            match_indexes += list(contains.head(4 - len(match_indexes)).index)
        matches = index.loc[match_indexes]
        for _, row in matches.head(4).iterrows():
            name = str(row["Metabolite Name"])
            hmdb_id = str(row["HMDB ID"])
            rows.append(html.A([
                html.I(className="fas fa-atom me-2"),
                html.Span(name),
                html.Span(hmdb_id, className="home-search-suggestion-meta"),
            ], href=f"/metabolite?id={quote_plus(hmdb_id)}", className="home-search-suggestion-row"))

    return rows


@callback(
    Output("home-search-redirect", "href", allow_duplicate=True),
    Input("home-search-btn", "n_clicks"),
    Input("home-search-input", "n_submit"),
    State("home-search-input", "value"),
    prevent_initial_call=True,
)
def home_search_go(_n_clicks, _n_submit, value):
    """Redirect home search to the full search page for free-text queries."""
    if not value or len(value.strip()) < 2:
        return no_update
    return _home_search_href(value)


# Export page_content for use in app/main.py page routing
page_content = layout
