"""
Redesigned Database Page (v2), stat dashboard, sidebar filters, styled DataTable,
CSV export, and dynamic organism donut chart.

All callback IDs prefixed with 'dbv2-'.
"""

import io
import json
import math

import pandas as pd
import plotly.graph_objects as go
from dash import (
    dcc,
    html,
    dash_table,
    Input,
    Output,
    State,
    callback,
    no_update,
    ctx,
)
import dash_bootstrap_components as dbc

from components.page_header import make_page_header
from app.config import Config
from app.services.mei_service import get_mei_db
from app.services.mdi_service import get_mdi_db
from app.services.mmi_service import get_mmi_db
from app.services.mdri_service import get_mdri_db
from app.services.mgwas_service import get_mgwas_db
from app.services.mgi_service import get_mgi_db
from functools import lru_cache

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

_cfg = Config()
_df_raw = None  # Lazy-loaded on first access


def _get_df():
    """Lazy-load the database CSV on first access."""
    global _df_raw
    if _df_raw is not None:
        return _df_raw
    try:
        from pathlib import Path as _P
        _rel = _P(__file__).parent.parent / "data" / "databases" / "release" / "coremetdb_mpi.csv"
        from app.services.csv_loader import load_optimized
        _df_raw = load_optimized(_rel if _rel.exists() else _cfg.MPI_DB_PATH)
    except FileNotFoundError:
        _df_raw = pd.DataFrame(columns=[
            "Species", "Metabolite Name", "HMDB ID", "SMILES",
            "Uniprot ID", "Protein Name", "Gene Name",
            "Pathway_ID", "Pathway_Name", "Evidence_Source",
        ])
    # Ensure v2 columns
    for _col, _default in [("Evidence_Source", "original"), ("Pathway_ID", ""), ("Pathway_Name", "")]:
        if _col not in _df_raw.columns:
            _df_raw[_col] = _default
    # Pre-build concatenated search column for fast vectorized search
    _df_raw["_search_text"] = (
        _df_raw["Species"].astype("string").fillna("").str.lower() + " " +
        _df_raw["Metabolite Name"].astype("string").fillna("").str.lower() + " " +
        _df_raw["HMDB ID"].astype("string").fillna("").str.lower() + " " +
        _df_raw["Uniprot ID"].astype("string").fillna("").str.lower() + " " +
        _df_raw["Protein Name"].astype("string").fillna("").str.lower() + " " +
        _df_raw["Gene Name"].astype("string").fillna("").str.lower() + " " +
        _df_raw["Pathway_Name"].astype("string").fillna("").str.lower() + " " +
        _df_raw["Evidence_Source"].astype("string").fillna("").str.lower()
    )
    return _df_raw


# Placeholder for pre-computed values (filled on first callback)
_df_raw_placeholder = pd.DataFrame(columns=[
    "Species", "Metabolite Name", "HMDB ID", "SMILES",
    "Uniprot ID", "Protein Name", "Gene Name",
    "Pathway_ID", "Pathway_Name", "Evidence_Source",
])

_df_combined = None  # Cached MPI + MEI combined DataFrame


_MISSING_TOKENS = {"", "nan", "none", "null", "<na>"}


def _clean_series(series: pd.Series) -> pd.Series:
    """Return a stripped string series with common missing tokens as empty strings."""
    cleaned = series.astype("string").fillna("").str.strip()
    return cleaned.mask(cleaned.str.lower().isin(_MISSING_TOKENS), "")


def _valid_option_values(values):
    """Clean a value collection and remove empty/missing sidebar options."""
    cleaned = []
    seen = set()
    for value in values:
        text = str(value).strip()
        if text.lower() in _MISSING_TOKENS or not text or text == "-":
            continue
        if text not in seen:
            cleaned.append(text)
            seen.add(text)
    return sorted(cleaned)


def _fill_display_blanks(df: pd.DataFrame) -> pd.DataFrame:
    """Replace remaining display blanks with a visible placeholder."""
    df = df.copy()
    for col in df.columns:
        blank = _clean_series(df[col]).eq("")
        df.loc[blank, col] = "-"
    return df


def _valid_hmdb_mask(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.match(r"^HMDB\d+$", na=False)


def _first_real_value_map(df: pd.DataFrame, key_col: str, value_col: str) -> dict:
    """Map IDs to the first informative value, skipping blanks and ID placeholders."""
    if key_col not in df.columns or value_col not in df.columns:
        return {}
    keys = _clean_series(df[key_col])
    values = _clean_series(df[value_col])
    valid = keys.ne("") & values.ne("") & ~values.str.casefold().eq(keys.str.casefold())
    if not valid.any():
        return {}
    mapped = pd.DataFrame({"key": keys[valid], "value": values[valid]})
    return mapped.drop_duplicates("key").set_index("key")["value"].to_dict()


def _refine_combined_mpi(df: pd.DataFrame) -> pd.DataFrame:
    """Fill display-critical MPI/MEI fields so Browse does not show empty cells."""
    from app.services.metabolite_names import refine_metabolite_names

    df = df.copy()
    for col in [
        "Species", "Metabolite Name", "HMDB ID", "Uniprot ID", "Protein Name",
        "Gene Name", "Pathway_ID", "Pathway_Name", "Evidence_Source", "EC_Number",
    ]:
        if col in df.columns:
            df[col] = _clean_series(df[col])

    hmdb_to_name = _first_real_value_map(df, "HMDB ID", "Metabolite Name")
    uniprot_to_protein = _first_real_value_map(df, "Uniprot ID", "Protein Name")
    uniprot_to_gene = _first_real_value_map(df, "Uniprot ID", "Gene Name")

    if hmdb_to_name:
        hmdb = _clean_series(df["HMDB ID"])
        met = _clean_series(df["Metabolite Name"])
        needs_name = met.eq("") | (hmdb.ne("") & met.str.casefold().eq(hmdb.str.casefold()))
        df.loc[needs_name, "Metabolite Name"] = hmdb[needs_name].map(hmdb_to_name).fillna(met[needs_name])

    if uniprot_to_protein:
        uniprot = _clean_series(df["Uniprot ID"])
        protein = _clean_series(df["Protein Name"])
        needs_name = protein.eq("") | (uniprot.ne("") & protein.str.casefold().eq(uniprot.str.casefold()))
        df.loc[needs_name, "Protein Name"] = uniprot[needs_name].map(uniprot_to_protein).fillna(protein[needs_name])

    if uniprot_to_gene:
        gene = _clean_series(df["Gene Name"])
        needs_gene = gene.eq("")
        df.loc[needs_gene, "Gene Name"] = _clean_series(df.loc[needs_gene, "Uniprot ID"]).map(uniprot_to_gene).fillna("")

    hmdb = _clean_series(df["HMDB ID"])
    uniprot = _clean_series(df["Uniprot ID"])
    ec_number = _clean_series(df["EC_Number"])

    df = refine_metabolite_names(df, "Metabolite Name", "HMDB ID")

    protein = _clean_series(df["Protein Name"])
    df.loc[protein.eq("") & ec_number.ne(""), "Protein Name"] = "Enzyme EC " + ec_number[protein.eq("") & ec_number.ne("")]
    protein = _clean_series(df["Protein Name"])
    df.loc[protein.eq("") & uniprot.ne(""), "Protein Name"] = "UniProt " + uniprot[protein.eq("") & uniprot.ne("")]
    same_as_uniprot = uniprot.ne("") & _clean_series(df["Protein Name"]).str.casefold().eq(uniprot.str.casefold())
    df.loc[same_as_uniprot, "Protein Name"] = "UniProt " + uniprot[same_as_uniprot]
    df.loc[_clean_series(df["Protein Name"]).eq(""), "Protein Name"] = "Protein/enzyme unavailable"

    fallback_cols = {
        "Species": "Unspecified organism",
        "Gene Name": "-",
        "EC_Number": "-",
        "Pathway_Name": "-",
        "Evidence_Source": "Unspecified",
    }
    for col, fallback in fallback_cols.items():
        blank = _clean_series(df[col]).eq("")
        df.loc[blank, col] = fallback

    return df


def _get_combined_mpi():
    """Return a unified MPI + MEI DataFrame (enzymes presented as proteins)."""
    global _df_combined
    if _df_combined is not None:
        return _df_combined

    mpi = _get_df().copy()
    mpi["EC_Number"] = ""

    try:
        mei = get_mei_db()
        if not mei.empty:
            mei_norm = mei.rename(columns={
                "Metabolite_Name": "Metabolite Name",
                "HMDB_ID": "HMDB ID",
                "Uniprot_ID": "Uniprot ID",
                "Enzyme_Name": "Protein Name",
                "Gene_Name": "Gene Name",
            })
            if "EC_Number" not in mei_norm.columns:
                mei_norm["EC_Number"] = ""
            if "Evidence_Source" not in mei_norm.columns:
                mei_norm["Evidence_Source"] = "KEGG"
            # Keep only the columns that exist in MPI + EC_Number
            keep = ["Species", "Metabolite Name", "HMDB ID", "SMILES",
                    "Uniprot ID", "Protein Name", "Gene Name",
                    "Pathway_ID", "Pathway_Name", "Evidence_Source", "EC_Number"]
            for c in keep:
                if c not in mei_norm.columns:
                    mei_norm[c] = ""
            mpi = pd.concat([mpi[keep], mei_norm[keep]], ignore_index=True)
    except Exception:
        pass

    mpi = _refine_combined_mpi(mpi.fillna(""))
    # Pre-build search column
    mpi["_search_text"] = (
        mpi["Species"].str.lower() + " " +
        mpi["Metabolite Name"].str.lower() + " " +
        mpi["HMDB ID"].str.lower() + " " +
        mpi["Uniprot ID"].str.lower() + " " +
        mpi["Protein Name"].str.lower() + " " +
        mpi["Gene Name"].str.lower() + " " +
        mpi["Pathway_Name"].str.lower() + " " +
        mpi["Evidence_Source"].str.lower() + " " +
        mpi["EC_Number"].str.lower()
    )
    _df_combined = mpi
    return _df_combined

# Filter options are computed lazily
_ALL_ORGANISMS = []
_ALL_SOURCES = []

_ALL_PATHWAYS = []


def _ensure_filter_options():
    """Compute filter options on first use."""
    global _ALL_ORGANISMS, _ALL_SOURCES, _ALL_PATHWAYS
    if _ALL_ORGANISMS:
        return
    df = _get_df()
    _ALL_ORGANISMS = _valid_option_values(df["Species"].unique()) if "Species" in df.columns else []
    _ALL_SOURCES = _valid_option_values(df["Evidence_Source"].unique()) if "Evidence_Source" in df.columns else []
    pw_set = set()
    for pw in df["Pathway_Name"]:
        if pw:
            for token in str(pw).split(";"):
                token = token.strip()
                if token:
                    pw_set.add(token)
    _ALL_PATHWAYS = sorted(pw_set)


# Page-level constants
PAGE_SIZE = 20


@lru_cache(maxsize=1)
def _get_release_stats() -> dict:
    """Read lightweight release metadata without loading interaction CSVs."""
    stats_path = _cfg.DATA_DIR / "coremetdb_stats.json"
    try:
        return json.loads(stats_path.read_text())
    except Exception:
        return {}


def _stat_db(key: str) -> dict:
    return _get_release_stats().get("databases", {}).get(key, {})


def _stat_value(key: str, field: str, default: int = 0) -> int:
    value = _stat_db(key).get(field, default)
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return default


def _combined_mpi_stat(field: str) -> int:
    if field == "interactions":
        return _stat_value("MPI", field) + _stat_value("MEI", field)
    if field == "targets":
        return _stat_value("MPI", field) + _stat_value("MEI", field)
    return max(_stat_value("MPI", field), _stat_value("MEI", field))


def _fmt_stat(value: int) -> str:
    return f"{value:,}" if value else "0"


def _is_db_available(key: str) -> bool:
    return _stat_value(key, "interactions") > 0

# Evidence-source badge color mapping (CSS-safe inline colors)
SOURCE_COLORS = {
    "KEGG": {"bg": "rgba(49,130,206,0.12)", "fg": "#3182ce"},
    "Rhea": {"bg": "rgba(56,161,105,0.12)", "fg": "#38a169"},
    "original": {"bg": "rgba(26,54,93,0.08)", "fg": "#1a365d"},
    "BRENDA": {"bg": "rgba(128,90,213,0.12)", "fg": "#805ad5"},
    "HMDB": {"bg": "rgba(237,137,54,0.12)", "fg": "#ed8936"},
}

# Organism color map (matches design system)
ORG_COLORS = {
    "Homo sapiens": "#3182ce",
    "Mus musculus": "#e53e3e",
    "Rattus norvegicus": "#38a169",
    "Escherichia coli": "#d69e2e",
    "Bos taurus": "#b83280",
    "Pseudomonas aeruginosa": "#4299e1",
    "Arabidopsis thaliana": "#48bb78",
    "Saccharomyces cerevisiae": "#805ad5",
    "Drosophila melanogaster": "#ed8936",
    "Caenorhabditis elegans": "#319795",
}


# ---------------------------------------------------------------------------
# Helper: build the mini organism donut figure
# ---------------------------------------------------------------------------

def _make_organism_donut(df: pd.DataFrame, group_col: str = "Species") -> go.Figure:
    """Return a small Plotly donut chart of distribution (organism, disease, etc.)."""
    if group_col not in df.columns:
        # Fallback: try Disease_Name for MDI, then Category
        for fallback in ["Disease_Name", "Category", "Species"]:
            if fallback in df.columns:
                group_col = fallback
                break
        else:
            return go.Figure()
    
    group_values = _clean_series(df[group_col]).replace("", "Unspecified")
    counts = group_values.value_counts()
    labels = counts.index.tolist()
    values = counts.values.tolist()
    colors = [ORG_COLORS.get(lab, "#a0aec0") for lab in labels]

    fig = go.Figure(
        data=[
            go.Pie(
                labels=labels,
                values=values,
                hole=0.55,
                marker=dict(colors=colors, line=dict(color="#ffffff", width=1)),
                textposition="none",
                hovertemplate="%{label}<br>%{value:,} interactions (%{percent})<extra></extra>",
            )
        ]
    )
    fig.update_layout(
        showlegend=False,
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        height=160,
        width=160,
    )
    return fig


# ---------------------------------------------------------------------------
# Helper: build evidence-source badge markup
# ---------------------------------------------------------------------------

def _source_badge(source: str) -> str:
    """Return a Markdown-safe badge string for DataTable presentation cells."""
    # DataTable doesn't render HTML in data; we handle coloring via
    # conditional formatting.  Return the plain text.
    return str(source)


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def _make_stat_cards_row():
    """Four stat cards (dynamically updated via callback outputs)."""
    initial_total = _combined_mpi_stat("interactions")
    initial_metabolites = _combined_mpi_stat("metabolites")
    initial_targets = _combined_mpi_stat("targets")
    initial_organisms = _combined_mpi_stat("organisms")
    return dbc.Row(
        [
            dbc.Col(
                html.Div(
                    [
                        html.Div(html.I(className="fas fa-dna"), className="stat-icon"),
                        html.Div(
                            _fmt_stat(initial_total),
                            id="dbv2-stat-total",
                            className="stat-value",
                        ),
                        html.Div("Total Interactions", className="stat-label"),
                    ],
                    className="cm-stat-card accent-interaction",
                ),
                xs=6, md=3,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.Div(html.I(className="fas fa-flask"), className="stat-icon"),
                        html.Div(
                            _fmt_stat(initial_metabolites),
                            id="dbv2-stat-metabolites",
                            className="stat-value",
                        ),
                        html.Div("Metabolites", className="stat-label"),
                    ],
                    className="cm-stat-card accent-metabolite",
                ),
                xs=6, md=3,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.Div(html.I(className="fas fa-microscope"), className="stat-icon"),
                        html.Div(
                            _fmt_stat(initial_targets),
                            id="dbv2-stat-proteins",
                            className="stat-value",
                        ),
                        html.Div("Proteins", id="dbv2-stat-proteins-label", className="stat-label"),
                    ],
                    className="cm-stat-card accent-protein",
                ),
                xs=6, md=3,
            ),
            dbc.Col(
                html.Div(
                    [
                        html.Div(html.I(className="fas fa-globe"), className="stat-icon"),
                        html.Div(
                            _fmt_stat(initial_organisms),
                            id="dbv2-stat-organisms",
                            className="stat-value",
                        ),
                        html.Div("Organisms", id="dbv2-stat-organisms-label", className="stat-label"),
                    ],
                    className="cm-stat-card accent-accent",
                ),
                xs=6, md=3,
            ),
        ],
        className="g-3 mb-4",
    )


def _make_sidebar():
    """280px sidebar with organism checkboxes, pathway dropdown, evidence source checkboxes, reset."""
    return html.Div(
        [
            # ── Organism filter ──────────────────────────
            html.H6(
                [html.I(className="fas fa-globe me-2"), "Organisms"],
                id="dbv2-label-organisms",
                className="mb-2",
                style={"fontWeight": "600"},
            ),
            dbc.Checklist(
                id="dbv2-filter-organisms",
                options=[{"label": org, "value": org} for org in _ALL_ORGANISMS],
                value=[],
                inline=False,
                className="mb-3",
                style={"maxHeight": "200px", "overflowY": "auto", "fontSize": "0.85rem"},
            ),

            html.Hr(),

            # ── Pathway dropdown ──────────────────────────
            html.H6(
                [html.I(className="fas fa-route me-2"), "Pathway"],
                id="dbv2-label-pathway",
                className="mb-2",
                style={"fontWeight": "600"},
            ),
            dcc.Dropdown(
                id="dbv2-filter-pathway",
                options=[{"label": pw, "value": pw} for pw in _ALL_PATHWAYS],
                placeholder="Select pathway...",
                multi=False,
                clearable=True,
                style={"fontSize": "0.85rem"},
                className="mb-3",
            ),

            html.Hr(),

            # ── Evidence source checkboxes ────────────────
            html.H6(
                [html.I(className="fas fa-tag me-2"), "Evidence Source"],
                id="dbv2-label-sources",
                className="mb-2",
                style={"fontWeight": "600"},
            ),
            dbc.Checklist(
                id="dbv2-filter-sources",
                options=[{"label": s, "value": s} for s in _ALL_SOURCES],
                value=[],
                inline=False,
                className="mb-3",
                style={"fontSize": "0.85rem"},
            ),

            html.Hr(),

            # ── Mini donut chart ─────────────────────────
            html.H6(
                [html.I(className="fas fa-chart-pie me-2"), "Distribution"],
                id="dbv2-label-donut",
                className="mb-2",
                style={"fontWeight": "600"},
            ),
            dcc.Graph(
                id="dbv2-organism-donut",
                config={"displayModeBar": False},
                style={"height": "160px"},
            ),

            html.Hr(),

            # ── Reset button ─────────────────────────────
            dbc.Button(
                [html.I(className="fas fa-undo me-2"), "Reset Filters"],
                id="dbv2-btn-reset",
                className="cm-btn-secondary w-100",
                size="sm",
            ),
        ],
        className="cm-sidebar",
        style={
            "width": "280px",
            "minWidth": "280px",
            "padding": "16px",
            "background": "var(--cm-surface)",
            "borderRight": "1px solid var(--cm-border)",
            "borderRadius": "8px",
            "boxShadow": "var(--shadow-sm)",
            "overflowY": "auto",
            "maxHeight": "calc(100vh - 200px)",
        },
    )


def _make_table_section():
    """Global search bar, DataTable, pagination info, and CSV download."""
    # The callback populates rows once /database is opened. Keeping layout
    # construction data-free prevents app startup from loading large CSVs.
    initial_records = []
    return html.Div(
        [
            # ── Global search bar ─────────────────────────
            dbc.InputGroup(
                [
                    dbc.InputGroupText(html.I(className="fas fa-search")),
                    dbc.Input(
                        id="dbv2-search",
                        placeholder="Search by metabolite, protein, HMDB ID, UniProt ID, gene...",
                        type="text",
                        debounce=True,
                        style={"fontSize": "0.9rem"},
                    ),
                ],
                className="mb-3",
            ),

            # ── Action row (download + result count) ──────
            dbc.Row(
                [
                    dbc.Col(
                        html.Span(id="dbv2-result-count", style={"fontSize": "0.85rem", "color": "var(--cm-text-secondary)"}),
                        width="auto",
                    ),
                    dbc.Col(
                        dbc.Button(
                            [html.I(className="fas fa-download me-2"), "Download CSV"],
                            id="dbv2-btn-download",
                            className="cm-btn-primary",
                            size="sm",
                        ),
                        width="auto",
                        className="ms-auto",
                    ),
                ],
                className="mb-2 align-items-center",
            ),
            dcc.Download(id="dbv2-download-csv"),

            # ── DataTable ─────────────────────────────────
            dash_table.DataTable(
                id="dbv2-table",
                columns=[
                    {"name": "Species", "id": "Species"},
                    {"name": "Metabolite", "id": "Metabolite Name", "presentation": "markdown"},
                    {"name": "HMDB ID", "id": "HMDB ID", "presentation": "markdown"},
                    {"name": "UniProt ID", "id": "Uniprot ID", "presentation": "markdown"},
                    {"name": "Protein", "id": "Protein Name"},
                    {"name": "Gene", "id": "Gene Name"},
                    {"name": "Pathway", "id": "Pathway_Name"},
                    {"name": "Source", "id": "Evidence_Source"},
                ],
                data=initial_records,
                page_size=PAGE_SIZE,
                page_current=0,
                page_action="custom",
                sort_action="custom",
                sort_mode="multi",
                style_table={
                    "overflowX": "auto",
                    "border": "1px solid var(--cm-border)",
                    "borderRadius": "8px",
                },
                style_header={
                    "backgroundColor": "#1a365d",
                    "color": "#f7fafc",
                    "fontWeight": "600",
                    "fontSize": "0.82rem",
                    "position": "sticky",
                    "top": 0,
                    "zIndex": 10,
                    "textAlign": "left",
                    "padding": "10px 12px",
                },
                style_cell={
                    "fontSize": "0.82rem",
                    "padding": "8px 12px",
                    "textAlign": "left",
                    "maxWidth": "220px",
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                    "whiteSpace": "nowrap",
                    "fontFamily": "Arial, Helvetica, sans-serif",
                    "lineHeight": "1.4",
                    "verticalAlign": "middle",
                    "height": "auto",
                },
                style_data_conditional=[
                    # Alternating row colors
                    {
                        "if": {"row_index": "odd"},
                        "backgroundColor": "#f7fafc",
                    },
                    {
                        "if": {"row_index": "even"},
                        "backgroundColor": "#ffffff",
                    },
                    # Evidence source colored badges via conditional styling
                    {
                        "if": {
                            "filter_query": '{Evidence_Source} = "KEGG"',
                            "column_id": "Evidence_Source",
                        },
                        "backgroundColor": "rgba(49,130,206,0.12)",
                        "color": "#3182ce",
                        "fontWeight": "600",
                    },
                    {
                        "if": {
                            "filter_query": '{Evidence_Source} = "Rhea"',
                            "column_id": "Evidence_Source",
                        },
                        "backgroundColor": "rgba(56,161,105,0.12)",
                        "color": "#38a169",
                        "fontWeight": "600",
                    },
                    {
                        "if": {
                            "filter_query": '{Evidence_Source} = "original"',
                            "column_id": "Evidence_Source",
                        },
                        "backgroundColor": "rgba(26,54,93,0.08)",
                        "color": "#1a365d",
                        "fontWeight": "600",
                    },
                    {
                        "if": {
                            "filter_query": '{Evidence_Source} = "BRENDA"',
                            "column_id": "Evidence_Source",
                        },
                        "backgroundColor": "rgba(128,90,213,0.12)",
                        "color": "#805ad5",
                        "fontWeight": "600",
                    },
                    {
                        "if": {
                            "filter_query": '{Evidence_Source} = "HMDB"',
                            "column_id": "Evidence_Source",
                        },
                        "backgroundColor": "rgba(237,137,54,0.12)",
                        "color": "#ed8936",
                        "fontWeight": "600",
                    },
                ],
                style_cell_conditional=[
                    # Narrow columns
                    {"if": {"column_id": "Species"}, "width": "80px", "minWidth": "80px"},
                    {"if": {"column_id": "Rank"}, "width": "65px", "minWidth": "65px"},
                    {"if": {"column_id": "Evidence_Source"}, "width": "75px", "minWidth": "75px"},
                    {"if": {"column_id": "Source"}, "width": "75px", "minWidth": "75px"},
                    {"if": {"column_id": "Evidence_Level"}, "width": "80px", "minWidth": "80px"},
                    {"if": {"column_id": "Evidence"}, "width": "80px", "minWidth": "80px"},
                    {"if": {"column_id": "PMID"}, "width": "75px", "minWidth": "75px"},
                    {"if": {"column_id": "EC_Number"}, "width": "80px", "minWidth": "80px"},
                    {"if": {"column_id": "Taxonomy_ID"}, "width": "80px", "minWidth": "80px"},
                    {"if": {"column_id": "HMDB ID"}, "width": "120px", "minWidth": "120px"},
                    {"if": {"column_id": "HMDB_ID"}, "width": "120px", "minWidth": "120px"},
                    {"if": {"column_id": "KEGG_ID"}, "width": "80px", "minWidth": "80px"},
                    {"if": {"column_id": "Uniprot ID"}, "width": "90px", "minWidth": "90px"},
                    {"if": {"column_id": "DrugBank_ID"}, "width": "100px", "minWidth": "100px"},
                    {"if": {"column_id": "Disease_ID"}, "width": "90px", "minWidth": "90px"},
                    {"if": {"column_id": "Association_Type"}, "width": "100px", "minWidth": "100px"},
                    {"if": {"column_id": "Interaction_Type"}, "width": "100px", "minWidth": "100px"},
                    {"if": {"column_id": "Avg_Network_Score"}, "width": "100px", "minWidth": "100px"},
                    # Wider columns
                    {"if": {"column_id": "Metabolite Name"}, "width": "150px", "minWidth": "130px"},
                    {"if": {"column_id": "Metabolite_Name"}, "width": "150px", "minWidth": "130px"},
                    {"if": {"column_id": "Protein Name"}, "width": "160px", "minWidth": "130px"},
                    {"if": {"column_id": "Gene Name"}, "width": "90px", "minWidth": "80px"},
                    {"if": {"column_id": "Pathway_Name"}, "width": "180px", "minWidth": "140px"},
                    {"if": {"column_id": "Disease_Name"}, "width": "170px", "minWidth": "140px"},
                    {"if": {"column_id": "Microbe_Name"}, "width": "160px", "minWidth": "130px"},
                    {"if": {"column_id": "Drug_Name"}, "width": "140px", "minWidth": "110px"},
                    {"if": {"column_id": "Category"}, "width": "120px", "minWidth": "100px"},
                    {"if": {"column_id": "Tissue"}, "width": "100px", "minWidth": "80px"},
                    {"if": {"column_id": "Organism"}, "width": "100px", "minWidth": "80px"},
                    {"if": {"column_id": "Relationship_Type"}, "width": "100px", "minWidth": "80px"},
                    {"if": {"column_id": "Description"}, "width": "200px", "minWidth": "160px"},
                ],
                style_as_list_view=True,
                css=[{"selector": ".dash-spreadsheet", "rule": "font-family: Arial, Helvetica, sans-serif;"}],
            ),
        ],
        style={"flex": 1, "minWidth": 0},
    )


# ---------------------------------------------------------------------------
# LAYOUT
# ---------------------------------------------------------------------------

layout = html.Div(
    [
        html.Div(
            [
                make_page_header(
                    title="CoreMet Database",
                    subtitle="Browse and filter curated metabolite interaction databases with pathway annotations and evidence sources.",
                    breadcrumb_items=[("Home", "/home"), ("Database", None)],
                ),

                # ── Interaction Type Selector ────────────────────────
                html.Div(
                    [
                        dbc.RadioItems(
                            id="dbv2-interaction-type",
                            options=[
                                {"label": html.Span([
                                    html.I(className="fas fa-dna me-1"),
                                    "Metabolite–Protein (MPI)",
                                    dbc.Badge(_fmt_stat(_combined_mpi_stat("interactions")), className="ms-2", color="primary", pill=True),
                                ]), "value": "mpi"},
                                {"label": html.Span([
                                    html.I(className="fas fa-heartbeat me-1"),
                                    "Metabolite–Disease (MDI)",
                                    dbc.Badge(_fmt_stat(_stat_value("MDI", "interactions")), className="ms-2", color="danger", pill=True),
                                ]), "value": "mdi"},
                                {"label": html.Span([
                                    html.I(className="fas fa-bacteria me-1") if _is_db_available("MMI") else html.I(className="fas fa-bacterium me-1"),
                                    "Metabolite–Microbe (MMI)",
                                    dbc.Badge(
                                        _fmt_stat(_stat_value("MMI", "interactions")) if _is_db_available("MMI") else "Coming Soon",
                                        className="ms-2",
                                        color="warning" if _is_db_available("MMI") else "secondary",
                                        pill=True,
                                    ),
                                ]), "value": "mmi", "disabled": not _is_db_available("MMI")},
                                {"label": html.Span([
                                    html.I(className="fas fa-pills me-1"),
                                    "Metabolite–Drug (MDrI)",
                                    dbc.Badge(
                                        _fmt_stat(_stat_value("MDrI", "interactions")) if _is_db_available("MDrI") else "Coming Soon",
                                        className="ms-2",
                                        color="info" if _is_db_available("MDrI") else "secondary",
                                        pill=True,
                                    ),
                                ]), "value": "mdri", "disabled": not _is_db_available("MDrI")},
                                {"label": html.Span([
                                    html.I(className="fas fa-map-marker-alt me-1"),
                                    "Metabolite–SNP (mGWAS)",
                                    dbc.Badge(
                                        _fmt_stat(_stat_value("mGWAS", "interactions")) if _is_db_available("mGWAS") else "Coming Soon",
                                        className="ms-2",
                                        color="success" if _is_db_available("mGWAS") else "secondary",
                                        pill=True,
                                    ),
                                ]), "value": "mgwas", "disabled": not _is_db_available("mGWAS")},
                                {"label": html.Span([
                                    html.I(className="fas fa-dna me-1"),
                                    "Metabolite–Gene (MGI)",
                                    dbc.Badge(
                                        _fmt_stat(_stat_value("MGI", "interactions")) if _is_db_available("MGI") else "Coming Soon",
                                        className="ms-2",
                                        color="warning" if _is_db_available("MGI") else "secondary",
                                        pill=True,
                                    ),
                                ]), "value": "mgi", "disabled": not _is_db_available("MGI")},
                            ],
                            value="mpi",
                            inline=True,
                            className="mb-0",
                            inputClassName="me-1",
                        ),
                    ],
                    className="cm-card mb-3",
                    style={"padding": "12px 20px"},
                ),

                # Static summary, updates via callback
                html.Div(
                    [
                        html.Span("Database summary: ", style={"fontWeight": "600"}),
                        html.Span(
                            id="dbv2-summary-text",
                            children=(
                                f"{_fmt_stat(_combined_mpi_stat('interactions'))} interactions across "
                                f"{_fmt_stat(_combined_mpi_stat('metabolites'))} metabolites and "
                                f"{_fmt_stat(_combined_mpi_stat('targets'))} proteins/enzymes "
                                f"in {_fmt_stat(_combined_mpi_stat('organisms'))} organisms."
                            ),
                            style={"color": "var(--cm-text-secondary)"},
                        ),
                    ],
                    className="cm-card mb-3",
                    style={"fontSize": "0.9rem"},
                ),

                # Stat cards
                _make_stat_cards_row(),

                # Sidebar + Table in responsive Bootstrap row so content aligns nicely
                dbc.Row(
                    [
                        dbc.Col(
                            _make_sidebar(),
                            xs=12,
                            md=4,
                            lg=3,
                            className="mb-3 mb-md-0",
                        ),
                        dbc.Col(
                            _make_table_section(),
                            xs=12,
                            md=8,
                            lg=9,
                        ),
                    ],
                    className="g-3 align-items-start",
                ),
            ],
            className="cm-page-container",
        ),

        # Hidden store for filtered data (JSON)
        dcc.Store(id="dbv2-filtered-store"),
    ],
)


# ---------------------------------------------------------------------------
# Helpers used by callbacks
# ---------------------------------------------------------------------------

def _apply_filters(organisms, pathway, sources, search_text, interaction_type="mpi"):
    """Return a filtered view of the raw DataFrame using vectorized ops."""
    if interaction_type == "mdi":
        df = get_mdi_db()
    elif interaction_type == "mmi":
        df = get_mmi_db()
    elif interaction_type == "mdri":
        df = get_mdri_db()
    elif interaction_type == "mgwas":
        df = get_mgwas_db()
    elif interaction_type == "mgi":
        df = get_mgi_db()
    else:
        df = _get_combined_mpi()

    if df.empty:
        return df

    mask = pd.Series(True, index=df.index)

    # Organism / Category / Tissue filter
    if organisms:
        if interaction_type == "mmi" and "Organism" in df.columns:
            mask &= df["Organism"].isin(organisms)
        elif interaction_type == "mdri" and "Interaction_Type" in df.columns:
            mask &= df["Interaction_Type"].isin(organisms)
        elif interaction_type == "mgwas" and "Chromosome" in df.columns:
            mask &= df["Chromosome"].isin(organisms)
        elif interaction_type == "mgi" and "Organism" in df.columns:
            mask &= df["Organism"].isin(organisms)
        elif "Species" in df.columns:
            mask &= df["Species"].isin(organisms)
        elif "Category" in df.columns:
            mask &= df["Category"].isin(organisms)
        elif "Tissue" in df.columns:
            mask &= df["Tissue"].isin(organisms)

    # Evidence / Relationship source filter
    if sources:
        if interaction_type == "mmi" and "Relationship_Type" in df.columns:
            mask &= df["Relationship_Type"].isin(sources)
        elif interaction_type == "mdri" and "Source" in df.columns:
            mask &= df["Source"].isin(sources)
        elif interaction_type == "mgwas" and "Source" in df.columns:
            mask &= df["Source"].isin(sources)
        elif interaction_type == "mgi" and "Source" in df.columns:
            mask &= df["Source"].isin(sources)
        elif "Evidence_Source" in df.columns:
            mask &= df["Evidence_Source"].isin(sources)
        elif "Evidence_Level" in df.columns:
            mask &= df["Evidence_Level"].isin(sources)

    # Pathway / Disease / Tissue / Drug filter
    if pathway:
        if interaction_type == "mdi" and "Disease_Name" in df.columns:
            mask &= df["Disease_Name"].str.contains(pathway, case=False, na=False)
        elif interaction_type == "mmi" and "Tissue" in df.columns:
            mask &= df["Tissue"].str.contains(pathway, case=False, na=False)
        elif interaction_type == "mdri" and "Drug_Name" in df.columns:
            mask &= df["Drug_Name"].str.contains(pathway, case=False, na=False)
        elif interaction_type == "mgwas" and "Mapped_Gene" in df.columns:
            mask &= df["Mapped_Gene"].str.contains(pathway, case=False, na=False)
        elif interaction_type == "mgi" and "Gene_Symbol" in df.columns:
            mask &= df["Gene_Symbol"].str.contains(pathway, case=False, na=False)
        else:
            pw_col = "Pathway_Name" if "Pathway_Name" in df.columns else "Pathway Name"
            if pw_col in df.columns:
                mask &= df[pw_col].str.contains(pathway, case=False, na=False)

    # Global search
    if search_text:
        search_lower = search_text.lower()
        if "_search_text" in df.columns:
            mask &= df["_search_text"].str.contains(search_lower, na=False, regex=False)
        else:
            # Build ad-hoc search across all string columns
            search_mask = pd.Series(False, index=df.index)
            for col in df.columns:
                if df[col].dtype == object:
                    search_mask |= df[col].astype(str).str.lower().str.contains(search_lower, na=False, regex=False)
            mask &= search_mask

    return df.loc[mask]


def _add_markdown_links(df: pd.DataFrame, interaction_type="mpi") -> pd.DataFrame:
    """Convert IDs and metabolite names to markdown hyperlinks for the DataTable."""
    df = df.copy()

    # Determine column names
    hmdb_col = "HMDB ID" if "HMDB ID" in df.columns else "HMDB_ID"
    met_col = "Metabolite Name" if "Metabolite Name" in df.columns else "Metabolite_Name"
    uni_col = "Uniprot ID" if "Uniprot ID" in df.columns else "Uniprot_ID"

    # Metabolite name → link to /metabolite?id=HMDB... or /metabolite?name=...
    if met_col in df.columns:
        met_strs = df[met_col].astype(str)
        hmdb_strs = df[hmdb_col].astype(str) if hmdb_col in df.columns else pd.Series("", index=df.index)
        has_hmdb = _valid_hmdb_mask(hmdb_strs)
        has_hmdb &= ~hmdb_strs.str.lower().isin(_MISSING_TOKENS | {"-"})
        # With HMDB ID: link by id
        df.loc[has_hmdb, met_col] = (
            "[" + met_strs[has_hmdb] + "](/metabolite?id=" + hmdb_strs[has_hmdb] + ")"
        )
        # Without HMDB ID: link by name
        no_hmdb = ~has_hmdb & (met_strs.str.len() > 0) & ~met_strs.str.lower().isin(_MISSING_TOKENS | {"-"})
        df.loc[no_hmdb, met_col] = (
            "[" + met_strs[no_hmdb] + "](/metabolite?name=" + met_strs[no_hmdb] + ")"
        )

    if hmdb_col in df.columns:
        hmdb_strs = df[hmdb_col].astype(str)
        hmdb_mask = _valid_hmdb_mask(hmdb_strs)
        hmdb_mask &= ~hmdb_strs.str.lower().isin(_MISSING_TOKENS | {"-"})
        # Skip if already turned into metabolite link (would double-encode)
        hmdb_mask &= ~df[hmdb_col].astype(str).str.startswith("[")
        df.loc[hmdb_mask, hmdb_col] = (
            "[" + df.loc[hmdb_mask, hmdb_col] + "](https://hmdb.ca/metabolites/" + df.loc[hmdb_mask, hmdb_col] + ")"
        )
    if uni_col in df.columns:
        uni_strs = df[uni_col].astype(str)
        uni_mask = uni_strs.str.len() > 0
        uni_mask &= ~uni_strs.str.lower().isin(_MISSING_TOKENS | {"-"})
        df.loc[uni_mask, uni_col] = (
            "[" + df.loc[uni_mask, uni_col] + "](https://www.uniprot.org/uniprot/" + df.loc[uni_mask, uni_col] + ")"
        )
    # DrugBank IDs
    if "DrugBank_ID" in df.columns:
        db_mask = df["DrugBank_ID"].astype(str).str.startswith("DB")
        df.loc[db_mask, "DrugBank_ID"] = (
            "[" + df.loc[db_mask, "DrugBank_ID"] + "](https://go.drugbank.com/drugs/" + df.loc[db_mask, "DrugBank_ID"] + ")"
        )
    # rsIDs → dbSNP links
    if "rsID" in df.columns:
        rs_mask = df["rsID"].astype(str).str.startswith("rs")
        df.loc[rs_mask, "rsID"] = (
            "[" + df.loc[rs_mask, "rsID"] + "](https://www.ncbi.nlm.nih.gov/snp/" + df.loc[rs_mask, "rsID"] + ")"
        )
    # Gene IDs → NCBI Gene links (MGI)
    if "Gene_ID" in df.columns and interaction_type == "mgi":
        gid_mask = df["Gene_ID"].astype(str).str.match(r"^\d+$", na=False)
        df.loc[gid_mask, "Gene_ID"] = (
            "[" + df.loc[gid_mask, "Gene_ID"] + "](https://www.ncbi.nlm.nih.gov/gene/" + df.loc[gid_mask, "Gene_ID"] + ")"
        )
    # PMIDs → PubMed links
    if "PMID" in df.columns:
        pmid_mask = df["PMID"].astype(str).str.match(r"^\d+$", na=False)
        df.loc[pmid_mask, "PMID"] = (
            "[" + df.loc[pmid_mask, "PMID"].astype(str) +
            "](https://pubmed.ncbi.nlm.nih.gov/" +
            df.loc[pmid_mask, "PMID"].astype(str) + ")"
        )
    return df


# ---------------------------------------------------------------------------
# CALLBACKS
# ---------------------------------------------------------------------------

@callback(
    Output("dbv2-filtered-store", "data"),
    Output("dbv2-stat-total", "children"),
    Output("dbv2-stat-metabolites", "children"),
    Output("dbv2-stat-proteins", "children"),
    Output("dbv2-stat-organisms", "children"),
    Output("dbv2-organism-donut", "figure"),
    Output("dbv2-table", "data"),
    Output("dbv2-table", "columns"),
    Output("dbv2-table", "page_current"),
    Output("dbv2-table", "page_count"),
    Output("dbv2-result-count", "children"),
    Output("dbv2-summary-text", "children"),
    Output("dbv2-stat-proteins-label", "children"),
    Output("dbv2-stat-organisms-label", "children"),
    Input("dbv2-interaction-type", "value"),
    Input("dbv2-filter-organisms", "value"),
    Input("dbv2-filter-pathway", "value"),
    Input("dbv2-filter-sources", "value"),
    Input("dbv2-search", "value"),
    Input("dbv2-table", "page_current"),
    Input("dbv2-table", "page_size"),
    Input("dbv2-table", "sort_by"),
)
def update_database_view(interaction_type, organisms, pathway, sources, search_text, page_current, page_size, sort_by):
    """Master callback: filter data, update stats, donut chart, and paginated table."""
    if not interaction_type:
        interaction_type = "mpi"

    df_filtered = _apply_filters(organisms, pathway, sources, search_text, interaction_type)

    # Stats, column names differ between interaction types
    total = f"{len(df_filtered):,}"
    hmdb_col = "HMDB ID" if "HMDB ID" in df_filtered.columns else "HMDB_ID"
    uni_col = "Uniprot ID" if "Uniprot ID" in df_filtered.columns else "Uniprot_ID"
    n_metabolites = f"{df_filtered[hmdb_col].nunique():,}" if hmdb_col in df_filtered.columns else "0"

    if interaction_type == "mdi":
        n_organisms = f"{df_filtered['Category'].nunique():,}" if "Category" in df_filtered.columns else "0"
        n_proteins = f"{df_filtered['Disease_Name'].nunique():,}" if "Disease_Name" in df_filtered.columns else "0"
    elif interaction_type == "mmi":
        n_organisms = f"{df_filtered['Organism'].nunique():,}" if "Organism" in df_filtered.columns else "0"
        n_proteins = f"{df_filtered['Microbe_Name'].nunique():,}" if "Microbe_Name" in df_filtered.columns else "0"
    elif interaction_type == "mdri":
        n_organisms = f"{df_filtered['Interaction_Type'].nunique():,}" if "Interaction_Type" in df_filtered.columns else "0"
        n_proteins = f"{df_filtered['Drug_Name'].nunique():,}" if "Drug_Name" in df_filtered.columns else "0"
    elif interaction_type == "mgwas":
        n_organisms = f"{df_filtered['Chromosome'].nunique():,}" if "Chromosome" in df_filtered.columns else "0"
        n_proteins = f"{df_filtered['rsID'].nunique():,}" if "rsID" in df_filtered.columns else "0"
    elif interaction_type == "mgi":
        n_organisms = f"{df_filtered['Organism'].nunique():,}" if "Organism" in df_filtered.columns else "0"
        n_proteins = f"{df_filtered['Gene_Symbol'].nunique():,}" if "Gene_Symbol" in df_filtered.columns else "0"
    else:
        n_organisms = f"{df_filtered['Species'].nunique():,}" if "Species" in df_filtered.columns else "0"
        n_proteins = f"{df_filtered[uni_col].nunique():,}" if uni_col in df_filtered.columns else "0"

    # Donut chart
    if interaction_type == "mdi":
        donut_fig = _make_organism_donut(df_filtered, group_col="Disease_Name")
    elif interaction_type == "mmi":
        donut_fig = _make_organism_donut(df_filtered, group_col="Microbe_Name")
    elif interaction_type == "mdri":
        donut_fig = _make_organism_donut(df_filtered, group_col="Drug_Name")
    elif interaction_type == "mgwas":
        donut_fig = _make_organism_donut(df_filtered, group_col="Chromosome")
    elif interaction_type == "mgi":
        donut_fig = _make_organism_donut(df_filtered, group_col="Organism")
    else:
        donut_fig = _make_organism_donut(df_filtered, group_col="Species")

    # Server-side sorting
    if sort_by and len(sort_by):
        col_ids = [s["column_id"] for s in sort_by]
        ascending = [s["direction"] == "asc" for s in sort_by]
        # Only sort by columns that exist in the filtered data
        valid_cols = [c for c in col_ids if c in df_filtered.columns]
        if valid_cols:
            valid_asc = [ascending[i] for i, c in enumerate(col_ids) if c in df_filtered.columns]
            df_filtered = df_filtered.sort_values(valid_cols, ascending=valid_asc)

    # Server-side pagination: only send the current page to the browser
    if page_size is None:
        page_size = PAGE_SIZE
    if page_current is None:
        page_current = 0

    # Reset to page 0 when filters change (not when page/sort itself changes)
    triggered = ctx.triggered_id
    if triggered and triggered != "dbv2-table":
        page_current = 0

    start = page_current * page_size
    end = start + page_size

    # Dynamic columns based on interaction type
    if interaction_type == "mdi":
        display_cols = [
            "Metabolite_Name", "HMDB_ID", "Disease_Name", "Disease_ID",
            "Category", "Association_Type", "Evidence_Level", "Avg_Network_Score", "Source",
        ]
        table_columns = [
            {"name": "Metabolite", "id": "Metabolite_Name", "presentation": "markdown"},
            {"name": "HMDB ID", "id": "HMDB_ID", "presentation": "markdown"},
            {"name": "Disease", "id": "Disease_Name"},
            {"name": "Disease ID", "id": "Disease_ID"},
            {"name": "Category", "id": "Category"},
            {"name": "Association", "id": "Association_Type"},
            {"name": "Evidence", "id": "Evidence_Level"},
            {"name": "Network Score", "id": "Avg_Network_Score"},
            {"name": "Source", "id": "Source"},
        ]
        label = "MDI"
    elif interaction_type == "mmi":
        display_cols = [
            "Metabolite_Name", "HMDB_ID", "KEGG_ID",
            "Microbe_Name", "Taxonomy_ID", "Rank",
            "Relationship_Type", "Tissue", "Organism", "PMID",
        ]
        table_columns = [
            {"name": "Metabolite", "id": "Metabolite_Name", "presentation": "markdown"},
            {"name": "HMDB ID", "id": "HMDB_ID", "presentation": "markdown"},
            {"name": "KEGG", "id": "KEGG_ID"},
            {"name": "Microbe", "id": "Microbe_Name"},
            {"name": "Tax. ID", "id": "Taxonomy_ID"},
            {"name": "Rank", "id": "Rank"},
            {"name": "Relationship", "id": "Relationship_Type"},
            {"name": "Tissue", "id": "Tissue"},
            {"name": "Organism", "id": "Organism"},
            {"name": "PMID", "id": "PMID", "presentation": "markdown"},
        ]
        label = "MMI"
    elif interaction_type == "mdri":
        display_cols = [
            "Metabolite_Name", "HMDB_ID", "Drug_Name", "DrugBank_ID",
            "Interaction_Type", "Tissue", "Evidence_Level", "Description", "PMID", "Source",
        ]
        table_columns = [
            {"name": "Metabolite", "id": "Metabolite_Name", "presentation": "markdown"},
            {"name": "HMDB ID", "id": "HMDB_ID", "presentation": "markdown"},
            {"name": "Drug", "id": "Drug_Name"},
            {"name": "DrugBank ID", "id": "DrugBank_ID", "presentation": "markdown"},
            {"name": "Interaction", "id": "Interaction_Type"},
            {"name": "Tissue/Cell Line", "id": "Tissue"},
            {"name": "Evidence", "id": "Evidence_Level"},
            {"name": "Description", "id": "Description"},
            {"name": "PMID", "id": "PMID", "presentation": "markdown"},
            {"name": "Source", "id": "Source"},
        ]
        label = "MDrI"
    elif interaction_type == "mgwas":
        display_cols = [
            "Metabolite_Name", "HMDB_ID", "rsID", "Chromosome", "Position",
            "Mapped_Gene", "P_Value", "Beta", "Trait", "Source", "PMID",
        ]
        table_columns = [
            {"name": "Metabolite", "id": "Metabolite_Name", "presentation": "markdown"},
            {"name": "HMDB ID", "id": "HMDB_ID", "presentation": "markdown"},
            {"name": "rsID", "id": "rsID", "presentation": "markdown"},
            {"name": "Chr", "id": "Chromosome"},
            {"name": "Position", "id": "Position"},
            {"name": "Gene", "id": "Mapped_Gene"},
            {"name": "P-Value", "id": "P_Value"},
            {"name": "Beta", "id": "Beta"},
            {"name": "Trait", "id": "Trait"},
            {"name": "Source", "id": "Source"},
            {"name": "PMID", "id": "PMID", "presentation": "markdown"},
        ]
        label = "mGWAS"
    elif interaction_type == "mgi":
        display_cols = [
            "Metabolite_Name", "HMDB_ID", "Gene_Symbol", "Gene_ID",
            "Organism", "Interaction_Type", "Interaction_Actions", "Source", "PMID",
        ]
        table_columns = [
            {"name": "Metabolite", "id": "Metabolite_Name", "presentation": "markdown"},
            {"name": "HMDB ID", "id": "HMDB_ID", "presentation": "markdown"},
            {"name": "Gene", "id": "Gene_Symbol"},
            {"name": "Gene ID", "id": "Gene_ID"},
            {"name": "Organism", "id": "Organism"},
            {"name": "Interaction", "id": "Interaction_Type"},
            {"name": "Actions", "id": "Interaction_Actions"},
            {"name": "Source", "id": "Source"},
            {"name": "PMID", "id": "PMID", "presentation": "markdown"},
        ]
        label = "MGI"
    else:
        display_cols = [
            "Species", "Metabolite Name", "HMDB ID", "Uniprot ID",
            "Protein Name", "Gene Name", "EC_Number", "Pathway_Name", "Evidence_Source",
        ]
        table_columns = [
            {"name": "Species", "id": "Species"},
            {"name": "Metabolite", "id": "Metabolite Name", "presentation": "markdown"},
            {"name": "HMDB ID", "id": "HMDB ID", "presentation": "markdown"},
            {"name": "UniProt ID", "id": "Uniprot ID", "presentation": "markdown"},
            {"name": "Protein/Enzyme", "id": "Protein Name"},
            {"name": "Gene", "id": "Gene Name"},
            {"name": "EC Number", "id": "EC_Number"},
            {"name": "Pathway", "id": "Pathway_Name"},
            {"name": "Source", "id": "Evidence_Source"},
        ]
        label = "MPI"

    # Only select columns that exist
    display_cols = [c for c in display_cols if c in df_filtered.columns]
    df_page = _fill_display_blanks(_add_markdown_links(df_filtered[display_cols].iloc[start:end], interaction_type))
    table_records = df_page.to_dict("records")

    result_text = f"Showing {min(start + 1, len(df_filtered)):,}–{min(end, len(df_filtered)):,} of {len(df_filtered):,} {label} associations"

    # Compute total page count for the pagination controls
    page_count = max(1, math.ceil(len(df_filtered) / page_size))

    # Store the filtered DataFrame index list + interaction type for CSV download
    store_data = {"indices": df_filtered.index.tolist(), "type": interaction_type}

    # Summary text
    if interaction_type == "mdi":
        summary = (
            f"The MDI database contains {len(df_filtered):,} metabolite–disease associations "
            f"across {n_proteins} diseases in {n_organisms} categories, curated from literature "
            f"and CoreMet disease case studies."
        )
    elif interaction_type == "mmi":
        summary = (
            f"The MMI database contains {len(df_filtered):,} metabolite–microbe interactions "
            f"across {n_proteins} microbes in {n_organisms} tissue types, curated from "
            f"gutMGene and AGORA2 sources."
        ) if len(df_filtered) > 0 else (
            "The MMI database is awaiting data population. "
            "Download gutMGene and/or AGORA2 data, then run scripts/build_mmi_database.py."
        )
    elif interaction_type == "mdri":
        summary = (
            f"The MDrI database contains {len(df_filtered):,} metabolite\u2013drug interactions "
            f"involving {n_proteins} drugs across {n_organisms} interaction types, curated from "
            f"HMDB, CTD, and published pharmacokinetic/pharmacodynamic literature."
        )
    elif interaction_type == "mgwas":
        summary = (
            f"The mGWAS database contains {len(df_filtered):,} metabolite\u2013SNP associations "
            f"involving {n_proteins} SNPs across {n_organisms} chromosomes, curated from "
            f"the GWAS Catalog (genome-wide significant, p < 5\u00d710\u207b\u2078)."
        )
    elif interaction_type == "mgi":
        summary = (
            f"The MGI database contains {len(df_filtered):,} metabolite–gene interactions "
            f"involving {n_proteins} genes across {n_organisms} organisms, curated from "
            f"the Comparative Toxicogenomics Database (CTD)."
        )
    else:
        summary = (
            f"The MPI database contains {len(df_filtered):,} metabolite–protein/enzyme interactions "
            f"from {n_organisms} organisms, curated from KEGG, Rhea/UniProt, and HMDB."
        )

    # Stat card labels
    if interaction_type == "mdi":
        protein_label = "Diseases"
        organism_label = "Categories"
    elif interaction_type == "mmi":
        protein_label = "Microbes"
        organism_label = "Organisms"
    elif interaction_type == "mdri":
        protein_label = "Drugs"
        organism_label = "Interaction Types"
    elif interaction_type == "mgwas":
        protein_label = "SNPs"
        organism_label = "Chromosomes"
    elif interaction_type == "mgi":
        protein_label = "Genes"
        organism_label = "Organisms"
    else:
        protein_label = "Proteins/Enzymes"
        organism_label = "Organisms"

    return (
        store_data,
        total,
        n_metabolites,
        n_proteins,
        n_organisms,
        donut_fig,
        table_records,
        table_columns,
        page_current,
        page_count,
        result_text,
        summary,
        protein_label,
        organism_label,
    )


@callback(
    Output("dbv2-filter-organisms", "value"),
    Output("dbv2-filter-pathway", "value"),
    Output("dbv2-filter-sources", "value"),
    Output("dbv2-search", "value"),
    Input("dbv2-btn-reset", "n_clicks"),
    prevent_initial_call=True,
)
def reset_filters(n_clicks):
    """Clear all filter controls."""
    if n_clicks:
        return [], None, [], ""
    return no_update, no_update, no_update, no_update


@callback(
    Output("dbv2-filter-organisms", "options"),
    Output("dbv2-filter-pathway", "options"),
    Output("dbv2-filter-sources", "options"),
    Output("dbv2-label-organisms", "children"),
    Output("dbv2-label-pathway", "children"),
    Output("dbv2-label-sources", "children"),
    Output("dbv2-label-donut", "children"),
    Output("dbv2-filter-pathway", "placeholder"),
    Output("dbv2-search", "placeholder"),
    Input("dbv2-interaction-type", "value"),
)
def update_sidebar_options(interaction_type):
    """Recompute sidebar filter options and labels when interaction type changes."""
    if interaction_type == "mdi":
        df = get_mdi_db()
    elif interaction_type == "mmi":
        df = get_mmi_db()
    elif interaction_type == "mdri":
        df = get_mdri_db()
    elif interaction_type == "mgwas":
        df = get_mgwas_db()
    elif interaction_type == "mgi":
        df = get_mgi_db()
    else:
        df = _get_combined_mpi()

    # Organisms: MDI uses Category, MMI uses Organism, MDrI uses Interaction_Type, mGWAS uses Chromosome, others use Species
    if interaction_type == "mmi" and "Organism" in df.columns:
        orgs = _valid_option_values(df["Organism"].unique())
    elif interaction_type == "mdri" and "Interaction_Type" in df.columns:
        orgs = _valid_option_values(df["Interaction_Type"].unique())
    elif interaction_type == "mgwas" and "Chromosome" in df.columns:
        orgs = sorted(_valid_option_values(df["Chromosome"].unique()), key=lambda x: (int(x) if x.isdigit() else 99, x))
    elif interaction_type == "mgi" and "Organism" in df.columns:
        orgs = _valid_option_values(df["Organism"].unique())
    elif "Species" in df.columns:
        orgs = _valid_option_values(df["Species"].unique())
    elif "Category" in df.columns:
        orgs = _valid_option_values(df["Category"].unique())
    else:
        orgs = []

    # Evidence / Relationship filter: MMI uses Relationship_Type, MDrI uses Source, MDI/others use Evidence
    if interaction_type == "mmi" and "Relationship_Type" in df.columns:
        sources = _valid_option_values(df["Relationship_Type"].unique())
    elif interaction_type == "mdri" and "Source" in df.columns:
        sources = _valid_option_values(df["Source"].unique())
    elif interaction_type == "mgwas" and "Source" in df.columns:
        sources = _valid_option_values(df["Source"].unique())
    elif interaction_type == "mgi" and "Source" in df.columns:
        sources = _valid_option_values(df["Source"].unique())
    elif "Evidence_Source" in df.columns:
        sources = _valid_option_values(df["Evidence_Source"].unique())
    elif "Evidence_Level" in df.columns:
        sources = _valid_option_values(df["Evidence_Level"].unique())
    else:
        sources = []

    # Pathways, Diseases, or Tissues
    if interaction_type == "mdi" and "Disease_Name" in df.columns:
        pathways = _valid_option_values(df["Disease_Name"].unique())
    elif interaction_type == "mmi" and "Tissue" in df.columns:
        pathways = _valid_option_values(df["Tissue"].unique())
    elif interaction_type == "mdri" and "Drug_Name" in df.columns:
        pathways = _valid_option_values(df["Drug_Name"].unique())
    elif interaction_type == "mgwas" and "Mapped_Gene" in df.columns:
        pathways = _valid_option_values(df["Mapped_Gene"].unique())
    elif interaction_type == "mgi" and "Gene_Symbol" in df.columns:
        pathways = _valid_option_values(df["Gene_Symbol"].unique())
    else:
        pw_col = "Pathway_Name" if "Pathway_Name" in df.columns else "Pathway Name"
        pw_set = set()
        if pw_col in df.columns:
            for pw in df[pw_col]:
                if pw:
                    for token in str(pw).split(";"):
                        token = token.strip()
                        if token and token != "-":
                            pw_set.add(token)
        pathways = sorted(pw_set)

    # Dynamic sidebar labels
    if interaction_type == "mdi":
        org_label = [html.I(className="fas fa-th-large me-2"), "Disease Categories"]
        pw_label = [html.I(className="fas fa-heartbeat me-2"), "Disease"]
        src_label = [html.I(className="fas fa-layer-group me-2"), "Evidence Level"]
        donut_label = [html.I(className="fas fa-chart-pie me-2"), "Disease Distribution"]
        pw_placeholder = "Select disease..."
        search_placeholder = "Search by metabolite, HMDB ID, disease, category..."
    elif interaction_type == "mmi":
        org_label = [html.I(className="fas fa-globe me-2"), "Organisms"]
        pw_label = [html.I(className="fas fa-lungs me-2"), "Tissue"]
        src_label = [html.I(className="fas fa-link me-2"), "Relationship Type"]
        donut_label = [html.I(className="fas fa-chart-pie me-2"), "Microbe Distribution"]
        pw_placeholder = "Select tissue..."
        search_placeholder = "Search by metabolite, HMDB ID, microbe, tissue..."
    elif interaction_type == "mdri":
        org_label = [html.I(className="fas fa-exchange-alt me-2"), "Interaction Type"]
        pw_label = [html.I(className="fas fa-pills me-2"), "Drug"]
        src_label = [html.I(className="fas fa-database me-2"), "Data Source"]
        donut_label = [html.I(className="fas fa-chart-pie me-2"), "Drug Distribution"]
        pw_placeholder = "Select drug..."
        search_placeholder = "Search by metabolite, drug, tissue, DrugBank ID..."
    elif interaction_type == "mgwas":
        org_label = [html.I(className="fas fa-dna me-2"), "Chromosome"]
        pw_label = [html.I(className="fas fa-map-marker-alt me-2"), "Gene"]
        src_label = [html.I(className="fas fa-database me-2"), "Data Source"]
        donut_label = [html.I(className="fas fa-chart-pie me-2"), "Chromosome Distribution"]
        pw_placeholder = "Select gene..."
        search_placeholder = "Search by metabolite, rsID, gene, chromosome..."
    elif interaction_type == "mgi":
        org_label = [html.I(className="fas fa-globe me-2"), "Organism"]
        pw_label = [html.I(className="fas fa-dna me-2"), "Gene"]
        src_label = [html.I(className="fas fa-database me-2"), "Data Source"]
        donut_label = [html.I(className="fas fa-chart-pie me-2"), "Organism Distribution"]
        pw_placeholder = "Select gene..."
        search_placeholder = "Search by metabolite, HMDB ID, gene symbol, organism..."
    else:
        org_label = [html.I(className="fas fa-globe me-2"), "Organisms"]
        pw_label = [html.I(className="fas fa-route me-2"), "Pathway"]
        src_label = [html.I(className="fas fa-tag me-2"), "Evidence Source"]
        donut_label = [html.I(className="fas fa-chart-pie me-2"), "Organism Distribution"]
        pw_placeholder = "Select pathway..."
        search_placeholder = "Search by metabolite, protein, HMDB ID, UniProt ID, gene..."

    return (
        [{"label": o, "value": o} for o in orgs],
        [{"label": p, "value": p} for p in pathways],
        [{"label": s, "value": s} for s in sources],
        org_label,
        pw_label,
        src_label,
        donut_label,
        pw_placeholder,
        search_placeholder,
    )


@callback(
    Output("dbv2-download-csv", "data"),
    Input("dbv2-btn-download", "n_clicks"),
    State("dbv2-filtered-store", "data"),
    prevent_initial_call=True,
)
def download_filtered_csv(n_clicks, store_data):
    """Generate and send a CSV file of the currently filtered data."""
    if not n_clicks:
        return no_update

    # Unpack store: new format is dict with indices + type
    if isinstance(store_data, dict):
        stored_indices = store_data.get("indices")
        interaction_type = store_data.get("type", "mpi")
    else:
        stored_indices = store_data
        interaction_type = "mpi"

    if interaction_type == "mdi":
        df_full = get_mdi_db()
        export_cols = [
            "Metabolite_Name", "HMDB_ID", "SMILES", "Disease_Name",
            "Disease_ID", "MeSH_ID", "Category", "Association_Type",
            "Evidence_Level", "Avg_Network_Score", "Source",
        ]
        fname = "cormet_mdi_database_filtered.csv"
    elif interaction_type == "mmi":
        df_full = get_mmi_db()
        export_cols = [
            "Metabolite_Name", "HMDB_ID", "KEGG_ID", "PubChem_CID", "ChEBI_ID", "SMILES",
            "Microbe_Name", "Taxonomy_ID", "Rank",
            "Relationship_Type", "Tissue", "Organism",
            "Experimental_Method", "PMID", "Source",
        ]
        fname = "cormet_mmi_database_filtered.csv"
    elif interaction_type == "mdri":
        df_full = get_mdri_db()
        export_cols = [
            "Metabolite_Name", "HMDB_ID", "SMILES", "Drug_Name", "DrugBank_ID",
            "Interaction_Type", "Tissue", "Cell_Location", "Biospecimen",
            "Evidence_Level", "Description", "PMID", "Source",
        ]
        fname = "cormet_mdri_database_filtered.csv"
    elif interaction_type == "mgwas":
        df_full = get_mgwas_db()
        export_cols = [
            "HMDB_ID", "Metabolite_Name", "SMILES",
            "rsID", "Chromosome", "Position", "Mapped_Gene",
            "P_Value", "Beta", "Trait", "Source", "PMID",
        ]
        fname = "cormet_mgwas_database_filtered.csv"
    elif interaction_type == "mgi":
        df_full = get_mgi_db()
        export_cols = [
            "HMDB_ID", "Metabolite_Name", "SMILES",
            "Gene_ID", "Gene_Symbol", "Organism",
            "Interaction_Type", "Interaction_Actions", "Source", "PMID",
        ]
        fname = "cormet_mgi_database_filtered.csv"
    else:
        df_full = _get_combined_mpi()
        export_cols = [
            "Species", "Metabolite Name", "HMDB ID",
            "Uniprot ID", "Protein Name", "Gene Name",
            "EC_Number", "Pathway_ID", "Pathway_Name", "Evidence_Source",
        ]
        fname = "cormet_mpi_database_filtered.csv"

    if stored_indices is not None:
        df_out = df_full.loc[stored_indices]
    else:
        df_out = df_full

    available = [c for c in export_cols if c in df_out.columns]
    csv_string = df_out[available].to_csv(index=False)

    return dcc.send_string(csv_string, filename=fname)


# Export for main.py routing (navbar/footer handled globally)
page_content = layout
