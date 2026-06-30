"""
Module Factory, generates browse and predict pages from configuration dicts.

Each interaction type (MPI, MDI, MMI, MDrI, MGI, mGWAS) provides a config dict
specifying its database loader, columns, filters, visualizations, and stat cards.
The factory produces Dash layout components and registers callbacks.
"""

import logging
import json
from functools import lru_cache
from pathlib import Path

import pandas as pd
import dash_bootstrap_components as dbc
from dash import dcc, html, dash_table, Input, Output, State, callback, no_update
import plotly.express as px
import plotly.graph_objects as go

from components.page_header import make_page_header

logger = logging.getLogger(__name__)

_MISSING_LINK_VALUES = {"", "-", "nan", "none", "null", "<na>"}
_DISPLAY_MISSING = "N/A"


def _is_valid_external_value(value, url_template: str) -> bool:
    text = str(value or "").strip()
    if text.lower() in _MISSING_LINK_VALUES:
        return False
    template = str(url_template).lower()
    upper = text.upper()
    lower = text.lower()
    if "hmdb.ca/metabolites" in template:
        return upper.startswith("HMDB") and upper[4:].isdigit()
    if "uniprot.org/uniprot" in template:
        return bool(text) and " " not in text
    if "drugbank.com/drugs" in template:
        return upper.startswith("DB")
    if "ncbi.nlm.nih.gov/snp" in template:
        return lower.startswith("rs")
    if "ncbi.nlm.nih.gov/gene" in template:
        return text.isdigit()
    return True


def _format_external_link(value, url_template: str) -> str:
    text = str(value or "").strip()
    if not _is_valid_external_value(text, url_template):
        return _DISPLAY_MISSING if text.lower() in _MISSING_LINK_VALUES else text
    return f"[{text}]({url_template.format(value=text)})"


def _fill_module_display_values(df: pd.DataFrame) -> pd.DataFrame:
    """Make module tables readable without changing filtered source data."""
    out = df.copy()
    for col in out.columns:
        out[col] = out[col].astype("string").fillna("").str.strip()
        values = out[col]
        out.loc[values.str.lower().isin(_MISSING_LINK_VALUES), col] = _DISPLAY_MISSING

    protein_pairs = [
        ("Protein Name", "Uniprot ID"),
        ("Protein_Name", "Uniprot_ID"),
        ("Enzyme_Name", "Uniprot_ID"),
    ]
    for name_col, id_col in protein_pairs:
        if name_col not in out.columns or id_col not in out.columns:
            continue
        names = out[name_col].astype("string").fillna("").str.strip()
        ids = out[id_col].astype("string").fillna("").str.strip()
        missing_name = names.str.lower().isin(_MISSING_LINK_VALUES)
        id_as_name = ids.ne("") & names.str.casefold().eq(ids.str.casefold())
        mask = ids.ne("") & (missing_name | id_as_name)
        out.loc[mask, name_col] = "UniProt " + ids[mask]

    return out


# ── Reusable stat card ───────────────────────────────────────────────

def _stat_card(value, label, color="#1a365d", icon="fas fa-database"):
    return dbc.Col(
        html.Div([
            html.I(className=f"{icon} me-2", style={"color": color, "fontSize": "1.2rem"}),
            html.Div([
                html.Div(f"{value:,}" if isinstance(value, (int, float)) else str(value),
                         style={"fontSize": "1.5rem", "fontWeight": "700", "color": color}),
                html.Div(label, style={"fontSize": "0.75rem", "color": "#718096"}),
            ]),
        ], className="d-flex align-items-center p-3",
           style={"background": "#f7fafc", "borderRadius": "8px"}),
        xs=6, md=3,
    )


# ── Visualization builders ───────────────────────────────────────────

def _build_donut(df, column, title, color_map=None):
    """Build a Plotly donut chart from value counts."""
    counts = df[column].value_counts().head(10)
    fig = px.pie(values=counts.values, names=counts.index, hole=0.45,
                  title=title, color_discrete_sequence=px.colors.qualitative.Set2)
    fig.update_layout(margin=dict(t=40, b=10, l=10, r=10), height=280,
                       font=dict(family="Arial, Helvetica, sans-serif", size=10), showlegend=True,
                       legend=dict(font=dict(size=8)))
    fig.update_traces(textposition='inside', textinfo='percent+label',
                       textfont_size=8)
    return fig


def _build_bar_top(df, column, title, n=15, color="#0072B2"):
    """Build a horizontal bar chart of top N values."""
    counts = df[column].value_counts().head(n)
    fig = px.bar(x=counts.values, y=counts.index, orientation='h',
                  title=title, labels={'x': 'Count', 'y': ''},
                  color_discrete_sequence=[color])
    fig.update_layout(margin=dict(t=40, b=30, l=10, r=10), height=350,
                       font=dict(family="Arial, Helvetica, sans-serif", size=9), yaxis=dict(autorange='reversed'),
                       showlegend=False)
    return fig


def _build_histogram(df, column, title, color="#319795", nbins=20):
    """Build a histogram."""
    fig = px.histogram(df, x=column, nbins=nbins, title=title,
                        color_discrete_sequence=[color])
    fig.update_layout(margin=dict(t=40, b=30, l=10, r=10), height=280,
                       font=dict(family="Arial, Helvetica, sans-serif", size=9), showlegend=False)
    return fig


def _build_stacked_bar(df, x_col, color_col, title):
    """Build a stacked bar chart."""
    ct = pd.crosstab(df[x_col], df[color_col])
    fig = px.bar(ct, barmode='stack', title=title,
                  color_discrete_sequence=px.colors.qualitative.Set2)
    fig.update_layout(margin=dict(t=40, b=30, l=10, r=10), height=300,
                       font=dict(family="Arial, Helvetica, sans-serif", size=9), legend=dict(font=dict(size=8)))
    return fig


def _build_coverage_bar(df, column, title, color="#0072B2"):
    """Build a horizontal stacked bar showing % annotated vs unannotated."""
    total = len(df)
    annotated = df[column].notna().sum() if column in df.columns else 0
    pct = 100 * annotated / total if total > 0 else 0
    fig = go.Figure()
    fig.add_trace(go.Bar(y=[''], x=[pct], name=f'Annotated ({annotated:,})',
                          orientation='h', marker_color=color))
    fig.add_trace(go.Bar(y=[''], x=[100 - pct], name=f'Unannotated ({total - annotated:,})',
                          orientation='h', marker_color='#e2e8f0'))
    fig.update_layout(barmode='stack', title=title,
                       margin=dict(t=40, b=10, l=10, r=10), height=120,
                       font=dict(family="Arial, Helvetica, sans-serif", size=9), legend=dict(font=dict(size=8), orientation='h'),
                       xaxis=dict(title='%', range=[0, 100]))
    return fig


def _build_heatmap(df, row_col, col_col, title, n=10, color_scale='YlOrRd'):
    """Build a heatmap of top N × top N cross-tabulation."""
    top_rows = df[row_col].value_counts().head(n).index
    top_cols = df[col_col].value_counts().head(n).index
    sub = df[df[row_col].isin(top_rows) & df[col_col].isin(top_cols)]
    ct = pd.crosstab(sub[row_col], sub[col_col]).reindex(
        index=top_rows, columns=top_cols, fill_value=0)
    fig = px.imshow(ct, title=title, color_continuous_scale=color_scale,
                     aspect='auto')
    fig.update_layout(margin=dict(t=40, b=10, l=10, r=10), height=350,
                       font=dict(family="Arial, Helvetica, sans-serif", size=8))
    return fig


def _build_log_histogram(df, column, title, color="#319795", nbins=30):
    """Build a histogram of -log10 transformed values."""
    import numpy as np
    vals = pd.to_numeric(df[column], errors='coerce').dropna()
    vals = vals[vals > 0]
    log_vals = -np.log10(vals)
    fig = px.histogram(log_vals, nbins=nbins, title=title,
                        labels={'value': '-log10(p-value)', 'count': 'Count'},
                        color_discrete_sequence=[color])
    fig.update_layout(margin=dict(t=40, b=30, l=10, r=10), height=280,
                       font=dict(family="Arial, Helvetica, sans-serif", size=9), showlegend=False)
    return fig


def _build_counts_chart(counts, title, chart_type="bar_top", color="#0072B2"):
    """Build a chart from precomputed value counts."""
    series = pd.Series(counts or {}, dtype="int64")
    if series.empty:
        return None
    if chart_type == "donut":
        fig = px.pie(values=series.values, names=series.index, hole=0.45,
                     title=title, color_discrete_sequence=px.colors.qualitative.Set2)
        fig.update_layout(margin=dict(t=40, b=10, l=10, r=10), height=280,
                          font=dict(family="Arial, Helvetica, sans-serif", size=10),
                          showlegend=True, legend=dict(font=dict(size=8)))
        fig.update_traces(textposition='inside', textinfo='percent+label',
                          textfont_size=8)
        return fig
    fig = px.bar(x=series.values, y=series.index, orientation='h',
                 title=title, labels={'x': 'Count', 'y': ''},
                 color_discrete_sequence=[color])
    fig.update_layout(margin=dict(t=40, b=30, l=10, r=10), height=350,
                      font=dict(family="Arial, Helvetica, sans-serif", size=9),
                      yaxis=dict(autorange='reversed'), showlegend=False)
    return fig


VIZ_BUILDERS = {
    'donut': _build_donut,
    'bar_top': _build_bar_top,
    'histogram': _build_histogram,
    'stacked_bar': _build_stacked_bar,
    'coverage_bar': _build_coverage_bar,
    'heatmap': _build_heatmap,
    'log_histogram': _build_log_histogram,
}


# ── Browse Page Layout ───────────────────────────────────────────────

def create_browse_layout(config):
    """Generate a database browse page layout from a config dict.

    Args:
        config: dict with keys: type_key, type_name, short_name, icon, color,
                description, stat_cards, columns, filters, external_links
    """
    key = config['type_key']
    prefix = f"mod-{key}"

    # Header
    header = make_page_header(
        title=config['type_name'],
        subtitle=config['description'],
    )

    # Stat cards (populated by callback)
    stats_row = dbc.Row(id=f"{prefix}-stats", className="g-3 mb-3")

    # Visualization row (populated by callback)
    viz_row = html.Div(id=f"{prefix}-viz", className="mb-3")

    # Search + controls
    controls = dbc.Row([
        dbc.Col([
            dbc.InputGroup([
                dbc.InputGroupText(html.I(className="fas fa-search")),
                dbc.Input(id=f"{prefix}-search", type="text",
                          placeholder=f"Search {config['short_name']} database...",
                          debounce=True),
            ]),
        ], md=6),
        dbc.Col([
            html.Div([
                html.Span(id=f"{prefix}-count", className="me-3 text-muted"),
                dbc.Button([html.I(className="fas fa-download me-2"), "CSV"],
                           id=f"{prefix}-download-btn", size="sm",
                           outline=True, color="secondary"),
                dcc.Download(id=f"{prefix}-download"),
            ], className="d-flex align-items-center justify-content-end"),
        ], md=6),
    ], className="mb-3")

    # Filters sidebar, always present with evidence/confidence + type-specific
    filter_items = []

    # Evidence type filter (universal, from standardized provenance)
    filter_items.append(html.Div([
        html.Label("Evidence Type", className="fw-semibold mb-1",
                    style={"fontSize": "0.8rem"}),
        dcc.Dropdown(id=f"{prefix}-filter-evidence", multi=True,
                      placeholder="All evidence types"),
    ], className="mb-3"))

    # Confidence slider (universal)
    filter_items.append(html.Div([
        html.Label("Min Confidence", className="fw-semibold mb-1",
                    style={"fontSize": "0.8rem"}),
        dcc.Slider(id=f"{prefix}-filter-confidence", min=0, max=1, step=0.1,
                    value=0, marks={0: '0', 0.5: '0.5', 1: '1'},
                    tooltip={"placement": "bottom"}),
    ], className="mb-3"))

    # Type-specific filters
    if config.get('filters', {}).get('organism'):
        filter_items.append(html.Div([
            html.Label("Organism", className="fw-semibold mb-1",
                        style={"fontSize": "0.8rem"}),
            dcc.Dropdown(id=f"{prefix}-filter-organism", multi=True,
                          placeholder="All organisms"),
        ], className="mb-3"))

    if config.get('filters', {}).get('category'):
        filter_items.append(html.Div([
            html.Label("Category", className="fw-semibold mb-1",
                        style={"fontSize": "0.8rem"}),
            dcc.Dropdown(id=f"{prefix}-filter-category", multi=True,
                          placeholder="All categories"),
        ], className="mb-3"))

    # Extra biological filters (disease, tissue, chromosome, etc.)
    extra_filters = config.get('filters', {}).get('extra', [])
    for i, ef in enumerate(extra_filters):
        filter_items.append(html.Div([
            html.Label(ef['label'], className="fw-semibold mb-1",
                        style={"fontSize": "0.8rem"}),
            dcc.Dropdown(id=f"{prefix}-filter-extra-{i}", multi=True,
                          placeholder=f"All {ef['label'].lower()}s"),
        ], className="mb-3"))

    # Reset button
    filter_items.append(
        dbc.Button("Reset Filters", id=f"{prefix}-reset-filters",
                    size="sm", color="secondary", outline=True, className="w-100 mt-2"))

    # "Hide predicted" toggle, filters out Predicted evidence
    filter_items.insert(2, html.Div([
        dbc.Switch(
            id=f"{prefix}-hide-predicted",
            label="Show experimental only",
            value=False,
            style={"fontSize": "0.8rem"},
        ),
    ], className="mb-3"))

    sidebar = dbc.Col([
        html.Div([
            html.H6([html.I(className="fas fa-filter me-2"), "Filters"],
                      className="mb-3", style={"color": config['color']}),
        ] + filter_items,
                  style={"background": "#f7fafc", "borderRadius": "8px",
                         "padding": "16px"}),
    ], md=3, className="d-none d-md-block")
    table_col_width = 9

    # DataTable, use markdown for columns with external links, entity links, source links, or PMID
    display_cols = config['columns']['display']
    ext_link_cols = set(config.get('external_links', {}).keys())
    entity_link_cols = set(config.get('entity_links', {}).keys())
    source_link_cols = set(config.get('source_links', {}).keys())
    pmid_col = config.get('pmid_column')
    markdown_cols = ext_link_cols | entity_link_cols | source_link_cols
    if pmid_col and pmid_col in display_cols:
        markdown_cols.add(pmid_col)
    table_columns = []
    for c in display_cols:
        col_def = {"name": c, "id": c}
        if c in markdown_cols:
            col_def["type"] = "text"
            col_def["presentation"] = "markdown"
        table_columns.append(col_def)
    table = dash_table.DataTable(
        id=f"{prefix}-table",
        columns=table_columns,
        page_size=20,
        page_action="native",
        sort_action="native",
        sort_mode="multi",
        style_table={"overflowX": "auto"},
        style_header={
            "backgroundColor": config['color'],
            "color": "white",
            "fontWeight": "bold",
            "fontSize": "0.75rem",
            "textAlign": "left",
            "padding": "8px",
        },
        style_cell={
            "fontSize": "0.75rem",
            "padding": "6px 8px",
            "textAlign": "left",
            "maxWidth": "200px",
            "overflow": "hidden",
            "textOverflow": "ellipsis",
            "fontFamily": "Arial, Helvetica, sans-serif",
        },
        style_data_conditional=[
            {"if": {"row_index": "odd"}, "backgroundColor": "#f8fafc"},
        ],
    )

    # Assemble layout
    table_section = dbc.Row([
        sidebar,
        dbc.Col([table], md=table_col_width),
    ] if sidebar else [dbc.Col([table], md=12)], className="mb-4")

    # Predict link
    predict_banner = html.Div([
        dbc.Button([
            html.I(className="fas fa-magic me-2"),
            f"Predict Novel {config['short_name']} Interactions →"
        ], href=f"/{key}/predict", color="primary", outline=True, size="sm"),
    ], className="text-center mb-4")

    # Source references section
    refs = config.get('references', [])
    ref_section = html.Div()
    if refs:
        ref_items = []
        for r in refs:
            pmid_link = ""
            if r.get('pmid'):
                pmid_link = html.A(f"PMID: {r['pmid']}",
                                   href=f"https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/",
                                   target="_blank",
                                   style={"color": "#3182ce", "fontSize": "0.78rem", "marginLeft": "6px"})
            ref_items.append(html.Li([
                html.Span(f"{r['authors']} ", style={"fontWeight": "600", "fontSize": "0.8rem"}),
                html.Span(f"{r['title']}. ", style={"fontSize": "0.8rem"}),
                html.Span(f"{r['journal']} ({r['year']}). ", style={
                    "fontSize": "0.8rem", "fontStyle": "italic", "color": "#718096"}),
                pmid_link,
            ], style={"marginBottom": "6px"}))

        ref_section = html.Div([
            html.H6([html.I(className="fas fa-book me-2"), "Source References"],
                     style={"fontWeight": "600", "color": "#2d3748", "marginBottom": "12px"}),
            html.Ul(ref_items, style={"listStyleType": "none", "paddingLeft": "0"}),
        ], className="cm-card", style={"padding": "16px", "marginTop": "16px", "marginBottom": "24px"})

    layout = html.Div([
        header,
        dbc.Container([
            stats_row,
            viz_row,
            controls,
            table_section,
            predict_banner,
            ref_section,
        ], fluid=True, className="px-4"),
    ])

    return layout


# ── Predict Placeholder Layout ───────────────────────────────────────

def create_predict_layout(config):
    """Generate a prediction placeholder page."""
    key = config['type_key']
    return html.Div([
        make_page_header(
            title=f"Predict {config['type_name']}",
            subtitle="Powered by CoreMet-FM",
        ),
        dbc.Container([
            dbc.Alert([
                html.I(className="fas fa-flask me-2"),
                f"Prediction for {config['short_name']} interactions is powered by ",
                html.Strong("CoreMet-FM"),
                ", a multi-task heterogeneous graph neural network. ",
                "This module will be enabled upon publication of the CoreMet-FM method paper.",
            ], color="info", className="mt-4"),
            html.Div([
                dbc.Button("← Back to Browse", href=f"/{key}/browse",
                           color="secondary", outline=True),
            ], className="text-center mt-3"),
        ], fluid=True, className="px-4"),
    ])


# ── Callback Registration ────────────────────────────────────────────

def register_module_callbacks(app, config, db_loader_fn):
    """Register all callbacks for a module's browse page.

    Args:
        app: Dash app instance.
        config: module config dict.
        db_loader_fn: callable that returns a pandas DataFrame.
    """
    key = config['type_key']
    prefix = f"mod-{key}"
    display_cols = config['columns']['display']
    search_cols = config['columns'].get('search', display_cols)

    @lru_cache(maxsize=1)
    def _get_precomputed_summary():
        summary_path = config.get('precomputed_summary')
        if not summary_path:
            return {}
        path = Path(summary_path)
        if not path.is_absolute():
            path = Path(__file__).parent.parent.parent / path
        try:
            return json.loads(path.read_text())
        except Exception as exc:
            logger.warning("Precomputed summary unavailable for %s: %s", key, exc)
            return {}

    @lru_cache(maxsize=1)
    def _get_df():
        df = db_loader_fn()
        # Apply data transform if configured (e.g., add evidence_type)
        transform_fn = config.get('data_transform')
        if transform_fn and callable(transform_fn):
            df = transform_fn(df)
        try:
            from app.services.metabolite_names import refine_metabolite_names_for_known_columns
            df = refine_metabolite_names_for_known_columns(df)
        except Exception as exc:
            logger.debug("Metabolite-name refinement skipped for %s: %s", key, exc)
        # Ensure display columns exist
        for col in display_cols:
            if col not in df.columns:
                df[col] = ""
        return df

    # Stats callback
    @app.callback(
        Output(f"{prefix}-stats", "children"),
        Input(f"{prefix}-search", "id"),  # trigger on load
    )
    def update_stats(_):
        summary = _get_precomputed_summary()
        stat_values = summary.get("stat_cards", {})
        if stat_values:
            cards = []
            for sc in config.get('stat_cards', []):
                label = sc['label']
                val = stat_values.get(label, 0)
                cards.append(_stat_card(val, label, color=config['color'],
                                        icon=sc.get('icon', 'fas fa-database')))
            return cards

        df = _get_df()
        cards = []
        for sc in config.get('stat_cards', []):
            label = sc['label']
            col = sc.get('column')
            if col and col in df.columns:
                val = df[col].nunique()
            elif sc.get('key') == 'total':
                val = len(df)
            else:
                val = 0
            cards.append(_stat_card(val, label, color=config['color'],
                                     icon=sc.get('icon', 'fas fa-database')))
        return cards

    # Visualization callback
    @app.callback(
        Output(f"{prefix}-viz", "children"),
        Input(f"{prefix}-search", "id"),
    )
    def update_viz(_):
        summary = _get_precomputed_summary()
        chart_counts = summary.get("charts", {})
        if chart_counts:
            graphs = []
            for vc in config.get('viz', []):
                title = vc.get('title', '')
                counts = chart_counts.get(title)
                if counts is None:
                    continue
                fig = _build_counts_chart(
                    counts, title,
                    chart_type='donut' if vc.get('type') == 'donut' else 'bar_top',
                    color=config['color'],
                )
                if fig is not None:
                    graphs.append(dbc.Col(
                        dcc.Graph(figure=fig, config={"displayModeBar": False}),
                        md=4))
            if graphs:
                return dbc.Row(graphs, className="g-3")

        df = _get_df()
        viz_configs = config.get('viz', [])
        graphs = []
        for vc in viz_configs:
            vtype = vc['type']
            col = vc.get('column')
            title = vc.get('title', '')
            fig = None
            try:
                if vtype == 'donut' and col and col in df.columns:
                    fig = _build_donut(df, col, title)
                elif vtype == 'bar_top' and col and col in df.columns:
                    fig = _build_bar_top(df, col, title,
                                          n=vc.get('n', 15), color=config['color'])
                elif vtype == 'histogram' and col and col in df.columns:
                    fig = _build_histogram(df, col, title, color=config['color'])
                elif vtype == 'coverage_bar' and col and col in df.columns:
                    fig = _build_coverage_bar(df, col, title, color=config['color'])
                elif vtype == 'heatmap' and vc.get('row_col') in df.columns:
                    fig = _build_heatmap(df, vc['row_col'], vc['col_col'], title,
                                          n=vc.get('n', 10))
                elif vtype == 'log_histogram' and col and col in df.columns:
                    fig = _build_log_histogram(df, col, title, color=config['color'])
                elif vtype == 'stacked_bar' and vc.get('x_col') in df.columns:
                    fig = _build_stacked_bar(df, vc['x_col'], vc['color_col'], title)
            except Exception as e:
                logger.warning(f"Viz error ({vtype}): {e}")
                continue
            if fig is not None:
                graphs.append(dbc.Col(
                    dcc.Graph(figure=fig, config={"displayModeBar": False}),
                    md=4))
        return dbc.Row(graphs, className="g-3") if graphs else html.Div()

    # Table + search + filter callback
    extra_filters = config.get('filters', {}).get('extra', [])
    extra_inputs = [Input(f"{prefix}-filter-extra-{i}", "value") for i in range(len(extra_filters))]
    has_organism = config.get('filters', {}).get('organism', False)
    org_col = config.get('filters', {}).get('organism_column', 'Species')
    if has_organism:
        extra_inputs.append(Input(f"{prefix}-filter-organism", "value"))

    @app.callback(
        [Output(f"{prefix}-table", "data"),
         Output(f"{prefix}-count", "children")],
        [Input(f"{prefix}-search", "value"),
         Input(f"{prefix}-filter-evidence", "value"),
         Input(f"{prefix}-filter-confidence", "value"),
         Input(f"{prefix}-hide-predicted", "value")] + extra_inputs,
    )
    def update_table(search_term, evidence_filter, min_confidence, hide_predicted, *extra_vals):
        df = _get_df()
        # Text search
        if search_term:
            mask = df[search_cols].apply(
                lambda col: col.astype(str).str.contains(search_term, case=False,
                                                          na=False, regex=False)
            ).any(axis=1)
            df = df[mask]
        # Evidence type filter
        if evidence_filter and 'evidence_type' in df.columns:
            df = df[df['evidence_type'].isin(evidence_filter)]
        # Hide predicted
        if hide_predicted and 'evidence_type' in df.columns:
            df = df[df['evidence_type'] == 'Experimental']
        # Confidence filter
        if min_confidence and min_confidence > 0 and 'confidence' in df.columns:
            conf = pd.to_numeric(df['confidence'], errors='coerce').fillna(0)
            df = df[conf >= min_confidence]
        # Extra biological filters
        n_extra = len(extra_filters)
        for i, ef in enumerate(extra_filters):
            if i < len(extra_vals) and extra_vals[i]:
                col = ef['column']
                if col in df.columns:
                    df = df[df[col].astype(str).isin([str(v) for v in extra_vals[i]])]
        # Organism filter (appended after extra filters in extra_vals)
        if has_organism and len(extra_vals) > n_extra:
            org_val = extra_vals[n_extra]
            if org_val and org_col in df.columns:
                df = df[df[org_col].astype(str).isin([str(v) for v in org_val])]
        out = _fill_module_display_values(df[display_cols].head(5000))
        # Add external links as markdown
        ext_links = config.get('external_links', {})
        for col, url_template in ext_links.items():
            if col in out.columns:
                out[col] = out[col].apply(lambda v, _tpl=url_template: _format_external_link(v, _tpl))
        # Add source database links as markdown
        src_links = config.get('source_links', {})
        for col, url_map in src_links.items():
            if col in out.columns and col not in ext_links:
                def _fmt_source(v, _map=url_map):
                    if pd.isna(v) or not str(v).strip():
                        return str(v) if pd.notna(v) else ""
                    s = str(v).strip()
                    if s.lower() in _MISSING_LINK_VALUES:
                        return _DISPLAY_MISSING
                    url = _map.get(s)
                    if url:
                        return f"[{s}]({url})"
                    # Try partial match
                    for key, u in _map.items():
                        if key in s:
                            return f"[{s}]({u})"
                    return s
                out[col] = out[col].apply(_fmt_source)
        # Add PMID links to PubMed
        pmid_col = config.get('pmid_column')
        if pmid_col and pmid_col in out.columns:
            def _fmt_pmid(v):
                if pd.isna(v) or not str(v).strip():
                    return _DISPLAY_MISSING
                s = str(v).strip().rstrip('.0')
                # Handle multiple PMIDs separated by ; or ,
                first = s.split(';')[0].split(',')[0].strip()
                if first.isdigit():
                    return f"[{first}](https://pubmed.ncbi.nlm.nih.gov/{first}/)"
                return s
            out[pmid_col] = out[pmid_col].apply(_fmt_pmid)
        # Add entity detail links as markdown
        from urllib.parse import quote_plus as _qp
        entity_links = config.get('entity_links', {})
        for col, link_cfg in entity_links.items():
            if col not in out.columns or col in ext_links:
                continue  # skip if already handled by external_links
            route = link_cfg['route']
            param = link_cfg['param']
            id_col = link_cfg.get('id_col')  # use a different column for the param value
            if id_col and id_col in out.columns:
                out[col] = out.apply(
                    lambda r, _c=col, _r=route, _p=param, _ic=id_col: (
                        f"[{r[_c]}]({_r}?{_p}={_qp(str(r[_ic]))})"
                        if pd.notna(r[_c]) and str(r[_c]).strip().lower() not in _MISSING_LINK_VALUES
                        and pd.notna(r[_ic]) and str(r[_ic]).strip().lower() not in _MISSING_LINK_VALUES
                        else str(r[_c]) if pd.notna(r[_c]) and str(r[_c]).strip().lower() not in _MISSING_LINK_VALUES else _DISPLAY_MISSING
                    ), axis=1)
            else:
                out[col] = out[col].apply(
                    lambda v, _r=route, _p=param: (
                        f"[{v}]({_r}?{_p}={_qp(str(v))})"
                        if pd.notna(v) and str(v).strip().lower() not in _MISSING_LINK_VALUES else _DISPLAY_MISSING
                    ))
        count_text = f"{len(df):,} records"
        return out.to_dict("records"), count_text

    # Filter callbacks (organism)
    if config.get('filters', {}).get('organism'):
        org_col = config['filters'].get('organism_column', 'Species')

        @app.callback(
            Output(f"{prefix}-filter-organism", "options"),
            Input(f"{prefix}-search", "id"),
        )
        def populate_organisms(_):
            df = _get_df()
            if org_col in df.columns:
                orgs = sorted(df[org_col].dropna().unique())
                return [{"label": o, "value": o} for o in orgs]
            return []

    # Evidence type filter population
    @app.callback(
        Output(f"{prefix}-filter-evidence", "options"),
        Input(f"{prefix}-search", "id"),
    )
    def populate_evidence(_):
        df = _get_df()
        if 'evidence_type' in df.columns:
            types = sorted(df['evidence_type'].dropna().unique())
            return [{"label": t.capitalize(), "value": t} for t in types]
        return []

    # Extra filter population
    extra_filters_cfg = config.get('filters', {}).get('extra', [])
    for idx, ef_cfg in enumerate(extra_filters_cfg):
        _col = ef_cfg['column']
        _idx = idx

        @app.callback(
            Output(f"{prefix}-filter-extra-{_idx}", "options"),
            Input(f"{prefix}-search", "id"),
        )
        def populate_extra(_, __col=_col):
            df = _get_df()
            if __col in df.columns:
                vals = df[__col].dropna().astype(str).value_counts().head(100).index.tolist()
                return [{"label": v, "value": v} for v in vals]
            return []

    # Reset filters, include extra filters
    extra_reset_outputs = [Output(f"{prefix}-filter-extra-{i}", "value") for i in range(len(extra_filters_cfg))]
    if has_organism:
        extra_reset_outputs.append(Output(f"{prefix}-filter-organism", "value"))

    # Reset filters
    @app.callback(
        [Output(f"{prefix}-search", "value"),
         Output(f"{prefix}-filter-evidence", "value"),
         Output(f"{prefix}-filter-confidence", "value"),
         Output(f"{prefix}-hide-predicted", "value")] + extra_reset_outputs,
        Input(f"{prefix}-reset-filters", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_filters(_):
        n_resets = len(extra_filters_cfg) + (1 if has_organism else 0)
        return ("", None, 0, False) + tuple(None for _ in range(n_resets))

    # CSV download
    @app.callback(
        Output(f"{prefix}-download", "data"),
        Input(f"{prefix}-download-btn", "n_clicks"),
        State(f"{prefix}-search", "value"),
        prevent_initial_call=True,
    )
    def download_csv(n_clicks, search_term):
        if not n_clicks:
            return no_update
        df = _get_df()
        if search_term:
            mask = df[search_cols].apply(
                lambda col: col.astype(str).str.contains(search_term, case=False,
                                                          na=False, regex=False)
            ).any(axis=1)
            df = df[mask]
        return dcc.send_data_frame(df[display_cols].to_csv, f"cormet_{key}_export.csv",
                                    index=False)
