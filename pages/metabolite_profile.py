"""
Metabolite Profile Page, unified view of all interactions for a single metabolite.

Aggregates data from MPI, MEI, MDI, MMI, and MDrI databases to show a comprehensive
metabolite profile including structure, interactions, disease associations,
enzyme connections, and drug interactions.
"""

import io
from pathlib import Path

from dash import dcc, html, Input, Output, State, callback, dash_table, no_update, ctx
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go

from components.page_header import make_page_header


# ---------------------------------------------------------------------------
# Data loaders (lazy import to avoid circular deps)
# ---------------------------------------------------------------------------

def _search_all_databases(query: str) -> dict:
    """Search for a metabolite across all databases, return summary dict."""
    results = {"mpi": pd.DataFrame(), "mei": pd.DataFrame(),
               "mdi": pd.DataFrame(), "mmi": pd.DataFrame(),
               "mdri": pd.DataFrame()}

    if not query or len(query.strip()) < 2:
        return results

    q = query.strip()

    try:
        from app.services.mei_service import get_mei_db
        mei = get_mei_db()
        if not mei.empty:
            mask = pd.Series(False, index=mei.index)
            for col in ["Metabolite_Name", "HMDB_ID"]:
                if col in mei.columns:
                    mask |= mei[col].astype(str).str.lower().str.contains(q.lower(), na=False)
            results["mei"] = mei[mask]
    except Exception:
        pass

    try:
        from app.services.mdi_service import get_mdi_db
        mdi = get_mdi_db()
        if not mdi.empty:
            mask = pd.Series(False, index=mdi.index)
            for col in ["Metabolite_Name", "HMDB_ID"]:
                if col in mdi.columns:
                    mask |= mdi[col].astype(str).str.lower().str.contains(q.lower(), na=False)
            results["mdi"] = mdi[mask]
    except Exception:
        pass

    try:
        from app.services.mmi_service import get_mmi_db
        mmi = get_mmi_db()
        if not mmi.empty:
            mask = pd.Series(False, index=mmi.index)
            for col in ["Metabolite_Name", "HMDB_ID"]:
                if col in mmi.columns:
                    mask |= mmi[col].astype(str).str.lower().str.contains(q.lower(), na=False)
            results["mmi"] = mmi[mask]
    except Exception:
        pass

    try:
        from app.config import Config
        mpi_path = Path(Config.MPI_DB_PATH)
        if mpi_path.exists():
            mpi = pd.read_csv(mpi_path)
            mask = pd.Series(False, index=mpi.index)
            for col in ["Metabolite Name", "HMDB ID"]:
                if col in mpi.columns:
                    mask |= mpi[col].astype(str).str.lower().str.contains(q.lower(), na=False)
            results["mpi"] = mpi[mask]
    except Exception:
        pass

    try:
        from app.services.mdri_service import get_mdri_db
        mdri = get_mdri_db()
        if not mdri.empty:
            mask = pd.Series(False, index=mdri.index)
            for col in ["Metabolite_Name", "HMDB_ID"]:
                if col in mdri.columns:
                    mask |= mdri[col].astype(str).str.lower().str.contains(q.lower(), na=False)
            results["mdri"] = mdri[mask]
    except Exception:
        pass

    return results


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        make_page_header(
            "Metabolite Profile",
            "Search for any metabolite to see its interactions across all CoreMet databases, proteins, enzymes, diseases, microbes, and drugs.",
            [("Home", "/home"), ("Profile", None)],
        ),

        # ── Search Section ─────────────────────────────────────
        html.Div([
            html.H5([html.I(className="fas fa-search me-2"),
                     "Search Metabolite"], className="cm-card-title"),
            dbc.Row([
                dbc.Col([
                    dbc.Input(
                        id="profile-search",
                        placeholder="Enter metabolite name or HMDB ID (e.g., Glucose, HMDB0000122)",
                        type="text",
                        debounce=True,
                        className="mb-2",
                    ),
                ], md=8),
                dbc.Col([
                    dbc.Button(
                        [html.I(className="fas fa-search me-2"), "Search"],
                        id="profile-search-btn",
                        className="cm-btn-primary w-100",
                    ),
                ], md=2),
                dbc.Col([
                    dbc.Button(
                        [html.I(className="fas fa-flask me-2"), "Example"],
                        id="profile-example-btn",
                        className="cm-btn-secondary w-100",
                    ),
                ], md=2),
            ]),
            html.Small(
                "Searches across MPI, MEI, MDI, MMI, and MDrI databases simultaneously.",
                className="text-muted",
            ),
        ], className="cm-card mb-4"),

        # ── Results ────────────────────────────────────────────
        dcc.Loading(
            id="profile-loading",
            type="default",
            children=html.Div(id="profile-results"),
        ),

    ], className="cm-page-container"),
])


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@callback(
    Output("profile-results", "children"),
    Output("profile-search", "value"),
    Input("profile-search-btn", "n_clicks"),
    Input("profile-search", "n_submit"),
    Input("profile-example-btn", "n_clicks"),
    State("profile-search", "value"),
    prevent_initial_call=True,
)
def search_metabolite_profile(n_search, n_submit, n_example, query):
    """Handle search: look up metabolite across all databases."""
    trigger = ctx.triggered_id

    if trigger == "profile-example-btn":
        query = "Glucose"

    search_val = query if trigger != "profile-example-btn" else no_update

    if not query or len(query.strip()) < 2:
        return html.Div(
            [html.I(className="fas fa-info-circle me-2"),
             "Enter a metabolite name or HMDB ID to search."],
            className="cm-empty-state p-4",
        ), search_val if trigger == "profile-example-btn" else no_update

    results = _search_all_databases(query)

    # Extract basic info from first non-empty result
    met_name = query
    hmdb_id = ""
    smiles = ""
    for db_key in ["mpi", "mei", "mdi", "mmi", "mdri"]:
        df = results[db_key]
        if df.empty:
            continue
        name_col = "Metabolite Name" if "Metabolite Name" in df.columns else "Metabolite_Name"
        if name_col in df.columns:
            met_name = df[name_col].iloc[0]
        hmdb_col = "HMDB ID" if "HMDB ID" in df.columns else "HMDB_ID"
        if hmdb_col in df.columns and pd.notna(df[hmdb_col].iloc[0]):
            hmdb_id = str(df[hmdb_col].iloc[0])
        if "SMILES" in df.columns and pd.notna(df["SMILES"].iloc[0]):
            smiles = str(df["SMILES"].iloc[0])
        break

    total_hits = sum(len(results[k]) for k in results)

    if total_hits == 0:
        return html.Div(
            [html.I(className="fas fa-search me-2"),
             f'No results found for "{query}" across any database.'],
            className="cm-empty-state p-4",
        ), "Glucose" if trigger == "profile-example-btn" else no_update

    # ── Overview card ──────────────────────────────────
    detail_href = f"/metabolite?id={hmdb_id}" if hmdb_id else f"/metabolite?name={met_name}"
    overview_items = [
        html.H5([html.I(className="fas fa-atom me-2"),
                 dcc.Link(met_name, href=detail_href, style={"color": "inherit", "textDecoration": "underline"})],
                className="cm-card-title mb-3"),
        html.Small([
            html.I(className="fas fa-arrow-right me-1"),
            dcc.Link("View full metabolite detail page →", href=detail_href),
        ], className="d-block mb-2"),
    ]
    if hmdb_id:
        overview_items.append(html.P([
            html.Strong("HMDB ID: "),
            html.A(hmdb_id, href=f"https://hmdb.ca/metabolites/{hmdb_id}", target="_blank"),
        ]))
    if smiles:
        overview_items.append(html.P([
            html.Strong("SMILES: "),
            html.Code(smiles[:80] + ("..." if len(smiles) > 80 else ""),
                      style={"fontSize": "0.8rem"}),
        ]))

    # Interaction type badges
    badges = []
    type_labels = {"mpi": ("Protein", "primary"), "mei": ("Enzyme", "success"),
                   "mdi": ("Disease", "danger"), "mmi": ("Microbe", "warning"),
                   "mdri": ("Drug", "info")}
    for k, (label, color) in type_labels.items():
        n = len(results[k])
        if n > 0:
            badges.append(dbc.Badge(f"{n} {label}", color=color, className="me-2", pill=True))
        else:
            badges.append(dbc.Badge(f"0 {label}", color="secondary", className="me-2", pill=True))
    overview_items.append(html.Div(badges, className="mb-3"))

    # Distribution pie chart
    counts = {type_labels[k][0]: len(results[k]) for k in results if len(results[k]) > 0}
    if counts:
        pie_colors = {"Protein": "#3182ce", "Enzyme": "#38a169",
                      "Disease": "#e53e3e", "Microbe": "#d69e2e",
                      "Drug": "#319795"}
        fig = go.Figure(go.Pie(
            labels=list(counts.keys()),
            values=list(counts.values()),
            marker_colors=[pie_colors.get(k, "#a0aec0") for k in counts],
            hole=0.4,
            textinfo="label+percent",
        ))
        fig.update_layout(
            height=250,
            margin=dict(l=20, r=20, t=20, b=20),
            showlegend=False,
            font=dict(family="Arial", size=11),
        )
        overview_items.append(dcc.Graph(figure=fig, config={"displaylogo": False}))

    # ── Build section cards for each database type ──────────
    sections = []

    # MPI section
    if not results["mpi"].empty:
        mpi_df = results["mpi"]
        display_cols = [c for c in ["Species", "Protein Name", "Gene Name",
                                     "Pathway_Name", "Evidence_Source"] if c in mpi_df.columns]
        sections.append(_make_section_card(
            "Metabolite–Protein Interactions (MPI)",
            "fas fa-dna", "primary",
            mpi_df, display_cols,
        ))

    # MEI section
    if not results["mei"].empty:
        mei_df = results["mei"]
        display_cols = [c for c in ["Species", "EC_Number", "Enzyme_Name",
                                     "Gene_Name", "Pathway_Name", "Evidence_Source"] if c in mei_df.columns]
        sections.append(_make_section_card(
            "Metabolite–Enzyme Interactions (MEI)",
            "fas fa-flask", "success",
            mei_df, display_cols,
        ))

    # MDI section
    if not results["mdi"].empty:
        mdi_df = results["mdi"]
        display_cols = [c for c in ["Disease_Name", "Category", "Association_Type",
                                     "Evidence_Level", "Source"] if c in mdi_df.columns]
        sections.append(_make_section_card(
            "Metabolite–Disease Associations (MDI)",
            "fas fa-heartbeat", "danger",
            mdi_df, display_cols,
        ))

    # MMI section
    if not results["mmi"].empty:
        mmi_df = results["mmi"]
        display_cols = [c for c in ["Microbe_Name", "Taxonomy_ID", "Relationship_Type",
                                     "Tissue", "Evidence_Level", "Source"] if c in mmi_df.columns]
        sections.append(_make_section_card(
            "Metabolite–Microbe Interactions (MMI)",
            "fas fa-bacterium", "warning",
            mmi_df, display_cols,
        ))

    # MDrI section
    if not results["mdri"].empty:
        mdri_df = results["mdri"]
        display_cols = [c for c in ["Drug_Name", "DrugBank_ID", "Interaction_Type",
                                     "Evidence_Level", "Description", "Source"] if c in mdri_df.columns]
        sections.append(_make_section_card(
            "Metabolite–Drug Interactions (MDrI)",
            "fas fa-pills", "info",
            mdri_df, display_cols,
        ))

    return html.Div([
        html.Div(overview_items, className="cm-card mb-4"),
        *sections,
    ]), "Glucose" if trigger == "profile-example-btn" else no_update


def _make_section_card(title, icon, color, df, display_cols):
    """Create a result section card with a data table."""
    return html.Div([
        html.H5([
            html.I(className=f"{icon} me-2"),
            title,
            dbc.Badge(f"{len(df)}", color=color, className="ms-2", pill=True),
        ], className="cm-card-title"),
        dash_table.DataTable(
            data=df[display_cols].head(50).to_dict("records") if display_cols else [],
            columns=[{"name": c.replace("_", " "), "id": c} for c in display_cols],
            sort_action="native",
            page_size=10,
            style_cell={
                "textAlign": "left", "fontSize": "0.85rem",
                "padding": "6px 10px",
                "fontFamily": "Arial, Helvetica, sans-serif",
            },
            style_header={
                "fontWeight": "600", "backgroundColor": "#edf2f7",
                "fontFamily": "Arial, Helvetica, sans-serif",
            },
            style_data_conditional=[
                {"if": {"row_index": "odd"}, "backgroundColor": "#f7fafc"},
            ],
            style_table={"overflowX": "auto"},
        ),
        html.Small(
            f"Showing {min(50, len(df))} of {len(df)} records",
            className="text-muted mt-1",
        ) if len(df) > 50 else html.Span(),
    ], className="cm-card mb-4")


# Export for main.py routing
page_content = layout
