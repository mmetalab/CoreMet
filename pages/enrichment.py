"""
Enrichment Analysis page, Pathway, Disease, Microbe, Drug & Gene enrichment.
Upload predictions or select disease, run Fisher's exact + FDR.
Pathway enrichment: enriched KEGG pathways among predicted proteins.
Disease enrichment: enriched diseases among predicted metabolites (via MDI).
Microbe enrichment: enriched gut microbes among predicted metabolites (via MMI).
Drug enrichment: enriched drug associations among predicted metabolites (via MDrI).
Gene enrichment: enriched gene associations among predicted metabolites (via MGI).
"""

import base64
import io
import logging
from pathlib import Path

from dash import dcc, html, Input, Output, State, callback, dash_table, no_update, ctx
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go

from components.page_header import make_page_header

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"
DISEASE_PRED_DIR = DATA_DIR / "mpidatabase" / "disease_mpi"

DISEASE_OPTIONS = [
    {"label": "-- Select a disease --", "value": ""},
    {"label": "Hepatocellular Carcinoma (HCC)", "value": "hcc"},
    {"label": "Papillary Thyroid Cancer", "value": "thyroid_cancer"},
    {"label": "Schizophrenia", "value": "schizophrenia"},
    {"label": "Breast Cancer", "value": "breast_cancer"},
    {"label": "Alzheimer's Disease", "value": "alzheimers"},
    {"label": "Type 2 Diabetes", "value": "t2_diabetes"},
]

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        make_page_header(
            "Enrichment Analysis",
            "For any metabolite set, identify enriched KEGG pathways, diseases, microbes, drugs, "
            "and genes from the CoreMet interaction layers using Fisher's exact test with BH FDR correction.",
            [("Home", "/home"), ("Enrichment", None)],
        ),

        # ── Input section ──────────────────────────────────────
        html.Div([
            html.H5("Input Data", className="cm-card-title"),
            dbc.Row([
                dbc.Col([
                    html.H6([html.I(className="fas fa-cloud-upload-alt me-2"),
                             "Option A: Upload Predictions CSV"]),
                    html.P(
                        "Upload a CSV with columns: Metabolite, Protein, Uniprot_ID (or Uniprot ID), Prediction Score.",
                        className="text-muted", style={"fontSize": "0.85rem"},
                    ),
                    dcc.Upload(
                        id="enrich-upload",
                        children=html.Div([
                            html.I(className="fas fa-file-csv me-2"),
                            "Drag & drop predictions CSV, or click to browse",
                        ]),
                        className="cm-upload-zone",
                        style={"padding": "24px", "cursor": "pointer"},
                        multiple=False,
                    ),
                    html.Div(id="enrich-upload-status", className="mt-2"),
                ], md=6),
                dbc.Col([
                    html.H6([html.I(className="fas fa-disease me-2"),
                             "Option B: Select Disease"]),
                    html.P(
                        "Load pre-computed interaction predictions for a disease case study.",
                        className="text-muted", style={"fontSize": "0.85rem"},
                    ),
                    dbc.Select(
                        id="enrich-disease-select",
                        options=DISEASE_OPTIONS,
                        value="",
                        className="mb-3",
                    ),
                    html.Div(id="enrich-disease-status", className="mt-2"),
                ], md=6),
            ]),

            html.Hr(className="my-3"),

            # Example data button + parameters
            dbc.Row([
                dbc.Col([
                    dbc.Button(
                        [html.I(className="fas fa-flask me-2"),
                         "Load Example (T2 Diabetes)"],
                        id="enrich-load-example",
                        className="cm-btn-secondary",
                        size="sm",
                    ),
                    html.Small(
                        ", loads 300 predicted interactions for Type 2 Diabetes",
                        className="text-muted ms-2",
                    ),
                ], md=12, className="mb-3"),
            ]),

            # Parameters row
            dbc.Row([
                dbc.Col([
                    html.Label("FDR Threshold", className="fw-semibold"),
                    dbc.Input(id="enrich-fdr", type="number", value=0.25,
                              min=0.001, max=1, step=0.01),
                ], md=4),
                dbc.Col([
                    html.Label("Min Prediction Score", className="fw-semibold"),
                    dbc.Input(id="enrich-score", type="number", value=0.3,
                              min=0, max=1, step=0.05),
                ], md=4),
            ], className="mt-2 mb-3"),
            # Enrichment buttons, full width row
            html.Div([
                dbc.Button(
                    [html.I(className="fas fa-play me-2"),
                     "Pathway"],
                    id="enrich-run",
                    className="cm-btn-primary",
                    size="sm",
                    style={"flex": "1"},
                ),
                dbc.Button(
                    [html.I(className="fas fa-disease me-2"),
                     "Disease"],
                    id="enrich-run-disease",
                    color="danger",
                    outline=True,
                    size="sm",
                    style={"flex": "1"},
                ),
                dbc.Button(
                    [html.I(className="fas fa-bacterium me-2"),
                     "Microbe"],
                    id="enrich-run-microbe",
                    color="success",
                    outline=True,
                    size="sm",
                    style={"flex": "1"},
                ),
                dbc.Button(
                    [html.I(className="fas fa-pills me-2"),
                     "Drug"],
                    id="enrich-run-drug",
                    color="info",
                    outline=True,
                    size="sm",
                    style={"flex": "1"},
                ),
                dbc.Button(
                    [html.I(className="fas fa-dna me-2"),
                     "Gene"],
                    id="enrich-run-gene",
                    color="warning",
                    outline=True,
                    size="sm",
                    style={"flex": "1"},
                ),
            ], style={"display": "flex", "gap": "8px"}, className="mb-3"),

            # PDF Report download button
            html.Div([
                dbc.Button(
                    [html.I(className="fas fa-file-pdf me-2"),
                     "Download Unified PDF Report"],
                    id="enrich-pdf-btn",
                    className="cm-btn-secondary",
                    color="dark",
                    outline=True,
                    size="sm",
                ),
                html.Small(
                    ", runs all 5 enrichment analyses and generates a single PDF",
                    className="text-muted ms-2",
                ),
            ], className="mt-3"),
        ], className="cm-card mb-4"),

        # ── Data preview ────────────────────────────────────────
        html.Div(id="enrich-data-preview", className="mb-4"),

        # ── Stores ──────────────────────────────────────────────
        dcc.Store(id="enrich-data-store"),

        # ── Results (tabbed: pathway + disease) ────────────────
        dbc.Tabs([
            dbc.Tab(
                dcc.Loading(
                    id="enrich-loading",
                    type="default",
                    children=html.Div(id="enrich-results"),
                ),
                label="Pathway Enrichment",
                tab_id="tab-pathway",
                active_label_style={"fontWeight": "bold"},
            ),
            dbc.Tab(
                dcc.Loading(
                    id="enrich-disease-loading",
                    type="default",
                    children=html.Div(id="enrich-disease-results"),
                ),
                label="Disease Enrichment",
                tab_id="tab-disease",
                active_label_style={"fontWeight": "bold"},
                tab_class_name="ms-2",
            ),
            dbc.Tab(
                dcc.Loading(
                    id="enrich-microbe-loading",
                    type="default",
                    children=html.Div(id="enrich-microbe-results"),
                ),
                label="Microbe Enrichment",
                tab_id="tab-microbe",
                active_label_style={"fontWeight": "bold"},
                tab_class_name="ms-2",
            ),
            dbc.Tab(
                dcc.Loading(
                    id="enrich-drug-loading",
                    type="default",
                    children=html.Div(id="enrich-drug-results"),
                ),
                label="Drug Enrichment",
                tab_id="tab-drug",
                active_label_style={"fontWeight": "bold"},
                tab_class_name="ms-2",
            ),
            dbc.Tab(
                dcc.Loading(
                    id="enrich-gene-loading",
                    type="default",
                    children=html.Div(id="enrich-gene-results"),
                ),
                label="Gene Enrichment (MGI)",
                tab_id="tab-gene",
                active_label_style={"fontWeight": "bold"},
                tab_class_name="ms-2",
            ),
        ], id="enrich-tabs", active_tab="tab-pathway", className="mb-4"),
    ], className="cm-page-container"),

    dcc.Download(id="enrich-download"),
    dcc.Download(id="enrich-disease-download"),
    dcc.Download(id="enrich-microbe-download"),
    dcc.Download(id="enrich-drug-download"),
    dcc.Download(id="enrich-gene-download"),
    dcc.Download(id="enrich-pdf-download"),
])


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@callback(
    Output("enrich-data-store", "data"),
    Output("enrich-upload-status", "children"),
    Output("enrich-disease-status", "children"),
    Output("enrich-data-preview", "children"),
    Output("enrich-disease-select", "value"),
    Input("enrich-upload", "contents"),
    State("enrich-upload", "filename"),
    Input("enrich-disease-select", "value"),
    Input("enrich-load-example", "n_clicks"),
    prevent_initial_call=True,
)
def load_data(contents, filename, disease, example_clicks):
    trigger = ctx.triggered_id
    df = None

    upload_status = ""
    disease_status = ""
    disease_val = no_update

    if trigger == "enrich-upload" and contents:
        try:
            ct, cs = contents.split(',')
            decoded = base64.b64decode(cs)
            if filename and filename.endswith('.xlsx'):
                df = pd.read_excel(io.BytesIO(decoded))
            else:
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            upload_status = html.Div(
                [html.I(className="fas fa-check-circle me-1 text-success"),
                 f"Loaded {len(df)} rows from {filename}"],
                style={"color": "#38a169", "fontWeight": "500",
                       "fontSize": "0.85rem"},
            )
            disease_val = ""
        except Exception as e:
            upload_status = html.Div(
                [html.I(className="fas fa-exclamation-circle me-1"),
                 f"Error: {e}"],
                style={"color": "#e53e3e", "fontSize": "0.85rem"},
            )
            return no_update, upload_status, "", no_update, no_update

    elif trigger == "enrich-disease-select" and disease:
        pred_path = DISEASE_PRED_DIR / disease / "predictions.csv"
        if pred_path.exists():
            df = pd.read_csv(pred_path)
            label = next(
                (o["label"] for o in DISEASE_OPTIONS if o["value"] == disease),
                disease,
            )
            disease_status = html.Div(
                [html.I(className="fas fa-check-circle me-1 text-success"),
                 f"Loaded {len(df)} predictions for {label}"],
                style={"color": "#38a169", "fontWeight": "500",
                       "fontSize": "0.85rem"},
            )
        else:
            disease_status = html.Div(
                [html.I(className="fas fa-exclamation-circle me-1"),
                 f"No prediction data found for {disease}"],
                style={"color": "#e53e3e", "fontSize": "0.85rem"},
            )
            return no_update, "", disease_status, no_update, no_update

    elif trigger == "enrich-load-example":
        pred_path = DISEASE_PRED_DIR / "t2_diabetes" / "predictions.csv"
        if pred_path.exists():
            df = pd.read_csv(pred_path)
            disease_status = html.Div(
                [html.I(className="fas fa-check-circle me-1 text-success"),
                 f"Example loaded: {len(df)} predictions for Type 2 Diabetes"],
                style={"color": "#38a169", "fontWeight": "500",
                       "fontSize": "0.85rem"},
            )
            disease_val = "t2_diabetes"
        else:
            disease_status = html.Div(
                [html.I(className="fas fa-exclamation-circle me-1"),
                 "Example data not found"],
                style={"color": "#e53e3e", "fontSize": "0.85rem"},
            )
            return no_update, "", disease_status, no_update, no_update

    if df is None or df.empty:
        return no_update, upload_status, disease_status, no_update, disease_val

    # Build preview card
    preview = html.Div([
        html.H6(
            [html.I(className="fas fa-table me-2"),
             f"Data Preview \u2014 {len(df)} rows"],
            className="cm-card-title",
        ),
        dash_table.DataTable(
            data=df.head(5).round(3).to_dict('records'),
            columns=[{'name': c, 'id': c} for c in df.columns],
            style_cell={
                'textAlign': 'left', 'fontSize': '0.8rem',
                'padding': '4px 8px',
                'fontFamily': 'Arial, Helvetica, sans-serif',
            },
            style_header={
                'fontWeight': '600', 'backgroundColor': '#edf2f7',
                'fontFamily': 'Arial, Helvetica, sans-serif',
            },
            style_table={'overflowX': 'auto'},
        ),
    ], className="cm-card")

    return (df.to_json(orient='split'), upload_status, disease_status,
            preview, disease_val)


@callback(
    Output("enrich-results", "children"),
    Input("enrich-run", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def run_enrichment_callback(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return html.Div(
            [html.I(className="fas fa-info-circle me-2"),
             "Please load data first, then click 'Run Enrichment'."],
            className="cm-empty-state p-4",
        )

    fdr_thresh = fdr_thresh or 0.25
    score_thresh = score_thresh or 0.3

    df = pd.read_json(io.StringIO(data), orient='split')

    try:
        from app.services.enrichment_service import run_enrichment
        results = run_enrichment(df, "All", fdr_thresh, score_thresh)
    except Exception as e:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             f"Enrichment failed: {e}"],
            className="alert alert-danger",
        )

    if results.empty:
        return html.Div(
            [html.I(className="fas fa-search me-2"),
             "No enriched pathways found. Try lowering the FDR threshold "
             "or min score filter."],
            className="cm-empty-state p-4",
        )

    # Display results
    sig = results[results['FDR'] <= fdr_thresh]
    display_df = sig if not sig.empty else results.head(20)
    top_chart = display_df.head(20).sort_values('FDR', ascending=False)

    # ── Bar chart ──────────────────────────────────────────
    colors = [
        '#38a169' if fdr < 0.01 else
        '#3182ce' if fdr < 0.05 else
        '#d69e2e' if fdr < 0.1 else '#a0aec0'
        for fdr in top_chart['FDR']
    ]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=top_chart['Pathway_Name'].str[:55],
        x=top_chart['Fold_Enrichment'],
        orientation='h',
        marker_color=colors,
        text=[f"FDR={f:.2e}" for f in top_chart['FDR']],
        textposition='outside',
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Fold Enrichment: %{x:.2f}<br>"
            "%{text}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        height=max(350, len(top_chart) * 35),
        margin=dict(l=300, r=100, t=30, b=50),
        xaxis_title="Fold Enrichment",
        yaxis=dict(automargin=True),
        font=dict(family="Arial", size=11),
        plot_bgcolor="white",
    )
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

    # ── Summary stats ──────────────────────────────────────
    n_total = len(results)
    n_sig = len(results[results['FDR'] <= fdr_thresh])

    summary = dbc.Row([
        dbc.Col(html.Div([
            html.Div(str(n_total), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Total Pathways Tested", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(n_sig), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700",
                            "color": "#38a169" if n_sig > 0 else "#a0aec0"}),
            html.Div(f"Significant (FDR \u2264 {fdr_thresh})",
                     className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{results['Fold_Enrichment'].max():.1f}\u00d7",
                     className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Max Fold Enrichment", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{results['FDR'].min():.2e}", className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Best FDR", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
    ], className="mb-4")

    # ── Results table ──────────────────────────────────────
    display_cols = [c for c in results.columns if c != 'Significant']
    table = dash_table.DataTable(
        id="enrich-results-table",
        data=results.round(4).to_dict('records'),
        columns=[{'name': c, 'id': c} for c in display_cols],
        sort_action='native',
        filter_action='native',
        page_size=20,
        style_cell={
            'textAlign': 'left', 'fontSize': '0.85rem',
            'padding': '6px 10px',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_header={
            'fontWeight': '600', 'backgroundColor': '#1a365d',
            'color': 'white',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#f7fafc'},
            {
                'if': {
                    'filter_query': f'{{FDR}} <= {fdr_thresh}',
                    'column_id': 'FDR',
                },
                'color': '#38a169',
                'fontWeight': '600',
            },
        ],
        style_table={'overflowX': 'auto'},
    )

    # ── Color legend ───────────────────────────────────────
    legend = html.Div([
        html.Span("Bar colors: ",
                   style={"fontWeight": "600", "fontSize": "0.8rem"}),
        html.Span("\u25cf FDR < 0.01 ",
                   style={"color": "#38a169", "fontSize": "0.8rem"}),
        html.Span("\u25cf FDR < 0.05 ",
                   style={"color": "#3182ce", "fontSize": "0.8rem"}),
        html.Span("\u25cf FDR < 0.10 ",
                   style={"color": "#d69e2e", "fontSize": "0.8rem"}),
        html.Span("\u25cf FDR \u2265 0.10",
                   style={"color": "#a0aec0", "fontSize": "0.8rem"}),
    ], className="mb-3")

    return html.Div([
        html.Div([
            html.H5([html.I(className="fas fa-chart-bar me-2"),
                     "Enrichment Summary"], className="cm-card-title"),
            summary,
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-chart-bar me-2"),
                     "Enriched Pathways"], className="cm-card-title"),
            legend,
            dcc.Graph(id="enrich-bar-chart", figure=fig,
                      config={"displaylogo": False}),
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-table me-2"),
                     "Full Results Table"], className="cm-card-title"),
            table,
            html.Div(
                dbc.Button(
                    [html.I(className="fas fa-download me-2"),
                     "Download Results CSV"],
                    id="enrich-download-btn",
                    className="cm-btn-secondary mt-3",
                ),
                className="text-end",
            ),
        ], className="cm-card"),
    ])


@callback(
    Output("enrich-download", "data"),
    Input("enrich-download-btn", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def download_enrichment(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return no_update
    try:
        df = pd.read_json(io.StringIO(data), orient='split')
        from app.services.enrichment_service import run_enrichment
        results = run_enrichment(df, "All", fdr_thresh or 0.25,
                                 score_thresh or 0.3)
        if results.empty:
            return no_update
        return dict(content=results.to_csv(index=False),
                    filename="pathway_enrichment.csv")
    except Exception:
        return no_update


# ---------------------------------------------------------------------------
# Disease Enrichment callback
# ---------------------------------------------------------------------------

@callback(
    Output("enrich-disease-results", "children"),
    Input("enrich-run-disease", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def run_disease_enrichment_callback(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return html.Div(
            [html.I(className="fas fa-info-circle me-2"),
             "Please load data first, then click 'Disease Enrichment'."],
            className="cm-empty-state p-4",
        )

    fdr_thresh = fdr_thresh or 0.25
    score_thresh = score_thresh or 0.3

    df = pd.read_json(io.StringIO(data), orient='split')

    # Filter by min prediction score
    score_col = None
    for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
        if c in df.columns:
            score_col = c
            break
    if score_col:
        df = df[df[score_col] >= score_thresh]

    # Extract metabolite identifiers (names or HMDB IDs)
    metabolite_ids = set()
    for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break
    for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break

    if not metabolite_ids:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             "No metabolite identifiers found in your data. "
             "Expected columns: Metabolite, HMDB ID, or Metabolite_Name."],
            className="alert alert-warning",
        )

    try:
        from app.services.enrichment_service import run_disease_enrichment
        results = run_disease_enrichment(list(metabolite_ids), fdr_thresh)
    except Exception as e:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             f"Disease enrichment failed: {e}"],
            className="alert alert-danger",
        )

    if results.empty:
        return html.Div(
            [html.I(className="fas fa-search me-2"),
             "No enriched disease associations found. This may be due to "
             "limited overlap between your metabolites and the MDI database."],
            className="cm-empty-state p-4",
        )

    # Display results  ─────────────────────────────────────
    sig = results[results['FDR'] <= fdr_thresh]
    display_df = sig if not sig.empty else results.head(20)
    top_chart = display_df.head(15).sort_values('FDR', ascending=False)

    # ── Bar chart ──────────────────────────────────────────
    cat_colors = {
        'Cancer': '#e53e3e',
        'Neurodegenerative': '#805ad5',
        'Neuropsychiatric': '#3182ce',
        'Metabolic': '#d69e2e',
    }
    colors = [cat_colors.get(cat, '#a0aec0') for cat in top_chart['Category']]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=top_chart['Disease_Name'].str[:55],
        x=top_chart['Fold_Enrichment'],
        orientation='h',
        marker_color=colors,
        text=[f"FDR={f:.2e}" for f in top_chart['FDR']],
        textposition='outside',
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Fold Enrichment: %{x:.2f}<br>"
            "%{text}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        height=max(300, len(top_chart) * 45),
        margin=dict(l=250, r=100, t=30, b=50),
        xaxis_title="Fold Enrichment",
        yaxis=dict(automargin=True),
        font=dict(family="Arial", size=11),
        plot_bgcolor="white",
    )
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

    # ── Summary stats ──────────────────────────────────────
    n_total = len(results)
    n_sig = len(sig)

    summary = dbc.Row([
        dbc.Col(html.Div([
            html.Div(str(n_total), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Diseases Tested", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(n_sig), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700",
                            "color": "#e53e3e" if n_sig > 0 else "#a0aec0"}),
            html.Div(f"Significant (FDR \u2264 {fdr_thresh})",
                     className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{results['Fold_Enrichment'].max():.1f}\u00d7",
                     className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Max Fold Enrichment", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(len(metabolite_ids)), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Input Metabolites", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
    ], className="mb-4")

    # ── Results table ──────────────────────────────────────
    display_cols = [c for c in results.columns if c != 'Significant']
    table = dash_table.DataTable(
        id="enrich-disease-results-table",
        data=results.round(4).to_dict('records'),
        columns=[{'name': c.replace('_', ' '), 'id': c} for c in display_cols],
        sort_action='native',
        filter_action='native',
        page_size=20,
        style_cell={
            'textAlign': 'left', 'fontSize': '0.85rem',
            'padding': '6px 10px',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_header={
            'fontWeight': '600', 'backgroundColor': '#742a2a',
            'color': 'white',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#fff5f5'},
            {
                'if': {
                    'filter_query': f'{{FDR}} <= {fdr_thresh}',
                    'column_id': 'FDR',
                },
                'color': '#e53e3e',
                'fontWeight': '600',
            },
        ],
        style_table={'overflowX': 'auto'},
    )

    # ── Color legend (by disease category) ─────────────────
    legend = html.Div([
        html.Span("Bar colors by category: ",
                   style={"fontWeight": "600", "fontSize": "0.8rem"}),
        html.Span("\u25cf Cancer ",
                   style={"color": "#e53e3e", "fontSize": "0.8rem"}),
        html.Span("\u25cf Neurodegenerative ",
                   style={"color": "#805ad5", "fontSize": "0.8rem"}),
        html.Span("\u25cf Neuropsychiatric ",
                   style={"color": "#3182ce", "fontSize": "0.8rem"}),
        html.Span("\u25cf Metabolic ",
                   style={"color": "#d69e2e", "fontSize": "0.8rem"}),
    ], className="mb-3")

    return html.Div([
        html.Div([
            html.H5([html.I(className="fas fa-heartbeat me-2"),
                     "Disease Enrichment Summary"], className="cm-card-title"),
            summary,
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-chart-bar me-2"),
                     "Enriched Disease Associations"], className="cm-card-title"),
            legend,
            dcc.Graph(id="enrich-disease-bar-chart", figure=fig,
                      config={"displaylogo": False}),
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-table me-2"),
                     "Full Disease Results"], className="cm-card-title"),
            table,
            html.Div(
                dbc.Button(
                    [html.I(className="fas fa-download me-2"),
                     "Download Disease Results CSV"],
                    id="enrich-disease-download-btn",
                    className="cm-btn-secondary mt-3",
                    color="danger",
                    outline=True,
                ),
                className="text-end",
            ),
        ], className="cm-card"),
    ])


@callback(
    Output("enrich-disease-download", "data"),
    Input("enrich-disease-download-btn", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def download_disease_enrichment(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return no_update
    try:
        fdr_thresh = fdr_thresh or 0.25
        score_thresh = score_thresh or 0.3
        df = pd.read_json(io.StringIO(data), orient='split')

        # Filter by score
        score_col = None
        for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
            if c in df.columns:
                score_col = c
                break
        if score_col:
            df = df[df[score_col] >= score_thresh]

        # Extract metabolite identifiers
        metabolite_ids = set()
        for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break
        for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break

        if not metabolite_ids:
            return no_update

        from app.services.enrichment_service import run_disease_enrichment
        results = run_disease_enrichment(list(metabolite_ids), fdr_thresh)
        if results.empty:
            return no_update
        return dict(content=results.to_csv(index=False),
                    filename="disease_enrichment.csv")
    except Exception:
        return no_update


# ---------------------------------------------------------------------------
# Microbe Enrichment callback
# ---------------------------------------------------------------------------

@callback(
    Output("enrich-microbe-results", "children"),
    Input("enrich-run-microbe", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def run_microbe_enrichment_callback(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return html.Div(
            [html.I(className="fas fa-info-circle me-2"),
             "Please load data first, then click 'Microbe Enrichment'."],
            className="cm-empty-state p-4",
        )

    fdr_thresh = fdr_thresh or 0.25
    score_thresh = score_thresh or 0.3

    df = pd.read_json(io.StringIO(data), orient='split')

    # Filter by min prediction score
    score_col = None
    for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
        if c in df.columns:
            score_col = c
            break
    if score_col:
        df = df[df[score_col] >= score_thresh]

    # Extract metabolite identifiers
    metabolite_ids = set()
    for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break
    for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break

    if not metabolite_ids:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             "No metabolite identifiers found in your data. "
             "Expected columns: Metabolite, HMDB ID, or Metabolite_Name."],
            className="alert alert-warning",
        )

    try:
        from app.services.enrichment_service import run_microbe_enrichment
        results = run_microbe_enrichment(list(metabolite_ids), fdr_thresh)
    except Exception as e:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             f"Microbe enrichment failed: {e}"],
            className="alert alert-danger",
        )

    if results.empty:
        return html.Div(
            [html.I(className="fas fa-search me-2"),
             "No enriched microbe associations found. This may be due to "
             "limited overlap between your metabolites and the MMI database."],
            className="cm-empty-state p-4",
        )

    # Display results
    sig = results[results['FDR'] <= fdr_thresh]
    display_df = sig if not sig.empty else results.head(20)
    top_chart = display_df.head(20).sort_values('FDR', ascending=False)

    # ── Bar chart ──────────────────────────────────────────
    org_colors = {'human': '#3182ce', 'mouse': '#38a169'}
    colors = [org_colors.get(org, '#a0aec0') for org in top_chart['Organism']]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=top_chart['Microbe_Name'].str[:55],
        x=top_chart['Fold_Enrichment'],
        orientation='h',
        marker_color=colors,
        text=[f"FDR={f:.2e}" for f in top_chart['FDR']],
        textposition='outside',
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Fold Enrichment: %{x:.2f}<br>"
            "%{text}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        height=max(350, len(top_chart) * 35),
        margin=dict(l=300, r=100, t=30, b=50),
        xaxis_title="Fold Enrichment",
        yaxis=dict(automargin=True),
        font=dict(family="Arial", size=11),
        plot_bgcolor="white",
    )
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

    # ── Summary stats ──────────────────────────────────────
    n_total = len(results)
    n_sig = len(sig)

    summary = dbc.Row([
        dbc.Col(html.Div([
            html.Div(str(n_total), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Microbes Tested", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(n_sig), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700",
                            "color": "#38a169" if n_sig > 0 else "#a0aec0"}),
            html.Div(f"Significant (FDR \u2264 {fdr_thresh})",
                     className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{results['Fold_Enrichment'].max():.1f}\u00d7",
                     className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Max Fold Enrichment", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(len(metabolite_ids)), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Input Metabolites", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
    ], className="mb-4")

    # ── Results table ──────────────────────────────────────
    display_cols = [c for c in results.columns if c != 'Significant']
    table = dash_table.DataTable(
        id="enrich-microbe-results-table",
        data=results.round(4).to_dict('records'),
        columns=[{'name': c.replace('_', ' '), 'id': c} for c in display_cols],
        sort_action='native',
        filter_action='native',
        page_size=20,
        style_cell={
            'textAlign': 'left', 'fontSize': '0.85rem',
            'padding': '6px 10px',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_header={
            'fontWeight': '600', 'backgroundColor': '#276749',
            'color': 'white',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#f0fff4'},
            {
                'if': {
                    'filter_query': f'{{FDR}} <= {fdr_thresh}',
                    'column_id': 'FDR',
                },
                'color': '#38a169',
                'fontWeight': '600',
            },
        ],
        style_table={'overflowX': 'auto'},
    )

    # ── Color legend ───────────────────────────────────────
    legend = html.Div([
        html.Span("Bar colors by organism: ",
                   style={"fontWeight": "600", "fontSize": "0.8rem"}),
        html.Span("\u25cf Human ",
                   style={"color": "#3182ce", "fontSize": "0.8rem"}),
        html.Span("\u25cf Mouse ",
                   style={"color": "#38a169", "fontSize": "0.8rem"}),
    ], className="mb-3")

    return html.Div([
        html.Div([
            html.H5([html.I(className="fas fa-bacterium me-2"),
                     "Microbe Enrichment Summary"], className="cm-card-title"),
            summary,
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-chart-bar me-2"),
                     "Enriched Gut Microbes"], className="cm-card-title"),
            legend,
            dcc.Graph(id="enrich-microbe-bar-chart", figure=fig,
                      config={"displaylogo": False}),
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-table me-2"),
                     "Full Microbe Results"], className="cm-card-title"),
            table,
            html.Div(
                dbc.Button(
                    [html.I(className="fas fa-download me-2"),
                     "Download Microbe Results CSV"],
                    id="enrich-microbe-download-btn",
                    className="cm-btn-secondary mt-3",
                    color="success",
                    outline=True,
                ),
                className="text-end",
            ),
        ], className="cm-card"),
    ])


@callback(
    Output("enrich-microbe-download", "data"),
    Input("enrich-microbe-download-btn", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def download_microbe_enrichment(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return no_update
    try:
        fdr_thresh = fdr_thresh or 0.25
        score_thresh = score_thresh or 0.3
        df = pd.read_json(io.StringIO(data), orient='split')

        score_col = None
        for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
            if c in df.columns:
                score_col = c
                break
        if score_col:
            df = df[df[score_col] >= score_thresh]

        metabolite_ids = set()
        for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break
        for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break

        if not metabolite_ids:
            return no_update

        from app.services.enrichment_service import run_microbe_enrichment
        results = run_microbe_enrichment(list(metabolite_ids), fdr_thresh)
        if results.empty:
            return no_update
        return dict(content=results.to_csv(index=False),
                    filename="microbe_enrichment.csv")
    except Exception:
        return no_update


# ---------------------------------------------------------------------------
# Drug Enrichment callback
# ---------------------------------------------------------------------------

@callback(
    Output("enrich-drug-results", "children"),
    Input("enrich-run-drug", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def run_drug_enrichment_callback(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return html.Div(
            [html.I(className="fas fa-info-circle me-2"),
             "Please load data first, then click 'Drug Enrichment'."],
            className="cm-empty-state p-4",
        )

    fdr_thresh = fdr_thresh or 0.25
    score_thresh = score_thresh or 0.3

    df = pd.read_json(io.StringIO(data), orient='split')

    # Filter by min prediction score
    score_col = None
    for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
        if c in df.columns:
            score_col = c
            break
    if score_col:
        df = df[df[score_col] >= score_thresh]

    # Extract metabolite identifiers
    metabolite_ids = set()
    for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break
    for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break

    if not metabolite_ids:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             "No metabolite identifiers found in your data. "
             "Expected columns: Metabolite, HMDB ID, or Metabolite_Name."],
            className="alert alert-warning",
        )

    try:
        from app.services.enrichment_service import run_drug_enrichment
        results = run_drug_enrichment(list(metabolite_ids), fdr_thresh)
    except Exception as e:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             f"Drug enrichment failed: {e}"],
            className="alert alert-danger",
        )

    if results.empty:
        return html.Div(
            [html.I(className="fas fa-search me-2"),
             "No enriched drug associations found. This may be due to "
             "limited overlap between your metabolites and the MDrI database."],
            className="cm-empty-state p-4",
        )

    # Display results
    sig = results[results['FDR'] <= fdr_thresh]
    display_df = sig if not sig.empty else results.head(20)
    top_chart = display_df.head(20).sort_values('FDR', ascending=False)

    # ── Bar chart ──────────────────────────────────────────
    type_colors = {'Pharmacokinetic': '#319795', 'Pharmacodynamic': '#2b6cb0'}
    colors = [type_colors.get(t, '#a0aec0') for t in top_chart['Interaction_Type']]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=top_chart['Drug_Name'].str[:55],
        x=top_chart['Fold_Enrichment'],
        orientation='h',
        marker_color=colors,
        text=[f"FDR={f:.2e}" for f in top_chart['FDR']],
        textposition='outside',
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Fold Enrichment: %{x:.2f}<br>"
            "%{text}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        height=max(350, len(top_chart) * 35),
        margin=dict(l=300, r=100, t=30, b=50),
        xaxis_title="Fold Enrichment",
        yaxis=dict(automargin=True),
        font=dict(family="Arial", size=11),
        plot_bgcolor="white",
    )
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

    # ── Summary stats ──────────────────────────────────────
    n_total = len(results)
    n_sig = len(sig)

    summary = dbc.Row([
        dbc.Col(html.Div([
            html.Div(str(n_total), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Drugs Tested", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(n_sig), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700",
                            "color": "#319795" if n_sig > 0 else "#a0aec0"}),
            html.Div(f"Significant (FDR \u2264 {fdr_thresh})",
                     className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{results['Fold_Enrichment'].max():.1f}\u00d7",
                     className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Max Fold Enrichment", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(len(metabolite_ids)), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Input Metabolites", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
    ], className="mb-4")

    # ── Results table ──────────────────────────────────────
    display_cols = [c for c in results.columns if c != 'Significant']
    table = dash_table.DataTable(
        id="enrich-drug-results-table",
        data=results.round(4).to_dict('records'),
        columns=[{'name': c.replace('_', ' '), 'id': c} for c in display_cols],
        sort_action='native',
        filter_action='native',
        page_size=20,
        style_cell={
            'textAlign': 'left', 'fontSize': '0.85rem',
            'padding': '6px 10px',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_header={
            'fontWeight': '600', 'backgroundColor': '#285e61',
            'color': 'white',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#e6fffa'},
            {
                'if': {
                    'filter_query': f'{{FDR}} <= {fdr_thresh}',
                    'column_id': 'FDR',
                },
                'color': '#319795',
                'fontWeight': '600',
            },
        ],
        style_table={'overflowX': 'auto'},
    )

    # ── Color legend ───────────────────────────────────────
    legend = html.Div([
        html.Span("Bar colors by interaction type: ",
                   style={"fontWeight": "600", "fontSize": "0.8rem"}),
        html.Span("\u25cf Pharmacokinetic ",
                   style={"color": "#319795", "fontSize": "0.8rem"}),
        html.Span("\u25cf Pharmacodynamic ",
                   style={"color": "#2b6cb0", "fontSize": "0.8rem"}),
    ], className="mb-3")

    return html.Div([
        html.Div([
            html.H5([html.I(className="fas fa-pills me-2"),
                     "Drug Enrichment Summary"], className="cm-card-title"),
            summary,
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-chart-bar me-2"),
                     "Enriched Drug Associations"], className="cm-card-title"),
            legend,
            dcc.Graph(id="enrich-drug-bar-chart", figure=fig,
                      config={"displaylogo": False}),
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-table me-2"),
                     "Full Drug Results"], className="cm-card-title"),
            table,
            html.Div(
                dbc.Button(
                    [html.I(className="fas fa-download me-2"),
                     "Download Drug Results CSV"],
                    id="enrich-drug-download-btn",
                    className="cm-btn-secondary mt-3",
                    color="info",
                    outline=True,
                ),
                className="text-end",
            ),
        ], className="cm-card"),
    ])


@callback(
    Output("enrich-drug-download", "data"),
    Input("enrich-drug-download-btn", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def download_drug_enrichment(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return no_update
    try:
        fdr_thresh = fdr_thresh or 0.25
        score_thresh = score_thresh or 0.3
        df = pd.read_json(io.StringIO(data), orient='split')

        score_col = None
        for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
            if c in df.columns:
                score_col = c
                break
        if score_col:
            df = df[df[score_col] >= score_thresh]

        metabolite_ids = set()
        for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break
        for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break

        if not metabolite_ids:
            return no_update

        from app.services.enrichment_service import run_drug_enrichment
        results = run_drug_enrichment(list(metabolite_ids), fdr_thresh)
        if results.empty:
            return no_update
        return dict(content=results.to_csv(index=False),
                    filename="drug_enrichment.csv")
    except Exception:
        return no_update


# ---------------------------------------------------------------------------
# Gene Enrichment callback (MGI)
# ---------------------------------------------------------------------------

@callback(
    Output("enrich-gene-results", "children"),
    Input("enrich-run-gene", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def run_gene_enrichment_callback(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return html.Div(
            [html.I(className="fas fa-info-circle me-2"),
             "Please load data first, then click 'Gene Enrichment'."],
            className="cm-empty-state p-4",
        )

    fdr_thresh = fdr_thresh or 0.25
    score_thresh = score_thresh or 0.3

    df = pd.read_json(io.StringIO(data), orient='split')

    # Filter by min prediction score
    score_col = None
    for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
        if c in df.columns:
            score_col = c
            break
    if score_col:
        df = df[df[score_col] >= score_thresh]

    # Extract metabolite identifiers
    metabolite_ids = set()
    for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break
    for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break

    if not metabolite_ids:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             "No metabolite identifiers found in your data. "
             "Expected columns: Metabolite, HMDB ID, or Metabolite_Name."],
            className="alert alert-warning",
        )

    try:
        from app.services.enrichment_service import run_gene_enrichment
        results = run_gene_enrichment(list(metabolite_ids), fdr_thresh)
    except Exception as e:
        return html.Div(
            [html.I(className="fas fa-exclamation-triangle me-2"),
             f"Gene enrichment failed: {e}"],
            className="alert alert-danger",
        )

    if results.empty:
        return html.Div(
            [html.I(className="fas fa-search me-2"),
             "No enriched gene associations found. This may be due to "
             "limited overlap between your metabolites and the MGI database."],
            className="cm-empty-state p-4",
        )

    # Display results
    sig = results[results['FDR'] <= fdr_thresh]
    display_df = sig if not sig.empty else results.head(20)
    top_chart = display_df.head(20).sort_values('FDR', ascending=False)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=top_chart['Gene_Symbol'].str[:40],
        x=top_chart['Fold_Enrichment'],
        orientation='h',
        marker_color='#d69e2e',
        text=[f"FDR={f:.2e}" for f in top_chart['FDR']],
        textposition='outside',
        hovertemplate=(
            "<b>%{y}</b><br>"
            "Fold Enrichment: %{x:.2f}<br>"
            "%{text}<br>"
            "<extra></extra>"
        ),
    ))
    fig.update_layout(
        height=max(350, len(top_chart) * 35),
        margin=dict(l=200, r=100, t=30, b=50),
        xaxis_title="Fold Enrichment",
        yaxis=dict(automargin=True),
        font=dict(family="Arial", size=11),
        plot_bgcolor="white",
    )
    fig.add_vline(x=1, line_dash="dash", line_color="gray", opacity=0.5)

    n_total = len(results)
    n_sig = len(sig)

    summary = dbc.Row([
        dbc.Col(html.Div([
            html.Div(str(n_total), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Genes Tested", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(n_sig), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700",
                            "color": "#d69e2e" if n_sig > 0 else "#a0aec0"}),
            html.Div(f"Significant (FDR \u2264 {fdr_thresh})",
                     className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{results['Fold_Enrichment'].max():.1f}\u00d7",
                     className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Max Fold Enrichment", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(str(len(metabolite_ids)), className="stat-value",
                     style={"fontSize": "1.5rem", "fontWeight": "700"}),
            html.Div("Input Metabolites", className="stat-label text-muted",
                     style={"fontSize": "0.8rem"}),
        ], className="text-center"), md=3),
    ], className="mb-4")

    display_cols = [c for c in results.columns if c != 'Significant']
    table = dash_table.DataTable(
        id="enrich-gene-results-table",
        data=results.round(4).to_dict('records'),
        columns=[{'name': c.replace('_', ' '), 'id': c} for c in display_cols],
        sort_action='native',
        filter_action='native',
        page_size=20,
        style_cell={
            'textAlign': 'left', 'fontSize': '0.85rem',
            'padding': '6px 10px',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_header={
            'fontWeight': '600', 'backgroundColor': '#744210',
            'color': 'white',
            'fontFamily': 'Arial, Helvetica, sans-serif',
        },
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#fffff0'},
            {
                'if': {
                    'filter_query': f'{{FDR}} <= {fdr_thresh}',
                    'column_id': 'FDR',
                },
                'color': '#d69e2e',
                'fontWeight': '600',
            },
        ],
        style_table={'overflowX': 'auto'},
    )

    return html.Div([
        html.Div([
            html.H5([html.I(className="fas fa-dna me-2"),
                     "Gene Enrichment Summary (MGI)"], className="cm-card-title"),
            summary,
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-chart-bar me-2"),
                     "Enriched Gene Associations"], className="cm-card-title"),
            dcc.Graph(id="enrich-gene-bar-chart", figure=fig,
                      config={"displaylogo": False}),
        ], className="cm-card mb-4"),

        html.Div([
            html.H5([html.I(className="fas fa-table me-2"),
                     "Full Gene Results"], className="cm-card-title"),
            table,
            html.Div(
                dbc.Button(
                    [html.I(className="fas fa-download me-2"),
                     "Download Gene Results CSV"],
                    id="enrich-gene-download-btn",
                    className="cm-btn-secondary mt-3",
                    color="warning",
                    outline=True,
                ),
                className="text-end",
            ),
        ], className="cm-card"),
    ])


@callback(
    Output("enrich-gene-download", "data"),
    Input("enrich-gene-download-btn", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def download_gene_enrichment(n, data, fdr_thresh, score_thresh):
    if not n or not data:
        return no_update
    try:
        fdr_thresh = fdr_thresh or 0.25
        score_thresh = score_thresh or 0.3
        df = pd.read_json(io.StringIO(data), orient='split')

        score_col = None
        for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
            if c in df.columns:
                score_col = c
                break
        if score_col:
            df = df[df[score_col] >= score_thresh]

        metabolite_ids = set()
        for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break
        for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
            if c in df.columns:
                metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
                break

        if not metabolite_ids:
            return no_update

        from app.services.enrichment_service import run_gene_enrichment
        results = run_gene_enrichment(list(metabolite_ids), fdr_thresh)
        if results.empty:
            return no_update
        return dict(content=results.to_csv(index=False),
                    filename="gene_enrichment.csv")
    except Exception:
        return no_update


# ---------------------------------------------------------------------------
# Unified PDF Report callback
# ---------------------------------------------------------------------------

@callback(
    Output("enrich-pdf-download", "data"),
    Input("enrich-pdf-btn", "n_clicks"),
    State("enrich-data-store", "data"),
    State("enrich-fdr", "value"),
    State("enrich-score", "value"),
    prevent_initial_call=True,
)
def generate_pdf_report(n, data, fdr_thresh, score_thresh):
    """Run all 4 enrichment analyses and generate a unified PDF report."""
    if not n or not data:
        return no_update

    import base64 as b64

    fdr_thresh = fdr_thresh or 0.25
    score_thresh = score_thresh or 0.3
    df = pd.read_json(io.StringIO(data), orient='split')

    # Filter by min prediction score
    score_col = None
    for c in ['Prediction Score', 'prediction_score', 'Score', 'score']:
        if c in df.columns:
            score_col = c
            break
    if score_col:
        df = df[df[score_col] >= score_thresh]

    # Extract metabolite + protein identifiers
    metabolite_ids = set()
    for c in ['Metabolite', 'Metabolite Name', 'Metabolite_Name', 'metabolite']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break
    for c in ['HMDB ID', 'HMDB_ID', 'hmdb_id']:
        if c in df.columns:
            metabolite_ids.update(df[c].dropna().astype(str).str.strip().tolist())
            break

    # Run all enrichments
    pathway_results = None
    disease_results = None
    microbe_results = None
    drug_results = None
    gene_results = None

    try:
        from app.services.enrichment_service import run_enrichment
        pathway_results = run_enrichment(df, fdr_threshold=fdr_thresh, score_threshold=score_thresh)
    except Exception:
        pass

    if metabolite_ids:
        try:
            from app.services.enrichment_service import run_disease_enrichment
            disease_results = run_disease_enrichment(list(metabolite_ids), fdr_thresh)
        except Exception:
            pass
        try:
            from app.services.enrichment_service import run_microbe_enrichment
            microbe_results = run_microbe_enrichment(list(metabolite_ids), fdr_thresh)
        except Exception:
            pass
        try:
            from app.services.enrichment_service import run_drug_enrichment
            drug_results = run_drug_enrichment(list(metabolite_ids), fdr_thresh)
        except Exception:
            pass
        try:
            from app.services.enrichment_service import run_gene_enrichment
            gene_results = run_gene_enrichment(list(metabolite_ids), fdr_thresh)
        except Exception:
            pass

    # Generate PDF
    try:
        from app.services.pdf_report_service import generate_enrichment_pdf
        pdf_bytes = generate_enrichment_pdf(
            pathway_results=pathway_results,
            disease_results=disease_results,
            microbe_results=microbe_results,
            drug_results=drug_results,
            input_summary={
                "n_metabolites": len(metabolite_ids),
                "n_proteins": len(df),
                "source": "User upload / Disease selection",
            },
            fdr_threshold=fdr_thresh,
        )
        encoded = b64.b64encode(pdf_bytes).decode()
        return dict(
            content=encoded,
            filename="cormet_enrichment_report.pdf",
            base64=True,
        )
    except Exception as e:
        logger.error(f"PDF generation failed: {e}")
        return no_update


# Export for main.py routing (navbar/footer handled globally)
page_content = layout