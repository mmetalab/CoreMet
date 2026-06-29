"""
Predict page, 3-step wizard: Input → Configure → Results
Supports 4 interaction models: MPI, MDI, MMI, MDrI with model selection.
"""

import base64
import io
import json
import logging

from dash import dcc, html, Input, Output, State, callback, dash_table, no_update, ctx
import dash_bootstrap_components as dbc
import pandas as pd

from components.page_header import make_page_header

logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────

GENOME_DICT = {
    'Homo sapiens': "Homo_sapiens",
    'Mus musculus': "Mus_musculus",
    'Rattus norvegicus': "Rattus_norvegicus",
    'Escherichia coli': "Escherichia_coli",
    'Bos taurus': "Bos_taurus",
    'Pseudomonas aeruginosa': "Pseudomonas_aeruginosa",
    'Arabidopsis thaliana': "Arabidopsis_thaliana",
    'Saccharomyces cerevisiae': "Saccharomyces_cerevisiae",
    'Drosophila melanogaster': "Drosophila_melanogaster",
    'Caenorhabditis elegans': "Caenorhabditis_elegans",
}

# Model type definitions
MODEL_TYPES = {
    'mpi': {
        'label': 'Metabolite–Protein (MPI)',
        'icon': 'fas fa-dna',
        'color': '#3182ce',
        'target_label': 'Proteins',
        'target_icon': 'fas fa-microscope',
        'description': 'Predict novel metabolite–protein interactions (model trained on 38,061 curated edges)',
        'placeholder': 'Paste proteins, CSV: UniprotID,Protein Name,Gene Name,Organism,Sequence\nOr one ID/name per line (P04217 or ALB)',
        'met_placeholder': 'Paste metabolites, CSV: Name,HMDB_ID,SMILES\nOr one per line: HMDB0000122 or Glucose',
        'example_met_file': 'example_metabolites_hcc.csv',
        'example_target_file': 'example_proteins_hcc.csv',
    },
    'mdi': {
        'label': 'Metabolite–Disease (MDI)',
        'icon': 'fas fa-disease',
        'color': '#e53e3e',
        'target_label': 'Diseases',
        'target_icon': 'fas fa-heartbeat',
        'description': 'Predict novel metabolite–disease associations (model trained on 82,882 curated edges)',
        'placeholder': 'Paste diseases, CSV: Name,MeSH_ID,Category\nOr one per line: MESH:D006973 or Hypertension',
        'met_placeholder': 'Paste metabolites, CSV: Name,HMDB_ID,SMILES\nOr one per line: HMDB0001877 or Valproic acid',
        'example_met_file': 'example_metabolites_mdi.csv',
        'example_target_file': 'example_diseases.csv',
    },
    'mmi': {
        'label': 'Metabolite–Microbe (MMI)',
        'icon': 'fas fa-bacterium',
        'color': '#38a169',
        'target_label': 'Microbes',
        'target_icon': 'fas fa-bacterium',
        'description': 'Predict novel metabolite–microbe interactions (model trained on 77,605 curated edges)',
        'placeholder': 'Paste microbes, CSV: Name,Taxonomy_ID,Rank\nOr one per line: 511145 or Escherichia coli',
        'met_placeholder': 'Paste metabolites, CSV: Name,HMDB_ID,SMILES\nOr one per line: HMDB ID or name',
        'example_met_file': 'example_metabolites_mmi.csv',
        'example_target_file': 'example_microbes.csv',
    },
    'mdri': {
        'label': 'Metabolite–Drug (MDrI)',
        'icon': 'fas fa-pills',
        'color': '#805ad5',
        'target_label': 'Drugs',
        'target_icon': 'fas fa-pills',
        'description': 'Predict metabolite–drug interactions using GraphSAGE (3,500 curated)',
        'placeholder': 'Paste drugs, CSV: Name,DrugBank_ID,SMILES\nOr one per line: DB00331 or Metformin',
        'met_placeholder': 'Paste metabolites, CSV: Name,HMDB_ID,SMILES\nOr one per line: HMDB ID or name',
        'example_met_file': 'example_metabolites_mdri.csv',
        'example_target_file': 'example_drugs.csv',
    },
    'mgi': {
        'label': 'Metabolite–Gene (MGI)',
        'icon': 'fas fa-dna',
        'color': '#d69e2e',
        'target_label': 'Genes',
        'target_icon': 'fas fa-dna',
        'description': 'Predict metabolite–gene interactions using GraphSAGE (1,658,745 from CTD)',
        'placeholder': 'Paste genes, one per line: CYP1A2 or TP53\nOr CSV: Gene_Symbol,Organism',
        'met_placeholder': 'Paste metabolites, CSV: Name,HMDB_ID,SMILES\nOr one per line: HMDB ID or name',
        'example_met_file': 'example_metabolites_mgi.csv',
        'example_target_file': 'example_genes.csv',
    },
    'mgwas': {
        'label': 'Metabolite–SNP (mGWAS)',
        'icon': 'fas fa-map-marker-alt',
        'color': '#319795',
        'target_label': 'SNPs',
        'target_icon': 'fas fa-map-marker-alt',
        'description': 'Predict metabolite–SNP associations using GraphSAGE (44,344 from GWAS Catalog)',
        'placeholder': 'Paste SNPs, one per line: rs1260326 or rs174547',
        'met_placeholder': 'Paste metabolites, CSV: Name,HMDB_ID,SMILES\nOr one per line: HMDB ID or name',
        'example_met_file': 'example_metabolites_mgwas.csv',
        'example_target_file': 'example_snps.csv',
    },
}


# ── Helper functions ──────────────────────────────────────────────────

def make_stepper(active_step=1):
    """Create the stepper indicator."""
    steps = [("1", "Input Data"), ("2", "Configure"), ("3", "Results")]
    children = []
    for i, (num, label) in enumerate(steps):
        step_num = i + 1
        if step_num < active_step:
            cls = "completed"
        elif step_num == active_step:
            cls = "active"
        else:
            cls = ""
        children.append(
            html.Div(
                [
                    html.Div(
                        html.I(className="fas fa-check") if cls == "completed" else num,
                        className="cm-stepper-circle",
                    ),
                    html.Span(label, className="cm-stepper-label"),
                ],
                className=f"cm-stepper-step {cls}",
            )
        )
        if i < len(steps) - 1:
            conn_cls = "completed" if step_num < active_step else ""
            children.append(html.Div(className=f"cm-stepper-connector {conn_cls}"))
    return html.Div(children, className="cm-stepper")


def make_upload_card(title, upload_id, textarea_id, preview_id, status_id,
                     icon_class, download_id=None, example_btn_id=None,
                     placeholder=None):
    """Create an upload card with drag-and-drop, paste area, and inline example load."""
    children = [
        html.H5([html.I(className=f"{icon_class} me-2"), title], className="cm-card-title"),
        dcc.Upload(
            id=upload_id,
            children=html.Div(
                [
                    html.Div(html.I(className="fas fa-cloud-upload-alt"), className="upload-icon"),
                    html.Div("Drag & drop CSV or XLSX, or click to browse", className="upload-text"),
                    html.Div("Supported: .csv, .xlsx", className="upload-hint"),
                ],
            ),
            className="cm-upload-zone mb-3",
            multiple=False,
        ),
    ]
    # Inline buttons: example load + download
    btn_row = []
    if example_btn_id:
        btn_row.append(
            dbc.Button(
                [html.I(className="fas fa-flask me-1"), f"Load Example {title}"],
                id=example_btn_id,
                className="cm-btn-secondary me-2",
                size="sm",
            )
        )
    if download_id:
        btn_row.append(
            dbc.Button(
                [html.I(className="fas fa-download me-1"), "Download Example CSV"],
                id=download_id,
                className="cm-btn-secondary",
                size="sm",
            )
        )
    if btn_row:
        children.append(html.Div(btn_row, className="d-flex mb-3"))

    children.extend([
        html.Div("- or paste data manually -", className="text-center text-muted mb-2",
                  style={"fontSize": "0.8rem"}),
        dbc.Textarea(
            id=textarea_id,
            placeholder=placeholder or f"Paste {title.lower()} data (one per line, or comma-separated)...",
            rows=3,
            className="mb-3",
            style={"fontSize": "0.85rem"},
        ),
        html.Div(id=preview_id),
        html.Div(id=status_id, className="mt-2"),
    ])
    return html.Div(children, className="cm-card")


# ── Layout ────────────────────────────────────────────────────────────

layout = html.Div(
    [
        html.Div(
            [
                make_page_header(
                    "Predict Interactions (beta)",
                    "An optional research feature built on the CoreMet graph: upload metabolites "
                    "and a graph-neural-network model scores candidate interactions across types. "
                    "Predictions are hypotheses for follow-up, the curated database remains the core resource.",
                    [("Home", "/home"), ("Predict", None)],
                ),

                # Action buttons
                dbc.Row(
                    [
                        dbc.Col([
                            dbc.Button(
                                [html.I(className="fas fa-question-circle me-2"), "View Tutorial"],
                                href="/help",
                                className="cm-btn-secondary",
                            ),
                        ]),
                    ],
                    className="mb-4",
                ),

                # Stepper
                html.Div(id="predict-stepper", children=make_stepper(1)),

                # Hidden stores
                dcc.Store(id="predict-model-type", data="mpi"),
                dcc.Store(id="predict-metabolite-store"),
                dcc.Store(id="predict-protein-store"),   # legacy name, holds *target* data
                dcc.Store(id="predict-results-store"),
                dcc.Store(id="predict-current-step", data=1),
                dcc.Store(id="predict-job-id"),
                dcc.Download(id="predict-example-download"),
                # Hidden dummy for removed global example button
                html.Div(id="predict-load-example", style={"display": "none"}),

                # ═══ Step 1: Input ═══
                html.Div(
                    [
                        # Model type selector
                        html.Div(
                            [
                                html.H5(
                                    [html.I(className="fas fa-project-diagram me-2"),
                                     "Select Interaction Model"],
                                    className="cm-card-title",
                                ),
                                html.P(
                                    "Choose the type of interaction to predict.",
                                    style={"color": "#718096", "fontSize": "0.9rem",
                                           "marginBottom": "16px"},
                                ),
                                dbc.RadioItems(
                                    id="predict-model-selector",
                                    options=[
                                        {
                                            "label": html.Span([
                                                html.I(className=f"{v['icon']} me-2",
                                                       style={"color": v['color'], "width": "20px"}),
                                                html.Span(v['label'], style={"fontWeight": "500"}),
                                                html.Span(f", {v['description']}",
                                                          style={"color": "#a0aec0",
                                                                 "fontSize": "0.85rem"}),
                                            ]),
                                            "value": k,
                                        }
                                        for k, v in MODEL_TYPES.items()
                                    ],
                                    value="mpi",
                                    className="cm-radio-group",
                                    inline=False,
                                    labelClassName="d-flex align-items-center py-2",
                                ),
                            ],
                            className="cm-card mb-4",
                        ),

                        # Upload cards
                        dbc.Row(
                            [
                                dbc.Col(
                                    html.Div(
                                        id="predict-met-card-container",
                                        children=make_upload_card(
                                            "Metabolites",
                                            "predict-met-upload",
                                            "predict-met-textarea",
                                            "predict-met-preview",
                                            "predict-met-status",
                                            "fas fa-flask",
                                            download_id="predict-download-met-example",
                                            example_btn_id="predict-load-met-example",
                                            placeholder=MODEL_TYPES['mpi']['met_placeholder'],
                                        ),
                                    ),
                                    md=6,
                                ),
                                dbc.Col(
                                    html.Div(
                                        id="predict-target-card-container",
                                        children=make_upload_card(
                                            "Proteins",
                                            "predict-prot-upload",
                                            "predict-prot-textarea",
                                            "predict-prot-preview",
                                            "predict-prot-status",
                                            "fas fa-microscope",
                                            download_id="predict-download-prot-example",
                                            example_btn_id="predict-load-target-example",
                                            placeholder=MODEL_TYPES['mpi']['placeholder'],
                                        ),
                                    ),
                                    md=6,
                                ),
                            ],
                            className="g-3",
                        ),

                        # Target info note
                        html.Div(
                            id="predict-target-note",
                            children=html.Div(style={"display": "none"}),
                        ),

                        html.Div(
                            dbc.Button(
                                [html.Span("Next Step"),
                                 html.I(className="fas fa-arrow-right ms-2")],
                                id="predict-next-1",
                                className="cm-btn-primary mt-4",
                                disabled=True,
                            ),
                            className="text-end",
                        ),
                    ],
                    id="predict-step-1",
                ),

                # ═══ Step 2: Configure ═══
                html.Div(
                    [
                        html.Div(
                            [
                                # Organism selector (MPI only)
                                html.Div(
                                    [
                                        html.H5("Select Organism", className="cm-card-title"),
                                        dbc.Select(
                                            id="predict-organism",
                                            options=[{"label": k, "value": k} for k in GENOME_DICT],
                                            value="Homo sapiens",
                                            className="mb-3",
                                        ),
                                    ],
                                    id="predict-organism-section",
                                ),
                                # Foundation model toggle
                                html.Div(
                                    [
                                        dbc.Checkbox(
                                            id="predict-use-foundation",
                                            label="Use Foundation Model (CoreMet-FM), unified cross-type embeddings",
                                            value=False,
                                            className="mb-2",
                                        ),
                                        html.Small(
                                            "Foundation model learns shared metabolite embeddings across all 6 interaction types. "
                                            "Best for metabolites with diverse interaction profiles.",
                                            className="text-muted",
                                        ),
                                    ],
                                    id="predict-foundation-section",
                                    className="mb-3 p-2",
                                    style={"background": "#f7fafc", "borderRadius": "8px"},
                                ),
                                # Model info summary
                                html.Div(id="predict-config-summary"),
                                # Feature visualization
                                html.H5("Feature Space Visualization",
                                         className="cm-card-title mt-4"),
                                dcc.Graph(
                                    id="predict-umap-plot",
                                    config={"displaylogo": False},
                                    style={"height": "400px"},
                                ),
                            ],
                            className="cm-card",
                        ),
                        html.Div(
                            [
                                dbc.Button(
                                    [html.I(className="fas fa-arrow-left me-2"), "Back"],
                                    id="predict-back-2",
                                    className="cm-btn-secondary me-2",
                                ),
                                dbc.Button(
                                    [html.I(className="fas fa-play me-2"), "Run Prediction"],
                                    id="predict-run",
                                    className="cm-btn-primary",
                                ),
                            ],
                            className="text-end mt-4",
                        ),
                    ],
                    id="predict-step-2",
                    style={"display": "none"},
                ),

                # ═══ Step 3: Results ═══
                html.Div(
                    [
                        # Job info bar
                        html.Div(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col([
                                            html.Span("Job ID: ",
                                                       style={"fontWeight": "600"}),
                                            html.Span(id="predict-job-id-display"),
                                        ], width="auto"),
                                        dbc.Col([
                                            html.Span(id="predict-status-badge"),
                                        ], width="auto"),
                                        dbc.Col(
                                            html.Span(id="predict-count-display"),
                                            width="auto",
                                        ),
                                        dbc.Col(
                                            html.Span(id="predict-model-badge"),
                                            width="auto",
                                        ),
                                    ],
                                    className="align-items-center",
                                ),
                            ],
                            className="cm-card mb-3",
                            style={"padding": "12px 24px"},
                        ),

                        # Summary stats badges
                        html.Div(id="predict-summary-badges", className="mb-3"),

                        # Results table
                        html.Div(
                            [
                                html.H5("Prediction Results", className="cm-card-title"),
                                dcc.Loading(
                                    id="predict-results-loading",
                                    type="circle",
                                    color="#1a365d",
                                    children=html.Div(id="predict-results-table"),
                                ),
                            ],
                            className="cm-card mb-3",
                        ),

                        # Export bar
                        html.Div(
                            [
                                html.H5("Export", className="cm-card-title"),
                                dbc.ButtonGroup(
                                    [
                                        dbc.Button(
                                            [html.I(className="fas fa-file-csv me-1"), "CSV"],
                                            id="predict-export-csv",
                                            className="cm-btn-secondary",
                                        ),
                                        dbc.Button(
                                            [html.I(className="fas fa-file-code me-1"), "JSON"],
                                            id="predict-export-json",
                                            className="cm-btn-secondary",
                                        ),
                                        dbc.Button(
                                            [html.I(className="fas fa-project-diagram me-1"),
                                             "GraphML"],
                                            id="predict-export-graphml",
                                            className="cm-btn-secondary",
                                        ),
                                        dbc.Button(
                                            [html.I(className="fas fa-share-alt me-1"), "SIF"],
                                            id="predict-export-sif",
                                            className="cm-btn-secondary",
                                        ),
                                    ],
                                ),
                                dcc.Download(id="predict-download"),
                            ],
                            className="cm-card mb-3",
                        ),

                        # New prediction
                        html.Div(
                            dbc.Button(
                                [html.I(className="fas fa-redo me-2"), "New Prediction"],
                                id="predict-new",
                                className="cm-btn-accent",
                            ),
                            className="text-center mt-3",
                        ),
                    ],
                    id="predict-step-3",
                    style={"display": "none"},
                ),
            ],
            className="cm-page-container",
        ),
    ]
)


# ══════════════════════════════════════════════════════════════════════
# Callbacks
# ══════════════════════════════════════════════════════════════════════

def parse_upload(contents, filename):
    """Parse uploaded file contents."""
    if contents is None:
        return None
    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        if 'csv' in filename.lower():
            return pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename.lower():
            return pd.read_excel(io.BytesIO(decoded))
    except Exception as e:
        print(f"Upload parse error: {e}")
    return None


# ── Model type syncing ────────────────────────────────────────────────

@callback(
    Output("predict-model-type", "data"),
    Output("predict-target-note", "children"),
    Output("predict-organism-section", "style"),
    Output("predict-config-summary", "children"),
    Input("predict-model-selector", "value"),
)
def update_model_type(model_key):
    """Sync model type store and update dynamic UI elements."""
    info = MODEL_TYPES.get(model_key, MODEL_TYPES['mpi'])

    # Target note: show for non-MPI models
    if model_key == 'mpi':
        note = html.Div(style={"display": "none"})
        organism_style = {}
    else:
        note = html.Div(
            [
                html.I(className="fas fa-info-circle me-2",
                       style={"color": "#3182ce"}),
                html.Span(
                    f"Target ({info['target_label'].lower()}) list is optional, "
                    f"leave empty to predict against all known "
                    f"{info['target_label'].lower()} in the database.",
                    style={"fontSize": "0.85rem", "color": "#718096"},
                ),
            ],
            className="mt-2 px-2",
        )
        organism_style = {"display": "none"}

    # Config summary for Step 2
    summary = html.Div(
        [
            html.H5(
                [html.I(className=f"{info['icon']} me-2"), "Model Configuration"],
                className="cm-card-title",
            ),
            dbc.Row([
                dbc.Col([
                    html.Div("Interaction Type",
                             style={"fontSize": "0.8rem", "color": "#a0aec0"}),
                    html.Div(info['label'], style={"fontWeight": "600"}),
                ], sm=4),
                dbc.Col([
                    html.Div("Architecture",
                             style={"fontSize": "0.8rem", "color": "#a0aec0"}),
                    html.Div("GraphSAGE + MLP", style={"fontWeight": "600"}),
                ], sm=4),
                dbc.Col([
                    html.Div("Source",
                             style={"fontSize": "0.8rem", "color": "#a0aec0"}),
                    html.Div(f"Metabolites → {info['target_label']}",
                             style={"fontWeight": "600"}),
                ], sm=4),
            ], className="mt-2"),
        ],
        className="mb-3",
    )

    return model_key, note, organism_style, summary


# ── Dynamic target card ──────────────────────────────────────────────

@callback(
    Output("predict-target-card-container", "children"),
    Output("predict-met-card-container", "children"),
    Input("predict-model-selector", "value"),
)
def update_upload_cards(model_key):
    """Update both upload cards based on selected model type."""
    info = MODEL_TYPES.get(model_key, MODEL_TYPES['mpi'])

    # Metabolite card, always present, placeholder varies by model
    met_card = make_upload_card(
        "Metabolites",
        "predict-met-upload",
        "predict-met-textarea",
        "predict-met-preview",
        "predict-met-status",
        "fas fa-flask",
        download_id="predict-download-met-example",
        example_btn_id="predict-load-met-example",
        placeholder=info.get('met_placeholder', 'Paste metabolites, one per line'),
    )

    # Target card, label and icons change by model
    target_card = make_upload_card(
        info['target_label'],
        "predict-prot-upload",
        "predict-prot-textarea",
        "predict-prot-preview",
        "predict-prot-status",
        info['target_icon'],
        download_id="predict-download-prot-example",
        example_btn_id="predict-load-target-example",
        placeholder=info.get('placeholder'),
    )
    return target_card, met_card


# ── Metabolite input ─────────────────────────────────────────────────

@callback(
    Output("predict-metabolite-store", "data"),
    Output("predict-met-preview", "children"),
    Output("predict-met-status", "children"),
    Input("predict-met-upload", "contents"),
    State("predict-met-upload", "filename"),
    Input("predict-met-textarea", "value"),
    Input("predict-load-met-example", "n_clicks"),
    State("predict-model-type", "data"),
    prevent_initial_call=True,
)
def handle_metabolite_input(contents, filename, textarea_val, example_clicks,
                            model_type):
    trigger = ctx.triggered_id
    df = None
    model_type = model_type or "mpi"
    info = MODEL_TYPES.get(model_type, MODEL_TYPES['mpi'])

    if trigger == "predict-load-met-example" and example_clicks:
        from pathlib import Path
        example_file = info.get('example_met_file', 'example_metabolites_hcc.csv')
        example_path = (Path(__file__).parent.parent
                        / "data" / "examples" / example_file)
        if not example_path.exists():
            example_path = (Path(__file__).parent.parent
                            / "data" / "example_metabolite_data.csv")
        if example_path.exists():
            df = pd.read_csv(example_path)
    elif trigger == "predict-load-example":
        # Legacy global button (hidden), fall back to HCC
        from pathlib import Path
        example_path = (Path(__file__).parent.parent
                        / "data" / "examples" / "example_metabolites_hcc.csv")
        if example_path.exists():
            df = pd.read_csv(example_path)
    elif trigger == "predict-met-upload" and contents:
        df = parse_upload(contents, filename)
    elif trigger == "predict-met-textarea" and textarea_val:
        rows = []
        for line in textarea_val.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 3:
                # Name, HMDB_ID, SMILES
                rows.append({'Metabolite Name': parts[0],
                             'HMDB ID': parts[1], 'SMILES': parts[2]})
            elif len(parts) == 2:
                # Could be Name,ID or ID,Name
                if parts[0].startswith('HMDB'):
                    rows.append({'HMDB ID': parts[0], 'Metabolite Name': parts[1]})
                elif parts[1].startswith('HMDB'):
                    rows.append({'Metabolite Name': parts[0], 'HMDB ID': parts[1]})
                else:
                    rows.append({'Metabolite Name': parts[0], 'HMDB ID': parts[1]})
            elif len(parts) == 1:
                val = parts[0]
                if val.startswith('HMDB'):
                    rows.append({'HMDB ID': val})
                else:
                    rows.append({'Metabolite Name': val})
        if rows:
            df = pd.DataFrame(rows)

    if df is None or df.empty:
        return no_update, no_update, no_update

    preview = dash_table.DataTable(
        data=df.head(3).to_dict('records'),
        columns=[{'name': c, 'id': c} for c in df.columns],
        style_cell={'textAlign': 'left', 'fontSize': '0.8rem', 'padding': '4px 8px',
                    'maxWidth': '180px', 'overflow': 'hidden',
                    'textOverflow': 'ellipsis',
                    'fontFamily': 'Arial, Helvetica, sans-serif'},
        style_header={'fontWeight': '600', 'backgroundColor': '#edf2f7',
                      'fontFamily': 'Arial, Helvetica, sans-serif'},
        tooltip_data=[
            {col: {'value': str(val), 'type': 'markdown'}
             for col, val in row.items()}
            for row in df.head(3).to_dict('records')
        ],
        tooltip_duration=None,
    )
    status = html.Div(
        [html.I(className="fas fa-check-circle me-1 text-success"),
         f"{len(df)} metabolites loaded"],
        style={"color": "#38a169", "fontWeight": "500"},
    )
    return df.to_json(orient='split'), preview, status


# ── Target input (proteins / diseases / microbes / drugs) ────────────

@callback(
    Output("predict-protein-store", "data"),
    Output("predict-prot-preview", "children"),
    Output("predict-prot-status", "children"),
    Input("predict-prot-upload", "contents"),
    State("predict-prot-upload", "filename"),
    Input("predict-prot-textarea", "value"),
    Input("predict-load-target-example", "n_clicks"),
    State("predict-model-type", "data"),
    prevent_initial_call=True,
)
def handle_target_input(contents, filename, textarea_val, example_clicks,
                        model_type):
    trigger = ctx.triggered_id
    df = None
    model_type = model_type or "mpi"
    info = MODEL_TYPES.get(model_type, MODEL_TYPES['mpi'])

    if trigger == "predict-load-target-example" and example_clicks:
        from pathlib import Path
        example_file = info.get('example_target_file')
        if example_file:
            example_path = (Path(__file__).parent.parent
                            / "data" / "examples" / example_file)
            if example_path.exists():
                df = pd.read_csv(example_path)
    elif trigger == "predict-load-example":
        # Legacy global button (hidden), fall back to MPI example
        from pathlib import Path
        if model_type == 'mpi':
            example_path = (Path(__file__).parent.parent
                            / "data" / "examples" / "example_proteins_hcc.csv")
            if example_path.exists():
                df = pd.read_csv(example_path)
    elif trigger == "predict-prot-upload" and contents:
        df = parse_upload(contents, filename)
    elif trigger == "predict-prot-textarea" and textarea_val:
        rows = []
        for line in textarea_val.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(',')]
            if model_type == 'mpi':
                if len(parts) >= 5:
                    rows.append({
                        'UniprotID': parts[0], 'Protein Name': parts[1],
                        'Gene Name': parts[2], 'Organism': parts[3],
                        'Sequence': ','.join(parts[4:]),
                    })
                elif len(parts) == 3:
                    rows.append({
                        'UniprotID': parts[0], 'Protein Name': parts[1],
                        'Gene Name': parts[2],
                    })
                elif len(parts) == 1:
                    val = parts[0]
                    if val.startswith(('P', 'Q', 'O', 'A')) and len(val) == 6:
                        rows.append({'UniprotID': val})
                    else:
                        rows.append({'Gene Name': val})
                else:
                    rows.append({'UniprotID': parts[0],
                                 'Protein Name': parts[1] if len(parts) > 1 else parts[0]})
            elif model_type == 'mdi':
                if len(parts) >= 3:
                    # Name,MeSH_ID,Category
                    rows.append({'Disease_Name': parts[0], 'Disease_ID': parts[1],
                                 'Category': parts[2]})
                elif len(parts) == 2:
                    if parts[0].startswith('MESH:'):
                        rows.append({'Disease_ID': parts[0], 'Disease_Name': parts[1]})
                    elif parts[1].startswith('MESH:'):
                        rows.append({'Disease_Name': parts[0], 'Disease_ID': parts[1]})
                    else:
                        rows.append({'Disease_Name': parts[0], 'Disease_ID': parts[1]})
                else:
                    val = parts[0]
                    if val.startswith('MESH:'):
                        rows.append({'Disease_ID': val})
                    else:
                        rows.append({'Disease_Name': val})
            elif model_type == 'mmi':
                if len(parts) >= 3:
                    # Name,Taxonomy_ID,Rank
                    rows.append({'Microbe_Name': parts[0], 'Taxonomy_ID': parts[1],
                                 'Rank': parts[2]})
                elif len(parts) == 2:
                    if parts[0].isdigit():
                        rows.append({'Taxonomy_ID': parts[0], 'Microbe_Name': parts[1]})
                    elif parts[1].isdigit():
                        rows.append({'Microbe_Name': parts[0], 'Taxonomy_ID': parts[1]})
                    else:
                        rows.append({'Microbe_Name': parts[0], 'Taxonomy_ID': parts[1]})
                else:
                    val = parts[0]
                    if val.isdigit():
                        rows.append({'Taxonomy_ID': val})
                    else:
                        rows.append({'Microbe_Name': val})
            elif model_type == 'mdri':
                if len(parts) >= 3:
                    # Name,DrugBank_ID,SMILES
                    rows.append({'Drug_Name': parts[0], 'DrugBank_ID': parts[1],
                                 'SMILES': ','.join(parts[2:])})
                elif len(parts) == 2:
                    if parts[0].startswith('DB'):
                        rows.append({'DrugBank_ID': parts[0], 'Drug_Name': parts[1]})
                    elif parts[1].startswith('DB'):
                        rows.append({'Drug_Name': parts[0], 'DrugBank_ID': parts[1]})
                    else:
                        rows.append({'Drug_Name': parts[0], 'DrugBank_ID': parts[1]})
                else:
                    val = parts[0]
                    if val.startswith('DB'):
                        rows.append({'DrugBank_ID': val})
                    else:
                        rows.append({'Drug_Name': val})
            elif model_type == 'mgi':
                if len(parts) >= 2:
                    # Gene_Symbol,Organism
                    rows.append({'Gene_Symbol': parts[0], 'Organism': parts[1]})
                else:
                    rows.append({'Gene_Symbol': parts[0]})
            elif model_type == 'mgwas':
                val = parts[0]
                rows.append({'rsID': val})
            else:
                rows.append({'ID': parts[0]})
        if rows:
            df = pd.DataFrame(rows)

    if df is None or df.empty:
        return no_update, no_update, no_update

    preview = dash_table.DataTable(
        data=df.head(3).to_dict('records'),
        columns=[{'name': c, 'id': c} for c in df.columns],
        style_cell={'textAlign': 'left', 'fontSize': '0.8rem', 'padding': '4px 8px',
                    'maxWidth': '150px', 'overflow': 'hidden',
                    'textOverflow': 'ellipsis',
                    'fontFamily': 'Arial, Helvetica, sans-serif'},
        style_header={'fontWeight': '600', 'backgroundColor': '#edf2f7',
                      'fontFamily': 'Arial, Helvetica, sans-serif'},
    )
    label = info['target_label'].lower()
    status = html.Div(
        [html.I(className="fas fa-check-circle me-1 text-success"),
         f"{len(df)} {label} loaded"],
        style={"color": "#38a169", "fontWeight": "500"},
    )
    return df.to_json(orient='split'), preview, status


# ── Enable Next button ────────────────────────────────────────────────

@callback(
    Output("predict-next-1", "disabled"),
    Input("predict-metabolite-store", "data"),
    Input("predict-protein-store", "data"),
    Input("predict-model-type", "data"),
)
def enable_next_button(met_data, prot_data, model_type):
    """For MPI: need both metabolites AND proteins.
    For others: only metabolites required."""
    if not met_data:
        return True
    if model_type == 'mpi' and not prot_data:
        return True
    return False


# ── Navigation ────────────────────────────────────────────────────────

@callback(
    Output("predict-step-1", "style"),
    Output("predict-step-2", "style"),
    Output("predict-step-3", "style"),
    Output("predict-stepper", "children"),
    Input("predict-next-1", "n_clicks"),
    Input("predict-back-2", "n_clicks"),
    Input("predict-run", "n_clicks"),
    Input("predict-new", "n_clicks"),
    prevent_initial_call=True,
)
def navigate_steps(next1, back2, run, new):
    trigger = ctx.triggered_id
    if trigger == "predict-next-1":
        return {"display": "none"}, {}, {"display": "none"}, make_stepper(2)
    elif trigger == "predict-back-2":
        return {}, {"display": "none"}, {"display": "none"}, make_stepper(1)
    elif trigger == "predict-run":
        return {"display": "none"}, {"display": "none"}, {}, make_stepper(3)
    elif trigger == "predict-new":
        return {}, {"display": "none"}, {"display": "none"}, make_stepper(1)
    return {}, {"display": "none"}, {"display": "none"}, make_stepper(1)


# ── Run prediction ────────────────────────────────────────────────────

@callback(
    Output("predict-results-store", "data"),
    Output("predict-results-table", "children"),
    Output("predict-job-id-display", "children"),
    Output("predict-status-badge", "children"),
    Output("predict-count-display", "children"),
    Output("predict-summary-badges", "children"),
    Output("predict-model-badge", "children"),
    Input("predict-run", "n_clicks"),
    State("predict-metabolite-store", "data"),
    State("predict-protein-store", "data"),
    State("predict-organism", "value"),
    State("predict-model-type", "data"),
    State("predict-use-foundation", "value"),
    prevent_initial_call=True,
)
def run_prediction(n_clicks, met_data, prot_data, organism, model_type, use_foundation):
    if not n_clicks or not met_data:
        return (no_update,) * 7

    model_type = model_type or "mpi"
    info = MODEL_TYPES.get(model_type, MODEL_TYPES['mpi'])

    df_met = pd.read_json(io.StringIO(met_data), orient='split')
    df_target = (pd.read_json(io.StringIO(prot_data), orient='split')
                 if prot_data else None)

    # Create job
    from app.services.job_service import create_job, update_job
    job_id = create_job(met_data, prot_data or "{}", organism)

    try:
        if model_type == 'mpi' and not use_foundation:
            results_df = _predict_mpi(df_met, df_target, organism)
        else:
            results_df = _predict_multi(model_type, df_met, df_target,
                                         use_foundation=use_foundation)
    except Exception as e:
        logger.error(f"Prediction error ({model_type}): {e}")
        import traceback
        traceback.print_exc()
        results_df = pd.DataFrame({
            'Source': [], 'Target': [],
            'Prediction Score': [], 'Existing': [],
        })

    update_job(job_id, "complete", results_df.to_json(orient='split'))

    # ── Build results table ───────────────────────────────────────────
    if model_type == 'mpi':
        preferred_order = [
            'Metabolite', 'Protein Name', 'EC_Number', 'Enzyme_Name',
            'Interaction_Type', 'Prediction Score', 'Existing',
            'similarity_score', 'Protein',
        ]
    elif model_type == 'mdi':
        preferred_order = [
            'Metabolite', 'Target', 'Disease_Name',
            'Prediction Score', 'Existing',
        ]
    elif model_type == 'mmi':
        preferred_order = [
            'Metabolite', 'Target', 'Microbe_Name',
            'Prediction Score', 'Existing',
        ]
    elif model_type == 'mdri':
        preferred_order = [
            'Metabolite', 'Target', 'Drug_Name',
            'Prediction Score', 'Existing',
        ]
    elif model_type == 'mgi':
        preferred_order = [
            'Metabolite', 'Target', 'Gene_Symbol', 'Organism',
            'Prediction Score', 'Existing',
        ]
    elif model_type == 'mgwas':
        preferred_order = [
            'Metabolite', 'rsID', 'Chromosome', 'Mapped_Gene',
            'Prediction Score', 'Existing',
        ]
    else:
        preferred_order = list(results_df.columns)

    display_cols = [c for c in preferred_order if c in results_df.columns]
    display_cols += [c for c in results_df.columns
                     if c not in display_cols and c != 'Is_Enzyme']
    results_df_display = results_df[display_cols] if display_cols else results_df

    col_name_map = {
        'EC_Number': 'EC Number',
        'Enzyme_Name': 'Enzyme',
        'Interaction_Type': 'Type',
        'similarity_score': 'Similarity',
        'Disease_Name': 'Disease',
        'Microbe_Name': 'Microbe',
        'Drug_Name': 'Drug',
    }

    table = dash_table.DataTable(
        data=results_df_display.to_dict('records'),
        columns=[{'name': col_name_map.get(c, c), 'id': c}
                 for c in results_df_display.columns],
        sort_action='native',
        filter_action='native',
        page_size=20,
        style_cell={'textAlign': 'left', 'fontSize': '0.85rem', 'padding': '8px',
                    'fontFamily': 'Arial, Helvetica, sans-serif'},
        style_header={'fontWeight': '600', 'backgroundColor': '#1a365d',
                      'color': 'white',
                      'fontFamily': 'Arial, Helvetica, sans-serif'},
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#f7fafc'},
        ],
    )

    # ── Summary badges ────────────────────────────────────────────────
    n_total = len(results_df)
    if 'Prediction Score' in results_df.columns:
        scores = pd.to_numeric(results_df['Prediction Score'], errors='coerce')
        n_high = int((scores >= 0.7).sum())
        n_med = int(((scores >= 0.4) & (scores < 0.7)).sum())
        n_low = int((scores < 0.4).sum())
    else:
        n_high = n_med = n_low = 0

    badge_items = [
        dbc.Col(html.Span(f"{n_total} interactions",
                          className="cm-badge cm-badge-source"), width="auto"),
        dbc.Col(html.Span(f"{n_high} high confidence",
                          className="cm-badge cm-badge-high"), width="auto"),
        dbc.Col(html.Span(f"{n_med} medium",
                          className="cm-badge cm-badge-medium"), width="auto"),
        dbc.Col(html.Span(f"{n_low} low",
                          className="cm-badge cm-badge-low"), width="auto"),
    ]
    badges = dbc.Row(badge_items, className="g-2")

    status_badge = html.Span("Complete", className="cm-badge cm-badge-high")
    model_badge = html.Span(
        info['label'],
        className="cm-badge",
        style={"backgroundColor": f"{info['color']}1a",
               "color": info['color']},
    )

    return (
        results_df.to_json(orient='split'),
        table,
        job_id,
        status_badge,
        f"{n_total} interactions",
        badges,
        model_badge,
    )


def _predict_mpi(df_met, df_prot, organism):
    """Run MPI prediction using existing PredictionService."""
    from app.services.prediction_service import PredictionService
    if not hasattr(_predict_mpi, '_service'):
        _predict_mpi._service = PredictionService()
    return _predict_mpi._service.predict_interactions(df_met, df_prot, organism)


def _predict_multi(model_type, df_met, df_target, use_foundation=False):
    """Run MDI/MMI/MDrI/MGI/mGWAS/MPI prediction using MultiPredictionService."""
    from app.services.multi_prediction_service import MultiPredictionService

    if not hasattr(_predict_multi, '_service'):
        _predict_multi._service = MultiPredictionService()

    svc = _predict_multi._service

    # Extract metabolite IDs
    met_ids = []
    for col in ['HMDB ID', 'HMDB_ID', 'hmdb_id', 'Metabolite ID']:
        if col in df_met.columns:
            met_ids = df_met[col].dropna().tolist()
            break
    if not met_ids and 'Metabolite Name' in df_met.columns:
        met_ids = df_met['Metabolite Name'].dropna().tolist()
    if not met_ids:
        met_ids = df_met.iloc[:, 0].dropna().tolist()

    # Extract target IDs (if provided)
    target_ids = None
    if df_target is not None and not df_target.empty:
        id_col_map = {
            'mdi': ['Disease_ID', 'Disease ID', 'MESH_ID'],
            'mmi': ['Taxonomy_ID', 'Taxonomy ID', 'Tax_ID'],
            'mdri': ['DrugBank_ID', 'DrugBank ID', 'Drug_ID'],
            'mgi': ['Gene_Symbol', 'Gene Symbol', 'Gene_ID'],
            'mgwas': ['rsID', 'rs_ID', 'SNP_ID', 'SNP'],
        }
        for col in id_col_map.get(model_type, []):
            if col in df_target.columns:
                target_ids = df_target[col].dropna().tolist()
                break
        if target_ids is None:
            target_ids = df_target.iloc[:, 0].dropna().tolist()

    # Run prediction
    results_df = svc.predict_with_metadata(model_type, met_ids, target_ids)
    return results_df


# ── Export ─────────────────────────────────────────────────────────────

@callback(
    Output("predict-download", "data"),
    Input("predict-export-csv", "n_clicks"),
    Input("predict-export-json", "n_clicks"),
    Input("predict-export-graphml", "n_clicks"),
    Input("predict-export-sif", "n_clicks"),
    State("predict-results-store", "data"),
    prevent_initial_call=True,
)
def handle_export(csv_clicks, json_clicks, graphml_clicks, sif_clicks,
                  results_data):
    if not results_data:
        return no_update
    trigger = ctx.triggered_id
    df = pd.read_json(io.StringIO(results_data), orient='split')

    from app.services.export_service import (export_csv, export_json,
                                             export_graphml, export_sif)

    if trigger == "predict-export-csv":
        return dict(content=export_csv(df), filename="cormet_predictions.csv")
    elif trigger == "predict-export-json":
        return dict(content=export_json(df),
                    filename="cormet_predictions.json")
    elif trigger == "predict-export-graphml":
        return dict(content=export_graphml(df),
                    filename="cormet_predictions.graphml")
    elif trigger == "predict-export-sif":
        return dict(content=export_sif(df), filename="cormet_predictions.sif")
    return no_update


# ── Example download ──────────────────────────────────────────────────

@callback(
    Output("predict-example-download", "data"),
    Input("predict-download-met-example", "n_clicks"),
    Input("predict-download-prot-example", "n_clicks"),
    State("predict-model-type", "data"),
    prevent_initial_call=True,
)
def download_example_data(met_clicks, prot_clicks, model_type):
    # Guard: ignore spurious fires when dynamic components are recreated
    if not met_clicks and not prot_clicks:
        return no_update

    trigger = ctx.triggered_id
    if trigger is None:
        return no_update

    model_type = model_type or "mpi"
    info = MODEL_TYPES.get(model_type, MODEL_TYPES['mpi'])
    from pathlib import Path
    base = Path(__file__).parent.parent / "data" / "examples"

    if trigger == "predict-download-met-example" and met_clicks:
        fname = info.get('example_met_file', 'example_metabolites_hcc.csv')
        path = base / fname
    elif trigger == "predict-download-prot-example" and prot_clicks:
        fname = info.get('example_target_file', 'example_proteins_hcc.csv')
        path = base / fname
    else:
        return no_update

    if not path.exists():
        return no_update
    return dcc.send_file(str(path), filename=fname)


# ── UMAP scatter (Step 2) ─────────────────────────────────────────────

def _morgan_fp(smiles_str, nbits=128):
    """Compute Morgan fingerprint as numpy array, or None if invalid."""
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    if not smiles_str or smiles_str == 'nan':
        return None
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits)
    arr = np.zeros(nbits)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def _aa_composition(seq_str, dim=128):
    """Compute amino-acid composition vector, padded to *dim*."""
    import numpy as np
    if not seq_str or seq_str == 'nan' or len(seq_str) < 5:
        return None
    aa = "ACDEFGHIKLMNPQRSTVWY"
    feat = np.array([seq_str.upper().count(a) / len(seq_str) for a in aa])
    return np.pad(feat, (0, dim - len(feat)))


_DISEASE_CATS = [
    'CTD_curated', 'Cancer', 'Cardiovascular', 'Metabolic',
    'Neurodegenerative', 'Hepatic', 'Respiratory', 'Neuropsychiatric',
    'Renal', 'Autoimmune', 'Infectious', 'Other',
]

def _disease_feature(category_str, dim=128):
    """One-hot category + small noise → padded to *dim*."""
    import numpy as np
    cat = str(category_str).strip() if category_str else 'Other'
    vec = np.zeros(len(_DISEASE_CATS))
    if cat in _DISEASE_CATS:
        vec[_DISEASE_CATS.index(cat)] = 1.0
    else:
        vec[_DISEASE_CATS.index('Other')] = 1.0
    noise = np.random.default_rng(42).normal(0, 0.01, len(vec))
    vec = vec + noise
    return np.pad(vec, (0, dim - len(vec)))


_MICROBE_RANKS = ['strain', 'species', 'genus', 'family', 'order',
                  'class', 'phylum', 'no rank', 'subspecies']

def _microbe_feature(rank_str, dim=128):
    """One-hot rank + small noise → padded to *dim*."""
    import numpy as np
    rank = str(rank_str).strip().lower() if rank_str else 'no rank'
    vec = np.zeros(len(_MICROBE_RANKS))
    if rank in _MICROBE_RANKS:
        vec[_MICROBE_RANKS.index(rank)] = 1.0
    else:
        vec[len(_MICROBE_RANKS) - 2] = 1.0  # 'no rank'
    noise = np.random.default_rng(42).normal(0, 0.01, len(vec))
    vec = vec + noise
    return np.pad(vec, (0, dim - len(vec)))


@callback(
    Output("predict-umap-plot", "figure"),
    Input("predict-next-1", "n_clicks"),
    State("predict-metabolite-store", "data"),
    State("predict-protein-store", "data"),
    State("predict-model-type", "data"),
    prevent_initial_call=True,
)
def render_umap_plot(n_clicks, met_data, prot_data, model_type):
    """Generate UMAP/PCA scatter of metabolite and target features."""
    import numpy as np
    import plotly.graph_objects as go

    model_type = model_type or "mpi"
    info = MODEL_TYPES.get(model_type, MODEL_TYPES['mpi'])

    fig = go.Figure()
    fig.update_layout(
        template="plotly_white",
        font=dict(family="Arial", size=11),
        xaxis_title="Component 1",
        yaxis_title="Component 2",
        height=400,
        margin=dict(l=50, r=20, t=30, b=50),
        legend=dict(orientation="h", yanchor="bottom", y=1.02,
                    xanchor="right", x=1),
    )

    if not met_data:
        fig.add_annotation(
            text="Upload data in Step 1 to see feature space",
            showarrow=False, font=dict(size=14, color="#718096"),
        )
        return fig

    try:
        df_met = pd.read_json(io.StringIO(met_data), orient='split')

        # ── Metabolite features (Morgan FP from SMILES) ──────────────
        met_fps, met_names = [], []
        for _, row in df_met.iterrows():
            fp = _morgan_fp(str(row.get('SMILES', '')))
            if fp is not None:
                met_fps.append(fp)
                met_names.append(
                    str(row.get('Metabolite Name',
                                row.get('HMDB ID', '')))[:20])

        # ── Target features (model-dependent) ────────────────────────
        target_fps, target_names = [], []
        if prot_data:
            df_target = pd.read_json(io.StringIO(prot_data), orient='split')

            if model_type == 'mpi':
                # Protein: AA composition from Sequence
                for _, row in df_target.iterrows():
                    fp = _aa_composition(str(row.get('Sequence', '')))
                    if fp is not None:
                        target_fps.append(fp)
                        target_names.append(
                            str(row.get('Protein Name',
                                        row.get('UniprotID', '')))[:20])

            elif model_type == 'mdi':
                # Disease: Category one-hot
                for _, row in df_target.iterrows():
                    fp = _disease_feature(row.get('Category', 'Other'))
                    target_fps.append(fp)
                    target_names.append(
                        str(row.get('Disease_Name',
                                    row.get('Disease Name',
                                            row.get('Disease_ID', ''))))[:20])

            elif model_type == 'mmi':
                # Microbe: Rank one-hot
                for _, row in df_target.iterrows():
                    fp = _microbe_feature(row.get('Rank', 'no rank'))
                    target_fps.append(fp)
                    target_names.append(
                        str(row.get('Microbe_Name',
                                    row.get('Microbe Name',
                                            row.get('Taxonomy_ID', ''))))[:20])

            elif model_type == 'mdri':
                # Drug: Morgan FP from SMILES
                for _, row in df_target.iterrows():
                    fp = _morgan_fp(str(row.get('SMILES', '')))
                    if fp is not None:
                        target_fps.append(fp)
                        target_names.append(
                            str(row.get('Drug_Name',
                                        row.get('Drug Name',
                                                row.get('DrugBank_ID', ''))))[:20])

            elif model_type == 'mgi':
                # Gene: random feature vector (no structural feature available)
                for _, row in df_target.iterrows():
                    gene = str(row.get('Gene_Symbol', row.get('Gene Symbol', '')))
                    fp = np.random.RandomState(hash(gene) % 2**31).randn(128).astype(np.float32) * 0.01
                    target_fps.append(fp)
                    target_names.append(gene[:20])

            elif model_type == 'mgwas':
                # SNP: random feature vector
                for _, row in df_target.iterrows():
                    snp = str(row.get('rsID', row.get('SNP', '')))
                    fp = np.random.RandomState(hash(snp) % 2**31).randn(128).astype(np.float32) * 0.01
                    target_fps.append(fp)
                    target_names.append(snp[:20])

        total = len(met_fps) + len(target_fps)
        if total < 3:
            fig.add_annotation(
                text="Need ≥ 3 entities with valid features (SMILES / Sequence / metadata) for visualization",
                showarrow=False, font=dict(size=12, color="#718096"),
            )
            return fig

        all_fps = np.array(met_fps + target_fps)

        try:
            import umap
            reducer = umap.UMAP(
                n_components=2, random_state=42,
                n_neighbors=min(15, len(all_fps) - 1),
            )
            coords = reducer.fit_transform(all_fps)
        except Exception:
            from sklearn.decomposition import PCA
            pca = PCA(n_components=2)
            coords = pca.fit_transform(all_fps)

        n_met = len(met_fps)
        fig.add_trace(go.Scatter(
            x=coords[:n_met, 0], y=coords[:n_met, 1],
            mode='markers+text', name='Metabolites',
            text=met_names, textposition='top center',
            textfont=dict(size=8),
            marker=dict(size=10, color='#3182ce', symbol='circle'),
        ))
        if target_fps:
            fig.add_trace(go.Scatter(
                x=coords[n_met:, 0], y=coords[n_met:, 1],
                mode='markers+text', name=info['target_label'],
                text=target_names, textposition='top center',
                textfont=dict(size=8),
                marker=dict(size=10, color=info['color'], symbol='square'),
            ))
    except Exception as e:
        fig.add_annotation(
            text=f"Feature visualization error: {str(e)[:60]}",
            showarrow=False,
        )

    return fig


# Navbar callbacks


# Export for main.py routing
page_content = layout
