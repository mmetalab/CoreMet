"""
Results page, bookmarkable view at /results/<job_id>
"""

from dash import dcc, html, Input, Output, callback, dash_table
import dash_bootstrap_components as dbc
import pandas as pd
import io

from components.page_header import make_page_header


layout = html.Div(
    [
        html.Div(
            [
                make_page_header(
                    "Prediction Results",
                    "View saved prediction results.",
                    [("Home", "/home"), ("Predict", "/predict"), ("Results", None)],
                ),
                dcc.Location(id="results-url", refresh=False),
                html.Div(id="results-content"),
            ],
            className="cm-page-container",
        ),
    ]
)


@callback(
    Output("results-content", "children"),
    Input("results-url", "pathname"),
)
def load_results(pathname):
    if not pathname or not pathname.startswith("/results/"):
        return html.Div(
            [
                html.Div(html.I(className="fas fa-search"), className="empty-icon"),
                html.Div("No job ID specified", className="empty-title"),
                html.Div("Use /results/<job_id> to view saved results.", className="empty-description"),
                dbc.Button("Go to Predict", href="/predict", className="cm-btn-primary"),
            ],
            className="cm-empty-state",
        )

    job_id = pathname.split("/results/")[-1]

    try:
        from app.services.job_service import get_job
        job = get_job(job_id)
    except Exception:
        job = None

    if job is None:
        return html.Div(
            [
                html.Div(html.I(className="fas fa-exclamation-triangle"), className="empty-icon"),
                html.Div("Job not found", className="empty-title"),
                html.Div(f"Job '{job_id}' was not found or has expired.", className="empty-description"),
                dbc.Button("New Prediction", href="/predict", className="cm-btn-primary"),
            ],
            className="cm-empty-state",
        )

    # Parse results
    if job['results']:
        df = pd.read_json(io.StringIO(job['results']), orient='split')
    else:
        df = pd.DataFrame()

    return html.Div(
        [
            # Job info
            html.Div(
                [
                    html.Span(f"Job ID: {job['job_id']}", style={"fontWeight": "600"}),
                    html.Span(f" | Status: {job['status']}", className="ms-3"),
                    html.Span(f" | Organism: {job['organism']}", className="ms-3"),
                    html.Span(f" | Created: {job['created_at']}", className="ms-3"),
                ],
                className="cm-card mb-3",
                style={"padding": "12px 24px", "fontSize": "0.9rem"},
            ),

            # Results table
            html.Div(
                [
                    html.H5(f"{len(df)} Predicted Interactions", className="cm-card-title"),
                    dash_table.DataTable(
                        data=df.to_dict('records') if not df.empty else [],
                        columns=[{'name': c, 'id': c} for c in df.columns] if not df.empty else [],
                        sort_action='native',
                        filter_action='native',
                        page_size=20,
                        style_cell={'textAlign': 'left', 'fontSize': '0.85rem', 'padding': '8px', 'fontFamily': 'Arial, Helvetica, sans-serif'},
                        style_header={'fontWeight': '600', 'backgroundColor': '#1a365d', 'color': 'white', 'fontFamily': 'Arial, Helvetica, sans-serif'},
                        style_data_conditional=[
                            {'if': {'row_index': 'odd'}, 'backgroundColor': '#f7fafc'},
                        ],
                    ),
                ],
                className="cm-card",
            ),
        ]
    )






# Export for main.py routing (navbar/footer handled globally)
page_content = layout
