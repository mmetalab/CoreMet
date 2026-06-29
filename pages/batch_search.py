"""
Batch Metabolite Search page, paste a list of HMDB IDs or metabolite names,
get a summary table showing interaction counts across all five databases,
with links to individual metabolite detail pages and CSV download.

All callback IDs prefixed with 'batch-'.
"""

import io

import pandas as pd
from dash import dcc, html, dash_table, Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc

from components.page_header import make_page_header

# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

_EXAMPLE_INPUT = (
    "HMDB0000001\n"
    "HMDB0000064\n"
    "HMDB0000122\n"
    "HMDB0000158\n"
    "HMDB0000167\n"
    "L-Glutamic acid\n"
    "Tryptophan"
)

page_content = html.Div([
    make_page_header(
        "Batch Metabolite Search",
        "Look up multiple metabolites at once across all five CoreMet databases.",
    ),
    dbc.Container([
        dbc.Row([
            # Left panel, input
            dbc.Col([
                dbc.Card([
                    dbc.CardBody([
                        html.H5([
                            html.I(className="fas fa-paste me-2"),
                            "Input Metabolites",
                        ], className="mb-3"),
                        html.P(
                            "Enter HMDB IDs or metabolite names, one per line.",
                            className="text-muted mb-2",
                            style={"fontSize": "0.85rem"},
                        ),
                        dbc.Textarea(
                            id="batch-input",
                            placeholder="HMDB0000001\nHMDB0000064\nL-Glutamic acid\n...",
                            value="",
                            rows=10,
                            style={"fontFamily": "monospace", "fontSize": "0.82rem"},
                        ),
                        html.Div([
                            dbc.Button(
                                [html.I(className="fas fa-search me-1"), "Search"],
                                id="batch-submit",
                                color="primary",
                                className="me-2 mt-2",
                            ),
                            dbc.Button(
                                [html.I(className="fas fa-undo me-1"), "Clear"],
                                id="batch-clear",
                                color="secondary",
                                outline=True,
                                className="me-2 mt-2",
                            ),
                            dbc.Button(
                                [html.I(className="fas fa-flask me-1"), "Load Example"],
                                id="batch-example",
                                color="info",
                                outline=True,
                                className="mt-2",
                                size="sm",
                            ),
                        ]),
                    ]),
                ], className="cm-card mb-3"),
            ], md=4),

            # Right panel, results
            dbc.Col([
                dcc.Loading(
                    html.Div(id="batch-results"),
                    type="circle",
                    color="#3182ce",
                ),
            ], md=8),
        ]),
    ], fluid=True, className="px-4 pb-5"),
])


# ---------------------------------------------------------------------------
# Callbacks
# ---------------------------------------------------------------------------

@callback(
    Output("batch-input", "value"),
    Input("batch-example", "n_clicks"),
    prevent_initial_call=True,
)
def load_example(n):
    return _EXAMPLE_INPUT


@callback(
    Output("batch-input", "value", allow_duplicate=True),
    Input("batch-clear", "n_clicks"),
    prevent_initial_call=True,
)
def clear_input(n):
    return ""


@callback(
    Output("batch-results", "children"),
    Input("batch-submit", "n_clicks"),
    State("batch-input", "value"),
    prevent_initial_call=True,
)
def run_batch_search(n_clicks, raw_text):
    if not raw_text or not raw_text.strip():
        return dbc.Alert("Please enter at least one metabolite.", color="warning")

    lines = [l.strip() for l in raw_text.strip().splitlines() if l.strip()]
    if len(lines) > 200:
        return dbc.Alert("Maximum 200 metabolites per batch.", color="danger")

    from pages.metabolite_detail import _lookup_metabolite

    rows = []
    for query in lines:
        is_hmdb = query.upper().startswith("HMDB")
        data = _lookup_metabolite(
            hmdb_id=query if is_hmdb else "",
            name="" if is_hmdb else query,
        )
        hmdb_id = data.get("hmdb_id", "")
        name = data.get("name", query)
        rows.append({
            "Query": query,
            "Metabolite": name if name and name != hmdb_id else query,
            "HMDB ID": hmdb_id,
            "MPI": len(data["mpi"]),
            "MEI": len(data["mei"]),
            "MDI": len(data["mdi"]),
            "MMI": len(data["mmi"]),
            "MDrI": len(data["mdri"]),
            "Total": (len(data["mpi"]) + len(data["mei"]) + len(data["mdi"])
                      + len(data["mmi"]) + len(data["mdri"])),
        })

    df = pd.DataFrame(rows)

    found = df[df["Total"] > 0]
    not_found = df[df["Total"] == 0]

    summary = html.Div([
        html.H5([
            html.I(className="fas fa-table me-2"),
            "Results",
            dbc.Badge(f"{len(found)}/{len(df)} found", className="ms-2",
                      pill=True, color="success" if len(found) > 0 else "secondary"),
        ], className="mb-3"),
    ])

    # Results DataTable
    table = dash_table.DataTable(
        id="batch-results-table",
        columns=[
            {"name": c, "id": c, "presentation": "markdown"}
            if c in ("Metabolite", "HMDB ID") else {"name": c, "id": c}
            for c in ["Query", "Metabolite", "HMDB ID", "MPI", "MEI", "MDI", "MMI", "MDrI", "Total"]
        ],
        data=_make_linked_rows(df),
        style_table={"overflowX": "auto"},
        style_header={
            "backgroundColor": "#2d3748", "color": "white",
            "fontWeight": "bold", "fontSize": "0.82rem",
        },
        style_cell={
            "textAlign": "center", "fontSize": "0.8rem",
            "padding": "6px 10px",
        },
        style_data_conditional=[
            {"if": {"filter_query": "{Total} = 0"},
             "backgroundColor": "#fff5f5", "color": "#c53030"},
            {"if": {"filter_query": "{Total} > 0"},
             "backgroundColor": "#f0fff4"},
        ],
        sort_action="native",
        page_action="native",
        page_size=50,
        markdown_options={"link_target": "_self"},
    )

    # CSV download
    csv_section = html.Div([
        dbc.Button(
            [html.I(className="fas fa-download me-1"), "Download CSV"],
            id="batch-download-btn",
            color="primary",
            size="sm",
            className="mt-2",
        ),
        dcc.Download(id="batch-download-csv"),
    ])

    # Not-found alert
    alert = html.Div()
    if len(not_found) > 0:
        alert = dbc.Alert(
            f"{len(not_found)} metabolite(s) not found: "
            + ", ".join(not_found["Query"].tolist()[:10])
            + ("…" if len(not_found) > 10 else ""),
            color="warning",
            className="mt-2",
            dismissable=True,
        )

    # Store raw CSV in hidden div for download
    csv_store = dcc.Store(id="batch-csv-store", data=df.to_csv(index=False))

    return html.Div([summary, table, csv_section, alert, csv_store])


def _make_linked_rows(df):
    """Add markdown links for Metabolite name and HMDB ID."""
    records = []
    for _, row in df.iterrows():
        r = row.to_dict()
        hmdb = r.get("HMDB ID", "")
        if hmdb:
            r["HMDB ID"] = f"[{hmdb}](/metabolite?id={hmdb})"
            name = r.get("Metabolite", "")
            if name:
                r["Metabolite"] = f"[{name}](/metabolite?id={hmdb})"
        records.append(r)
    return records


@callback(
    Output("batch-download-csv", "data"),
    Input("batch-download-btn", "n_clicks"),
    State("batch-csv-store", "data"),
    prevent_initial_call=True,
)
def download_csv(n, csv_str):
    if not csv_str:
        return no_update
    return dict(content=csv_str, filename="cormet_batch_results.csv")
