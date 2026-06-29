"""
API Documentation Page, interactive reference for CoreMet REST API endpoints.

Route: /api-docs
"""

from dash import html, dcc
import dash_bootstrap_components as dbc

from components.page_header import make_page_header

# ---------------------------------------------------------------------------
# Endpoint definitions (single source of truth)
# ---------------------------------------------------------------------------

BASE = "/api/v1"

ENDPOINTS = [
    {
        "method": "GET",
        "path": f"{BASE}/health",
        "summary": "Health check",
        "desc": "Returns server status and version information.",
        "params": [],
        "example_response": '{"status": "healthy", "version": "3.0.0"}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/species",
        "summary": "List supported species",
        "desc": "Returns all organism species available in the MPI database with interaction counts.",
        "params": [],
        "example_response": '[{"species": "Homo sapiens", "count": 8494}, ...]',
    },
    {
        "method": "POST",
        "path": f"{BASE}/predict",
        "summary": "Submit prediction job",
        "desc": "Submit an asynchronous interaction prediction job. Returns a job ID for polling results.",
        "params": [
            ("metabolites", "body (JSON)", "List of HMDB IDs or metabolite names", True),
            ("species", "body (JSON)", 'Target organism (default: "Homo sapiens")', False),
            ("model", "body (JSON)", 'Model type: mpi, mdi, mmi, mdri (default: "mpi")', False),
            ("top_k", "body (JSON)", "Number of top predictions to return (default: 50)", False),
        ],
        "example_response": '{"job_id": "abc123", "status": "queued"}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/results/<job_id>",
        "summary": "Get prediction results",
        "desc": "Retrieve results for a previously submitted prediction job.",
        "params": [
            ("job_id", "path", "Job identifier returned by /predict", True),
        ],
        "example_response": '{"job_id": "abc123", "status": "completed", "results": [...]}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/database/search",
        "summary": "Search MPI database",
        "desc": "Query the metabolite–protein interaction database with filters.",
        "params": [
            ("q", "query", "Search term (metabolite name or HMDB ID)", True),
            ("species", "query", 'Filter by organism (default: all)', False),
            ("limit", "query", "Maximum results (default: 100, max: 1000)", False),
        ],
        "example_response": '{"results": [...], "total": 42}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/mmi/stats",
        "summary": "MMI database statistics",
        "desc": "Returns summary statistics for the Metabolite–Microbe Interaction database.",
        "params": [],
        "example_response": '{"total": 83149, "metabolites": 1043, "microbes": 1262}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/mmi/search",
        "summary": "Search MMI database",
        "desc": "Search metabolite–microbe interactions by metabolite name, HMDB ID, or microbe name.",
        "params": [
            ("q", "query", "Search term", True),
            ("limit", "query", "Maximum results (default: 100)", False),
        ],
        "example_response": '{"results": [...], "total": 15}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/mdri/stats",
        "summary": "MDrI database statistics",
        "desc": "Returns summary statistics for the Metabolite–Drug Interaction database.",
        "params": [],
        "example_response": '{"total": 3500, "metabolites": 312, "drugs": 97}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/mdri/search",
        "summary": "Search MDrI database",
        "desc": "Search metabolite–drug interactions by metabolite name, HMDB ID, or drug name.",
        "params": [
            ("q", "query", "Search term", True),
            ("limit", "query", "Maximum results (default: 100)", False),
        ],
        "example_response": '{"results": [...], "total": 8}',
    },
    {
        "method": "GET",
        "path": f"{BASE}/export/metabolite",
        "summary": "Export metabolite data as CSV",
        "desc": "Download interaction data for a single metabolite as a CSV file.",
        "params": [
            ("id", "query", "HMDB ID (e.g. HMDB0000122)", False),
            ("name", "query", "Metabolite name (alternative to id)", False),
            ("db", "query", 'Database: mpi, mei, mdi, mmi, mdri, or all (default: "all")', False),
        ],
        "example_response": "(CSV file download)",
    },    {
        "method": "GET",
        "path": f"{BASE}/autocomplete",
        "summary": "Metabolite name autocomplete",
        "desc": "Return top matching metabolites for a partial name or HMDB ID prefix. Used by the navbar search dropdown.",
        "params": [
            ("q", "query", "Partial metabolite name or HMDB ID (min 2 chars)", True),
            ("limit", "query", "Maximum results (default: 10, max: 30)", False),
        ],
        "example_response": '[{"label": "L-Glutamic acid (HMDB0000148)", "name": "L-Glutamic acid", "hmdb_id": "HMDB0000148"}, ...]',
    },]

METHOD_COLORS = {"GET": "success", "POST": "primary", "PUT": "warning", "DELETE": "danger"}


def _build_endpoint_card(ep):
    """Build an accordion item for one API endpoint."""
    method = ep["method"]
    badge_color = METHOD_COLORS.get(method, "secondary")

    # Parameter table
    param_rows = []
    for p in ep.get("params", []):
        name, location, desc, required = p[0], p[1], p[2], p[3] if len(p) > 3 else False
        param_rows.append(html.Tr([
            html.Td(html.Code(name)),
            html.Td(location, style={"fontSize": "0.8rem", "color": "#718096"}),
            html.Td(desc, style={"fontSize": "0.85rem"}),
            html.Td(dbc.Badge("required", color="danger", pill=True) if required
                     else html.Span("optional", style={"color": "#a0aec0", "fontSize": "0.8rem"})),
        ]))

    param_section = []
    if param_rows:
        param_section = [
            html.H6("Parameters", className="mt-3 mb-2", style={"fontWeight": "600"}),
            dbc.Table([
                html.Thead(html.Tr([
                    html.Th("Name"), html.Th("In"), html.Th("Description"), html.Th(""),
                ])),
                html.Tbody(param_rows),
            ], bordered=True, size="sm", className="mb-3",
               style={"fontSize": "0.85rem"}),
        ]

    # Example response
    example_section = []
    if ep.get("example_response"):
        example_section = [
            html.H6("Example Response", className="mt-2 mb-1", style={"fontWeight": "600"}),
            html.Pre(
                html.Code(ep["example_response"]),
                style={"backgroundColor": "#f7fafc", "padding": "10px",
                       "borderRadius": "6px", "fontSize": "0.8rem",
                       "border": "1px solid #e2e8f0", "overflowX": "auto"},
            ),
        ]

    # Example curl
    curl_cmd = f"curl -s 'http://localhost:8080{ep['path']}'"
    if method == "POST":
        curl_cmd = (f'curl -X POST "http://localhost:8080{ep["path"]}" '
                    f'-H "Content-Type: application/json" '
                    f'-d \'{{"metabolites": ["HMDB0000122"]}}\'')

    return dbc.AccordionItem([
        html.P(ep["desc"], className="text-muted", style={"fontSize": "0.9rem"}),
        *param_section,
        *example_section,
        html.H6("Try it", className="mt-3 mb-1", style={"fontWeight": "600"}),
        html.Pre(
            html.Code(curl_cmd),
            style={"backgroundColor": "#2d3748", "color": "#e2e8f0",
                   "padding": "10px", "borderRadius": "6px",
                   "fontSize": "0.78rem", "overflowX": "auto"},
        ),
    ], title=f"{method}  {ep['path']} , {ep['summary']}")


# ---------------------------------------------------------------------------
# Page layout
# ---------------------------------------------------------------------------

_cards = [_build_endpoint_card(ep) for ep in ENDPOINTS]

page_content = html.Div([
    make_page_header(
        "API Documentation",
        "Programmatic access to CoreMet databases and prediction engine.",
        [("Home", "/home"), ("API Docs", None)],
    ),
    html.Div([
        # Overview
        html.Div([
            html.H5([
                html.I(className="fas fa-plug me-2", style={"color": "#3182ce"}),
                "REST API Overview",
            ], className="cm-card-title mb-3"),
            html.P([
                "CoreMet provides a REST API at ",
                html.Code("http://your-host/api/v1/"),
                " for programmatic access to all databases and the prediction engine. "
                "All endpoints return JSON (except CSV exports). "
                "Rate limit: ",
                html.Strong("100 requests/hour"),
                " per IP.",
            ], style={"fontSize": "0.9rem"}),
            dbc.Row([
                dbc.Col([
                    html.Div([
                    html.Div("11", style={"fontWeight": "700", "fontSize": "1.5rem",
                                              "color": "#3182ce"}),
                        html.Div("Endpoints", style={"fontSize": "0.8rem", "color": "#718096"}),
                    ], className="text-center"),
                ], width=3),
                dbc.Col([
                    html.Div([
                        html.Div("JSON", style={"fontWeight": "700", "fontSize": "1.5rem",
                                                "color": "#38a169"}),
                        html.Div("Response Format", style={"fontSize": "0.8rem", "color": "#718096"}),
                    ], className="text-center"),
                ], width=3),
                dbc.Col([
                    html.Div([
                        html.Div("No Auth", style={"fontWeight": "700", "fontSize": "1.5rem",
                                                   "color": "#d69e2e"}),
                        html.Div("Required", style={"fontSize": "0.8rem", "color": "#718096"}),
                    ], className="text-center"),
                ], width=3),
                dbc.Col([
                    html.Div([
                        html.Div("CORS", style={"fontWeight": "700", "fontSize": "1.5rem",
                                                "color": "#319795"}),
                        html.Div("Enabled", style={"fontSize": "0.8rem", "color": "#718096"}),
                    ], className="text-center"),
                ], width=3),
            ], className="mt-3"),
        ], className="cm-card mb-4"),

        # Endpoints accordion
        html.Div([
            html.H5([
                html.I(className="fas fa-list me-2", style={"color": "#38a169"}),
                "Endpoints",
            ], className="cm-card-title mb-3"),
            dbc.Accordion(_cards, start_collapsed=True),
        ], className="cm-card mb-4"),
    ]),
], className="cm-page-container")
