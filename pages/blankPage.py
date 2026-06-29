from dash import dcc, html
import dash_bootstrap_components as dbc

from components.page_header import make_page_header

layout = html.Div([
    make_page_header(
        "Page Not Found",
        "The page you are looking for does not exist.",
        [("Home", "/home")],
    ),
    html.Div([
        html.Div([
            html.I(className="fas fa-map-signs fa-3x mb-3",
                   style={"color": "#a0aec0"}),
            html.H3("404, Page Not Found", className="mb-2",
                     style={"color": "#2d3748"}),
            html.P(
                "This URL doesn't match any CoreMet page. "
                "Use the navigation bar above or try one of the links below.",
                className="text-muted mb-4",
            ),
            html.Div([
                dcc.Link(dbc.Button([html.I(className="fas fa-home me-2"), "Home"],
                                    color="primary", className="me-2"), href="/home"),
                dcc.Link(dbc.Button([html.I(className="fas fa-database me-2"), "Database"],
                                    color="secondary", outline=True, className="me-2"),
                         href="/database"),
                dcc.Link(dbc.Button([html.I(className="fas fa-search me-2"), "Search"],
                                    color="secondary", outline=True),
                         href="/database"),
            ]),
        ], className="text-center py-5"),
    ], className="cm-card mb-4"),
], className="cm-page-container")

