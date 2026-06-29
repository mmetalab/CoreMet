"""
Reusable navbar component, minimal fixed top bar.
CoreMet | Search | Explore | CoreMet-AI | Downloads | API | Help
"""

from dash import html, dcc
import dash_bootstrap_components as dbc


def make_navbar():
    """Create the main navigation bar, clean, journal-like."""
    navbar = dbc.Navbar(
        dbc.Container(
            [
                # Brand / logo
                dbc.NavbarBrand(
                    html.Span("CoreMet", style={"fontWeight": "700"}),
                    href="/home",
                    className="d-flex align-items-center",
                ),
                # Toggler for mobile
                dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
                # Collapsible nav links, minimal set
                dbc.Collapse(
                    dbc.Nav(
                        [
                            # Database-first ordering: browse/search/explore the data,
                            # then download/API/help; prediction (AI) is a secondary feature.
                            dbc.NavItem(dbc.NavLink("Browse", href="/database", id="nav-browse")),
                            dbc.NavItem(dbc.NavLink("Search", href="/search", id="nav-search")),
                            dbc.NavItem(dbc.NavLink("Explore", href="/explore", id="nav-explore")),
                            dbc.NavItem(dbc.NavLink("Downloads", href="/downloads", id="nav-downloads")),
                            dbc.NavItem(dbc.NavLink("API", href="/api-docs", id="nav-api")),
                            dbc.NavItem(dbc.NavLink("Help", href="/help", id="nav-help")),
                            dbc.NavItem(dbc.NavLink("AI (beta)", href="/coremetai", id="nav-ai")),
                        ],
                        navbar=True,
                        className="ms-auto d-flex align-items-center",
                    ),
                    id="navbar-collapse",
                    is_open=False,
                    navbar=True,
                ),
            ],
            fluid=True,
        ),
        className="cm-navbar",
        dark=True,
        fixed="top",
    )
    return navbar
