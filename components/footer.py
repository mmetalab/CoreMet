"""
Reusable footer component.
"""

from dash import html


def make_footer():
    """Create the page footer with citation, version, and links."""
    return html.Div(
        [
            html.Span("CoreMet"),
            html.Span(" | "),
            html.A("GitHub", href="https://github.com/mmetalab/CoreMet", target="_blank"),
            html.Span(" | "),
            html.A("Downloads", href="/downloads"),
            html.Span(" | "),
            html.A("Help", href="/help"),
        ],
        className="cm-footer",
    )
