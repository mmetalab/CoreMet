"""
Reusable page header component, title + subtitle + optional breadcrumb.
"""

from dash import html


def make_page_header(title, subtitle=None, breadcrumb_items=None):
    """
    Create a consistent page header.

    Args:
        title: Page title string
        subtitle: Optional subtitle/description
        breadcrumb_items: Optional list of (label, href) tuples for breadcrumb
    """
    children = []

    if breadcrumb_items:
        crumbs = []
        for label, href in breadcrumb_items:
            if href:
                crumbs.append(html.A(label, href=href))
                crumbs.append(html.Span(" / ", style={"margin": "0 4px"}))
            else:
                crumbs.append(html.Span(label))
        children.append(html.Div(crumbs, className="page-breadcrumb"))

    children.append(html.H1(title, className="page-title"))

    if subtitle:
        children.append(html.P(subtitle, className="page-subtitle"))

    return html.Div(children, className="cm-page-header")
