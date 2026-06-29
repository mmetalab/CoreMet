"""
Reusable stat card component, icon + number + label with accent stripe.
"""

from dash import html
import dash_bootstrap_components as dbc


def make_stat_card(icon_class, value, label, accent="primary"):
    """
    Create a stat card with left accent stripe.

    Args:
        icon_class: FontAwesome icon class (e.g. "fas fa-dna")
        value: The number to display
        label: Description label
        accent: One of "metabolite", "protein", "interaction", "accent", "primary"
    """
    card = html.Div(
        [
            html.Div(html.I(className=icon_class), className="stat-icon"),
            html.Div(f"{value:,}" if isinstance(value, (int, float)) else str(value), className="stat-value"),
            html.Div(label, className="stat-label"),
        ],
        className=f"cm-stat-card accent-{accent}",
    )
    return card
