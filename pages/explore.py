"""
Explore page, interactive discovery interface.
Three modes: browse by entity type, browse by interaction layer, network explorer.

Routes:
  /explore                 , landing page with three modes
  /explore/metabolites     , browse unique metabolites
  /explore/proteins        , browse unique proteins
  /explore/genes           , browse unique genes
  /explore/diseases        , browse unique diseases
  /explore/microbes        , browse unique microbes
  /explore/drugs           , browse unique drugs
  /explore/snps            , browse unique SNPs
"""

import math
from functools import lru_cache
from urllib.parse import quote_plus

from dash import html, dcc, dash_table, Input, Output, State, callback, no_update
import dash_bootstrap_components as dbc
import pandas as pd

from components.page_header import make_page_header

# ── Entity type colours (consistent across site) ─────────────────────
ENTITY_COLORS = {
    "metabolite": "#e27a3f",
    "protein":    "#3182ce",
    "gene":       "#d69e2e",
    "disease":    "#e53e3e",
    "microbe":    "#38a169",
    "drug":       "#805ad5",
    "snp":        "#319795",
}

# ── Interaction layers ────────────────────────────────────────────────
LAYERS = [
    ("MPI", "Metabolite–Protein", "Curated physical and enzymatic metabolite–protein interactions.",
     "/mpi", ENTITY_COLORS["protein"], "fas fa-cubes"),
    ("MEI", "Metabolite–Enzyme", "Enzyme-catalysed metabolic reactions with EC number annotations.",
     "/mpi", ENTITY_COLORS["protein"], "fas fa-flask"),
    ("MDI", "Metabolite–Disease", "Associations between metabolites and disease phenotypes.",
     "/mdi", ENTITY_COLORS["disease"], "fas fa-heartbeat"),
    ("MMI", "Metabolite–Microbe", "Production, consumption and transport of metabolites by gut microbes.",
     "/mmi", ENTITY_COLORS["microbe"], "fas fa-bacterium"),
    ("MDrI", "Metabolite–Drug", "Drug–metabolite pharmacokinetic and pharmacodynamic interactions.",
     "/mdri", ENTITY_COLORS["drug"], "fas fa-pills"),
    ("MGI", "Metabolite–Gene", "Genetic associations with metabolite concentrations, including eQTL-mediated links.",
     "/mgi", ENTITY_COLORS["gene"], "fas fa-dna"),
    ("mGWAS", "Metabolite–SNP", "Genome-wide association signals linking SNPs to metabolite concentrations.",
     "/mgwas", ENTITY_COLORS["snp"], "fas fa-map-marker-alt"),
]

# Entity browse configuration, maps URL suffix to lookup parameters
ENTITY_BROWSE_CONFIG = {
    "metabolites": {
        "label": "Metabolites", "singular": "Metabolite",
        "icon": "fas fa-atom", "color": ENTITY_COLORS["metabolite"],
        "desc": "Endogenous and exogenous small molecules in the knowledge graph.",
        "detail_route": "/metabolite", "id_param": "id",
    },
    "proteins": {
        "label": "Proteins", "singular": "Protein",
        "icon": "fas fa-cubes", "color": ENTITY_COLORS["protein"],
        "desc": "Enzymes, receptors, and transporters that interact with metabolites.",
        "detail_route": "/protein", "id_param": "name",
    },
    "genes": {
        "label": "Genes", "singular": "Gene",
        "icon": "fas fa-dna", "color": ENTITY_COLORS["gene"],
        "desc": "Gene symbols with metabolite-level associations across CoreMet layers.",
        "detail_route": "/gene", "id_param": "name",
    },
    "diseases": {
        "label": "Diseases", "singular": "Disease",
        "icon": "fas fa-heartbeat", "color": ENTITY_COLORS["disease"],
        "desc": "Disease phenotypes linked to metabolite changes.",
        "detail_route": "/disease-detail", "id_param": "name",
    },
    "microbes": {
        "label": "Microbes", "singular": "Microbe",
        "icon": "fas fa-bacterium", "color": ENTITY_COLORS["microbe"],
        "desc": "Gut and environmental micro-organisms connected to metabolites.",
        "detail_route": "/microbe", "id_param": "name",
    },
    "drugs": {
        "label": "Drugs", "singular": "Drug",
        "icon": "fas fa-pills", "color": ENTITY_COLORS["drug"],
        "desc": "Pharmaceuticals with known metabolite interactions.",
        "detail_route": "/drug", "id_param": "name",
    },
    "snps": {
        "label": "SNPs", "singular": "SNP",
        "icon": "fas fa-map-marker-alt", "color": ENTITY_COLORS["snp"],
        "desc": "Genetic variants linked to metabolite concentrations via mGWAS.",
        "detail_route": "/snp", "id_param": "name",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════

def _entity_card(name, icon, color, desc, href):
    return dbc.Col(
        html.A(
            html.Div([
                html.Div(html.I(className=icon), style={
                    "width": "44px", "height": "44px", "borderRadius": "10px",
                    "display": "flex", "alignItems": "center", "justifyContent": "center",
                    "fontSize": "1.1rem", "color": color, "background": f"{color}14",
                    "marginBottom": "12px",
                }),
                html.H5(name, style={"fontWeight": "600", "fontSize": "0.95rem",
                                     "color": "#1a202c", "marginBottom": "4px"}),
                html.P(desc, style={"fontSize": "0.8rem", "color": "#718096", "margin": "0"}),
            ], className="cm-card", style={"padding": "20px", "height": "100%"}),
            href=href, style={"textDecoration": "none"},
        ),
        xs=6, md=4, lg=3, className="mb-3",
    )


def _layer_card(short, full, desc, href, color, icon):
    return dbc.Col(
        html.A(
            html.Div([
                html.Div([
                    html.Div(html.I(className=icon), style={
                        "width": "40px", "height": "40px", "borderRadius": "8px",
                        "display": "flex", "alignItems": "center", "justifyContent": "center",
                        "fontSize": "1rem", "color": color, "background": f"{color}14",
                        "marginRight": "12px", "flexShrink": "0",
                    }),
                    html.Div([
                        html.Div([
                            html.Span(short, style={"fontWeight": "700", "fontSize": "0.9rem", "color": color}),
                            html.Span(f", {full}", style={"fontSize": "0.85rem", "color": "#4a5568"}),
                        ]),
                        html.P(desc, style={"fontSize": "0.78rem", "color": "#718096", "margin": "4px 0 0"}),
                    ]),
                ], style={"display": "flex", "alignItems": "flex-start"}),
            ], className="cm-card", style={"padding": "16px", "height": "100%"}),
            href=href, style={"textDecoration": "none"},
        ),
        md=6, className="mb-3",
    )


# ═══════════════════════════════════════════════════════════════════════
# Static landing layout
# ═══════════════════════════════════════════════════════════════════════

layout = html.Div([
    html.Div([
        make_page_header(
            title="Explore CoreMet",
            subtitle="Browse entities, interaction layers, or explore the network.",
            breadcrumb_items=[("Home", "/home"), ("Explore", None)],
        ),

        # ── Mode 1: Browse by Entity Type ───────────────────
        html.H3("Browse by entity type", style={
            "fontSize": "1.15rem", "fontWeight": "600", "color": "#1a202c",
            "marginBottom": "16px",
        }),

        dbc.Row([
            _entity_card("Metabolites", "fas fa-atom", ENTITY_COLORS["metabolite"],
                         "Endogenous and exogenous small molecules.", "/explore/metabolites"),
            _entity_card("Proteins", "fas fa-cubes", ENTITY_COLORS["protein"],
                         "Enzymes, receptors, and transporters.", "/explore/proteins"),
            _entity_card("Genes", "fas fa-dna", ENTITY_COLORS["gene"],
                         "Gene symbols with metabolite-level associations.", "/explore/genes"),
            _entity_card("Diseases", "fas fa-heartbeat", ENTITY_COLORS["disease"],
                         "Disease phenotypes linked to metabolite changes.", "/explore/diseases"),
            _entity_card("Microbes", "fas fa-bacterium", ENTITY_COLORS["microbe"],
                         "Gut and environmental micro-organisms.", "/explore/microbes"),
            _entity_card("Drugs", "fas fa-pills", ENTITY_COLORS["drug"],
                         "Pharmaceuticals with metabolite interactions.", "/explore/drugs"),
            _entity_card("SNPs", "fas fa-map-marker-alt", ENTITY_COLORS["snp"],
                         "Genetic variants linked to metabolite levels.", "/explore/snps"),
        ], className="g-3 mb-5"),

        # ── Mode 2: Browse by Interaction Layer ─────────────
        html.H3("Browse by interaction layer", style={
            "fontSize": "1.15rem", "fontWeight": "600", "color": "#1a202c",
            "marginBottom": "16px",
        }),

        dbc.Row([_layer_card(*layer) for layer in LAYERS], className="g-3 mb-5"),

        # ── Mode 3: Network Explorer ────────────────────────
        html.H3("Network explorer", style={
            "fontSize": "1.15rem", "fontWeight": "600", "color": "#1a202c",
            "marginBottom": "8px",
        }),
        html.P(
            "Launch the cross-layer graph explorer to traverse paths between any two entities.",
            style={"fontSize": "0.85rem", "color": "#718096", "marginBottom": "16px"},
        ),
        html.A(
            html.Div([
                html.Div([
                    html.I(className="fas fa-project-diagram",
                           style={"fontSize": "1.5rem", "color": "#3182ce", "marginRight": "16px"}),
                    html.Div([
                        html.Span("Open Network Explorer",
                                  style={"fontWeight": "600", "fontSize": "1rem", "color": "#1a202c"}),
                        html.P("Choose start and target nodes, allowed layers, maximum hops, "
                               "and confidence thresholds. View subgraphs and ranked paths.",
                               style={"fontSize": "0.85rem", "color": "#718096", "margin": "4px 0 0"}),
                    ]),
                    html.I(className="fas fa-arrow-right",
                           style={"marginLeft": "auto", "color": "#e2e8f0"}),
                ], style={"display": "flex", "alignItems": "center"}),
            ], className="cm-card", style={"padding": "24px"}),
            href="/network", style={"textDecoration": "none"},
        ),
    ], className="cm-page-container", style={"paddingBottom": "48px"}),
])

page_content = layout


# ═══════════════════════════════════════════════════════════════════════
# Dynamic entity browse pages  (called from display_page in main.py)
# ═══════════════════════════════════════════════════════════════════════

def _clean_text(value) -> str:
    """Normalize empty pandas values for display."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "<na>"}:
        return ""
    return text


@lru_cache(maxsize=8)
def _get_entity_list(entity_key: str) -> pd.DataFrame:
    """Return a DataFrame with columns [Name, ID, Interactions] for a given entity type."""
    rows = []

    if entity_key == "metabolites":
        by_hmdb = {}

        def _add_metabolites(df, hmdb_col, name_col):
            if df.empty or hmdb_col not in df.columns:
                return
            for hmdb, g in df.groupby(hmdb_col, observed=False):
                hmdb = _clean_text(hmdb)
                if not hmdb.startswith("HMDB"):
                    continue
                name = ""
                if name_col in g.columns:
                    names = [_clean_text(v) for v in g[name_col].dropna().unique()]
                    name = next((v for v in names if v), "")
                current = by_hmdb.setdefault(hmdb, {"Name": name, "ID": hmdb, "Interactions": 0})
                if not current["Name"] and name:
                    current["Name"] = name
                current["Interactions"] += len(g)

        try:
            from pages.metabolite_detail import _get_mpi_db
            _add_metabolites(_get_mpi_db(), "HMDB ID", "Metabolite Name")
        except Exception:
            pass
        try:
            from app.services.mei_service import get_mei_db
            _add_metabolites(get_mei_db(), "HMDB_ID", "Metabolite_Name")
        except Exception:
            pass
        try:
            from app.services.mdi_service import get_mdi_db
            _add_metabolites(get_mdi_db(), "HMDB_ID", "Metabolite_Name")
        except Exception:
            pass
        try:
            from app.services.mmi_service import get_mmi_db
            _add_metabolites(get_mmi_db(), "HMDB_ID", "Metabolite_Name")
        except Exception:
            pass
        try:
            from app.services.mdri_service import get_mdri_db
            _add_metabolites(get_mdri_db(), "HMDB_ID", "Metabolite_Name")
        except Exception:
            pass
        try:
            from app.services.mgi_service import get_mgi_db
            _add_metabolites(get_mgi_db(), "HMDB_ID", "Metabolite_Name")
        except Exception:
            pass
        try:
            from app.services.mgwas_service import get_mgwas_db
            _add_metabolites(get_mgwas_db(), "HMDB_ID", "Metabolite_Name")
        except Exception:
            pass

        rows = list(by_hmdb.values())

    elif entity_key == "proteins":
        try:
            from pages.metabolite_detail import _get_mpi_db
            mpi = _get_mpi_db()
            if not mpi.empty and "Protein Name" in mpi.columns:
                for name, g in mpi.groupby("Protein Name", observed=False):
                    name = _clean_text(name)
                    if not name:
                        continue
                    uid = _clean_text(g.iloc[0].get("Uniprot ID", ""))
                    rows.append({"Name": name, "ID": uid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "genes":
        try:
            from app.services.mgi_service import get_mgi_db
            df = get_mgi_db()
            if not df.empty and "Gene_Symbol" in df.columns:
                for name, g in df.groupby("Gene_Symbol", observed=False):
                    name = _clean_text(name)
                    if not name:
                        continue
                    gid = _clean_text(g.iloc[0].get("Gene_ID", ""))
                    rows.append({"Name": name, "ID": gid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "diseases":
        try:
            from app.services.mdi_service import get_mdi_db
            df = get_mdi_db()
            if not df.empty and "Disease_Name" in df.columns:
                for name, g in df.groupby("Disease_Name", observed=False):
                    name = _clean_text(name)
                    if not name:
                        continue
                    did = _clean_text(g.iloc[0].get("Disease_ID", g.iloc[0].get("MeSH_ID", "")))
                    rows.append({"Name": name, "ID": did, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "microbes":
        try:
            from app.services.mmi_service import get_mmi_db
            df = get_mmi_db()
            if not df.empty and "Microbe_Name" in df.columns:
                for name, g in df.groupby("Microbe_Name", observed=False):
                    name = _clean_text(name)
                    if not name:
                        continue
                    tid = _clean_text(g.iloc[0].get("Taxonomy_ID", ""))
                    rows.append({"Name": name, "ID": tid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "drugs":
        try:
            from app.services.mdri_service import get_mdri_db
            df = get_mdri_db()
            if not df.empty and "Drug_Name" in df.columns:
                for name, g in df.groupby("Drug_Name", observed=False):
                    name = _clean_text(name)
                    if not name:
                        continue
                    dbid = _clean_text(g.iloc[0].get("DrugBank_ID", ""))
                    rows.append({"Name": name, "ID": dbid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "snps":
        try:
            from app.services.mgwas_service import get_mgwas_db
            df = get_mgwas_db()
            if not df.empty and "rsID" in df.columns:
                for name, g in df.groupby("rsID", observed=False):
                    name = _clean_text(name)
                    if not name:
                        continue
                    rows.append({"Name": name, "ID": name, "Interactions": len(g)})
        except Exception:
            pass

    result = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["Name", "ID", "Interactions"])
    if not result.empty:
        result["Name"] = result["Name"].map(_clean_text)
        result["ID"] = result["ID"].map(_clean_text).replace("", "-")
        result = result[result["Name"] != ""]
        if entity_key == "metabolites":
            result = result[result["ID"].str.startswith("HMDB", na=False)]
    return result.sort_values("Interactions", ascending=False).reset_index(drop=True)


def build_entity_browse_page(entity_key: str) -> html.Div:
    """Build a paginated browse page that lists all unique entities of a given type."""
    cfg = ENTITY_BROWSE_CONFIG.get(entity_key)
    if not cfg:
        return html.Div("Unknown entity type.", className="text-center text-muted p-5")

    color = cfg["color"]

    return html.Div([
        html.Div([
            make_page_header(
                title=f"Browse {cfg['label']}",
                subtitle=cfg["desc"],
                breadcrumb_items=[("Home", "/home"), ("Explore", "/explore"), (cfg["label"], None)],
            ),

            # Summary
            html.Div([
                html.Div([
                    html.I(className=f"{cfg['icon']} me-2", style={"color": color, "fontSize": "1.1rem"}),
                    html.Span(id="explore-entity-total", style={
                        "fontSize": "1.5rem", "fontWeight": "700", "color": color,
                    }),
                    html.Span(f"  unique {cfg['label'].lower()} in CoreMet",
                              style={"fontSize": "0.9rem", "color": "#718096", "marginLeft": "8px"}),
                ], style={"display": "flex", "alignItems": "center"}),
            ], className="cm-card mb-4", style={"padding": "16px"}),

            dbc.InputGroup(
                [
                    dbc.InputGroupText(html.I(className="fas fa-search")),
                    dbc.Input(
                        id="explore-entity-search",
                        type="text",
                        debounce=True,
                        placeholder=f"Search {cfg['label'].lower()}...",
                    ),
                ],
                className="mb-3",
            ),

            html.Div(id="explore-entity-result-count", className="mb-2", style={
                "fontSize": "0.85rem", "color": "var(--cm-text-secondary)",
            }),

            # Data table
            dash_table.DataTable(
                id="explore-entity-table",
                data=[],
                columns=[
                    {"name": "Name", "id": "Name", "type": "text", "presentation": "markdown"},
                    {"name": "ID", "id": "ID"},
                    {"name": "Interactions", "id": "Interactions"},
                ],
                page_size=25,
                page_action="custom",
                page_current=0,
                sort_action="custom",
                sort_mode="single",
                style_table={"overflowX": "auto", "border": "1px solid #e2e8f0", "borderRadius": "8px"},
                style_header={
                    "backgroundColor": color, "color": "white",
                    "fontWeight": "600", "fontSize": "0.8rem", "padding": "8px 12px",
                },
                style_cell={
                    "fontSize": "0.85rem", "padding": "8px 12px",
                    "textAlign": "left", "fontFamily": "Arial, Helvetica, sans-serif",
                },
                style_data_conditional=[
                    {"if": {"row_index": "odd"}, "backgroundColor": "#fafbfc"},
                ],
            ),
            dcc.Store(id="explore-entity-key", data=entity_key),
        ], className="cm-page-container"),
    ])


@callback(
    Output("explore-entity-table", "data"),
    Output("explore-entity-table", "page_count"),
    Output("explore-entity-table", "page_current"),
    Output("explore-entity-total", "children"),
    Output("explore-entity-result-count", "children"),
    Input("explore-entity-key", "data"),
    Input("explore-entity-search", "value"),
    Input("explore-entity-table", "page_current"),
    Input("explore-entity-table", "page_size"),
    Input("explore-entity-table", "sort_by"),
)
def update_entity_browse_table(entity_key, search_text, page_current, page_size, sort_by):
    """Serve one clean page of entity browse records at a time."""
    if not entity_key:
        return [], 0, 0, "", ""

    df = _get_entity_list(entity_key).copy()
    total = len(df)

    if search_text:
        q = search_text.strip().lower()
        if q:
            mask = (
                df["Name"].astype(str).str.lower().str.contains(q, na=False, regex=False) |
                df["ID"].astype(str).str.lower().str.contains(q, na=False, regex=False)
            )
            df = df[mask]

    if sort_by:
        for sort in reversed(sort_by):
            col = sort.get("column_id")
            if col in df.columns:
                df = df.sort_values(col, ascending=sort.get("direction") == "asc")
    else:
        df = df.sort_values("Interactions", ascending=False)

    page_size = page_size or 25
    page_count = max(1, math.ceil(len(df) / page_size)) if len(df) else 0
    page_current = min(page_current or 0, max(page_count - 1, 0))
    start = page_current * page_size
    page_df = df.iloc[start:start + page_size].copy()

    cfg = ENTITY_BROWSE_CONFIG.get(entity_key, {})
    detail_route = cfg.get("detail_route", "/search")
    id_param = cfg.get("id_param", "q")
    if not page_df.empty:
        page_df["Name"] = page_df.apply(
            lambda r: (
                f"[{r['Name']}]({detail_route}?{id_param}="
                f"{quote_plus(str(r['ID'] if entity_key == 'metabolites' and r['ID'] != '-' else r['Name']))})"
            ),
            axis=1,
        )

    visible = len(df)
    if visible:
        result_text = f"Showing {start + 1:,}-{min(start + page_size, visible):,} of {visible:,} records"
    else:
        result_text = "No matching records"

    return page_df.to_dict("records"), page_count, page_current, f"{total:,}", result_text
