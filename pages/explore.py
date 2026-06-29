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

def _get_entity_list(entity_key: str) -> pd.DataFrame:
    """Return a DataFrame with columns [Name, ID, Interactions] for a given entity type."""
    rows = []

    if entity_key == "metabolites":
        try:
            from pages.metabolite_detail import _get_mpi_db
            mpi = _get_mpi_db()
            if not mpi.empty:
                for _, g in mpi.groupby("HMDB ID"):
                    r = g.iloc[0]
                    rows.append({"Name": str(r.get("Metabolite Name", "")),
                                 "ID": str(r.get("HMDB ID", "")),
                                 "Interactions": len(g)})
        except Exception:
            pass
        try:
            from app.services.mdi_service import get_mdi_db
            mdi = get_mdi_db()
            if not mdi.empty and "HMDB_ID" in mdi.columns:
                for hmdb, g in mdi.groupby("HMDB_ID"):
                    if any(r.get("ID") == hmdb for r in rows):
                        continue
                    r = g.iloc[0]
                    rows.append({"Name": str(r.get("Metabolite_Name", "")),
                                 "ID": str(hmdb), "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "proteins":
        try:
            from pages.metabolite_detail import _get_mpi_db
            mpi = _get_mpi_db()
            if not mpi.empty and "Protein Name" in mpi.columns:
                for name, g in mpi.groupby("Protein Name"):
                    if not str(name).strip():
                        continue
                    uid = str(g.iloc[0].get("Uniprot ID", ""))
                    rows.append({"Name": str(name), "ID": uid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "genes":
        try:
            from app.services.mgi_service import get_mgi_db
            df = get_mgi_db()
            if not df.empty and "Gene_Symbol" in df.columns:
                for name, g in df.groupby("Gene_Symbol"):
                    gid = str(g.iloc[0].get("Gene_ID", ""))
                    rows.append({"Name": str(name), "ID": gid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "diseases":
        try:
            from app.services.mdi_service import get_mdi_db
            df = get_mdi_db()
            if not df.empty and "Disease_Name" in df.columns:
                for name, g in df.groupby("Disease_Name"):
                    did = str(g.iloc[0].get("Disease_ID", g.iloc[0].get("MeSH_ID", "")))
                    rows.append({"Name": str(name), "ID": did, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "microbes":
        try:
            from app.services.mmi_service import get_mmi_db
            df = get_mmi_db()
            if not df.empty and "Microbe_Name" in df.columns:
                for name, g in df.groupby("Microbe_Name"):
                    tid = str(g.iloc[0].get("Taxonomy_ID", ""))
                    rows.append({"Name": str(name), "ID": tid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "drugs":
        try:
            from app.services.mdri_service import get_mdri_db
            df = get_mdri_db()
            if not df.empty and "Drug_Name" in df.columns:
                for name, g in df.groupby("Drug_Name"):
                    dbid = str(g.iloc[0].get("DrugBank_ID", ""))
                    rows.append({"Name": str(name), "ID": dbid, "Interactions": len(g)})
        except Exception:
            pass

    elif entity_key == "snps":
        try:
            from app.services.mgwas_service import get_mgwas_db
            df = get_mgwas_db()
            if not df.empty and "rsID" in df.columns:
                for name, g in df.groupby("rsID"):
                    rows.append({"Name": str(name), "ID": str(name), "Interactions": len(g)})
        except Exception:
            pass

    result = pd.DataFrame(rows) if rows else pd.DataFrame(columns=["Name", "ID", "Interactions"])
    return result.sort_values("Interactions", ascending=False).reset_index(drop=True)


def build_entity_browse_page(entity_key: str) -> html.Div:
    """Build a paginated browse page that lists all unique entities of a given type."""
    cfg = ENTITY_BROWSE_CONFIG.get(entity_key)
    if not cfg:
        return html.Div("Unknown entity type.", className="text-center text-muted p-5")

    df = _get_entity_list(entity_key)
    color = cfg["color"]
    detail_route = cfg["detail_route"]
    id_param = cfg["id_param"]

    # Build linked name column, markdown links to entity detail pages
    if not df.empty:
        df["Name"] = df.apply(
            lambda r: f"[{r['Name']}]({detail_route}?{id_param}={quote_plus(str(r['ID'] if entity_key == 'metabolites' else r['Name']))})",
            axis=1,
        )

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
                    html.Span(f"{len(df):,}", style={"fontSize": "1.5rem", "fontWeight": "700", "color": color}),
                    html.Span(f"  unique {cfg['label'].lower()} in CoreMet",
                              style={"fontSize": "0.9rem", "color": "#718096", "marginLeft": "8px"}),
                ], style={"display": "flex", "alignItems": "center"}),
            ], className="cm-card mb-4", style={"padding": "16px"}),

            # Data table
            dash_table.DataTable(
                id=f"explore-{entity_key}-table",
                data=df.to_dict("records") if not df.empty else [],
                columns=[
                    {"name": "Name", "id": "Name", "type": "text", "presentation": "markdown"},
                    {"name": "ID", "id": "ID"},
                    {"name": "Interactions", "id": "Interactions"},
                ],
                page_size=25,
                page_action="native",
                sort_action="native",
                sort_mode="single",
                filter_action="native",
                style_table={"overflowX": "auto"},
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
        ], className="cm-page-container"),
    ])
