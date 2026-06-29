"""
Generic Entity Detail Page, URL-addressable detail view for any non-metabolite entity.

Routes:
  /disease?name=Colorectal+cancer
  /gene?name=LDHA
  /protein?name=Albumin  or  /protein?id=P02768
  /drug?name=Metformin
  /microbe?name=Lactobacillus+rhamnosus
  /snp?name=rs1260326

Each page aggregates data from all 7 CoreMet databases and shows a consistent
tabbed layout: Overview | Network | Interactions | Evidence | Downloads.
"""

from urllib.parse import parse_qs, quote_plus
from pathlib import Path

from dash import dcc, html, dash_table
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import plotly.graph_objects as go
import pandas as pd

from components.page_header import make_page_header

# ---------------------------------------------------------------------------
# Entity type configuration
# ---------------------------------------------------------------------------
ENTITY_TYPES = {
    "disease": {
        "label": "Disease",
        "color": "#e53e3e",
        "icon": "fas fa-heartbeat",
        "external_url": "https://meshb.nlm.nih.gov/record/ui?ui={}",
        "external_label": "MeSH",
        "id_col": "Disease_ID",
        "name_col": "Disease_Name",
    },
    "gene": {
        "label": "Gene",
        "color": "#d69e2e",
        "icon": "fas fa-dna",
        "external_url": "https://www.ncbi.nlm.nih.gov/gene/{}",
        "external_label": "NCBI Gene",
        "id_col": "Gene_ID",
        "name_col": "Gene_Symbol",
    },
    "protein": {
        "label": "Protein",
        "color": "#3182ce",
        "icon": "fas fa-cubes",
        "external_url": "https://www.uniprot.org/uniprot/{}",
        "external_label": "UniProt",
        "id_col": "Uniprot_ID",
        "name_col": "Protein_Name",
    },
    "drug": {
        "label": "Drug",
        "color": "#805ad5",
        "icon": "fas fa-pills",
        "external_url": "https://go.drugbank.com/drugs/{}",
        "external_label": "DrugBank",
        "id_col": "DrugBank_ID",
        "name_col": "Drug_Name",
    },
    "microbe": {
        "label": "Microbe",
        "color": "#38a169",
        "icon": "fas fa-bacterium",
        "external_url": "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}",
        "external_label": "NCBI Taxonomy",
        "id_col": "Taxonomy_ID",
        "name_col": "Microbe_Name",
    },
    "snp": {
        "label": "SNP",
        "color": "#319795",
        "icon": "fas fa-map-marker-alt",
        "external_url": "https://www.ncbi.nlm.nih.gov/snp/{}",
        "external_label": "dbSNP",
        "id_col": "rsID",
        "name_col": "rsID",
    },
}

# Metabolite color for network nodes
# Consistent entity-type colors (matches CSS design system)
_ENTITY_COLORS = {
    "metabolite": "#e27a3f",
    "protein": "#3182ce",
    "enzyme": "#3182ce",
    "gene": "#d69e2e",
    "disease": "#e53e3e",
    "microbe": "#38a169",
    "drug": "#805ad5",
    "snp": "#319795",
}
_MET_COLOR = _ENTITY_COLORS["metabolite"]

_NODE_SHAPES = {
    "metabolite": "ellipse",
    "protein": "round-rectangle",
    "enzyme": "round-rectangle",
    "disease": "triangle",
    "microbe": "hexagon",
    "drug": "diamond",
    "snp": "star",
    "gene": "round-pentagon",
}

# Cytoscape stylesheet, publication-quality, Nature-style
_STYLESHEET = [
    # --- default node ---
    {"selector": "node", "style": {
        "label": "data(label)",
        "font-size": "9px",
        "font-family": "Arial, Helvetica, sans-serif",
        "font-weight": "500",
        "text-valign": "center",
        "text-halign": "center",
        "width": "data(size)",
        "height": "data(size)",
        "background-color": "data(color)",
        "background-opacity": 0.92,
        "shape": "data(shape)",
        "text-wrap": "ellipsis",
        "text-max-width": "80px",
        "border-width": 1.5,
        "border-color": "#cbd5e0",
        "text-outline-color": "#ffffff",
        "text-outline-width": 2,
        "text-outline-opacity": 0.9,
        "color": "#1a202c",
    }},
    # --- center node ---
    {"selector": "node.center", "style": {
        "width": 56, "height": 56,
        "font-size": "11px", "font-weight": "bold",
        "border-width": 3, "border-color": "#1a202c",
        "text-outline-width": 2.5,
        "z-index": 10,
    }},
    # --- metabolite ring ---
    {"selector": "node.metabolite", "style": {
        "width": 32, "height": 32, "font-size": "8px",
    }},
    # --- outer ring entities (2nd hop) ---
    {"selector": "node.outer", "style": {
        "width": 24, "height": 24, "font-size": "7px",
        "background-opacity": 0.8, "border-width": 1,
    }},
    # --- edges ---
    {"selector": "edge", "style": {
        "width": "data(weight)",
        "line-color": "data(color)",
        "target-arrow-color": "data(color)",
        "opacity": 0.55,
        "curve-style": "bezier",
    }},
    {"selector": "edge.hop1", "style": {
        "width": 2.0, "opacity": 0.65,
    }},
    {"selector": "edge.hop2", "style": {
        "width": 1.0, "opacity": 0.35, "line-style": "dashed",
    }},
    # --- selection ---
    {"selector": ":selected", "style": {
        "border-width": 3, "border-color": "#f6e05e",
        "background-opacity": 1,
    }},
    # --- hover ---
    {"selector": "node:active", "style": {
        "overlay-color": "#f6e05e", "overlay-opacity": 0.15,
    }},
]


# ---------------------------------------------------------------------------
# Cross-layer lookup for 2-hop network expansion
# ---------------------------------------------------------------------------

def _lookup_metabolite_neighbors(hmdb_ids: list, exclude_etype: str, max_per_layer: int = 8):
    """Given HMDB IDs, find non-metabolite entities connected to these metabolites.

    Returns list of dicts: {id, label, etype, color, shape, met_id, layer}
    """
    neighbors = []
    seen = set()

    # Column mappings: (db_service, db_key, id_col_met, name_col_target, id_col_target, etype, layer_label)
    _LAYER_DEFS = [
        ("app.services.mdi_service", "get_mdi_db", "HMDB_ID", "Disease_Name", "Disease_ID", "disease", "MDI"),
        ("app.services.mmi_service", "get_mmi_db", "HMDB_ID", "Microbe_Name", "Taxonomy_ID", "microbe", "MMI"),
        ("app.services.mdri_service", "get_mdri_db", "HMDB_ID", "Drug_Name", "DrugBank_ID", "drug", "MDrI"),
        ("app.services.mgwas_service", "get_mgwas_db", "HMDB_ID", "rsID", "rsID", "snp", "mGWAS"),
        ("app.services.mgi_service", "get_mgi_db", "HMDB_ID", "Gene_Symbol", "Gene_ID", "gene", "MGI"),
    ]

    # Also check MPI for proteins
    try:
        from pages.metabolite_detail import _get_mpi_db
        mpi = _get_mpi_db()
        if not mpi.empty:
            hc = "HMDB ID" if "HMDB ID" in mpi.columns else "HMDB_ID"
            nc = "Protein Name" if "Protein Name" in mpi.columns else "Protein_Name"
            ic = "Uniprot ID" if "Uniprot ID" in mpi.columns else "Uniprot_ID"
            if hc in mpi.columns and nc in mpi.columns and "protein" != exclude_etype:
                mask = mpi[hc].isin(hmdb_ids)
                sub = mpi[mask].drop_duplicates(subset=[nc]).head(max_per_layer)
                for _, row in sub.iterrows():
                    name = str(row.get(nc, ""))
                    nid = str(row.get(ic, name))
                    met_id = str(row.get(hc, ""))
                    if name and name != "nan" and name not in seen:
                        seen.add(name)
                        neighbors.append({
                            "id": f"prot-{name}", "label": name, "etype": "protein",
                            "color": _ENTITY_COLORS["protein"], "shape": _NODE_SHAPES["protein"],
                            "met_id": met_id, "layer": "MPI", "ext_id": nid,
                        })
    except Exception:
        pass

    # Check MEI for enzymes (also proteins)
    try:
        from app.services.mei_service import get_mei_db
        mei = get_mei_db()
        if not mei.empty and "protein" != exclude_etype:
            hc = "HMDB_ID"
            nc = "Enzyme_Name"
            ic = "Uniprot_ID"
            if hc in mei.columns and nc in mei.columns:
                mask = mei[hc].isin(hmdb_ids)
                sub = mei[mask].drop_duplicates(subset=[nc]).head(max_per_layer)
                for _, row in sub.iterrows():
                    name = str(row.get(nc, ""))
                    nid = str(row.get(ic, name))
                    met_id = str(row.get(hc, ""))
                    if name and name != "nan" and name not in seen:
                        seen.add(name)
                        neighbors.append({
                            "id": f"enz-{name}", "label": name, "etype": "protein",
                            "color": _ENTITY_COLORS["protein"], "shape": _NODE_SHAPES["enzyme"],
                            "met_id": met_id, "layer": "MEI", "ext_id": nid,
                        })
    except Exception:
        pass

    for mod_path, fn_name, hc, nc, ic, etype, layer in _LAYER_DEFS:
        if etype == exclude_etype:
            continue
        try:
            import importlib
            mod = importlib.import_module(mod_path)
            df = getattr(mod, fn_name)()
            if df.empty or hc not in df.columns or nc not in df.columns:
                continue
            mask = df[hc].isin(hmdb_ids)
            sub = df[mask].drop_duplicates(subset=[nc]).head(max_per_layer)
            for _, row in sub.iterrows():
                name = str(row.get(nc, ""))
                nid = str(row.get(ic, name)) if ic in df.columns else name
                met_id = str(row.get(hc, ""))
                if name and name != "nan" and name not in seen:
                    seen.add(name)
                    neighbors.append({
                        "id": f"{etype[:3]}-{name}", "label": name, "etype": etype,
                        "color": _ENTITY_COLORS.get(etype, "#a0aec0"),
                        "shape": _NODE_SHAPES.get(etype, "ellipse"),
                        "met_id": met_id, "layer": layer, "ext_id": nid,
                    })
        except Exception:
            continue

    return neighbors


# ---------------------------------------------------------------------------
# Data lookup: given an entity type + name/id, find all related metabolites
# ---------------------------------------------------------------------------

def _lookup_entity(etype: str, name: str = "", eid: str = "") -> dict:
    """Look up an entity across all databases, returning linked metabolites and edges.

    Returns dict with:
      name, id, etype, and DataFrames for each database where this entity appears.
    """
    cfg = ENTITY_TYPES[etype]
    query = name.strip() if name else eid.strip()
    if not query or len(query) < 2:
        return {"name": query, "id": eid, "etype": etype,
                "mpi": pd.DataFrame(), "mei": pd.DataFrame(), "mdi": pd.DataFrame(),
                "mmi": pd.DataFrame(), "mdri": pd.DataFrame(), "mgwas": pd.DataFrame(),
                "mgi": pd.DataFrame()}

    result = {"name": name or eid, "id": eid, "etype": etype,
              "mpi": pd.DataFrame(), "mei": pd.DataFrame(), "mdi": pd.DataFrame(),
              "mmi": pd.DataFrame(), "mdri": pd.DataFrame(), "mgwas": pd.DataFrame(),
              "mgi": pd.DataFrame()}

    q = query.lower()

    # --- Disease: search MDI ---
    if etype == "disease":
        try:
            from app.services.mdi_service import get_mdi_db
            df = get_mdi_db()
            if not df.empty:
                mask = pd.Series(False, index=df.index)
                for col in ["Disease_Name", "Disease_ID", "MeSH_ID"]:
                    if col in df.columns:
                        mask |= df[col].astype(str).str.lower() == q
                if not mask.any():
                    for col in ["Disease_Name"]:
                        if col in df.columns:
                            mask |= df[col].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mdi"] = df[mask]
                if not result["mdi"].empty:
                    row = result["mdi"].iloc[0]
                    result["name"] = str(row.get("Disease_Name", query))
                    result["id"] = str(row.get("Disease_ID", "") or row.get("MeSH_ID", ""))
        except Exception:
            pass

    # --- Gene: search MGI + MPI (Gene_Name) + MEI (Gene_Name) ---
    elif etype == "gene":
        try:
            from app.services.mgi_service import get_mgi_db
            df = get_mgi_db()
            if not df.empty:
                mask = pd.Series(False, index=df.index)
                for col in ["Gene_Symbol", "Gene_ID"]:
                    if col in df.columns:
                        mask |= df[col].astype(str).str.lower() == q
                if not mask.any():
                    if "Gene_Symbol" in df.columns:
                        mask |= df["Gene_Symbol"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mgi"] = df[mask]
                if not result["mgi"].empty:
                    row = result["mgi"].iloc[0]
                    result["name"] = str(row.get("Gene_Symbol", query))
                    result["id"] = str(row.get("Gene_ID", ""))
        except Exception:
            pass
        # Also search MPI for Gene Name
        try:
            from pages.metabolite_detail import _get_mpi_db
            mpi = _get_mpi_db()
            if not mpi.empty and "Gene Name" in mpi.columns:
                mask = mpi["Gene Name"].astype(str).str.lower() == q
                if not mask.any():
                    mask = mpi["Gene Name"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mpi"] = mpi[mask]
                if not result["name"] or result["name"] == query:
                    if not result["mpi"].empty:
                        result["name"] = str(result["mpi"].iloc[0].get("Gene Name", query))
        except Exception:
            pass
        # Also search MEI for Gene_Name
        try:
            from app.services.mei_service import get_mei_db
            mei = get_mei_db()
            if not mei.empty and "Gene_Name" in mei.columns:
                mask = mei["Gene_Name"].astype(str).str.lower() == q
                if not mask.any():
                    mask = mei["Gene_Name"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mei"] = mei[mask]
        except Exception:
            pass
        # Also search mGWAS for Mapped_Gene
        try:
            from app.services.mgwas_service import get_mgwas_db
            mgwas = get_mgwas_db()
            if not mgwas.empty and "Mapped_Gene" in mgwas.columns:
                mask = mgwas["Mapped_Gene"].astype(str).str.lower() == q
                if not mask.any():
                    mask = mgwas["Mapped_Gene"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mgwas"] = mgwas[mask]
        except Exception:
            pass

    # --- Protein: search MPI + MEI ---
    elif etype == "protein":
        try:
            from pages.metabolite_detail import _get_mpi_db
            mpi = _get_mpi_db()
            if not mpi.empty:
                mask = pd.Series(False, index=mpi.index)
                for col in ["Protein Name", "Uniprot ID"]:
                    if col in mpi.columns:
                        mask |= mpi[col].astype(str).str.lower() == q
                if not mask.any():
                    for col in ["Protein Name"]:
                        if col in mpi.columns:
                            mask |= mpi[col].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mpi"] = mpi[mask]
                if not result["mpi"].empty:
                    row = result["mpi"].iloc[0]
                    result["name"] = str(row.get("Protein Name", query))
                    result["id"] = str(row.get("Uniprot ID", ""))
        except Exception:
            pass
        try:
            from app.services.mei_service import get_mei_db
            mei = get_mei_db()
            if not mei.empty:
                mask = pd.Series(False, index=mei.index)
                for col in ["Enzyme_Name", "Uniprot_ID"]:
                    if col in mei.columns:
                        mask |= mei[col].astype(str).str.lower() == q
                if not mask.any():
                    for col in ["Enzyme_Name"]:
                        if col in mei.columns:
                            mask |= mei[col].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mei"] = mei[mask]
                if not result["mei"].empty and (not result["name"] or result["name"] == query):
                    row = result["mei"].iloc[0]
                    result["name"] = str(row.get("Enzyme_Name", query))
                    result["id"] = str(row.get("Uniprot_ID", ""))
        except Exception:
            pass

    # --- Drug: search MDrI ---
    elif etype == "drug":
        try:
            from app.services.mdri_service import get_mdri_db
            df = get_mdri_db()
            if not df.empty:
                mask = pd.Series(False, index=df.index)
                for col in ["Drug_Name", "DrugBank_ID"]:
                    if col in df.columns:
                        mask |= df[col].astype(str).str.lower() == q
                if not mask.any():
                    if "Drug_Name" in df.columns:
                        mask |= df["Drug_Name"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mdri"] = df[mask]
                if not result["mdri"].empty:
                    row = result["mdri"].iloc[0]
                    result["name"] = str(row.get("Drug_Name", query))
                    result["id"] = str(row.get("DrugBank_ID", ""))
        except Exception:
            pass

    # --- Microbe: search MMI ---
    elif etype == "microbe":
        try:
            from app.services.mmi_service import get_mmi_db
            df = get_mmi_db()
            if not df.empty:
                mask = pd.Series(False, index=df.index)
                for col in ["Microbe_Name", "Taxonomy_ID"]:
                    if col in df.columns:
                        mask |= df[col].astype(str).str.lower() == q
                if not mask.any():
                    if "Microbe_Name" in df.columns:
                        mask |= df["Microbe_Name"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mmi"] = df[mask]
                if not result["mmi"].empty:
                    row = result["mmi"].iloc[0]
                    result["name"] = str(row.get("Microbe_Name", query))
                    result["id"] = str(row.get("Taxonomy_ID", ""))
        except Exception:
            pass

    # --- SNP: search mGWAS ---
    elif etype == "snp":
        try:
            from app.services.mgwas_service import get_mgwas_db
            df = get_mgwas_db()
            if not df.empty:
                mask = pd.Series(False, index=df.index)
                for col in ["rsID"]:
                    if col in df.columns:
                        mask |= df[col].astype(str).str.lower() == q
                if not mask.any():
                    mask |= df["rsID"].astype(str).str.lower().str.contains(q, na=False, regex=False)
                result["mgwas"] = df[mask]
                if not result["mgwas"].empty:
                    row = result["mgwas"].iloc[0]
                    result["name"] = str(row.get("rsID", query))
                    result["id"] = str(row.get("rsID", ""))
        except Exception:
            pass

    return result


# ---------------------------------------------------------------------------
# Public builder: called from display_page
# ---------------------------------------------------------------------------

def build_page(etype: str, search: str = "") -> html.Div:
    """Build an entity detail page for any supported entity type."""
    if etype not in ENTITY_TYPES:
        return html.Div("Unknown entity type.", className="text-center text-muted p-5")

    cfg = ENTITY_TYPES[etype]
    params = parse_qs((search or "").lstrip("?"))
    name = params.get("name", [""])[0].strip()
    eid = params.get("id", [""])[0].strip()

    if not name and not eid:
        return html.Div([
            html.Div([
                make_page_header(
                    f"{cfg['label']} Detail",
                    f"View all interactions for a {cfg['label'].lower()} in CoreMet.",
                    [("Home", "/home"), ("Search", "/search"), (f"{cfg['label']} Detail", None)],
                ),
                _empty_state(cfg, f"No {cfg['label'].lower()} specified. Use a link from Search or provide ?name=... in the URL."),
            ], className="cm-page-container"),
        ])

    data = _lookup_entity(etype, name=name, eid=eid)
    total = sum(len(data[k]) for k in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"])

    if total == 0:
        q = name or eid
        return html.Div([
            html.Div([
                make_page_header(
                    f"{cfg['label']} Detail",
                    f"View all interactions for a {cfg['label'].lower()} in CoreMet.",
                    [("Home", "/home"), ("Search", "/search"), (f"{cfg['label']} Detail", None)],
                ),
                _empty_state(cfg, f'No records found for "{q}" in any CoreMet database.'),
            ], className="cm-page-container"),
        ])

    content = _build_detail_view(data)

    return html.Div([
        html.Div([
            make_page_header(
                f"{cfg['label']} Detail",
                f"Cross-layer interactions for {data['name']} in CoreMet.",
                [("Home", "/home"), ("Search", "/search"), (data["name"], None)],
            ),
            content,
        ], className="cm-page-container"),
    ])


# ---------------------------------------------------------------------------
# Empty state
# ---------------------------------------------------------------------------

def _empty_state(cfg: dict, message: str):
    return html.Div([
        html.Div([
            html.I(className=f"{cfg['icon']} fa-2x mb-3", style={"color": "#a0aec0"}),
            html.P(message, className="text-muted mb-3"),
            html.Div([
                dcc.Link(
                    dbc.Button([html.I(className="fas fa-search me-2"), "Go to Search"],
                               className="cm-btn-secondary mt-2"),
                    href="/search",
                ),
            ]),
        ], className="text-center p-5"),
    ], className="cm-card")


# ---------------------------------------------------------------------------
# Build the full detail view
# ---------------------------------------------------------------------------

def _build_detail_view(data: dict):
    """Assemble the entity detail page with tabbed layout."""
    etype = data["etype"]
    cfg = ENTITY_TYPES[etype]
    ename = data["name"]
    eid = data["id"]

    # Count linked metabolites per layer
    layer_counts = {}
    db_labels = {
        "mpi": ("Proteins (MPI)", "#3182ce", "fas fa-cubes"),
        "mei": ("Enzymes (MEI)", "#38a169", "fas fa-flask"),
        "mdi": ("Diseases (MDI)", "#e53e3e", "fas fa-heartbeat"),
        "mmi": ("Microbes (MMI)", "#38a169", "fas fa-bacterium"),
        "mdri": ("Drugs (MDrI)", "#805ad5", "fas fa-pills"),
        "mgwas": ("SNPs (mGWAS)", "#319795", "fas fa-map-marker-alt"),
        "mgi": ("Genes (MGI)", "#d69e2e", "fas fa-dna"),
    }

    # For non-metabolite entities, the primary interest is the linked metabolites
    met_col_map = {
        "mpi": "Metabolite Name",
        "mei": "Metabolite_Name",
        "mdi": "Metabolite_Name",
        "mmi": "Metabolite_Name",
        "mdri": "Metabolite_Name",
        "mgwas": "Metabolite_Name",
        "mgi": "Metabolite_Name",
    }

    counts = {}
    metabolite_set = set()
    total = 0
    for db_key in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"]:
        df = data[db_key]
        n = len(df)
        if n > 0:
            counts[db_key] = n
            total += n
            mc = met_col_map.get(db_key, "Metabolite_Name")
            if mc in df.columns:
                metabolite_set.update(df[mc].dropna().unique())

    n_metabolites = len(metabolite_set)

    # ═══════════════════════════════════════════════════════════
    # ENTITY HEADER
    # ═══════════════════════════════════════════════════════════
    # Look up CoreMet ID
    coremetdb_id = ""
    try:
        from app.services.entity_registry import lookup_id
        coremetdb_id = lookup_id(ename, etype) or ""
        if not coremetdb_id and eid:
            coremetdb_id = lookup_id(eid) or ""
    except Exception:
        pass

    ext_links = []
    if eid and eid != "nan" and cfg.get("external_url"):
        ext_links.append(html.A(
            cfg["external_label"],
            href=cfg["external_url"].format(eid),
            target="_blank", className="me-3",
            style={"fontSize": "0.8rem", "fontWeight": "500"},
        ))

    entity_header = html.Div([
        html.Div([
            html.H2(ename, style={
                "fontSize": "1.8rem", "fontWeight": "700", "color": "#1a202c",
                "marginBottom": "4px", "marginRight": "12px", "display": "inline",
            }),
            html.Span(cfg["label"], style={
                "display": "inline-block", "padding": "3px 12px", "borderRadius": "20px",
                "fontSize": "0.75rem", "fontWeight": "600", "color": cfg["color"],
                "background": f"{cfg['color']}14",
                "border": f"1px solid {cfg['color']}30",
                "verticalAlign": "middle",
            }),
        ]),
        html.Div([
            html.Span(coremetdb_id, style={
                "fontSize": "0.85rem", "color": "#3182ce", "fontWeight": "600",
                "fontFamily": "Arial, Helvetica, sans-serif", "marginRight": "12px",
            }) if coremetdb_id else html.Span(),
            html.Span(eid, style={
                "fontSize": "0.85rem", "color": "#718096", "fontFamily": "Arial, Helvetica, sans-serif",
            }) if eid and eid != "nan" else html.Span(),
        ], style={"marginBottom": "8px"}),
        html.Div(ext_links, style={"display": "flex", "alignItems": "center", "flexWrap": "wrap"}) if ext_links else html.Div(),
    ], style={"marginBottom": "24px"})

    # ═══════════════════════════════════════════════════════════
    # SUMMARY CARDS ROW
    # ═══════════════════════════════════════════════════════════
    summary_items = [
        ("Metabolites", n_metabolites, _MET_COLOR, "fas fa-atom"),
    ]
    for db_key, (label, color, icon) in db_labels.items():
        if db_key in counts and counts[db_key] > 0:
            summary_items.append((label, counts[db_key], color, icon))

    summary_cards = dbc.Row([
        dbc.Col(
            html.Div([
                html.Div([
                    html.I(className=icon, style={"color": color, "fontSize": "0.85rem", "marginRight": "6px"}),
                    html.Span(f"{cnt:,}", style={"fontSize": "1.2rem", "fontWeight": "700", "color": "#1a202c"}),
                ], style={"display": "flex", "alignItems": "center"}),
                html.Div(label, style={"fontSize": "0.72rem", "color": "#718096", "fontWeight": "500"}),
            ], style={
                "textAlign": "center", "padding": "12px 8px",
                "border": "1px solid #e2e8f0", "borderRadius": "8px",
                "borderLeft": f"3px solid {color}",
            }),
            xs=6, sm=4, md=True, className="mb-2",
        ) for label, cnt, color, icon in summary_items
    ] + [
        dbc.Col(
            html.Div([
                html.Div(f"{total:,}", style={"fontSize": "1.2rem", "fontWeight": "700", "color": "#1a365d"}),
                html.Div("Total Edges", style={"fontSize": "0.72rem", "color": "#718096", "fontWeight": "500"}),
            ], style={
                "textAlign": "center", "padding": "12px 8px",
                "border": "1px solid #e2e8f0", "borderRadius": "8px", "background": "#f8f9fa",
            }),
            xs=6, sm=4, md=True, className="mb-2",
        ),
    ], className="g-2 mb-4")

    # ═══════════════════════════════════════════════════════════
    # TABS
    # ═══════════════════════════════════════════════════════════

    # ── Overview ──
    overview_tab = dbc.Tab(
        html.Div([
            _build_overview(data, cfg, ename, eid, n_metabolites, counts, total),
        ], className="pt-3"),
        label="Overview", tab_id="tab-overview",
    )

    # ── Network ──
    network_tab = dbc.Tab(
        html.Div([
            _build_network(data, cfg, ename),
        ], className="pt-3"),
        label="Network", tab_id="tab-network",
    )

    # ── Interactions ──
    interactions_tab = dbc.Tab(
        html.Div([
            _build_interactions(data, cfg),
        ], className="pt-3"),
        label="Interactions", tab_id="tab-interactions",
    )

    # ── Evidence ──
    evidence_tab = dbc.Tab(
        html.Div([
            _build_evidence(data, counts, total),
        ], className="pt-3"),
        label="Evidence", tab_id="tab-evidence",
    )

    # ── Downloads ──
    downloads_tab = dbc.Tab(
        html.Div([
            _build_downloads(data, cfg, ename, eid),
        ], className="pt-3"),
        label="Downloads", tab_id="tab-downloads",
    )

    return html.Div([
        entity_header,
        summary_cards,
        dbc.Tabs([
            overview_tab, network_tab, interactions_tab, evidence_tab, downloads_tab,
        ], id="entity-detail-tabs", active_tab="tab-overview"),
    ])


# ---------------------------------------------------------------------------
# Overview Tab
# ---------------------------------------------------------------------------

def _build_overview(data, cfg, ename, eid, n_metabolites, counts, total):
    """Summary overview with description, layer distribution chart, top entities."""
    etype = data["etype"]

    # Summary text
    layers_present = [k for k in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"] if not data[k].empty]
    layer_names = {"mpi": "MPI", "mei": "MEI", "mdi": "MDI", "mmi": "MMI", "mdri": "MDrI", "mgwas": "mGWAS", "mgi": "MGI"}
    active_layers = ", ".join(layer_names[l] for l in layers_present)

    summary_text = (
        f"{ename} is a {cfg['label'].lower()} with {total:,} interactions across "
        f"{len(layers_present)} CoreMet layers ({active_layers}), "
        f"connecting to {n_metabolites:,} unique metabolites."
    )

    # Layer distribution bar chart
    bar_labels = []
    bar_values = []
    bar_colors = []
    db_label_map = {
        "mpi": ("MPI", "#3182ce"), "mei": ("MEI", "#38a169"), "mdi": ("MDI", "#e53e3e"),
        "mmi": ("MMI", "#38a169"), "mdri": ("MDrI", "#805ad5"),
        "mgwas": ("mGWAS", "#319795"), "mgi": ("MGI", "#d69e2e"),
    }
    for db_key in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"]:
        if db_key in counts:
            lbl, clr = db_label_map[db_key]
            bar_labels.append(lbl)
            bar_values.append(counts[db_key])
            bar_colors.append(clr)

    dist_chart = html.Div()
    if bar_labels:
        fig = go.Figure(go.Bar(
            x=bar_labels, y=bar_values,
            marker_color=bar_colors,
            text=bar_values, textposition="outside",
        ))
        fig.update_layout(
            height=280, margin=dict(l=20, r=20, t=10, b=40),
            xaxis_title="Interaction Layer", yaxis_title="Edges",
            font=dict(family="Arial, Helvetica", size=11),
        )
        dist_chart = html.Div([
            html.H5("Layer Distribution", className="cm-card-title mb-2"),
            dcc.Graph(figure=fig, config={"displaylogo": False, "displayModeBar": False}),
        ], className="cm-card mb-4")

    # Top connected metabolites
    all_mets = []
    for db_key in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"]:
        df = data[db_key]
        if df.empty:
            continue
        mc = "Metabolite Name" if "Metabolite Name" in df.columns else "Metabolite_Name"
        hc = "HMDB ID" if "HMDB ID" in df.columns else "HMDB_ID"
        if mc in df.columns:
            for _, row in df.iterrows():
                mname = str(row.get(mc, ""))
                hmdb = str(row.get(hc, "")) if hc in df.columns else ""
                if mname and mname != "nan":
                    all_mets.append((mname, hmdb))

    # Count metabolites
    met_counts = {}
    met_hmdb = {}
    for mname, hmdb in all_mets:
        met_counts[mname] = met_counts.get(mname, 0) + 1
        if hmdb and hmdb != "nan":
            met_hmdb[mname] = hmdb

    top_mets_section = html.Div()
    if met_counts:
        top = sorted(met_counts.items(), key=lambda x: -x[1])[:15]
        rows = []
        for mname, cnt in top:
            hmdb = met_hmdb.get(mname, "")
            rows.append(html.Tr([
                html.Td(
                    html.A(mname, href=f"/metabolite?id={hmdb}" if hmdb else f"/metabolite?name={quote_plus(mname)}",
                           style={"color": "#3182ce", "textDecoration": "none", "fontWeight": "500"}),
                ),
                html.Td(hmdb, style={"fontSize": "0.8rem", "color": "#718096", "fontFamily": "Arial, Helvetica, sans-serif"}),
                html.Td(f"{cnt}", style={"fontWeight": "600"}),
            ]))

        top_mets_section = html.Div([
            html.H5("Top Connected Metabolites", className="cm-card-title mb-2"),
            html.Table([
                html.Thead(html.Tr([
                    html.Th("Metabolite"), html.Th("HMDB ID"), html.Th("Edges"),
                ], style={"fontSize": "0.8rem", "color": "#718096"})),
                html.Tbody(rows),
            ], className="table table-sm", style={"fontSize": "0.85rem"}),
        ], className="cm-card mb-4")

    # Identifiers card
    id_card = html.Div([
        html.H5("Identifiers", className="cm-card-title mb-2"),
        html.Div([
            html.Div([
                html.Strong("Name: ", style={"fontSize": "0.85rem"}),
                html.Span(ename, style={"fontSize": "0.85rem"}),
            ], className="mb-1"),
            html.Div([
                html.Strong(f"{cfg['id_col']}: ", style={"fontSize": "0.85rem"}),
                html.A(eid, href=cfg["external_url"].format(eid), target="_blank",
                       style={"fontSize": "0.85rem"}) if eid and eid != "nan" else html.Span("-", style={"fontSize": "0.85rem"}),
            ], className="mb-1"),
            html.Div([
                html.Strong("Type: ", style={"fontSize": "0.85rem"}),
                html.Span(cfg["label"], style={"fontSize": "0.85rem", "color": cfg["color"], "fontWeight": "600"}),
            ]),
        ]),
    ], className="cm-card mb-4")

    return html.Div([
        html.Div([
            html.I(className=f"{cfg['icon']} me-2", style={"color": cfg["color"]}),
            html.Span(summary_text, style={"fontSize": "0.9rem", "color": "#4a5568"}),
        ], className="cm-card mb-4", style={"padding": "16px"}),
        dbc.Row([
            dbc.Col(id_card, md=4),
            dbc.Col(dist_chart, md=8),
        ]),
        top_mets_section,
    ])


# ---------------------------------------------------------------------------
# Network Tab
# ---------------------------------------------------------------------------

def _build_network(data, cfg, ename):
    """Build a cross-layer Cytoscape network centered on the entity.

    Structure (concentric):
      Center: the queried entity
      Ring 1: linked metabolites (1-hop via direct DB edges)
      Ring 2: entities linked to those metabolites in OTHER layers (2-hop cross-layer)

    For example, Disease "Obesity" → metabolites (via MDI) → genes, proteins,
    drugs, microbes, SNPs (via MPI/MEI/MMI/MDrI/MGI/mGWAS).
    """
    MAX_MET_NODES = 25      # metabolites (ring 1)
    MAX_OUTER_PER_LAYER = 8  # entities per cross-layer (ring 2)

    elements = []
    etype = data["etype"]
    center_color = cfg["color"]
    center_shape = _NODE_SHAPES.get(etype, "ellipse")

    # --- Center node ---
    elements.append({
        "data": {
            "id": "center", "label": ename[:30], "full_label": ename,
            "etype": etype, "color": center_color, "shape": center_shape,
            "size": 56, "href": "",
        },
        "classes": "center",
    })

    # --- Ring 1: collect linked metabolites ---
    seen_mets = set()
    met_nodes = []      # (label, hmdb_id, edge_color, layer)
    met_hmdb_set = set()

    for db_key, (layer_label, clr) in {
        "mpi": ("MPI", "#3182ce"), "mei": ("MEI", "#38a169"),
        "mdi": ("MDI", "#e53e3e"), "mmi": ("MMI", "#48bb78"),
        "mdri": ("MDrI", "#805ad5"), "mgwas": ("mGWAS", "#319795"),
        "mgi": ("MGI", "#d69e2e"),
    }.items():
        df = data[db_key]
        if df.empty:
            continue
        nc = "Metabolite Name" if "Metabolite Name" in df.columns else "Metabolite_Name"
        hc = "HMDB ID" if "HMDB ID" in df.columns else "HMDB_ID"
        if nc not in df.columns:
            continue
        subset = df.drop_duplicates(subset=[nc])
        for _, row in subset.iterrows():
            lbl = str(row.get(nc, ""))
            hmdb = str(row.get(hc, "")) if hc in df.columns else ""
            if not lbl or lbl == "nan" or lbl in seen_mets:
                continue
            seen_mets.add(lbl)
            met_nodes.append((lbl, hmdb, clr, layer_label))
            if hmdb and hmdb != "nan":
                met_hmdb_set.add(hmdb)

    # Sort by name, cap
    met_nodes = met_nodes[:MAX_MET_NODES]

    # Add metabolite nodes + edges to center
    for lbl, hmdb, clr, layer_label in met_nodes:
        node_key = f"met-{lbl}"
        href = f"/metabolite?id={hmdb}" if hmdb and hmdb != "nan" else ""
        elements.append({
            "data": {
                "id": node_key, "label": lbl[:22], "full_label": lbl,
                "etype": "metabolite", "color": _MET_COLOR,
                "shape": "ellipse", "size": 32, "href": href,
            },
            "classes": "metabolite",
        })
        elements.append({
            "data": {"source": "center", "target": node_key,
                     "color": clr, "weight": 2.0, "layer": layer_label},
            "classes": "hop1",
        })

    # --- Ring 2: cross-layer traversal from metabolites ---
    hmdb_list = [hmdb for _, hmdb, _, _ in met_nodes if hmdb and hmdb != "nan"]
    outer_nodes = []
    if hmdb_list:
        outer_nodes = _lookup_metabolite_neighbors(
            hmdb_list, exclude_etype=etype, max_per_layer=MAX_OUTER_PER_LAYER,
        )

    # Add outer nodes and edges (metabolite → outer entity)
    etypes_present = set()
    for nd in outer_nodes:
        elements.append({
            "data": {
                "id": nd["id"], "label": nd["label"][:18], "full_label": nd["label"],
                "etype": nd["etype"], "color": nd["color"],
                "shape": nd["shape"], "size": 24,
                "href": "",
            },
            "classes": "outer",
        })
        # Edge from the metabolite to this outer entity
        met_node_key = f"met-{_find_met_label_by_hmdb(met_nodes, nd['met_id'])}"
        if any(e.get("data", {}).get("id") == met_node_key for e in elements):
            elements.append({
                "data": {"source": met_node_key, "target": nd["id"],
                         "color": nd["color"], "weight": 1.0, "layer": nd["layer"]},
                "classes": "hop2",
            })
        etypes_present.add(nd["etype"])

    if len(elements) <= 1:
        return html.Div([
            html.Div([
                html.I(className="fas fa-project-diagram fa-2x mb-3", style={"color": "#a0aec0"}),
                html.P("No interaction data available for network visualization.",
                       style={"color": "#718096"}),
            ], className="text-center p-5"),
        ], className="cm-card mb-4")

    # --- Legend ---
    legend_items = [(cfg["label"], cfg["color"])]
    legend_items.append(("Metabolite", _MET_COLOR))
    for et in sorted(etypes_present):
        lbl = et.capitalize()
        if et == "snp":
            lbl = "SNP"
        legend_items.append((lbl, _ENTITY_COLORS.get(et, "#a0aec0")))

    legend = html.Div([
        html.Span([
            html.Span(style={
                "display": "inline-block", "width": "10px", "height": "10px",
                "borderRadius": _legend_shape(et_lbl), "backgroundColor": c,
                "marginRight": "5px", "verticalAlign": "middle",
            }),
            html.Span(et_lbl, style={
                "fontSize": "0.72rem", "marginRight": "14px", "color": "#4a5568",
                "fontFamily": "Arial, Helvetica, sans-serif",
            }),
        ], style={"display": "inline-flex", "alignItems": "center"})
        for et_lbl, c in legend_items
    ], style={"marginBottom": "8px", "lineHeight": "1.8"})

    # --- Stats line ---
    n_nodes = sum(1 for e in elements if "source" not in e.get("data", {}))
    n_edges = sum(1 for e in elements if "source" in e.get("data", {}))
    stats_text = f"{n_nodes} nodes, {n_edges} edges"
    if outer_nodes:
        stats_text += f" (cross-layer: {len(outer_nodes)} entities via {len(met_nodes)} metabolites)"

    return html.Div([
        html.Div([
            html.H5([
                html.I(className="fas fa-project-diagram me-2", style={"color": "#4a5568"}),
                "Cross-Layer Interaction Network",
            ], style={"fontSize": "1rem", "fontWeight": "600", "marginBottom": "4px",
                      "fontFamily": "Arial, Helvetica, sans-serif"}),
            html.Div(stats_text, style={"fontSize": "0.75rem", "color": "#a0aec0", "marginBottom": "8px"}),
        ]),
        legend,
        html.Div([
            cyto.Cytoscape(
                id="entity-cytoscape",
                elements=elements,
                layout={"name": "preset"},
                stylesheet=_STYLESHEET,
                style={"width": "100%", "height": "560px",
                       "minHeight": "560px", "backgroundColor": "#fdfdfe"},
                responsive=False,
                autoRefreshLayout=False,
            ),
        ], style={
            "border": "1px solid #e2e8f0", "borderRadius": "8px", "overflow": "hidden",
            "position": "relative",
        }),
        html.Div([
            html.Small([
                html.Strong("Center: "), f"{ename}  |  ",
                html.Strong("Ring 1: "), f"metabolites ({len(met_nodes)})  |  ",
                html.Strong("Ring 2: "), f"cross-layer entities ({len(outer_nodes)})",
            ], style={"fontSize": "0.72rem", "color": "#a0aec0",
                      "fontFamily": "Arial, Helvetica, sans-serif"}),
        ], className="mt-2"),
        html.Small(
            "Drag to pan, scroll to zoom. Click a node to highlight connections.",
            className="text-muted d-block", style={"fontSize": "0.7rem"},
        ),
    ], className="cm-card mb-4")


def _find_met_label_by_hmdb(met_nodes, hmdb_id):
    """Find metabolite label for a given HMDB ID from the met_nodes list."""
    for lbl, hmdb, _, _ in met_nodes:
        if hmdb == hmdb_id:
            return lbl
    return hmdb_id


def _legend_shape(label):
    """Return CSS border-radius approximating the node shape for legend display."""
    shapes = {"Metabolite": "50%", "Protein": "3px", "Enzyme": "3px",
              "Gene": "3px", "Disease": "0", "Microbe": "2px",
              "Drug": "2px", "SNP": "0"}
    return shapes.get(label, "50%")


# ---------------------------------------------------------------------------
# Interactions Tab
# ---------------------------------------------------------------------------

def _build_interactions(data, cfg):
    """Build interaction tables grouped by layer, with clickable entity links."""
    sections = []

    # Column → entity detail link mapping
    _ENTITY_LINK_MAP = {
        "Metabolite Name":  {"route": "/metabolite", "param": "id", "id_col": "HMDB ID"},
        "Metabolite_Name":  {"route": "/metabolite", "param": "id", "id_col": "HMDB_ID"},
        "Disease_Name":     {"route": "/disease-detail", "param": "name"},
        "Gene_Symbol":      {"route": "/gene", "param": "name"},
        "Gene Name":        {"route": "/gene", "param": "name"},
        "Gene_Name":        {"route": "/gene", "param": "name"},
        "Mapped_Gene":      {"route": "/gene", "param": "name"},
        "Protein Name":     {"route": "/protein", "param": "name"},
        "Protein_Name":     {"route": "/protein", "param": "name"},
        "Enzyme_Name":      {"route": "/protein", "param": "name"},
        "Drug_Name":        {"route": "/drug", "param": "name"},
        "Microbe_Name":     {"route": "/microbe", "param": "name"},
        "rsID":             {"route": "/snp", "param": "name"},
    }

    _TABLE_CONFIGS = {
        "mpi": {
            "title": "Metabolite–Protein Interactions (MPI)",
            "icon": "fas fa-cubes", "color": "#3182ce",
            "cols": ["Metabolite Name", "HMDB ID", "Protein Name", "Uniprot ID",
                     "Gene Name", "Species", "Pathway_Name", "evidence_type",
                     "Evidence_Source", "confidence"],
        },
        "mei": {
            "title": "Metabolite–Enzyme Interactions (MEI)",
            "icon": "fas fa-flask", "color": "#38a169",
            "cols": ["Metabolite_Name", "HMDB_ID", "Enzyme_Name", "Uniprot_ID",
                     "EC_Number", "Gene_Name", "Species", "Pathway_Name",
                     "evidence_type", "Evidence_Source"],
        },
        "mdi": {
            "title": "Metabolite–Disease Associations (MDI)",
            "icon": "fas fa-heartbeat", "color": "#e53e3e",
            "cols": ["Metabolite_Name", "HMDB_ID", "Disease_Name",
                     "Category", "Association_Type", "evidence_type",
                     "confidence", "Source", "pmid"],
        },
        "mmi": {
            "title": "Metabolite–Microbe Interactions (MMI)",
            "icon": "fas fa-bacterium", "color": "#38a169",
            "cols": ["Metabolite_Name", "HMDB_ID", "Microbe_Name", "Rank",
                     "Relationship_Type", "Tissue", "evidence_type",
                     "Evidence_Level", "Source", "PMID"],
        },
        "mdri": {
            "title": "Metabolite–Drug Interactions (MDrI)",
            "icon": "fas fa-pills", "color": "#805ad5",
            "cols": ["Metabolite_Name", "HMDB_ID", "Drug_Name", "DrugBank_ID",
                     "Interaction_Type", "evidence_type", "Evidence_Level",
                     "Source", "PMID"],
        },
        "mgwas": {
            "title": "Metabolite–SNP Associations (mGWAS)",
            "icon": "fas fa-map-marker-alt", "color": "#319795",
            "cols": ["Metabolite_Name", "HMDB_ID", "rsID", "Chromosome",
                     "Mapped_Gene", "P_Value", "Beta", "Trait",
                     "evidence_type", "Source", "PMID"],
        },
        "mgi": {
            "title": "Metabolite–Gene Interactions (MGI)",
            "icon": "fas fa-dna", "color": "#d69e2e",
            "cols": ["Metabolite_Name", "HMDB_ID", "Gene_Symbol", "Gene_ID",
                     "Organism", "Interaction_Type", "evidence_type",
                     "Source", "PMID"],
        },
    }

    # Source database → verifiable URL mapping
    _SOURCE_URLS = {
        "CTD": "https://ctdbase.org/",
        "KEGG": "https://www.genome.jp/kegg/",
        "Rhea": "https://www.rhea-db.org/",
        "HMDB": "https://hmdb.ca/",
        "gutMGene": "http://bio-annotation.cn/gutmgene/",
        "AGORA2": "https://www.vmh.life/#microbe",
        "DrugBank_cross_ref": "https://go.drugbank.com/",
        "DrugBank_enzyme_bridge": "https://go.drugbank.com/",
        "CoreMet_curated": "https://coremetdb.org/",
        "GWAS_Catalog": "https://www.ebi.ac.uk/gwas/",
        "Shin2014": "https://pubmed.ncbi.nlm.nih.gov/24816252/",
        "Yin2024": "https://pubmed.ncbi.nlm.nih.gov/38120091/",
        "original": "https://www.genome.jp/kegg/",
    }

    for db_key, tcfg in _TABLE_CONFIGS.items():
        df = data[db_key]
        if df.empty:
            continue

        # Filter to available columns
        cols = [c for c in tcfg["cols"] if c in df.columns]
        if not cols:
            continue

        display_df = df[cols].head(500).copy()

        # Convert entity name columns to markdown links
        for col in cols:
            lnk = _ENTITY_LINK_MAP.get(col)
            if not lnk:
                continue
            route = lnk["route"]
            param = lnk["param"]
            id_col = lnk.get("id_col")
            if id_col and id_col in df.columns:
                display_df[col] = df.head(500).apply(
                    lambda r, _c=col, _r=route, _p=param, _ic=id_col: (
                        f"[{r[_c]}]({_r}?{_p}={quote_plus(str(r[_ic]))})"
                        if pd.notna(r[_c]) and str(r[_c]).strip() and pd.notna(r.get(_ic, "")) and str(r.get(_ic, "")).strip()
                        else str(r[_c]) if pd.notna(r[_c]) else ""
                    ), axis=1)
            else:
                display_df[col] = display_df[col].apply(
                    lambda v, _r=route, _p=param: (
                        f"[{v}]({_r}?{_p}={quote_plus(str(v))})"
                        if pd.notna(v) and str(v).strip() and str(v) != "nan" else str(v)
                    ))

        # Convert Source columns to verifiable links
        for col in ["Source", "Evidence_Source", "source"]:
            if col in display_df.columns:
                def _fmt_src(v, _urls=_SOURCE_URLS):
                    if pd.isna(v) or not str(v).strip():
                        return ""
                    s = str(v).strip()
                    url = _urls.get(s)
                    if url:
                        return f"[{s}]({url})"
                    for k, u in _urls.items():
                        if k in s:
                            return f"[{s}]({u})"
                    return s
                display_df[col] = display_df[col].apply(_fmt_src)

        # Convert PMID columns to PubMed links
        for col in ["PMID", "pmid"]:
            if col in display_df.columns:
                def _fmt_pmid(v):
                    if pd.isna(v) or not str(v).strip():
                        return ""
                    s = str(v).strip().rstrip('.0')
                    first = s.split(';')[0].split(',')[0].strip()
                    if first.isdigit():
                        return f"[{first}](https://pubmed.ncbi.nlm.nih.gov/{first}/)"
                    return s
                display_df[col] = display_df[col].apply(_fmt_pmid)

        sections.append(html.Div([
            html.Div([
                html.Span(style={
                    "width": "8px", "height": "8px", "borderRadius": "50%",
                    "background": tcfg["color"], "display": "inline-block", "marginRight": "8px",
                }),
                html.I(className=f"{tcfg['icon']} me-2", style={"color": tcfg["color"]}),
                html.Span(tcfg["title"], style={"fontWeight": "600", "fontSize": "0.95rem"}),
                dbc.Badge(f"{len(df):,}", className="ms-2", pill=True, color="secondary"),
            ], style={"marginBottom": "12px"}),
            dash_table.DataTable(
                data=display_df.to_dict("records"),
                columns=[{"name": c.replace("_", " "), "id": c, "presentation": "markdown"}
                         for c in cols],
                page_size=15,
                sort_action="native",
                filter_action="native",
                export_format="csv",
                style_table={"overflowX": "auto"},
                style_cell={
                    "fontSize": "0.8rem", "padding": "6px 10px",
                    "textAlign": "left", "whiteSpace": "normal",
                    "maxWidth": "200px", "overflow": "hidden", "textOverflow": "ellipsis",
                    "fontFamily": "Arial, Helvetica, sans-serif",
                },
                style_header={
                    "fontWeight": "600", "fontSize": "0.75rem",
                    "color": "#4a5568", "backgroundColor": "#f7fafc",
                },
                style_data_conditional=[
                    {"if": {"row_index": "odd"}, "backgroundColor": "#fafbfc"},
                ],
            ),
        ], className="cm-card mb-4"))

    if not sections:
        sections = [html.Div("No interactions found.", className="text-center text-muted p-5")]

    return html.Div(sections)


# ---------------------------------------------------------------------------
# Evidence Tab
# ---------------------------------------------------------------------------

def _build_evidence(data, counts, total):
    """Build evidence provenance section with distribution charts."""
    charts = []

    # Source database distribution
    source_counts = {}
    evidence_counts = {}
    pmid_set = set()

    for db_key in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"]:
        df = data[db_key]
        if df.empty:
            continue
        # Sources
        for col in ["Source", "Evidence_Source", "source_db"]:
            if col in df.columns:
                for v in df[col].dropna():
                    v = str(v).strip()
                    if v and v != "nan":
                        source_counts[v] = source_counts.get(v, 0) + 1
        # Evidence levels
        for col in ["Evidence_Level", "evidence_type"]:
            if col in df.columns:
                for v in df[col].dropna():
                    v = str(v).strip()
                    if v and v != "nan":
                        evidence_counts[v] = evidence_counts.get(v, 0) + 1
        # PMIDs
        if "PMID" in df.columns:
            for v in df["PMID"].dropna():
                v = str(v).strip()
                if v and v != "nan" and v.isdigit():
                    pmid_set.add(v)

    # Summary stats cards
    charts.append(dbc.Row([
        dbc.Col(html.Div([
            html.Div(f"{total:,}", style={"fontSize": "1.3rem", "fontWeight": "700", "color": "#1a365d"}),
            html.Div("Total Edges", style={"fontSize": "0.75rem", "color": "#718096"}),
        ], className="cm-card text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{len(source_counts):,}", style={"fontSize": "1.3rem", "fontWeight": "700", "color": "#1a365d"}),
            html.Div("Source Databases", style={"fontSize": "0.75rem", "color": "#718096"}),
        ], className="cm-card text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{len(pmid_set):,}", style={"fontSize": "1.3rem", "fontWeight": "700", "color": "#1a365d"}),
            html.Div("PubMed References", style={"fontSize": "0.75rem", "color": "#718096"}),
        ], className="cm-card text-center"), md=3),
        dbc.Col(html.Div([
            html.Div(f"{len(evidence_counts):,}", style={"fontSize": "1.3rem", "fontWeight": "700", "color": "#1a365d"}),
            html.Div("Evidence Types", style={"fontSize": "0.75rem", "color": "#718096"}),
        ], className="cm-card text-center"), md=3),
    ], className="g-3 mb-4"))

    # Source distribution chart
    if source_counts:
        top_src = sorted(source_counts.items(), key=lambda x: -x[1])[:15]
        fig = go.Figure(go.Bar(
            y=[s[0][:40] for s in top_src][::-1],
            x=[s[1] for s in top_src][::-1],
            orientation="h", marker_color="#3182ce",
            text=[s[1] for s in top_src][::-1], textposition="outside",
        ))
        fig.update_layout(
            height=max(200, len(top_src) * 26 + 60),
            margin=dict(l=10, r=30, t=10, b=10),
            xaxis_title="Edges", font=dict(family="Arial, Helvetica", size=11),
        )
        charts.append(html.Div([
            html.H5("Source Database Distribution", className="cm-card-title mb-2"),
            dcc.Graph(figure=fig, config={"displaylogo": False, "displayModeBar": False}),
        ], className="cm-card mb-4"))

    # Evidence type chart
    if evidence_counts:
        top_ev = sorted(evidence_counts.items(), key=lambda x: -x[1])[:10]
        fig2 = go.Figure(go.Pie(
            labels=[e[0] for e in top_ev],
            values=[e[1] for e in top_ev],
            hole=0.4,
        ))
        fig2.update_layout(
            height=300, margin=dict(l=10, r=10, t=10, b=10),
            font=dict(family="Arial, Helvetica", size=11),
        )
        charts.append(html.Div([
            html.H5("Evidence Type Distribution", className="cm-card-title mb-2"),
            dcc.Graph(figure=fig2, config={"displaylogo": False, "displayModeBar": False}),
        ], className="cm-card mb-4"))

    if not charts:
        charts = [html.Div("No evidence data available.", className="text-center text-muted p-5")]

    return html.Div(charts)


# ---------------------------------------------------------------------------
# Downloads Tab
# ---------------------------------------------------------------------------

def _build_downloads(data, cfg, ename, eid):
    """Build export section with CSV download buttons and API endpoint reference."""
    buttons = []

    for db_key, label in [
        ("mpi", "MPI"), ("mei", "MEI"), ("mdi", "MDI"), ("mmi", "MMI"),
        ("mdri", "MDrI"), ("mgwas", "mGWAS"), ("mgi", "MGI"),
    ]:
        if not data[db_key].empty:
            buttons.append(html.Div([
                html.Div([
                    html.Span(label, style={"fontWeight": "600"}),
                    html.Span(f" ({len(data[db_key]):,} rows)", style={"fontSize": "0.8rem", "color": "#718096"}),
                ], style={"marginBottom": "4px"}),
                html.Div("CSV format", style={"fontSize": "0.75rem", "color": "#a0aec0", "marginBottom": "8px"}),
                dbc.Button([
                    html.I(className="fas fa-download me-2"),
                    f"Download {label}",
                ], size="sm", outline=True, color="primary", className="w-100"),
            ], className="cm-card", style={"padding": "16px", "marginBottom": "12px"}))

    api_ref = html.Div([
        html.H5("API Reference", className="cm-card-title mb-2"),
        html.P("Retrieve this entity's data programmatically:", style={"fontSize": "0.85rem", "color": "#718096"}),
        html.Pre(
            f"GET /api/v1/entity/{data['etype']}?name={quote_plus(ename)}",
            style={
                "backgroundColor": "#f7fafc", "padding": "12px", "borderRadius": "6px",
                "fontSize": "0.8rem", "overflowX": "auto",
            },
        ),
        html.Pre(
            f'curl "https://cormet-db.org/api/v1/entity/{data["etype"]}?name={quote_plus(ename)}"',
            style={
                "backgroundColor": "#f7fafc", "padding": "12px", "borderRadius": "6px",
                "fontSize": "0.8rem", "overflowX": "auto",
            },
        ),
    ], className="cm-card mb-4")

    return html.Div([
        html.H5("Export Data", className="cm-card-title mb-3"),
        dbc.Row([dbc.Col(b, md=4) for b in buttons]),
        api_ref,
    ])
