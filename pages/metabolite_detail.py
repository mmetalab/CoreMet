"""
Metabolite Detail Page, URL-addressable detail view for a single metabolite.

Route: /metabolite?id=HMDB0000122  (or /metabolite?name=Glucose)

Aggregates data from MPI, MEI, MDI, MMI, and MDrI databases to show a
comprehensive metabolite page with chemical info, interaction tables,
summary stats, and external links.
"""

from pathlib import Path
from urllib.parse import parse_qs

from dash import dcc, html, dash_table
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import pandas as pd
import plotly.graph_objects as go

from components.page_header import make_page_header


# ---------------------------------------------------------------------------
# RDKit chemical property computation
# ---------------------------------------------------------------------------

def _smiles_to_img_src(smiles: str, size: int = 250) -> str:
    """Render SMILES to a base64 PNG data-URI using RDKit."""
    if not smiles:
        return ""
    try:
        import io, base64
        from rdkit import Chem
        from rdkit.Chem import Draw
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        img = Draw.MolToImage(mol, size=(size, size))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return f"data:image/png;base64,{base64.b64encode(buf.getvalue()).decode()}"
    except Exception:
        return ""


def _compute_chem_properties(smiles: str) -> dict:
    """Compute molecular properties from SMILES using RDKit.

    Returns dict with molecular_formula, exact_mass, mol_weight, logp, hbd, hba,
    rotatable_bonds, tpsa, num_rings, num_heavy_atoms, or empty dict on failure.
    """
    if not smiles:
        return {}
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        return {
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
            "mol_weight": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
        }
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# Colour palette (matches design system)
# ---------------------------------------------------------------------------
# Cached MPI database loader (avoids re-reading CSV on every lookup)
# ---------------------------------------------------------------------------
_mpi_db_cache = None


def _get_mpi_db():
    """Load MPI database once (cached globally like MEI/MMI/MDI/MDrI services)."""
    global _mpi_db_cache
    if _mpi_db_cache is None:
        try:
            from app.config import Config
            mpi_path = Path(Config().MPI_DB_PATH)
            if mpi_path.exists():
                _mpi_db_cache = pd.read_csv(mpi_path, dtype=str).fillna("")
            else:
                _mpi_db_cache = pd.DataFrame()
        except Exception:
            _mpi_db_cache = pd.DataFrame()
    return _mpi_db_cache


# ---------------------------------------------------------------------------
COLORS = {
    "protein":  "#3182ce",
    "enzyme":   "#3182ce",
    "disease":  "#e53e3e",
    "microbe":  "#38a169",
    "drug":     "#805ad5",
    "snp":      "#319795",
    "gene":     "#d69e2e",
}

# Node shape mapping for Cytoscape
NODE_SHAPES = {
    "metabolite": "ellipse",
    "protein":    "round-rectangle",
    "enzyme":     "round-rectangle",
    "disease":    "triangle",
    "microbe":    "hexagon",
    "drug":       "diamond",
    "snp":        "star",
    "gene":       "round-pentagon",
}

# Cytoscape stylesheet, publication-quality, Nature-style
MINI_NETWORK_STYLESHEET = [
    {"selector": "node", "style": {
        "label": "data(label)", "font-size": "9px",
        "font-family": "Arial, Helvetica, sans-serif", "font-weight": "500",
        "text-valign": "center", "text-halign": "center",
        "width": "data(size)", "height": "data(size)",
        "background-color": "data(color)", "background-opacity": 0.92,
        "shape": "data(shape)",
        "text-wrap": "ellipsis", "text-max-width": "80px",
        "border-width": 1.5, "border-color": "#cbd5e0",
        "text-outline-color": "#ffffff", "text-outline-width": 2,
        "text-outline-opacity": 0.9, "color": "#1a202c",
    }},
    {"selector": "node.center", "style": {
        "width": 56, "height": 56, "font-size": "11px", "font-weight": "bold",
        "background-color": "#f6ad55", "border-width": 3, "border-color": "#dd6b20",
        "z-index": 10,
    }},
    {"selector": "edge", "style": {
        "width": 1.8, "line-color": "data(color)", "opacity": 0.55,
        "curve-style": "bezier",
    }},
    {"selector": ":selected", "style": {
        "border-width": 3, "border-color": "#f6e05e", "background-opacity": 1,
    }},
    {"selector": "node:active", "style": {
        "overlay-color": "#f6e05e", "overlay-opacity": 0.15,
    }},
]


# ---------------------------------------------------------------------------
# Cytoscape mini-network builder
# ---------------------------------------------------------------------------

def _build_mini_network(data: dict, met_name: str) -> html.Div:
    """Build an interactive Cytoscape ego-network for the metabolite.

    Shows the metabolite at the center with connections to proteins, enzymes,
    diseases, microbes, drugs, SNPs, and genes across all CoreMet layers.
    """
    MAX_NEIGHBORS = 20  # per category

    elements = []
    # Center node = metabolite
    elements.append({
        "data": {"id": "center", "label": met_name[:28],
                 "full_label": met_name,
                 "color": "#f6ad55", "shape": "ellipse", "size": 56},
        "classes": "center",
    })

    seen = set()
    type_counts = {}

    def _add_nodes(df, id_col, label_col, ntype, color, edge_class):
        if df.empty:
            return
        subset = df.drop_duplicates(subset=[label_col] if label_col in df.columns else None)
        count = 0
        for _, row in subset.head(MAX_NEIGHBORS).iterrows():
            nid = str(row.get(id_col, row.get(label_col, "")))
            lbl = str(row.get(label_col, nid))
            if not lbl or lbl == "nan" or lbl in seen:
                continue
            seen.add(lbl)
            count += 1
            node_key = f"{ntype}-{lbl}"
            elements.append({
                "data": {"id": node_key, "label": lbl[:22],
                         "full_label": lbl,
                         "color": color,
                         "shape": NODE_SHAPES.get(ntype, "ellipse"), "size": 30},
            })
            elements.append({
                "data": {"source": "center", "target": node_key, "color": color},
                "classes": edge_class,
            })
        if count > 0:
            type_counts[ntype] = type_counts.get(ntype, 0) + count

    _add_nodes(data["mpi"], "Uniprot ID", "Protein Name", "protein",
               COLORS["protein"], "mpi-edge")
    _add_nodes(data["mei"], "Uniprot_ID", "Enzyme_Name", "enzyme",
               COLORS["enzyme"], "mei-edge")
    _add_nodes(data["mdi"], "Disease_ID", "Disease_Name", "disease",
               COLORS["disease"], "mdi-edge")
    _add_nodes(data["mmi"], "Taxonomy_ID", "Microbe_Name", "microbe",
               COLORS["microbe"], "mmi-edge")
    _add_nodes(data["mdri"], "DrugBank_ID", "Drug_Name", "drug",
               COLORS["drug"], "mdri-edge")
    _add_nodes(data["mgwas"], "rsID", "rsID", "snp",
               COLORS["snp"], "mgwas-edge")
    _add_nodes(data["mgi"], "Gene_ID", "Gene_Symbol", "gene",
               COLORS["gene"], "mgi-edge")

    if len(elements) <= 1:
        return html.Div()

    # Legend, only show types that are present
    legend_items = []
    label_map = {"protein": "Protein", "enzyme": "Enzyme", "disease": "Disease",
                 "microbe": "Microbe", "drug": "Drug", "snp": "SNP", "gene": "Gene"}
    for ntype in ["protein", "enzyme", "disease", "microbe", "drug", "snp", "gene"]:
        if ntype in type_counts:
            legend_items.append((label_map[ntype], COLORS[ntype], type_counts[ntype]))

    legend = html.Div([
        html.Span([
            html.Span(style={"display": "inline-block", "width": "10px", "height": "10px",
                              "borderRadius": "50%", "backgroundColor": c,
                              "marginRight": "5px", "verticalAlign": "middle"}),
            html.Span(f"{lbl} ({cnt})", style={
                "fontSize": "0.72rem", "marginRight": "14px", "color": "#4a5568",
                "fontFamily": "Arial, Helvetica, sans-serif",
            }),
        ], style={"display": "inline-flex", "alignItems": "center"})
        for lbl, c, cnt in legend_items
    ], style={"marginBottom": "8px", "lineHeight": "1.8"})

    n_nodes = sum(1 for e in elements if "source" not in e.get("data", {}))
    n_edges = sum(1 for e in elements if "source" in e.get("data", {}))

    return html.Div([
        html.Div([
            html.H5([
                html.I(className="fas fa-project-diagram me-2", style={"color": "#4a5568"}),
                "Interaction Network",
            ], style={"fontSize": "1rem", "fontWeight": "600", "marginBottom": "4px",
                      "fontFamily": "Arial, Helvetica, sans-serif"}),
            html.Div(f"{n_nodes} nodes, {n_edges} edges across {len(type_counts)} entity types",
                     style={"fontSize": "0.75rem", "color": "#a0aec0", "marginBottom": "8px"}),
        ]),
        legend,
        html.Div([
            cyto.Cytoscape(
                id="metdetail-cytoscape",
                elements=elements,
                layout={"name": "preset"},
                stylesheet=MINI_NETWORK_STYLESHEET,
                style={"width": "100%", "height": "560px",
                       "minHeight": "560px", "backgroundColor": "#fdfdfe"},
                responsive=False,
                autoRefreshLayout=False,
            ),
        ], style={"border": "1px solid #e2e8f0", "borderRadius": "8px", "overflow": "hidden",
                  "position": "relative"}),
        html.Small(
            f"Showing up to {MAX_NEIGHBORS} neighbors per category. "
            "Drag to pan, scroll to zoom, click a node to highlight.",
            className="text-muted mt-1 d-block", style={"fontSize": "0.7rem"},
        ),
    ], className="cm-card mb-4")


# ---------------------------------------------------------------------------
# Pathway enrichment builder
# ---------------------------------------------------------------------------

def _build_pathway_enrichment(data: dict) -> html.Div:
    """Show enriched KEGG pathways from MPI + MEI protein/enzyme sets."""
    pathway_counts = {}

    # Collect from MPI
    if not data["mpi"].empty and "Pathway_Name" in data["mpi"].columns:
        for pw in data["mpi"]["Pathway_Name"].dropna():
            for p in str(pw).split(";"):
                p = p.strip()
                if p and p != "nan":
                    pathway_counts[p] = pathway_counts.get(p, 0) + 1

    # Collect from MEI
    if not data["mei"].empty and "Pathway_Name" in data["mei"].columns:
        for pw in data["mei"]["Pathway_Name"].dropna():
            for p in str(pw).split(";"):
                p = p.strip()
                if p and p != "nan":
                    pathway_counts[p] = pathway_counts.get(p, 0) + 1

    if not pathway_counts:
        return html.Div()

    # Sort by count descending and take top 15
    top = sorted(pathway_counts.items(), key=lambda x: -x[1])[:15]
    names = [t[0][:50] for t in top]
    vals = [t[1] for t in top]

    fig = go.Figure(go.Bar(
        y=names[::-1], x=vals[::-1], orientation="h",
        marker_color=COLORS["enzyme"],
        text=vals[::-1], textposition="outside",
    ))
    fig.update_layout(
        height=max(200, len(top) * 28 + 60),
        margin=dict(l=10, r=30, t=10, b=10),
        xaxis_title="Number of Interacting Proteins/Enzymes",
        yaxis=dict(tickfont=dict(size=10)),
        font=dict(family="Arial", size=11),
    )

    return html.Div([
        html.H5([
            html.I(className="fas fa-route me-2", style={"color": COLORS["enzyme"]}),
            "Pathway Enrichment",
            dbc.Badge(f"{len(pathway_counts)} pathways", className="ms-2",
                      pill=True, color="success"),
        ], className="cm-card-title mb-3"),
        dcc.Graph(figure=fig, config={"displaylogo": False, "displayModeBar": False}),
        html.Small(
            f"Top {len(top)} pathways by number of associated proteins/enzymes in MPI + MEI.",
            className="text-muted d-block mt-1",
        ),
    ], className="cm-card mb-4")


# ---------------------------------------------------------------------------
# PMID linkifier
# ---------------------------------------------------------------------------

def _linkify_pmids(df: pd.DataFrame) -> pd.DataFrame:
    """Convert PMID columns to PubMed markdown links."""
    df = df.copy()
    if "PMID" in df.columns:
        mask = df["PMID"].astype(str).str.match(r"^\d+$", na=False)
        df.loc[mask, "PMID"] = (
            "[" + df.loc[mask, "PMID"].astype(str) +
            "](https://pubmed.ncbi.nlm.nih.gov/" +
            df.loc[mask, "PMID"].astype(str) + ")"
        )
    if "Source" in df.columns:
        src_mask = df["Source"].astype(str).str.match(r"^PMID:\d+$", na=False)
        if src_mask.any():
            pmid_nums = df.loc[src_mask, "Source"].str.replace("PMID:", "", regex=False)
            df.loc[src_mask, "Source"] = (
                "[" + df.loc[src_mask, "Source"] +
                "](https://pubmed.ncbi.nlm.nih.gov/" + pmid_nums + ")"
            )
    return df


def _linkify_entities(df: pd.DataFrame) -> pd.DataFrame:
    """Convert entity ID/name columns to external database markdown links."""
    df = df.copy()
    _LINK_MAP = {
        "Uniprot ID":   "https://www.uniprot.org/uniprot/{}",
        "Uniprot_ID":   "https://www.uniprot.org/uniprot/{}",
        "DrugBank_ID":  "https://go.drugbank.com/drugs/{}",
        "MeSH_ID":      "https://meshb.nlm.nih.gov/record/ui?ui={}",
        "Taxonomy_ID":  "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={}",
        "EC_Number":    "https://enzyme.expasy.org/EC/{}",
        "Gene_ID":      "https://www.ncbi.nlm.nih.gov/gene/{}",
        "rsID":         "https://www.ncbi.nlm.nih.gov/snp/{}",
    }
    for col, url_tpl in _LINK_MAP.items():
        if col not in df.columns:
            continue
        vals = df[col].astype(str)
        mask = vals.str.strip().ne("") & vals.ne("nan")
        if not mask.any():
            continue
        df.loc[mask, col] = (
            "[" + vals[mask] + "](" + vals[mask].map(lambda v: url_tpl.format(v.strip())) + ")"
        )
    return df


# ---------------------------------------------------------------------------
# PubChem / ChEBI cross-reference mapping (offline CSV from build_id_mapping.py)
# ---------------------------------------------------------------------------

_xref_cache: dict | None = None

def _get_xref_mapping() -> dict:
    """Load HMDB→PubChem/ChEBI mapping CSV. Returns dict[hmdb_id] → {pubchem_cid, chebi_id}."""
    global _xref_cache
    if _xref_cache is not None:
        return _xref_cache
    import csv
    csv_path = Path(__file__).parent / "data" / "mappings" / "hmdb_xref_mapping.csv"
    if not csv_path.exists():
        # Also try project root
        csv_path = Path(__file__).parent.parent / "data" / "mappings" / "hmdb_xref_mapping.csv"
    if not csv_path.exists():
        _xref_cache = {}
        return _xref_cache
    _xref_cache = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            hmdb = row.get("hmdb_id", "").strip()
            if hmdb:
                _xref_cache[hmdb] = {
                    "pubchem_cid": row.get("pubchem_cid", "").strip(),
                    "chebi_id": row.get("chebi_id", "").strip(),
                }
    return _xref_cache


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def _lookup_metabolite(hmdb_id: str = "", name: str = "") -> dict:
    """Look up a single metabolite across all 7 databases.

    Returns dict with keys: name, hmdb_id, smiles, kegg_id, pubchem_cid,
    chebi_id, and DataFrames for mpi, mei, mdi, mmi, mdri, mgwas, mgi.
    """
    result = {
        "name": name or hmdb_id,
        "hmdb_id": hmdb_id,
        "smiles": "",
        "kegg_id": "",
        "pubchem_cid": "",
        "chebi_id": "",
        "mpi": pd.DataFrame(),
        "mei": pd.DataFrame(),
        "mdi": pd.DataFrame(),
        "mmi": pd.DataFrame(),
        "mdri": pd.DataFrame(),
        "mgwas": pd.DataFrame(),
        "mgi": pd.DataFrame(),
    }

    query = hmdb_id.strip() if hmdb_id else name.strip()
    if not query or len(query) < 2:
        return result

    # --- MPI ---
    try:
        mpi = _get_mpi_db()
        if not mpi.empty:
            mask = pd.Series(False, index=mpi.index)
            if hmdb_id:
                if "HMDB ID" in mpi.columns:
                    mask |= mpi["HMDB ID"].astype(str).str.upper() == hmdb_id.upper()
            if name:
                if "Metabolite Name" in mpi.columns:
                    mask |= mpi["Metabolite Name"].astype(str).str.lower() == name.lower()
            result["mpi"] = mpi[mask]
            # Extract info
            if not result["mpi"].empty:
                row = result["mpi"].iloc[0]
                if not result["name"] or result["name"] == hmdb_id:
                    result["name"] = str(row.get("Metabolite Name", ""))
                if not result["hmdb_id"] and "HMDB ID" in mpi.columns:
                    result["hmdb_id"] = str(row["HMDB ID"]) if pd.notna(row["HMDB ID"]) else ""
                if not result["smiles"] and "SMILES" in mpi.columns and pd.notna(row.get("SMILES")):
                    result["smiles"] = str(row["SMILES"])
    except Exception:
        pass

    # --- MEI ---
    try:
        from app.services.mei_service import get_mei_db
        mei = get_mei_db()
        if not mei.empty:
            mask = pd.Series(False, index=mei.index)
            if hmdb_id and "HMDB_ID" in mei.columns:
                mask |= mei["HMDB_ID"].astype(str).str.upper() == hmdb_id.upper()
            if name and "Metabolite_Name" in mei.columns:
                mask |= mei["Metabolite_Name"].astype(str).str.lower() == name.lower()
            result["mei"] = mei[mask]
            if not result["mei"].empty:
                row = result["mei"].iloc[0]
                if not result["name"] or result["name"] == hmdb_id:
                    result["name"] = str(row.get("Metabolite_Name", ""))
                if not result["hmdb_id"] and pd.notna(row.get("HMDB_ID")):
                    result["hmdb_id"] = str(row["HMDB_ID"])
                if not result["smiles"] and pd.notna(row.get("SMILES")):
                    result["smiles"] = str(row["SMILES"])
                if not result["kegg_id"] and pd.notna(row.get("KEGG_Compound")):
                    result["kegg_id"] = str(row["KEGG_Compound"])
    except Exception:
        pass

    # --- MDI ---
    try:
        from app.services.mdi_service import get_mdi_db
        mdi = get_mdi_db()
        if not mdi.empty:
            mask = pd.Series(False, index=mdi.index)
            if hmdb_id and "HMDB_ID" in mdi.columns:
                mask |= mdi["HMDB_ID"].astype(str).str.upper() == hmdb_id.upper()
            if name and "Metabolite_Name" in mdi.columns:
                mask |= mdi["Metabolite_Name"].astype(str).str.lower() == name.lower()
            result["mdi"] = mdi[mask]
            if not result["mdi"].empty:
                row = result["mdi"].iloc[0]
                if not result["name"] or result["name"] == hmdb_id:
                    result["name"] = str(row.get("Metabolite_Name", ""))
                if not result["hmdb_id"] and pd.notna(row.get("HMDB_ID")):
                    result["hmdb_id"] = str(row["HMDB_ID"])
                if not result["smiles"] and pd.notna(row.get("SMILES")):
                    result["smiles"] = str(row["SMILES"])
    except Exception:
        pass

    # --- MMI ---
    try:
        from app.services.mmi_service import get_mmi_db
        mmi = get_mmi_db()
        if not mmi.empty:
            mask = pd.Series(False, index=mmi.index)
            if hmdb_id and "HMDB_ID" in mmi.columns:
                mask |= mmi["HMDB_ID"].astype(str).str.upper() == hmdb_id.upper()
            if name and "Metabolite_Name" in mmi.columns:
                mask |= mmi["Metabolite_Name"].astype(str).str.lower() == name.lower()
            result["mmi"] = mmi[mask]
            if not result["mmi"].empty:
                row = result["mmi"].iloc[0]
                if not result["name"] or result["name"] == hmdb_id:
                    result["name"] = str(row.get("Metabolite_Name", ""))
                if not result["hmdb_id"] and pd.notna(row.get("HMDB_ID")):
                    result["hmdb_id"] = str(row["HMDB_ID"])
                if not result["smiles"] and pd.notna(row.get("SMILES")):
                    result["smiles"] = str(row["SMILES"])
                if not result["kegg_id"] and pd.notna(row.get("KEGG_ID")):
                    result["kegg_id"] = str(row["KEGG_ID"])
                if not result["pubchem_cid"] and pd.notna(row.get("PubChem_CID")):
                    result["pubchem_cid"] = str(row["PubChem_CID"])
                if not result["chebi_id"] and pd.notna(row.get("ChEBI_ID")):
                    result["chebi_id"] = str(row["ChEBI_ID"])
    except Exception:
        pass

    # --- MDrI ---
    try:
        from app.services.mdri_service import get_mdri_db
        mdri = get_mdri_db()
        if not mdri.empty:
            mask = pd.Series(False, index=mdri.index)
            if hmdb_id and "HMDB_ID" in mdri.columns:
                mask |= mdri["HMDB_ID"].astype(str).str.upper() == hmdb_id.upper()
            if name and "Metabolite_Name" in mdri.columns:
                mask |= mdri["Metabolite_Name"].astype(str).str.lower() == name.lower()
            result["mdri"] = mdri[mask]
            if not result["mdri"].empty:
                row = result["mdri"].iloc[0]
                if not result["name"] or result["name"] == hmdb_id:
                    result["name"] = str(row.get("Metabolite_Name", ""))
                if not result["hmdb_id"] and pd.notna(row.get("HMDB_ID")):
                    result["hmdb_id"] = str(row["HMDB_ID"])
                if not result["smiles"] and pd.notna(row.get("SMILES")):
                    result["smiles"] = str(row["SMILES"])
    except Exception:
        pass

    # --- mGWAS ---
    try:
        from app.services.mgwas_service import get_mgwas_db
        mgwas = get_mgwas_db()
        if not mgwas.empty:
            mask = pd.Series(False, index=mgwas.index)
            if hmdb_id and "HMDB_ID" in mgwas.columns:
                mask |= mgwas["HMDB_ID"].astype(str).str.upper() == hmdb_id.upper()
            if name and "Metabolite_Name" in mgwas.columns:
                mask |= mgwas["Metabolite_Name"].astype(str).str.lower() == name.lower()
            result["mgwas"] = mgwas[mask]
    except Exception:
        pass

    # --- MGI ---
    try:
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if not mgi.empty:
            mask = pd.Series(False, index=mgi.index)
            if hmdb_id and "HMDB_ID" in mgi.columns:
                mask |= mgi["HMDB_ID"].astype(str).str.upper() == hmdb_id.upper()
            if name and "Metabolite_Name" in mgi.columns:
                mask |= mgi["Metabolite_Name"].astype(str).str.lower() == name.lower()
            result["mgi"] = mgi[mask]
    except Exception:
        pass

    # --- Xref mapping fallback (PubChem / ChEBI from offline CSV) ---
    hid = result.get("hmdb_id", "")
    if hid and (not result["pubchem_cid"] or not result["chebi_id"]):
        xref = _get_xref_mapping().get(hid, {})
        if not result["pubchem_cid"] and xref.get("pubchem_cid"):
            result["pubchem_cid"] = xref["pubchem_cid"]
        if not result["chebi_id"] and xref.get("chebi_id"):
            result["chebi_id"] = xref["chebi_id"]

    return result


# ---------------------------------------------------------------------------
# Public builder: called by display_page with the URL search string
# ---------------------------------------------------------------------------

def build_page(search: str = "") -> html.Div:
    """Build the complete metabolite detail page from URL search params.

    Called from app/main.py display_page() when pathname == '/metabolite'.
    """
    params = parse_qs((search or "").lstrip("?"))
    hmdb_id = params.get("id", [""])[0].strip()
    name = params.get("name", [""])[0].strip()

    og_tags = []  # Open Graph meta tags for SEO / sharing

    if not hmdb_id and not name:
        content = _empty_state(
            "No metabolite specified. Use a link from the Database page "
            "or provide ?id=HMDB0000122 in the URL."
        )
    else:
        data = _lookup_metabolite(hmdb_id=hmdb_id, name=name)
        total_hits = sum(len(data[k]) for k in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"])
        if total_hits == 0:
            query = hmdb_id or name
            content = _empty_state(
                f'No records found for "{query}" in any CoreMet database.'
            )
        else:
            content = _build_detail_view(data)
            # Build OG meta tags
            met_name = data["name"]
            mid = data["hmdb_id"]
            og_desc = (f"{met_name}, {total_hits:,} interactions across "
                       f"CoreMet databases (MPI, MEI, MDI, MMI, MDrI, mGWAS, MGI).")
            # Note: Dash html.Meta does not support the 'property' attribute
            # needed for true OG tags. Using 'name' attribute for SEO description.
            # For full OG support, set app.index_string with OG meta in <head>.
            og_tags = [
                html.Meta(name="description", content=og_desc),
            ]

    return html.Div([
        *og_tags,
        html.Div([
            make_page_header(
                "Metabolite Detail",
                "Comprehensive view of a metabolite's interactions across all 7 CoreMet databases.",
                [("Home", "/home"), ("Database", "/database"), ("Metabolite Detail", None)],
            ),
            content,
        ], className="cm-page-container"),
    ])


# ---------------------------------------------------------------------------
# Empty state helper
# ---------------------------------------------------------------------------

def _empty_state(message: str):
    suggestions = [
        ("Glucose", "HMDB0000122"),
        ("L-Tryptophan", "HMDB0000929"),
        ("Ornithine", "HMDB0000214"),
        ("Glutamic acid", "HMDB0000148"),
    ]
    links = [
        dcc.Link(name, href=f"/metabolite?id={hid}",
                 style={"color": "#3182ce", "marginRight": "16px"})
        for name, hid in suggestions
    ]
    return html.Div([
        html.Div([
            html.I(className="fas fa-search fa-2x mb-3", style={"color": "#a0aec0"}),
            html.P(message, className="text-muted mb-3"),
            html.P("Try one of these metabolites:", className="text-muted mb-2",
                   style={"fontSize": "0.85rem"}),
            html.Div(links, className="mb-3"),
            dcc.Link(
                dbc.Button([html.I(className="fas fa-arrow-left me-2"), "Back to Database"],
                           className="cm-btn-secondary mt-2"),
                href="/database",
            ),
        ], className="text-center p-5"),
    ], className="cm-card")


# ---------------------------------------------------------------------------
# Build the full detail view
# ---------------------------------------------------------------------------

def _build_detail_view(data: dict):
    """Assemble the full metabolite detail page with tabbed layout.

    Structure:
    - Entity header (name, ID, type badge, aliases, quick external links)
    - Summary cards row (counts per layer)
    - Tabs: Overview | Network | Interactions | Evidence | Downloads
    """
    met_name = data["name"]
    hmdb_id = data["hmdb_id"]
    smiles = data["smiles"]
    kegg_id = data["kegg_id"]
    pubchem_cid = data["pubchem_cid"]
    chebi_id = data["chebi_id"]

    counts = {
        "Protein": len(data["mpi"]),
        "Enzyme": len(data["mei"]),
        "Disease": len(data["mdi"]),
        "Microbe": len(data["mmi"]),
        "Drug": len(data["mdri"]),
        "SNP": len(data["mgwas"]),
        "Gene": len(data["mgi"]),
    }
    total = sum(counts.values())
    chem_props = _compute_chem_properties(smiles)

    # ═══════════════════════════════════════════════════════════
    # ENTITY HEADER
    # ═══════════════════════════════════════════════════════════
    # Look up CoreMet ID
    coremetdb_id = ""
    try:
        from app.services.entity_registry import lookup_id
        coremetdb_id = lookup_id(met_name, "metabolite") or ""
        if not coremetdb_id and hmdb_id:
            coremetdb_id = lookup_id(hmdb_id) or ""
    except Exception:
        pass

    ext_links = []
    if hmdb_id:
        ext_links.append(html.A("HMDB", href=f"https://hmdb.ca/metabolites/{hmdb_id}",
                                target="_blank", className="me-3",
                                style={"fontSize": "0.8rem", "fontWeight": "500"}))
    if kegg_id and kegg_id != "nan":
        ext_links.append(html.A("KEGG", href=f"https://www.genome.jp/entry/{kegg_id}",
                                target="_blank", className="me-3",
                                style={"fontSize": "0.8rem", "fontWeight": "500"}))
    if pubchem_cid and pubchem_cid != "nan":
        ext_links.append(html.A("PubChem", href=f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}",
                                target="_blank", className="me-3",
                                style={"fontSize": "0.8rem", "fontWeight": "500"}))
    if chebi_id and chebi_id != "nan":
        ext_links.append(html.A("ChEBI", href=f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}",
                                target="_blank", className="me-3",
                                style={"fontSize": "0.8rem", "fontWeight": "500"}))

    entity_header = html.Div([
        html.Div([
            html.H2(met_name, style={
                "fontSize": "1.8rem", "fontWeight": "700", "color": "#1a202c",
                "marginBottom": "4px", "marginRight": "12px", "display": "inline",
            }),
            html.Span("Metabolite", style={
                "display": "inline-block", "padding": "3px 12px", "borderRadius": "20px",
                "fontSize": "0.75rem", "fontWeight": "600", "color": "#e27a3f",
                "background": "rgba(226, 122, 63, 0.12)",
                "border": "1px solid rgba(226, 122, 63, 0.25)",
                "verticalAlign": "middle",
            }),
        ]),
        html.Div([
            html.Span(coremetdb_id, style={
                "fontSize": "0.85rem", "color": "#3182ce", "fontWeight": "600",
                "fontFamily": "Arial, Helvetica, sans-serif", "marginRight": "12px",
            }) if coremetdb_id else html.Span(),
            html.Span(hmdb_id, style={
                "fontSize": "0.85rem", "color": "#718096", "fontFamily": "Arial, Helvetica, sans-serif",
            }) if hmdb_id else html.Span(),
        ], style={"marginBottom": "8px"}),
        html.Div(ext_links, style={"display": "flex", "alignItems": "center", "flexWrap": "wrap"}),
    ], style={"marginBottom": "24px"})

    # ═══════════════════════════════════════════════════════════
    # SUMMARY CARDS ROW
    # ═══════════════════════════════════════════════════════════
    layer_info = [
        ("Proteins", counts["Protein"], "#3182ce", "fas fa-cubes"),
        ("Genes", counts["Gene"], "#d69e2e", "fas fa-dna"),
        ("Diseases", counts["Disease"], "#e53e3e", "fas fa-heartbeat"),
        ("Microbes", counts["Microbe"], "#38a169", "fas fa-bacterium"),
        ("Drugs", counts["Drug"], "#805ad5", "fas fa-pills"),
        ("SNPs", counts["SNP"], "#319795", "fas fa-map-marker-alt"),
        ("Enzymes", counts["Enzyme"], "#38a169", "fas fa-flask"),
    ]

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
        ) for label, cnt, color, icon in layer_info if cnt > 0
    ] + [
        dbc.Col(
            html.Div([
                html.Div(f"{total:,}", style={
                    "fontSize": "1.2rem", "fontWeight": "700", "color": "#1a365d",
                }),
                html.Div("Total", style={"fontSize": "0.72rem", "color": "#718096", "fontWeight": "500"}),
            ], style={
                "textAlign": "center", "padding": "12px 8px",
                "border": "1px solid #e2e8f0", "borderRadius": "8px",
                "background": "#f8f9fa",
            }),
            xs=6, sm=4, md=True, className="mb-2",
        ),
    ], className="g-2 mb-4")

    # ═══════════════════════════════════════════════════════════
    # TAB CONTENT BUILDERS
    # ═══════════════════════════════════════════════════════════

    # ── Overview Tab ──────────────────────────────────────────
    overview_children = [
        _build_chemical_info(met_name, hmdb_id, smiles, kegg_id, pubchem_cid, chebi_id, chem_props),
        _build_summary_text(met_name, hmdb_id, counts, total, chem_props),
        _build_summary_card(counts, total),
    ]
    pw = _build_pathway_enrichment(data)
    if pw.children:
        overview_children.append(pw)

    overview_tab = dbc.Tab(
        html.Div(overview_children, className="pt-3"),
        label="Overview",
        tab_id="tab-overview",
    )

    # ── Network Tab ───────────────────────────────────────────
    net = _build_mini_network(data, met_name)
    network_tab = dbc.Tab(
        html.Div([net] if net.children else [
            html.Div("No interaction data to display in the network.",
                     className="text-center text-muted p-5"),
        ], className="pt-3"),
        label="Network",
        tab_id="tab-network",
    )

    # ── Interactions Tab ──────────────────────────────────────
    interaction_sections = []

    if not data["mpi"].empty:
        cols = [c for c in ["Species", "Protein Name", "Gene Name", "Uniprot ID",
                            "Pathway_Name", "Evidence_Source"] if c in data["mpi"].columns]
        mpi_display = _linkify_entities(data["mpi"])
        interaction_sections.append(_build_interaction_table(
            "Metabolite–Protein Interactions (MPI)", "fas fa-cubes", COLORS["protein"],
            mpi_display, cols, len(data["mpi"]), markdown_cols=["Uniprot ID"],
        ))

    if not data["mei"].empty:
        cols = [c for c in ["Species", "Enzyme_Name", "Gene_Name", "Uniprot_ID",
                            "EC_Number", "Pathway_Name", "Evidence_Source"] if c in data["mei"].columns]
        mei_display = _linkify_entities(data["mei"])
        interaction_sections.append(_build_interaction_table(
            "Metabolite–Enzyme Interactions (MEI)", "fas fa-flask", COLORS["enzyme"],
            mei_display, cols, len(data["mei"]), markdown_cols=["Uniprot_ID", "EC_Number"],
        ))

    if not data["mdi"].empty:
        cols = [c for c in ["Disease_Name", "MeSH_ID", "Category",
                            "Association_Type", "Evidence_Level", "Source"] if c in data["mdi"].columns]
        mdi_display = _linkify_entities(data["mdi"])
        interaction_sections.append(_build_interaction_table(
            "Metabolite–Disease Associations (MDI)", "fas fa-heartbeat", COLORS["disease"],
            mdi_display, cols, len(data["mdi"]), markdown_cols=["MeSH_ID"],
        ))

    if not data["mmi"].empty:
        cols = [c for c in ["Microbe_Name", "Taxonomy_ID", "Rank",
                            "Relationship_Type", "Tissue", "Organism",
                            "Evidence_Level", "Source", "PMID"] if c in data["mmi"].columns]
        mmi_display = _linkify_pmids(data["mmi"])
        mmi_display = _linkify_entities(mmi_display)
        interaction_sections.append(_build_interaction_table(
            "Metabolite–Microbe Interactions (MMI)", "fas fa-bacterium", COLORS["microbe"],
            mmi_display, cols, len(data["mmi"]), markdown_cols=["PMID", "Source", "Taxonomy_ID"],
        ))

    if not data["mdri"].empty:
        cols = [c for c in ["Drug_Name", "DrugBank_ID", "Interaction_Type",
                            "Evidence_Level", "Description", "Source", "PMID"] if c in data["mdri"].columns]
        mdri_display = _linkify_pmids(data["mdri"])
        mdri_display = _linkify_entities(mdri_display)
        interaction_sections.append(_build_interaction_table(
            "Metabolite–Drug Interactions (MDrI)", "fas fa-pills", COLORS["drug"],
            mdri_display, cols, len(data["mdri"]), markdown_cols=["PMID", "Source", "DrugBank_ID"],
        ))

    if not data["mgwas"].empty:
        cols = [c for c in ["rsID", "Chromosome", "Position", "Mapped_Gene",
                            "P_Value", "Beta", "Trait", "Source", "PMID"] if c in data["mgwas"].columns]
        mgwas_display = _linkify_pmids(data["mgwas"])
        mgwas_display = _linkify_entities(mgwas_display)
        interaction_sections.append(_build_interaction_table(
            "Metabolite–SNP Associations (mGWAS)", "fas fa-map-marker-alt", COLORS["snp"],
            mgwas_display, cols, len(data["mgwas"]), markdown_cols=["PMID", "rsID"],
        ))

    if not data["mgi"].empty:
        cols = [c for c in ["Gene_Symbol", "Gene_ID", "Organism",
                            "Interaction_Type", "Interaction_Actions",
                            "Source", "PMID"] if c in data["mgi"].columns]
        mgi_display = _linkify_pmids(data["mgi"])
        interaction_sections.append(_build_interaction_table(
            "Metabolite–Gene Interactions (MGI)", "fas fa-dna", COLORS["gene"],
            mgi_display, cols, len(data["mgi"]), markdown_cols=["PMID"],
        ))

    if not interaction_sections:
        interaction_sections = [
            html.Div("No interactions found.", className="text-center text-muted p-5"),
        ]

    interactions_tab = dbc.Tab(
        html.Div(interaction_sections, className="pt-3"),
        label="Interactions",
        tab_id="tab-interactions",
    )

    # ── Evidence Tab ──────────────────────────────────────────
    evidence_tab = dbc.Tab(
        html.Div([_build_evidence_tab(data, counts, total)], className="pt-3"),
        label="Evidence",
        tab_id="tab-evidence",
    )

    # ── Downloads Tab ─────────────────────────────────────────
    downloads_tab = dbc.Tab(
        html.Div([_build_export_section(data, met_name, hmdb_id)], className="pt-3"),
        label="Downloads",
        tab_id="tab-downloads",
    )

    # ═══════════════════════════════════════════════════════════
    # ASSEMBLE PAGE
    # ═══════════════════════════════════════════════════════════
    return html.Div([
        entity_header,
        summary_cards,
        dbc.Tabs([
            overview_tab,
            network_tab,
            interactions_tab,
            evidence_tab,
            downloads_tab,
        ], id="entity-detail-tabs", active_tab="tab-overview"),
    ])


# ---------------------------------------------------------------------------
# Sub-builders
# ---------------------------------------------------------------------------

def _build_chemical_info(met_name, hmdb_id, smiles, kegg_id, pubchem_cid, chebi_id, chem_props):
    """Chemical information card with identifiers, computed properties, and structure."""
    items = []

    # Title row
    items.append(html.H4([
        html.I(className="fas fa-atom me-2", style={"color": COLORS["protein"]}),
        met_name,
    ], className="cm-card-title mb-3"))

    # ── Identifiers section ─────────────────────────────────
    id_rows = []
    _id = lambda label, val, url=None: dbc.Row([
        dbc.Col(html.Strong(label), width=4, lg=3),
        dbc.Col(
            html.A(val, href=url, target="_blank", rel="noopener noreferrer")
            if url else html.Span(val),
            width=8, lg=9,
        ),
    ], className="mb-1", style={"fontSize": "0.85rem"})

    if hmdb_id:
        id_rows.append(_id("HMDB ID", hmdb_id,
                           f"https://hmdb.ca/metabolites/{hmdb_id}"))
    if kegg_id and kegg_id != "nan":
        id_rows.append(_id("KEGG ID", kegg_id,
                           f"https://www.genome.jp/entry/{kegg_id}"))
    if pubchem_cid and pubchem_cid != "nan":
        id_rows.append(_id("PubChem CID", pubchem_cid,
                           f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}"))
    if chebi_id and chebi_id != "nan":
        id_rows.append(_id("ChEBI ID", chebi_id,
                           f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}"))
    if smiles:
        id_rows.append(dbc.Row([
            dbc.Col(html.Strong("SMILES"), width=4, lg=3),
            dbc.Col(html.Code(
                smiles[:120] + ("..." if len(smiles) > 120 else ""),
                style={"fontSize": "0.76rem", "wordBreak": "break-all"},
            ), width=8, lg=9),
        ], className="mb-1", style={"fontSize": "0.85rem"}))

    # ── Computed molecular properties ───────────────────────
    prop_cards = []
    if chem_props:
        _prop_items = [
            ("Molecular Formula", chem_props.get("molecular_formula", "")),
            ("Molecular Weight",  f"{chem_props['mol_weight']} Da" if "mol_weight" in chem_props else ""),
            ("Exact Mass",        f"{chem_props['exact_mass']} Da" if "exact_mass" in chem_props else ""),
            ("LogP",              str(chem_props.get("logp", ""))),
            ("H-Bond Donors",     str(chem_props.get("hbd", ""))),
            ("H-Bond Acceptors",  str(chem_props.get("hba", ""))),
            ("Rotatable Bonds",   str(chem_props.get("rotatable_bonds", ""))),
            ("TPSA",              f"{chem_props['tpsa']} \u00C5\u00B2" if "tpsa" in chem_props else ""),
            ("Rings",             str(chem_props.get("num_rings", ""))),
            ("Aromatic Rings",    str(chem_props.get("num_aromatic_rings", ""))),
            ("Heavy Atoms",       str(chem_props.get("num_heavy_atoms", ""))),
            ("Fraction Csp3",     str(chem_props.get("fraction_csp3", ""))),
        ]
        for label, val in _prop_items:
            if val:
                prop_cards.append(
                    dbc.Col(
                        html.Div([
                            html.Div(val, style={"fontWeight": "600", "fontSize": "0.9rem",
                                                 "color": "#1a365d"}),
                            html.Div(label, style={"fontSize": "0.72rem", "color": "#718096"}),
                        ], style={"textAlign": "center", "padding": "8px 4px",
                                  "border": "1px solid #e2e8f0", "borderRadius": "6px",
                                  "backgroundColor": "#f7fafc"}),
                        xs=6, sm=4, md=3, lg=2, className="mb-2",
                    )
                )

    # Lipinski Rule of 5 check
    lipinski_badge = html.Span()
    if chem_props:
        violations = 0
        if chem_props.get("mol_weight", 0) > 500: violations += 1
        if chem_props.get("logp", 0) > 5: violations += 1
        if chem_props.get("hbd", 0) > 5: violations += 1
        if chem_props.get("hba", 0) > 10: violations += 1
        if violations == 0:
            lipinski_badge = dbc.Badge("Lipinski Ro5: Pass", color="success",
                                       className="me-2", pill=True)
        else:
            lipinski_badge = dbc.Badge(f"Lipinski Ro5: {violations} violation(s)",
                                       color="warning", className="me-2", pill=True)

    # ── External links ──────────────────────────────────────
    ext_links = []
    if hmdb_id:
        ext_links.append(dbc.Button(
            [html.I(className="fas fa-external-link-alt me-1"), "HMDB"],
            href=f"https://hmdb.ca/metabolites/{hmdb_id}",
            target="_blank", rel="noopener noreferrer",
            color="primary", outline=True, size="sm", className="me-2 mb-1",
        ))
    if kegg_id and kegg_id != "nan":
        ext_links.append(dbc.Button(
            [html.I(className="fas fa-external-link-alt me-1"), "KEGG"],
            href=f"https://www.genome.jp/entry/{kegg_id}",
            target="_blank", rel="noopener noreferrer",
            color="success", outline=True, size="sm", className="me-2 mb-1",
        ))
    if pubchem_cid and pubchem_cid != "nan":
        ext_links.append(dbc.Button(
            [html.I(className="fas fa-external-link-alt me-1"), "PubChem"],
            href=f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_cid}",
            target="_blank", rel="noopener noreferrer",
            color="warning", outline=True, size="sm", className="me-2 mb-1",
        ))
    if chebi_id and chebi_id != "nan":
        ext_links.append(dbc.Button(
            [html.I(className="fas fa-external-link-alt me-1"), "ChEBI"],
            href=f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={chebi_id}",
            target="_blank", rel="noopener noreferrer",
            color="info", outline=True, size="sm", className="me-2 mb-1",
        ))

    # ── 2D structure image (RDKit-rendered from SMILES) ─────
    structure_col = html.Div()
    img_src = _smiles_to_img_src(smiles)
    if img_src:
        structure_col = html.Div([
            html.Img(
                src=img_src,
                alt=f"2D molecular structure of {met_name}",
                style={"maxWidth": "200px", "maxHeight": "200px", "border": "1px solid #e2e8f0",
                       "borderRadius": "8px", "padding": "8px", "backgroundColor": "white"},
                title=f"2D structure of {met_name}",
            ),
            html.Div("2D Structure", style={"fontSize": "0.75rem", "color": "#718096",
                                             "marginTop": "4px"}),
        ], className="text-center mb-3")

    # ── Assemble card ───────────────────────────────────────
    card_children = [
        *items,
        # Identifiers + structure side by side
        dbc.Row([
            dbc.Col(id_rows, md=7),
            dbc.Col([structure_col], md=5,
                    className="d-flex align-items-center justify-content-center"),
        ], className="mb-3"),
    ]

    # Molecular properties grid
    if prop_cards:
        card_children.append(html.Hr(style={"margin": "12px 0"}))
        card_children.append(html.H6([
            html.I(className="fas fa-vial me-2"),
            "Molecular Properties",
            html.Span(" "),
            lipinski_badge,
        ], style={"fontSize": "0.9rem", "fontWeight": "600", "marginBottom": "10px"}))
        card_children.append(dbc.Row(prop_cards))

    # External links
    if ext_links:
        card_children.append(html.Hr(style={"margin": "12px 0"}))
        card_children.append(html.Div(ext_links))

    return html.Div(card_children, className="cm-card mb-4")


def _build_summary_text(met_name, hmdb_id, counts, total, chem_props):
    """Build a natural-language summary with biochemical context first, then interactions."""
    # Section 1: Biochemical Context
    bio_parts = [f"**{met_name}**"]
    if hmdb_id:
        bio_parts[0] += f" ({hmdb_id})"

    if chem_props:
        formula = chem_props.get("molecular_formula", "")
        mw = chem_props.get("mol_weight")
        logp = chem_props.get("logp")
        hbd = chem_props.get("hbd")
        hba = chem_props.get("hba")
        tpsa = chem_props.get("tpsa")
        rings = chem_props.get("num_rings", 0)

        if formula and mw:
            bio_parts.append(f"is a small molecule with molecular formula **{formula}** "
                             f"(MW {mw} Da).")
        elif formula:
            bio_parts.append(f"has the molecular formula **{formula}**.")
        else:
            bio_parts.append("is a metabolite catalogued in CoreMet.")

        # Physicochemical profile
        props = []
        if logp is not None:
            props.append(f"LogP={logp:.1f}")
        if tpsa is not None:
            props.append(f"TPSA={tpsa:.1f} \u00C5\u00B2")
        if hbd is not None and hba is not None:
            props.append(f"{hbd} H-bond donors, {hba} acceptors")
        if rings:
            props.append(f"{rings} ring(s)")
        if props:
            bio_parts.append(f"Physicochemical profile: {', '.join(props)}.")

        # Drug-likeness
        if logp is not None and tpsa is not None:
            if logp <= 5 and tpsa <= 140:
                bio_parts.append("Passes Lipinski Rule of 5 (favorable oral bioavailability).")
    else:
        bio_parts.append("is a metabolite catalogued in CoreMet.")

    # Section 2: Interaction Summary
    interaction_parts = []
    type_labels = [
        ("Protein", "protein"), ("Enzyme", "enzyme"), ("Disease", "disease"),
        ("Microbe", "microbe"), ("Drug", "drug"), ("SNP", "SNP (mGWAS)"),
        ("Gene", "gene (MGI)"),
    ]
    for key, label in type_labels:
        if counts.get(key, 0):
            interaction_parts.append(f"**{counts[key]:,}** {label}")

    if interaction_parts:
        bio_parts.append(f"\n\n**Interactions:** "
                         f"{', '.join(interaction_parts)} "
                         f"({total:,} total records across {sum(1 for k,v in counts.items() if v)} database types).")

    text = " ".join(bio_parts)

    return html.Div([
        html.H5([
            html.I(className="fas fa-file-alt me-2"),
            "Biochemical Summary",
        ], className="cm-card-title mb-2"),
        dcc.Markdown(text, style={"fontSize": "0.88rem", "lineHeight": "1.6",
                                   "color": "#2d3748"}),
    ], className="cm-card mb-4")


def _build_summary_card(counts: dict, total: int):
    """Interaction summary with badges and pie chart."""
    color_map = {
        "Protein": COLORS["protein"],
        "Enzyme": COLORS["enzyme"],
        "Disease": COLORS["disease"],
        "Microbe": COLORS["microbe"],
        "Drug": COLORS["drug"],
        "SNP": COLORS["snp"],
        "Gene": COLORS["gene"],
    }
    badge_color_map = {
        "Protein": "primary",
        "Enzyme": "success",
        "Disease": "danger",
        "Microbe": "warning",
        "Drug": "info",
        "SNP": "secondary",
        "Gene": "dark",
    }

    badges = []
    for label, n in counts.items():
        c = badge_color_map.get(label, "secondary")
        badges.append(dbc.Badge(f"{n:,} {label}", color=c if n > 0 else "secondary",
                                className="me-2 mb-1", pill=True,
                                style={"fontSize": "0.85rem"}))

    # Pie chart
    non_zero = {k: v for k, v in counts.items() if v > 0}
    chart = html.Div()
    if non_zero:
        fig = go.Figure(go.Pie(
            labels=list(non_zero.keys()),
            values=list(non_zero.values()),
            marker_colors=[color_map.get(k, "#a0aec0") for k in non_zero],
            hole=0.45,
            textinfo="label+value",
            textfont_size=11,
        ))
        fig.update_layout(
            height=260, margin=dict(l=20, r=20, t=20, b=20),
            showlegend=False,
            font=dict(family="Arial", size=11),
        )
        chart = dcc.Graph(figure=fig, config={"displaylogo": False, "displayModeBar": False})

    return html.Div([
        html.H5([
            html.I(className="fas fa-chart-pie me-2"),
            f"Interaction Summary, {total:,} total",
        ], className="cm-card-title mb-3"),
        html.Div(badges, className="mb-3"),
        chart,
    ], className="cm-card mb-4")


def _build_interaction_table(title, icon, color, df, display_cols, count,
                             markdown_cols=None):
    """Build a collapsible interaction section with DataTable."""
    table_id = f"metdetail-table-{title.split('(')[-1].strip(')')}"
    markdown_cols = set(markdown_cols or [])

    columns = []
    for c in display_cols:
        col_def = {"name": c.replace("_", " "), "id": c}
        if c in markdown_cols:
            col_def["presentation"] = "markdown"
        columns.append(col_def)

    return html.Div([
        html.H5([
            html.I(className=f"{icon} me-2", style={"color": color}),
            title,
            dbc.Badge(f"{count:,}", className="ms-2", pill=True,
                      style={"backgroundColor": color, "color": "white"}),
        ], className="cm-card-title mb-3"),
        dash_table.DataTable(
            data=df[display_cols].head(200).to_dict("records") if display_cols else [],
            columns=columns,
            sort_action="native",
            filter_action="native",
            page_size=10,
            style_cell={
                "textAlign": "left",
                "fontSize": "0.82rem",
                "padding": "6px 10px",
                "fontFamily": "Arial, Helvetica, sans-serif",
                "maxWidth": "250px",
                "overflow": "hidden",
                "textOverflow": "ellipsis",
                "lineHeight": "1.4",
                "verticalAlign": "middle",
            },
            style_header={
                "fontWeight": "600",
                "backgroundColor": "#edf2f7",
                "fontFamily": "Arial, Helvetica, sans-serif",
                "fontSize": "0.82rem",
            },
            style_data_conditional=[
                {"if": {"row_index": "odd"}, "backgroundColor": "#f7fafc"},
            ],
            style_table={"overflowX": "auto"},
            style_filter={
                "fontSize": "0.78rem",
                "padding": "2px 4px",
            },
        ),
        html.Small(
            f"Showing {min(200, count)} of {count:,} records. Use column filters to narrow results.",
            className="text-muted mt-1 d-block",
        ) if count > 200 else html.Span(),
    ], className="cm-card mb-4")


def _build_evidence_tab(data: dict, counts: dict, total: int) -> html.Div:
    """Build the Evidence tab with provenance charts and source stats."""
    # Collect evidence-level and source distributions across all databases
    evidence_levels = {}
    source_dbs = {}
    pmid_count = 0
    total_refs = 0

    for key in ["mpi", "mei", "mdi", "mmi", "mdri", "mgwas", "mgi"]:
        df = data.get(key, pd.DataFrame())
        if df.empty:
            continue
        # Evidence level
        for col in ["Evidence_Level", "Evidence_Source", "evidence_level"]:
            if col in df.columns:
                for val in df[col].dropna():
                    val = str(val).strip()
                    if val and val != "nan":
                        evidence_levels[val] = evidence_levels.get(val, 0) + 1
        # Source database
        for col in ["Source", "source", "Database"]:
            if col in df.columns:
                for val in df[col].dropna():
                    val = str(val).strip()
                    if val and val != "nan":
                        # Strip PMID prefix for counting
                        if val.startswith("PMID:"):
                            pmid_count += 1
                            total_refs += 1
                        else:
                            source_dbs[val] = source_dbs.get(val, 0) + 1
                            total_refs += 1
        # Count PMIDs
        if "PMID" in df.columns:
            valid_pmids = df["PMID"].astype(str).str.match(r"^\d+$", na=False)
            pmid_count += valid_pmids.sum()

    sections = []

    # Summary stats
    sections.append(html.Div([
        dbc.Row([
            dbc.Col(html.Div([
                html.Div(f"{total:,}", style={"fontSize": "1.5rem", "fontWeight": "700", "color": "#1a365d"}),
                html.Div("Total Interactions", style={"fontSize": "0.75rem", "color": "#718096"}),
            ], style={"textAlign": "center", "padding": "16px"}), md=3),
            dbc.Col(html.Div([
                html.Div(f"{pmid_count:,}", style={"fontSize": "1.5rem", "fontWeight": "700", "color": "#38a169"}),
                html.Div("PubMed References", style={"fontSize": "0.75rem", "color": "#718096"}),
            ], style={"textAlign": "center", "padding": "16px"}), md=3),
            dbc.Col(html.Div([
                html.Div(f"{len(source_dbs):,}", style={"fontSize": "1.5rem", "fontWeight": "700", "color": "#3182ce"}),
                html.Div("Source Databases", style={"fontSize": "0.75rem", "color": "#718096"}),
            ], style={"textAlign": "center", "padding": "16px"}), md=3),
            dbc.Col(html.Div([
                html.Div(f"{sum(1 for v in counts.values() if v > 0)}", style={
                    "fontSize": "1.5rem", "fontWeight": "700", "color": "#805ad5"}),
                html.Div("Active Layers", style={"fontSize": "0.75rem", "color": "#718096"}),
            ], style={"textAlign": "center", "padding": "16px"}), md=3),
        ]),
    ], className="cm-card mb-4"))

    # Evidence type distribution chart
    if evidence_levels:
        ev_sorted = sorted(evidence_levels.items(), key=lambda x: -x[1])[:10]
        fig_ev = go.Figure(go.Bar(
            x=[e[0] for e in ev_sorted],
            y=[e[1] for e in ev_sorted],
            marker_color="#3182ce",
        ))
        fig_ev.update_layout(
            height=250, margin=dict(l=10, r=10, t=10, b=40),
            xaxis_title="Evidence Type", yaxis_title="Count",
            font=dict(family="Inter, Arial", size=11),
        )
        sections.append(html.Div([
            html.H5([
                html.I(className="fas fa-layer-group me-2", style={"color": "#3182ce"}),
                "Evidence Type Distribution",
            ], className="cm-card-title mb-3"),
            dcc.Graph(figure=fig_ev, config={"displaylogo": False, "displayModeBar": False}),
        ], className="cm-card mb-4"))

    # Source database distribution chart
    if source_dbs:
        src_sorted = sorted(source_dbs.items(), key=lambda x: -x[1])[:10]
        fig_src = go.Figure(go.Bar(
            x=[s[0] for s in src_sorted],
            y=[s[1] for s in src_sorted],
            marker_color="#38a169",
        ))
        fig_src.update_layout(
            height=250, margin=dict(l=10, r=10, t=10, b=40),
            xaxis_title="Source Database", yaxis_title="Count",
            font=dict(family="Inter, Arial", size=11),
        )
        sections.append(html.Div([
            html.H5([
                html.I(className="fas fa-database me-2", style={"color": "#38a169"}),
                "Source Database Distribution",
            ], className="cm-card-title mb-3"),
            dcc.Graph(figure=fig_src, config={"displaylogo": False, "displayModeBar": False}),
        ], className="cm-card mb-4"))

    # Layer distribution chart (interaction type breakdown)
    non_zero_counts = {k: v for k, v in counts.items() if v > 0}
    if non_zero_counts:
        color_map = {
            "Protein": COLORS["protein"], "Enzyme": COLORS["enzyme"],
            "Disease": COLORS["disease"], "Microbe": COLORS["microbe"],
            "Drug": COLORS["drug"], "SNP": COLORS["snp"], "Gene": COLORS["gene"],
        }
        fig_layer = go.Figure(go.Pie(
            labels=list(non_zero_counts.keys()),
            values=list(non_zero_counts.values()),
            marker_colors=[color_map.get(k, "#a0aec0") for k in non_zero_counts],
            hole=0.4, textinfo="label+value", textfont_size=11,
        ))
        fig_layer.update_layout(
            height=260, margin=dict(l=20, r=20, t=20, b=20),
            showlegend=False, font=dict(family="Inter, Arial", size=11),
        )
        sections.append(html.Div([
            html.H5([
                html.I(className="fas fa-chart-pie me-2", style={"color": "#805ad5"}),
                "Layer Distribution",
            ], className="cm-card-title mb-3"),
            dcc.Graph(figure=fig_layer, config={"displaylogo": False, "displayModeBar": False}),
        ], className="cm-card mb-4"))

    if not sections:
        sections = [html.Div("No evidence data available.", className="text-center text-muted p-5")]

    return html.Div(sections)


def _build_export_section(data: dict, met_name: str, hmdb_id: str) -> html.Div:
    """Build CSV export buttons + API endpoint info for the Downloads tab."""
    # Build query base: prefer HMDB ID, fall back to name
    q = f"id={hmdb_id}" if hmdb_id else f"name={met_name}"

    export_items = []
    total_rows = 0
    for key, label in [("mpi", "MPI"), ("mei", "MEI"), ("mdi", "MDI"),
                        ("mmi", "MMI"), ("mdri", "MDrI"),
                        ("mgwas", "mGWAS"), ("mgi", "MGI")]:
        df = data.get(key, pd.DataFrame())
        if df.empty:
            continue
        total_rows += len(df)
        export_items.append(
            html.A(
                dbc.Button(
                    [html.I(className="fas fa-download me-1"), f"{label} ({len(df):,})"],
                    color="primary", outline=True, size="sm", className="me-2 mb-1",
                ),
                href=f"/api/v1/export/metabolite?{q}&db={key}",
                target="_blank",
            )
        )

    if total_rows:
        export_items.insert(0,
            html.A(
                dbc.Button(
                    [html.I(className="fas fa-file-csv me-1"),
                     f"All Interactions ({total_rows:,})"],
                    color="primary", size="sm", className="me-2 mb-1",
                ),
                href=f"/api/v1/export/metabolite?{q}&db=all",
                target="_blank",
            )
        )

    if not export_items:
        return html.Div()

    # API endpoint reference
    api_ref = html.Div([
        html.Hr(style={"margin": "16px 0"}),
        html.H6([
            html.I(className="fas fa-code me-2"),
            "API Endpoint",
        ], style={"fontSize": "0.9rem", "fontWeight": "600", "marginBottom": "8px"}),
        html.Pre(
            html.Code(f"curl -s 'http://localhost:8080/api/v1/export/metabolite?{q}&db=all'"),
            style={"backgroundColor": "#2d3748", "color": "#e2e8f0",
                   "padding": "10px", "borderRadius": "6px",
                   "fontSize": "0.78rem", "overflowX": "auto"},
        ),
    ])

    return html.Div([
        html.H5([
            html.I(className="fas fa-download me-2"),
            "Export Data",
        ], className="cm-card-title mb-3"),
        html.Div(export_items, className="d-flex flex-wrap"),
        html.Small(
            "Download interaction data as CSV files for further analysis.",
            className="text-muted d-block mt-2",
        ),
        api_ref,
    ], className="cm-card mb-4")


# ---------------------------------------------------------------------------
# Module-level export (must be after all function definitions)
# ---------------------------------------------------------------------------
page_content = build_page("")
