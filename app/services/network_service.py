"""
Network Service — builds multi-type Cytoscape elements for a query node.

Supports searching by metabolite, protein, disease, microbe, drug, gene, or SNP
across MPI, MDI, MMI, MDrI, MGI, and mGWAS databases. MEI (enzyme) data is
merged into MPI since enzymes are proteins.

Colour scheme (edges):
  - Blue   (#3182ce) : Metabolite–Protein  (MPI)
  - Red    (#e53e3e) : Metabolite–Disease  (MDI)
  - Purple (#805ad5) : Metabolite–Microbe  (MMI)
  - Teal   (#319795) : Metabolite–Drug     (MDrI)
  - Gold   (#d69e2e) : Metabolite–Gene     (MGI)
  - Violet (#805ad5) : Metabolite–SNP      (mGWAS)
"""

from __future__ import annotations

import logging
import threading
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# ── Colour palette ────────────────────────────────────────────────────────
TYPE_COLOURS = {
    "MPI": "#3182ce",
    "MDI": "#e53e3e",
    "MMI": "#805ad5",
    "MDrI": "#319795",
    "MGI": "#d69e2e",
    "mGWAS": "#9f7aea",
}

NODE_COLOURS = {
    "metabolite": "#ed8936",   # orange
    "protein":    "#3182ce",   # blue
    "disease":    "#e53e3e",   # red
    "microbe":    "#805ad5",   # purple
    "drug":       "#319795",   # teal
    "gene":       "#d69e2e",   # gold
    "snp":        "#9f7aea",   # violet
}

# Module-level cache for MPI DB (only loaded once)
_mpi_cache: pd.DataFrame | None = None
_mpi_lock = threading.Lock()


def _get_mpi_db() -> pd.DataFrame:
    """Lazy-load and cache the MPI database (thread-safe)."""
    global _mpi_cache
    if _mpi_cache is not None:
        return _mpi_cache
    with _mpi_lock:
        if _mpi_cache is not None:
            return _mpi_cache
        try:
            from app.config import Config
            _cfg = Config()
            _mpi_cache = pd.read_csv(_cfg.MPI_DB_PATH)
            logger.info(f"MPI database loaded for network service: {len(_mpi_cache):,} records")
        except FileNotFoundError:
            logger.warning("MPI database not found for network service")
            _mpi_cache = pd.DataFrame()
        except Exception:
            _mpi_cache = pd.DataFrame()
    return _mpi_cache


# ── Helper: normalise metabolite query ─────────────────────────────────────

def _match_any(df: pd.DataFrame, query: str, cols: list[str]) -> pd.DataFrame:
    """Return rows where any of *cols* match *query* (case-insensitive substring)."""
    if df.empty:
        return df
    q = query.strip().lower()
    mask = pd.Series(False, index=df.index)
    for col in cols:
        if col in df.columns:
            mask |= df[col].astype(str).str.lower().str.contains(q, na=False)
    return df[mask]


# ── Public API ─────────────────────────────────────────────────────────────

def build_network_elements(
    query: str,
    *,
    query_type: str = "metabolite",
    include_types: Optional[list[str]] = None,
    organism_filter: Optional[str] = None,
    max_nodes: int = 300,
) -> dict:
    """
    Build Cytoscape elements + summary stats for *query*.

    Parameters
    ----------
    query : str
        Search term (metabolite name, protein name/gene, disease name, or microbe name).
    query_type : str
        One of "metabolite", "protein", "disease", "microbe", "drug".
    include_types : list[str] | None
        Subset of ["MPI", "MDI", "MMI", "MDrI"]. ``None`` means all.
    organism_filter : str | None
        Filter rows to this organism (species).
    max_nodes : int
        Hard cap on total nodes to keep the graph responsive.

    Returns
    -------
    dict with keys:
        elements   – list of Cytoscape element dicts
        stats      – {type_counts, total_nodes, total_edges, metabolite_label}
    """
    if not query or len(query.strip()) < 2:
        return {"elements": [], "stats": {}, "type_data": {}}

    types = set(include_types or ["MPI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"])

    # ── Gather matching rows from each DB ──────────────────────
    hits: dict[str, pd.DataFrame] = {}

    if "MPI" in types:
        try:
            mpi = _get_mpi_db()
            if mpi is not None and not mpi.empty:
                if query_type == "metabolite":
                    subset = _match_any(mpi, query, ["Metabolite Name", "HMDB ID"])
                elif query_type == "protein":
                    subset = _match_any(mpi, query, ["Protein Name", "Gene Name", "Uniprot ID"])
                else:
                    subset = _match_any(mpi, query, ["Metabolite Name", "HMDB ID", "Protein Name", "Gene Name"])
                if organism_filter and "Species" in subset.columns:
                    subset = subset[
                        subset["Species"].str.lower() == organism_filter.lower()
                    ]
                hits["MPI"] = subset.head(max_nodes)
        except Exception as exc:
            logger.debug("MPI lookup failed: %s", exc)

        # Also search MEI and merge results as MPI (enzymes are proteins)
        try:
            from app.services.mei_service import get_mei_db
            mei = get_mei_db()
            if not mei.empty:
                if query_type == "metabolite":
                    mei_sub = _match_any(mei, query, ["Metabolite_Name", "HMDB_ID"])
                elif query_type == "protein":
                    mei_sub = _match_any(mei, query, ["Enzyme_Name", "Gene_Name", "Uniprot_ID", "EC_Number"])
                else:
                    mei_sub = _match_any(mei, query, ["Metabolite_Name", "HMDB_ID", "Enzyme_Name"])
                if organism_filter and "Species" in mei_sub.columns:
                    mei_sub = mei_sub[
                        mei_sub["Species"].str.lower() == organism_filter.lower()
                    ]
                if not mei_sub.empty:
                    hits["MEI_as_MPI"] = mei_sub.head(max_nodes)
        except Exception as exc:
            logger.debug("MEI lookup failed: %s", exc)

    if "MDI" in types:
        try:
            from app.services.mdi_service import get_mdi_db
            mdi = get_mdi_db()
            if not mdi.empty:
                if query_type == "disease":
                    subset = _match_any(mdi, query, ["Disease_Name", "Disease_ID"])
                elif query_type == "metabolite":
                    subset = _match_any(mdi, query, ["Metabolite_Name", "HMDB_ID"])
                else:
                    subset = _match_any(mdi, query, ["Metabolite_Name", "HMDB_ID", "Disease_Name"])
                hits["MDI"] = subset.head(max_nodes)
        except Exception as exc:
            logger.debug("MDI lookup failed: %s", exc)

    if "MMI" in types:
        try:
            from app.services.mmi_service import get_mmi_db
            mmi = get_mmi_db()
            if not mmi.empty:
                if query_type == "microbe":
                    subset = _match_any(mmi, query, ["Microbe_Name"])
                elif query_type == "metabolite":
                    subset = _match_any(mmi, query, ["Metabolite_Name", "HMDB_ID"])
                else:
                    subset = _match_any(mmi, query, ["Metabolite_Name", "HMDB_ID", "Microbe_Name"])
                if organism_filter and "Organism" in subset.columns:
                    subset = subset[
                        subset["Organism"].str.lower() == organism_filter.lower()
                    ]
                hits["MMI"] = subset.head(max_nodes)
        except Exception as exc:
            logger.debug("MMI lookup failed: %s", exc)
    if "MDrI" in types:
        try:
            from app.services.mdri_service import get_mdri_db
            mdri = get_mdri_db()
            if not mdri.empty:
                if query_type == "drug":
                    subset = _match_any(mdri, query, ["Drug_Name", "DrugBank_ID"])
                elif query_type == "metabolite":
                    subset = _match_any(mdri, query, ["Metabolite_Name", "HMDB_ID"])
                else:
                    subset = _match_any(mdri, query, ["Metabolite_Name", "HMDB_ID", "Drug_Name", "DrugBank_ID"])
                hits["MDrI"] = subset.head(max_nodes)
        except Exception as exc:
            logger.debug("MDrI lookup failed: %s", exc)

    if "MGI" in types:
        try:
            from app.services.mgi_service import get_mgi_db
            mgi = get_mgi_db()
            if not mgi.empty:
                if query_type == "gene":
                    subset = _match_any(mgi, query, ["Gene_Symbol", "Gene_ID"])
                elif query_type == "metabolite":
                    subset = _match_any(mgi, query, ["Metabolite_Name", "HMDB_ID"])
                else:
                    subset = _match_any(mgi, query, ["Metabolite_Name", "HMDB_ID", "Gene_Symbol"])
                if organism_filter and "Organism" in subset.columns:
                    subset = subset[
                        subset["Organism"].str.lower() == organism_filter.lower()
                    ]
                hits["MGI"] = subset.head(max_nodes)
        except Exception as exc:
            logger.debug("MGI lookup failed: %s", exc)

    if "mGWAS" in types:
        try:
            from app.services.mgwas_service import get_mgwas_db
            mgwas = get_mgwas_db()
            if not mgwas.empty:
                if query_type == "snp":
                    subset = _match_any(mgwas, query, ["rsID"])
                elif query_type == "metabolite":
                    subset = _match_any(mgwas, query, ["Metabolite_Name", "HMDB_ID"])
                else:
                    subset = _match_any(mgwas, query, ["Metabolite_Name", "HMDB_ID", "rsID", "Mapped_Gene"])
                hits["mGWAS"] = subset.head(max_nodes)
        except Exception as exc:
            logger.debug("mGWAS lookup failed: %s", exc)

    # ── Build elements ─────────────────────────────────────────
    nodes: dict[str, dict] = {}   # id → element dict
    edges: list[dict] = []
    edge_seen: set[str] = set()   # track (source, target, type) to dedup
    type_counts = {t: 0 for t in ["MPI", "MDI", "MMI", "MDrI", "MGI", "mGWAS"]}

    # Determine the centre node label and type based on query
    centre_label = query.strip()
    centre_type = query_type
    q_icon = {"metabolite": "fas fa-atom", "protein": "fas fa-dna",
              "disease": "fas fa-heartbeat", "microbe": "fas fa-bacterium",
              "drug": "fas fa-pills", "gene": "fas fa-dna",
              "snp": "fas fa-map-marker-alt"}.get(query_type, "fas fa-atom")

    # Try to find a better label from the data
    if query_type == "metabolite":
        for t, df in hits.items():
            if df.empty:
                continue
            name_col = "Metabolite Name" if t == "MPI" else "Metabolite_Name"
            if name_col in df.columns:
                names = df[name_col].dropna().unique()
                if len(names):
                    centre_label = str(names[0])
                    break
    elif query_type == "protein":
        for t, df in hits.items():
            if df.empty:
                continue
            name_col = "Protein Name" if t == "MPI" else "Enzyme_Name"
            if name_col in df.columns:
                names = df[name_col].dropna().unique()
                if len(names):
                    centre_label = str(names[0])
                    break
    elif query_type == "drug":
        for t, df in hits.items():
            if df.empty:
                continue
            for col in ["Drug_Name", "Drug Name"]:
                if col in df.columns:
                    names = df[col].dropna().unique()
                    if len(names):
                        centre_label = str(names[0])
                        break
            if centre_label != query.strip():
                break

    centre_id = f"{centre_type[:3]}:{centre_label}"
    nodes[centre_id] = {
        "data": {"id": centre_id, "label": centre_label, "node_type": centre_type,
                 "degree": 0},
        "classes": f"{centre_type} centre",
    }

    def _add_node(nid: str, label: str, node_type: str):
        if nid not in nodes:
            nodes[nid] = {
                "data": {"id": nid, "label": label, "node_type": node_type,
                         "degree": 0},
                "classes": node_type,
            }
        nodes[nid]["data"]["degree"] += 1
        nodes[centre_id]["data"]["degree"] += 1

    def _add_edge(source_nid: str, target_nid: str, itype: str, label: str, css_class: str):
        """Add a deduplicated edge."""
        edge_key = f"{source_nid}|{target_nid}|{itype}"
        if edge_key in edge_seen:
            return
        edge_seen.add(edge_key)
        edges.append({
            "data": {
                "source": source_nid, "target": target_nid,
                "interaction_type": itype,
                "tooltip": f"{itype}: {label}",
            },
            "classes": css_class,
        })
        type_counts[itype] = type_counts.get(itype, 0) + 1

    # MPI → build metabolite–protein edges
    mpi_df = hits.get("MPI", pd.DataFrame())
    if not mpi_df.empty:
        for _, row in mpi_df.iterrows():
            met_name = str(row.get("Metabolite Name", ""))
            pname = str(row.get("Protein Name", ""))
            uid = str(row.get("Uniprot ID", ""))
            hmdb = str(row.get("HMDB ID", ""))

            met_nid = f"met:{met_name}" if met_name else f"met:{hmdb}"
            prot_nid = f"prot:{uid}" if uid else f"prot:{pname}"

            if query_type == "metabolite":
                # Centre is metabolite, add protein satellite
                _add_node(prot_nid, pname[:30] if pname else uid, "protein")
                _add_edge(centre_id, prot_nid, "MPI", f"{met_name} ↔ {pname}", "mpi-edge")
            elif query_type == "protein":
                # Centre is protein, add metabolite satellite
                _add_node(met_nid, met_name[:30] if met_name else hmdb, "metabolite")
                _add_edge(centre_id, met_nid, "MPI", f"{pname} ↔ {met_name}", "mpi-edge")
            else:
                # Generic: add both as satellites connected to each other
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_node(prot_nid, pname[:30], "protein")
                _add_edge(met_nid, prot_nid, "MPI", f"{met_name} ↔ {pname}", "mpi-edge")

    # MEI merged into MPI — enzymes shown as proteins
    mei_df = hits.get("MEI_as_MPI", pd.DataFrame())
    if not mei_df.empty:
        for _, row in mei_df.iterrows():
            met_name = str(row.get("Metabolite_Name", ""))
            ename = str(row.get("Enzyme_Name", ""))
            ec = str(row.get("EC_Number", ""))
            uid = str(row.get("Uniprot_ID", ""))
            label = ename if ename and ename != "nan" else ec

            met_nid = f"met:{met_name}"
            prot_nid = f"prot:{uid}" if uid and uid != "nan" else f"prot:{label}"

            if query_type == "metabolite":
                _add_node(prot_nid, label[:30], "protein")
                _add_edge(centre_id, prot_nid, "MPI", f"{met_name} ↔ {label}", "mpi-edge")
            elif query_type == "protein":
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_edge(centre_id, met_nid, "MPI", f"{label} ↔ {met_name}", "mpi-edge")

    # MDI → disease nodes
    mdi_df = hits.get("MDI", pd.DataFrame())
    if not mdi_df.empty:
        for _, row in mdi_df.iterrows():
            dname = str(row.get("Disease_Name", ""))
            met_name = str(row.get("Metabolite_Name", ""))

            dis_nid = f"dis:{dname}"
            met_nid = f"met:{met_name}"

            if query_type == "disease":
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_edge(centre_id, met_nid, "MDI", f"{dname} ↔ {met_name}", "mdi-edge")
            elif query_type == "metabolite":
                _add_node(dis_nid, dname, "disease")
                _add_edge(centre_id, dis_nid, "MDI", f"{met_name} ↔ {dname}", "mdi-edge")
            else:
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_node(dis_nid, dname, "disease")
                _add_edge(met_nid, dis_nid, "MDI", f"{met_name} ↔ {dname}", "mdi-edge")

    # MMI → microbe nodes
    mmi_df = hits.get("MMI", pd.DataFrame())
    if not mmi_df.empty:
        for _, row in mmi_df.iterrows():
            mname = str(row.get("Microbe_Name", ""))
            met_name = str(row.get("Metabolite_Name", ""))

            mic_nid = f"mic:{mname}"
            met_nid = f"met:{met_name}"

            if query_type == "microbe":
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_edge(centre_id, met_nid, "MMI", f"{mname} ↔ {met_name}", "mmi-edge")
            elif query_type == "metabolite":
                _add_node(mic_nid, mname[:30], "microbe")
                _add_edge(centre_id, mic_nid, "MMI", f"{met_name} ↔ {mname}", "mmi-edge")
            else:
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_node(mic_nid, mname[:30], "microbe")
                _add_edge(met_nid, mic_nid, "MMI", f"{met_name} ↔ {mname}", "mmi-edge")

    # MDrI → drug nodes
    mdri_df = hits.get("MDrI", pd.DataFrame())
    if not mdri_df.empty:
        for _, row in mdri_df.iterrows():
            drug_name = str(row.get("Drug_Name", ""))
            drug_dbid = str(row.get("DrugBank_ID", ""))
            met_name = str(row.get("Metabolite_Name", ""))

            drug_nid = f"drug:{drug_dbid}" if drug_dbid and drug_dbid != "nan" else f"drug:{drug_name}"
            met_nid = f"met:{met_name}"

            if query_type == "drug":
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_edge(centre_id, met_nid, "MDrI", f"{drug_name} ↔ {met_name}", "mdri-edge")
            elif query_type == "metabolite":
                _add_node(drug_nid, drug_name[:30], "drug")
                _add_edge(centre_id, drug_nid, "MDrI", f"{met_name} ↔ {drug_name}", "mdri-edge")
            else:
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_node(drug_nid, drug_name[:30], "drug")
                _add_edge(met_nid, drug_nid, "MDrI", f"{met_name} ↔ {drug_name}", "mdri-edge")

    # MGI → gene nodes
    mgi_df = hits.get("MGI", pd.DataFrame())
    if not mgi_df.empty:
        for _, row in mgi_df.iterrows():
            gene_sym = str(row.get("Gene_Symbol", ""))
            gene_id = str(row.get("Gene_ID", ""))
            met_name = str(row.get("Metabolite_Name", ""))

            gene_nid = f"gene:{gene_sym}" if gene_sym and gene_sym != "nan" else f"gene:{gene_id}"
            met_nid = f"met:{met_name}"

            if query_type == "gene":
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_edge(centre_id, met_nid, "MGI", f"{gene_sym} ↔ {met_name}", "mgi-edge")
            elif query_type == "metabolite":
                _add_node(gene_nid, gene_sym[:30], "gene")
                _add_edge(centre_id, gene_nid, "MGI", f"{met_name} ↔ {gene_sym}", "mgi-edge")
            else:
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_node(gene_nid, gene_sym[:30], "gene")
                _add_edge(met_nid, gene_nid, "MGI", f"{met_name} ↔ {gene_sym}", "mgi-edge")

    # mGWAS → SNP nodes
    mgwas_df = hits.get("mGWAS", pd.DataFrame())
    if not mgwas_df.empty:
        for _, row in mgwas_df.iterrows():
            rsid = str(row.get("rsID", ""))
            met_name = str(row.get("Metabolite_Name", ""))

            snp_nid = f"snp:{rsid}"
            met_nid = f"met:{met_name}"

            if query_type == "snp":
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_edge(centre_id, met_nid, "mGWAS", f"{rsid} ↔ {met_name}", "mgwas-edge")
            elif query_type == "metabolite":
                _add_node(snp_nid, rsid, "snp")
                _add_edge(centre_id, snp_nid, "mGWAS", f"{met_name} ↔ {rsid}", "mgwas-edge")
            else:
                _add_node(met_nid, met_name[:30], "metabolite")
                _add_node(snp_nid, rsid, "snp")
                _add_edge(met_nid, snp_nid, "mGWAS", f"{met_name} ↔ {rsid}", "mgwas-edge")

    # ── Enforce max-nodes cap (trim least-connected first) ─────
    if len(nodes) > max_nodes:
        sorted_ids = sorted(
            (nid for nid in nodes if nid != centre_id),
            key=lambda nid: nodes[nid]["data"]["degree"],
        )
        keep = {centre_id} | set(sorted_ids[-(max_nodes - 1):])
        nodes = {k: v for k, v in nodes.items() if k in keep}
        edges = [e for e in edges
                 if e["data"]["source"] in keep and e["data"]["target"] in keep]

    elements = list(nodes.values()) + edges

    stats = {
        "total_nodes": len(nodes),
        "total_edges": len(edges),
        "type_counts": type_counts,
        "query_label": centre_label,
        "query_type": query_type,
    }

    return {"elements": elements, "stats": stats, "type_data": type_counts}


def get_available_organisms() -> list[str]:
    """Return sorted list of all unique organisms across databases."""
    orgs: set[str] = set()
    try:
        mpi = _get_mpi_db()
        if mpi is not None and "Species" in mpi.columns:
            orgs.update(mpi["Species"].dropna().unique())
    except Exception:
        pass
    try:
        from app.services.mei_service import get_mei_db
        mei = get_mei_db()
        if not mei.empty and "Species" in mei.columns:
            orgs.update(mei["Species"].dropna().unique())
    except Exception:
        pass
    try:
        from app.services.mmi_service import get_mmi_db
        mmi = get_mmi_db()
        if not mmi.empty and "Organism" in mmi.columns:
            orgs.update(mmi["Organism"].dropna().unique())
    except Exception:
        pass
    try:
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if not mgi.empty and "Organism" in mgi.columns:
            orgs.update(mgi["Organism"].dropna().unique())
    except Exception:
        pass
    return sorted(orgs)
