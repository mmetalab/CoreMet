"""
Help & Documentation page, comprehensive guide for CoreMet.
Tabbed: Getting Started, Interaction Types, Evidence & Scores, Identifiers, FAQ, Contact.

Route: /help
"""

from dash import html, dcc, Input, Output, State, callback
import dash_bootstrap_components as dbc

from components.page_header import make_page_header


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _section(title, children):
    return html.Div([
        html.H4(title, style={"fontWeight": "600", "fontSize": "1.1rem", "color": "#1a202c",
                               "marginBottom": "12px", "marginTop": "24px"}),
        *children,
    ])


def _para(text):
    return html.P(text, style={"fontSize": "0.88rem", "color": "#4a5568", "lineHeight": "1.65", "marginBottom": "12px"})


def _def_row(term, definition, color="#1a202c"):
    return html.Div([
        html.Dt(term, style={"fontWeight": "600", "fontSize": "0.88rem", "color": color, "marginBottom": "2px"}),
        html.Dd(definition, style={"fontSize": "0.85rem", "color": "#718096", "marginBottom": "12px", "marginLeft": "0"}),
    ])


# Entity colors for reference
_ENTITY_COLORS = {
    "Metabolite": "#e27a3f", "Protein": "#3182ce", "Gene": "#d69e2e",
    "Disease": "#e53e3e", "Microbe": "#38a169", "Drug": "#805ad5", "SNP": "#319795",
}


# --------------------------------------------------------------------------
# LAYOUT
# --------------------------------------------------------------------------

layout = html.Div([
    html.Div([
        make_page_header(
            title="Help & Documentation",
            subtitle="Learn how to use CoreMet effectively.",
            breadcrumb_items=[("Home", "/home"), ("Help", None)],
        ),

        dbc.Tabs([
            # ══════════════════════════════════════════════════════
            # TAB 1: Getting Started
            # ══════════════════════════════════════════════════════
            dbc.Tab(
                html.Div([
                    _section("Quick start", [
                        _para("CoreMet is a metabolite-centered knowledge graph that connects "
                              "metabolites to proteins, genes, diseases, microbes, drugs, and SNPs. "
                              "The fastest way to use it:"),
                        html.Ol([
                            html.Li([html.Strong("Search"), ": Type any entity name in the homepage search bar."]),
                            html.Li([html.Strong("View entity page"), ": See a summary of all connected biological layers."]),
                            html.Li([html.Strong("Explore the network"), ": Open the Network tab to see an interactive graph."]),
                            html.Li([html.Strong("Filter and export"), ": Use the Interactions tab to filter and download data."]),
                        ], style={"fontSize": "0.88rem", "color": "#4a5568", "lineHeight": "1.8"}),
                    ]),
                    _section("Supported entity types", [
                        html.Div([
                            html.Div([
                                html.Span(style={
                                    "width": "10px", "height": "10px", "borderRadius": "50%",
                                    "background": color, "display": "inline-block", "marginRight": "8px",
                                }),
                                html.Span(name, style={"fontWeight": "500", "fontSize": "0.88rem"}),
                            ], style={"display": "flex", "alignItems": "center", "padding": "6px 0"})
                            for name, color in _ENTITY_COLORS.items()
                        ]),
                        _para("These colors are used consistently across the entire site, in search results, "
                              "graphs, filters, legends, and statistics."),
                    ]),
                    _section("Recommended workflow", [
                        _para(
                            "1. Land on the homepage → 2. Type a metabolite name (e.g., \"butyrate\") → "
                            "3. Select from autocomplete → 4. View the entity page with cross-layer summary → "
                            "5. Open the Network tab → 6. Filter to high-confidence gene and disease edges → "
                            "7. Export the evidence table."
                        ),
                    ]),
                ], className="pt-3"),
                label="Getting Started",
                tab_id="tab-start",
            ),

            # ══════════════════════════════════════════════════════
            # TAB 2: Interaction Types
            # ══════════════════════════════════════════════════════
            dbc.Tab(
                html.Div([
                    _section("Seven interaction layers", [
                        _para("CoreMet organizes its knowledge graph into seven distinct interaction types. "
                              "Each represents a different biological relationship class."),
                    ]),
                    _def_row("MPI: Metabolite–Protein Interaction",
                             "Physical binding between a metabolite and a protein. Sources include HMDB, "
                             "UniProt, STITCH, and KEGG. Includes both direct binding and transport.",
                             "#3182ce"),
                    _def_row("MEI: Metabolite–Enzyme Interaction",
                             "Enzymatic reactions where a metabolite is a substrate, product, or cofactor. "
                             "Annotated with EC numbers from KEGG and Rhea.",
                             "#3182ce"),
                    _def_row("MDI: Metabolite–Disease Interaction",
                             "Associations between altered metabolite levels and disease phenotypes. "
                             "Curated from HMDB, DisGeNET, CTD, and literature.",
                             "#e53e3e"),
                    _def_row("MMI: Metabolite–Microbe Interaction",
                             "Metabolite production, consumption, or biotransformation by gut microbes. "
                             "Sources include gutMGene, KEGG, and HMDB.",
                             "#38a169"),
                    _def_row("MDrI: Metabolite–Drug Interaction",
                             "Pharmacokinetic and pharmacodynamic interactions between drugs and endogenous "
                             "metabolites. Sources include DrugBank, HMDB, CTD.",
                             "#805ad5"),
                    _def_row("MGI: Metabolite–Gene Interaction",
                             "Genetic associations including mQTL, eQTL-mediated links, and pathway-derived "
                             "gene–metabolite connections.",
                             "#d69e2e"),
                    _def_row("mGWAS: Metabolite–SNP Interaction",
                             "Genome-wide significant associations between SNP genotypes and metabolite "
                             "concentrations, derived from published mGWAS studies.",
                             "#319795"),
                ], className="pt-3"),
                label="Interaction Types",
                tab_id="tab-types",
            ),

            # ══════════════════════════════════════════════════════
            # TAB 3: Evidence & Confidence
            # ══════════════════════════════════════════════════════
            dbc.Tab(
                html.Div([
                    _section("Evidence hierarchy", [
                        _para("Each interaction in CoreMet is annotated with an evidence type that "
                              "reflects how the association was established:"),
                        html.Dl([
                            _def_row("Experimental", "Direct biochemical or biophysical measurement (highest confidence)."),
                            _def_row("Literature", "Curated from published research, supported by PubMed references."),
                            _def_row("Database", "Imported from established public databases (HMDB, KEGG, DrugBank, etc.)."),
                            _def_row("Computational", "Predicted by machine learning or pathway inference (lowest confidence)."),
                        ]),
                    ]),
                    _section("Confidence scores", [
                        _para("Interactions are scored on a 0–1 scale reflecting the strength of evidence:"),
                        html.Div([
                            html.Div(style={
                                "background": "linear-gradient(90deg, #a0aec0, #d69e2e, #38a169)",
                                "height": "20px", "borderRadius": "10px", "marginBottom": "4px",
                            }),
                            dbc.Row([
                                dbc.Col(html.Span("0.0, Low", style={"fontSize": "0.75rem", "color": "#a0aec0"})),
                                dbc.Col(html.Span("0.4, Medium", style={"fontSize": "0.75rem", "color": "#d69e2e"}), className="text-center"),
                                dbc.Col(html.Span("1.0, High", style={"fontSize": "0.75rem", "color": "#38a169"}), className="text-end"),
                            ]),
                        ], style={"maxWidth": "400px", "marginBottom": "16px"}),
                        _para("Score < 0.4: Unlikely or weakly supported. "
                              "Score 0.4–0.7: Possible, may warrant validation. "
                              "Score > 0.7: Likely interaction with strong evidence."),
                    ]),
                    _section("PubMed references", [
                        _para("85.8% of interactions in CoreMet are backed by at least one PubMed ID (PMID). "
                              "PMIDs are shown in the Evidence tab of each entity page and link directly to "
                              "PubMed abstracts."),
                    ]),
                ], className="pt-3"),
                label="Evidence & Scores",
                tab_id="tab-evidence",
            ),

            # ══════════════════════════════════════════════════════
            # TAB 4: Identifiers
            # ══════════════════════════════════════════════════════
            dbc.Tab(
                html.Div([
                    _section("Identifier mapping", [
                        _para("CoreMet uses standardized identifiers for all entity types:"),
                        html.Dl([
                            _def_row("Metabolites", "HMDB IDs (e.g., HMDB0000122). Also accepts KEGG compound IDs, PubChem CIDs, and common names."),
                            _def_row("Proteins", "UniProt accessions (e.g., P00533). Gene names are also supported."),
                            _def_row("Genes", "HGNC gene symbols (e.g., LDHA). Entrez Gene IDs also accepted."),
                            _def_row("Diseases", "Disease names normalized to MeSH or DOID ontologies."),
                            _def_row("Microbes", "NCBI Taxonomy names (e.g., Lactobacillus rhamnosus)."),
                            _def_row("Drugs", "DrugBank IDs or generic drug names."),
                            _def_row("SNPs", "dbSNP rs-numbers (e.g., rs1260326)."),
                        ]),
                    ]),
                    _section("Cross-references", [
                        _para("Entity pages include links to external databases: HMDB, UniProt, KEGG, PubChem, "
                              "DrugBank, dbSNP, NCBI Gene, MeSH, and more."),
                    ]),
                ], className="pt-3"),
                label="Identifiers",
                tab_id="tab-ids",
            ),

            # ══════════════════════════════════════════════════════
            # TAB 5: FAQ
            # ══════════════════════════════════════════════════════
            dbc.Tab(
                html.Div([
                    _section("Frequently asked questions", []),
                    html.Div([
                        dbc.Accordion([
                            dbc.AccordionItem([
                                _para("CoreMet is free and open to all users. No login or registration required."),
                            ], title="Is CoreMet free to use?"),
                            dbc.AccordionItem([
                                _para("Data is updated quarterly. The current release is v1.0, updated in June 2026."),
                            ], title="How often is the database updated?"),
                            dbc.AccordionItem([
                                _para("MPI refers to physical metabolite–protein binding, while MEI refers to "
                                      "enzyme-catalyzed reactions. MPI captures transport and signaling; MEI captures "
                                      "metabolic transformations with EC number annotations."),
                            ], title="What is the difference between MPI and MEI?"),
                            dbc.AccordionItem([
                                _para("MGI captures gene-level associations (e.g., gene expression changes correlated "
                                      "with metabolite levels), while mGWAS captures SNP-level associations from "
                                      "genome-wide studies. MGI is broader; mGWAS is more specific."),
                            ], title="What is the difference between MGI and mGWAS?"),
                            dbc.AccordionItem([
                                _para("Yes. All datasets are available on the Downloads page in CSV format. "
                                      "You can also access data programmatically via the REST API."),
                            ], title="Can I download the raw data?"),
                            dbc.AccordionItem([
                                _para("Use the Network tab on any entity page, or launch the Network Explorer "
                                      "from the Explore section. Both support interactive Cytoscape.js graphs "
                                      "with filtering by layer, confidence, and species."),
                            ], title="How do I visualize the interaction network?"),
                            dbc.AccordionItem([
                                _para("All datasets are released under the CC BY 4.0 license. You may use, "
                                      "redistribute, and adapt them for any purpose, provided you cite the source."),
                            ], title="What is the license?"),
                        ], start_collapsed=True),
                    ]),
                ], className="pt-3"),
                label="FAQ",
                tab_id="tab-faq",
            ),

            # ══════════════════════════════════════════════════════
            # TAB 6: Contact
            # ══════════════════════════════════════════════════════
            dbc.Tab(
                html.Div([
                    _section("Contact", [
                        _para("For questions, bug reports, or collaboration inquiries, "
                              "please open an issue on GitHub or contact the development team:"),
                        html.Ul([
                            html.Li("Cheng Wang, Developer & Lead Researcher"),
                            html.Li([html.A("GitHub Repository",
                                            href="https://github.com/mmetalab/CoreMet", target="_blank")]),
                        ], style={"fontSize": "0.88rem", "color": "#4a5568"}),
                    ]),
                ], className="pt-3"),
                label="Contact",
                tab_id="tab-contact",
            ),
        ],
            id="help-tabs",
            active_tab="tab-start",
        ),

    ], className="cm-page-container"),
])


# Export for main.py routing
page_content = layout
