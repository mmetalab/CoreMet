"""
PDF Report Service — generates unified enrichment reports as downloadable PDFs.

Combines pathway, disease, microbe, and drug enrichment results into a
single well-formatted PDF with tables, summary stats, and colour coding.
"""

import io
import logging
from datetime import datetime

import pandas as pd
from fpdf import FPDF

logger = logging.getLogger(__name__)

# Brand colours
BRAND_NAVY = (26, 54, 93)      # #1a365d
BRAND_BLUE = (49, 130, 206)    # #3182ce
BRAND_RED = (229, 62, 62)      # #e53e3e
BRAND_GREEN = (56, 161, 105)   # #38a169
BRAND_TEAL = (49, 151, 149)    # #319795
BRAND_PURPLE = (128, 90, 213)  # #805ad5
GRAY_200 = (226, 232, 240)
WHITE = (255, 255, 255)


class EnrichmentPDF(FPDF):
    """Custom PDF class with CoreMet branding."""

    def header(self):
        self.set_font("Helvetica", "B", 14)
        self.set_text_color(*BRAND_NAVY)
        self.cell(0, 8, "CoreMet — Enrichment Analysis Report", ln=True, align="C")
        self.set_font("Helvetica", "", 8)
        self.set_text_color(100, 100, 100)
        self.cell(0, 5, f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}", ln=True, align="C")
        self.ln(4)
        # Divider line
        self.set_draw_color(*BRAND_NAVY)
        self.set_line_width(0.5)
        self.line(10, self.get_y(), self.w - 10, self.get_y())
        self.ln(6)

    def footer(self):
        self.set_y(-15)
        self.set_font("Helvetica", "I", 7)
        self.set_text_color(130, 130, 130)
        self.cell(0, 10, f"CoreMet Enrichment Report  |  Page {self.page_no()}/{{nb}}",
                  align="C")

    def section_title(self, title, icon_color=BRAND_NAVY):
        """Render a coloured section header."""
        self.set_font("Helvetica", "B", 12)
        self.set_text_color(*icon_color)
        # Colored bar
        self.set_fill_color(*icon_color)
        self.rect(10, self.get_y(), 3, 7, "F")
        self.set_x(16)
        self.cell(0, 7, title, ln=True)
        self.set_text_color(0, 0, 0)
        self.ln(2)

    def summary_stat(self, label, value, color=BRAND_NAVY):
        """Inline stat: bold value + label."""
        self.set_font("Helvetica", "B", 10)
        self.set_text_color(*color)
        self.cell(20, 6, str(value), align="R")
        self.set_font("Helvetica", "", 9)
        self.set_text_color(80, 80, 80)
        self.cell(50, 6, f"  {label}")

    def render_table(self, df, col_widths=None, header_color=BRAND_NAVY, max_rows=30):
        """Render a DataFrame as a styled PDF table."""
        if df.empty:
            self.set_font("Helvetica", "I", 9)
            self.set_text_color(120, 120, 120)
            self.cell(0, 8, "No results.", ln=True)
            return

        cols = list(df.columns)
        if col_widths is None:
            avail = self.w - 20  # 10 mm margin each side
            col_widths = [avail / len(cols)] * len(cols)

        # Header
        self.set_font("Helvetica", "B", 7)
        self.set_fill_color(*header_color)
        self.set_text_color(*WHITE)
        for i, c in enumerate(cols):
            self.cell(col_widths[i], 6, str(c).replace("_", " ")[:20], border=1,
                      fill=True, align="C")
        self.ln()

        # Rows
        self.set_font("Helvetica", "", 7)
        self.set_text_color(0, 0, 0)
        for row_idx, (_, row) in enumerate(df.head(max_rows).iterrows()):
            if row_idx % 2 == 1:
                self.set_fill_color(245, 248, 252)
                fill = True
            else:
                self.set_fill_color(*WHITE)
                fill = True
            for i, c in enumerate(cols):
                val = row[c]
                if isinstance(val, float):
                    txt = f"{val:.2e}" if abs(val) < 0.001 else f"{val:.3f}"
                else:
                    txt = str(val)[:28]
                self.cell(col_widths[i], 5, txt, border=1, fill=fill, align="C")
            self.ln()

        if len(df) > max_rows:
            self.set_font("Helvetica", "I", 7)
            self.set_text_color(120, 120, 120)
            self.cell(0, 5, f"... showing {max_rows} of {len(df)} rows", ln=True)


def generate_enrichment_pdf(
    pathway_results: pd.DataFrame | None = None,
    disease_results: pd.DataFrame | None = None,
    microbe_results: pd.DataFrame | None = None,
    drug_results: pd.DataFrame | None = None,
    input_summary: dict | None = None,
    fdr_threshold: float = 0.25,
) -> bytes:
    """
    Generate a unified enrichment PDF report.

    Parameters
    ----------
    pathway_results : DataFrame or None
    disease_results : DataFrame or None
    microbe_results : DataFrame or None
    drug_results : DataFrame or None
    input_summary : dict with keys like 'n_metabolites', 'n_proteins', 'source'
    fdr_threshold : float

    Returns
    -------
    bytes : PDF content
    """
    pdf = EnrichmentPDF(orientation="L", format="A4")
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=18)
    pdf.add_page()

    # ── Input summary ──────────────────────────────────────
    info = input_summary or {}
    pdf.section_title("Analysis Overview")
    pdf.set_font("Helvetica", "", 9)
    if info.get("source"):
        pdf.cell(0, 5, f"Data source: {info['source']}", ln=True)
    if info.get("n_metabolites"):
        pdf.cell(0, 5, f"Input metabolites: {info['n_metabolites']}", ln=True)
    if info.get("n_proteins"):
        pdf.cell(0, 5, f"Input proteins: {info['n_proteins']}", ln=True)
    pdf.cell(0, 5, f"FDR threshold: {fdr_threshold}", ln=True)
    pdf.ln(4)

    # ── Quick overview table ───────────────────────────────
    overview_data = []
    for label, df, color_name in [
        ("Pathway", pathway_results, "blue"),
        ("Disease", disease_results, "red"),
        ("Microbe", microbe_results, "green"),
        ("Drug", drug_results, "teal"),
    ]:
        if df is not None and not df.empty:
            sig = len(df[df['FDR'] <= fdr_threshold]) if 'FDR' in df.columns else 0
            overview_data.append({
                "Enrichment Type": label,
                "Total Tested": len(df),
                "Significant": sig,
                "Top Fold": df['Fold_Enrichment'].max() if 'Fold_Enrichment' in df.columns else 0,
            })
    if overview_data:
        ov_df = pd.DataFrame(overview_data)
        pdf.render_table(ov_df, col_widths=[50, 35, 35, 35], header_color=BRAND_NAVY)
        pdf.ln(6)

    # ── Pathway section ────────────────────────────────────
    if pathway_results is not None and not pathway_results.empty:
        pdf.add_page()
        pdf.section_title("Pathway Enrichment (KEGG)", BRAND_BLUE)
        sig = pathway_results[pathway_results['FDR'] <= fdr_threshold] if 'FDR' in pathway_results.columns else pathway_results
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 5, f"{len(sig)} significant of {len(pathway_results)} tested", ln=True)
        pdf.ln(3)
        show_cols = [c for c in ['Pathway', 'Pathway_Name', 'Fold_Enrichment', 'P_value', 'FDR',
                                  'Protein_Count'] if c in pathway_results.columns]
        widths = [30, 80, 30, 35, 35, 25][:len(show_cols)]
        pdf.render_table(pathway_results[show_cols], col_widths=widths, header_color=BRAND_BLUE)

    # ── Disease section ────────────────────────────────────
    if disease_results is not None and not disease_results.empty:
        pdf.add_page()
        pdf.section_title("Disease Enrichment (MDI)", BRAND_RED)
        sig = disease_results[disease_results['FDR'] <= fdr_threshold] if 'FDR' in disease_results.columns else disease_results
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 5, f"{len(sig)} significant of {len(disease_results)} tested", ln=True)
        pdf.ln(3)
        show_cols = [c for c in ['Disease_Name', 'Category', 'Fold_Enrichment', 'P_value', 'FDR',
                                  'Metabolite_Count'] if c in disease_results.columns]
        widths = [60, 50, 30, 35, 35, 30][:len(show_cols)]
        pdf.render_table(disease_results[show_cols], col_widths=widths, header_color=BRAND_RED)

    # ── Microbe section ────────────────────────────────────
    if microbe_results is not None and not microbe_results.empty:
        pdf.add_page()
        pdf.section_title("Microbe Enrichment (MMI)", BRAND_GREEN)
        sig = microbe_results[microbe_results['FDR'] <= fdr_threshold] if 'FDR' in microbe_results.columns else microbe_results
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 5, f"{len(sig)} significant of {len(microbe_results)} tested", ln=True)
        pdf.ln(3)
        show_cols = [c for c in ['Microbe_Name', 'Organism', 'Fold_Enrichment', 'P_value', 'FDR',
                                  'Metabolite_Count'] if c in microbe_results.columns]
        widths = [60, 40, 30, 35, 35, 30][:len(show_cols)]
        pdf.render_table(microbe_results[show_cols], col_widths=widths, header_color=BRAND_GREEN)

    # ── Drug section ───────────────────────────────────────
    if drug_results is not None and not drug_results.empty:
        pdf.add_page()
        pdf.section_title("Drug Enrichment (MDrI)", BRAND_TEAL)
        sig = drug_results[drug_results['FDR'] <= fdr_threshold] if 'FDR' in drug_results.columns else drug_results
        pdf.set_font("Helvetica", "", 9)
        pdf.cell(0, 5, f"{len(sig)} significant of {len(drug_results)} tested", ln=True)
        pdf.ln(3)
        show_cols = [c for c in ['Drug_Name', 'DrugBank_ID', 'Interaction_Type',
                                  'Fold_Enrichment', 'P_value', 'FDR'] if c in drug_results.columns]
        widths = [55, 30, 40, 30, 35, 35][:len(show_cols)]
        pdf.render_table(drug_results[show_cols], col_widths=widths, header_color=BRAND_TEAL)

    # ── Methods footnote ───────────────────────────────────
    pdf.add_page()
    pdf.section_title("Methods")
    pdf.set_font("Helvetica", "", 8)
    methods_text = (
        "Enrichment analysis was performed using Fisher's exact test for each term "
        "(pathway/disease/microbe/drug). P-values were corrected for multiple testing "
        "using the Benjamini-Hochberg false discovery rate (FDR) method. "
        "Fold enrichment is calculated as the ratio of observed overlap to expected "
        "overlap under the null hypothesis. Background sets are derived from the "
        "CoreMet database: KEGG pathways for pathway enrichment, MDI for disease, "
        "MMI for microbe, and MDrI for drug enrichment."
    )
    pdf.multi_cell(0, 4, methods_text)
    pdf.ln(4)
    pdf.set_font("Helvetica", "I", 7)
    pdf.cell(0, 4, "CoreMet — Comprehensive Metabolite Interaction Database", ln=True)
    pdf.cell(0, 4, f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    return pdf.output()
