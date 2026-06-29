#!/usr/bin/env python3
"""
generate_database_figures.py — CoreMet (NAR Database issue) figures.

Database-centric figure set, all numbers read from data/coremetdb_stats.json and
data/coremetdb_use_cases.json (no hard-coded statistics). Nature/Cell style:
Arial, Type42, 300 DPI, PDF + PNG.

  Figure 1  Architecture & harmonization pipeline (schematic)
  Figure 2  Database content & coverage (counts / entities / overlap / cross-refs)
  Figure 3  Cross-domain mechanism tracing (HCC / metformin / butyrate / GCKR)
  Figure 4  Web interface — generated from live screenshots after the redesign

Output: Submission_NAR/Figures/Figure{1,2,3}_*.{pdf,png}

    conda run -n mpi-vgae python scripts/generate_database_figures.py
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

ROOT = Path(__file__).resolve().parent.parent
STATS = json.loads((ROOT / "data" / "coremetdb_stats.json").read_text())
USECASES = json.loads((ROOT / "data" / "coremetdb_use_cases.json").read_text())
OUT = ROOT.parent / "Submission_NAR" / "Figures"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 300

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "lines.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False,
})

# Entity-type colors (match web theme.css)
EC = {
    "metabolite": "#e27a3f", "protein": "#3182ce", "enzyme": "#00a3c4",
    "disease": "#e53e3e", "microbe": "#38a169", "drug": "#805ad5",
    "gene": "#d69e2e", "SNP": "#319795",
}
# Interaction-type color = its target entity color
TC = {"MPI": EC["protein"], "MEI": EC["enzyme"], "MDI": EC["disease"],
      "MMI": EC["microbe"], "MDrI": EC["drug"], "MGI": EC["gene"], "mGWAS": EC["SNP"]}
INK = "#1a202c"


def _save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=DPI, bbox_inches="tight",
                    facecolor="white", pad_inches=0.05)
    plt.close(fig)
    print("wrote", OUT / f"{name}.png")


def _panel(ax, label):
    ax.text(-0.06, 1.04, label, fontsize=10, fontweight="bold", transform=ax.transAxes)


def _box(ax, x, y, w, h, text, fc, ec=None, fs=6.5, tc="white", lw=0.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc=fc, ec=ec or fc, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, color=tc, zorder=3, fontweight="medium")


def _arrow(ax, p1, p2, color="#718096", lw=0.8, style="-|>"):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=7,
                                 color=color, lw=lw, zorder=1,
                                 shrinkA=2, shrinkB=2))


# ───────────────────────── Figure 1: architecture ──────────────────────────
def figure1():
    fig = plt.figure(figsize=(7.2, 5.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    ax.text(2, 97, "A  Data sources", fontsize=9, fontweight="bold")
    sources = ["KEGG / Rhea", "CTD / HMDB", "gutMGene / AGORA2", "DrugBank", "GWAS Catalog"]
    for i, s in enumerate(sources):
        _box(ax, 2, 86 - i * 5.6, 20, 4.6, s, "#4a5568", fs=6)

    ax.text(2, 50, "B  Harmonization", fontsize=9, fontweight="bold")
    steps = ["Entity normalization\n(HMDB key · InChIKey dedup)",
             "Edge integration\n(dedup · directionality · curated MPI)",
             "Evidence standardization\n(5-level · confidence 0–1)",
             "Cross-reference mapping\n(7 identifier systems)"]
    for i, s in enumerate(steps):
        _box(ax, 2, 40 - i * 9.5, 24, 7.6, s, "#2c7a7b", fs=5.8)
        if i < 3:
            _arrow(ax, (14, 40 - i * 9.5), (14, 40 - (i + 1) * 9.5 + 7.6))

    # interaction types column
    types = list(STATS["databases"].items())
    ax.text(34, 97, "Seven interaction types", fontsize=8, fontweight="bold")
    for i, (k, d) in enumerate(types):
        y = 90 - i * 6.0
        _box(ax, 34, y, 30, 5.0, f"{k}  ·  {d['interactions']:,}", TC[k], fs=6)
        _arrow(ax, (22, 70), (34, y + 2.5), color="#cbd5e0", lw=0.5)

    # knowledge graph hub (C)
    ax.text(70, 97, "C  Knowledge graph", fontsize=9, fontweight="bold")
    cx, cy = 84, 74
    ax.add_patch(Circle((cx, cy), 3.6, fc=EC["metabolite"], ec="white", lw=0.8, zorder=4))
    ax.text(cx, cy, "MET", ha="center", va="center", fontsize=6, color="white", zorder=5, fontweight="bold")
    spokes = [("protein", "Protein"), ("enzyme", "Enzyme"), ("disease", "Disease"),
              ("microbe", "Microbe"), ("drug", "Drug"), ("gene", "Gene"), ("SNP", "SNP")]
    for i, (ekey, lab) in enumerate(spokes):
        ang = math.radians(90 + i * (360 / 7))
        nx, ny = cx + 11 * math.cos(ang), cy + 11 * math.sin(ang)
        _arrow(ax, (cx, cy), (nx, ny), color=EC[ekey], lw=0.8, style="-")
        ax.add_patch(Circle((nx, ny), 2.6, fc=EC[ekey], ec="white", lw=0.6, zorder=4))
        ax.text(nx, ny, lab[:3], ha="center", va="center", fontsize=5, color="white", zorder=5)
    ax.text(cx, cy - 17, f"{STATS['totals']['interactions']:,} edges",
            ha="center", fontsize=6.5, style="italic", color=INK)

    # web modules (D)
    ax.text(70, 46, "D  Web interface", fontsize=9, fontweight="bold")
    mods = ["Evidence-aware browsers (×7)", "Cross-layer metabolite profile",
            "Network exploration", "Bulk + REST download"]
    for i, m in enumerate(mods):
        _box(ax, 70, 37 - i * 8.3, 28, 6.4, m, "#1a365d", fs=6)
    _arrow(ax, (84, 57), (84, 43.4), color="#a0aec0")
    _save(fig, "Figure1_architecture")


# ───────────────────────── Figure 2: content & coverage ────────────────────
def figure2():
    db = STATS["databases"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6))
    fig.subplots_adjust(hspace=0.42, wspace=0.34, left=0.13, right=0.97, top=0.93, bottom=0.1)

    # A interaction counts (log)
    ax = axes[0, 0]; _panel(ax, "A")
    order = sorted(db.items(), key=lambda kv: kv[1]["interactions"])
    names = [k for k, _ in order]; vals = [v["interactions"] for _, v in order]
    ax.barh(names, vals, color=[TC[k] for k in names], height=0.66)
    ax.set_xscale("log"); ax.set_xlabel("Interactions (log scale)")
    ax.set_title("Interaction counts by type", fontsize=7.5)
    for i, v in enumerate(vals):
        ax.text(v * 1.15, i, f"{v:,}", va="center", fontsize=5.4)
    ax.set_xlim(right=vals[-1] * 6)

    # B entity counts
    ax = axes[0, 1]; _panel(ax, "B")
    ents = [("Genes", db["MGI"]["targets"], "gene"), ("SNPs", db["mGWAS"]["targets"], "SNP"),
            ("Metabolites", STATS["totals"]["unique_metabolites_union"], "metabolite"),
            ("Proteins", db["MPI"]["targets"], "protein"), ("Diseases", db["MDI"]["targets"], "disease"),
            ("EC numbers", db["MEI"]["targets"], "enzyme"), ("Drugs", db["MDrI"]["targets"], "drug"),
            ("Microbes", db["MMI"]["targets"], "microbe")]
    ents.sort(key=lambda t: t[1])
    ax.barh([e[0] for e in ents], [e[1] for e in ents],
            color=[EC[e[2]] for e in ents], height=0.66)
    ax.set_xscale("log"); ax.set_xlabel("Unique entities (log scale)")
    ax.set_title("Entity coverage", fontsize=7.5)
    for i, e in enumerate(ents):
        ax.text(e[1] * 1.15, i, f"{e[1]:,}", va="center", fontsize=5.4)
    ax.set_xlim(right=max(e[1] for e in ents) * 6)

    # C cross-type overlap
    ax = axes[1, 0]; _panel(ax, "C")
    ov = STATS["cross_type_overlap"]["metabolites_by_db_count"]
    ks = [str(i) for i in range(1, 8)]
    vs = [ov.get(k, 0) for k in ks]
    cmap = plt.cm.viridis(np.linspace(0.15, 0.9, 7))
    ax.bar(ks, vs, color=cmap, width=0.74)
    ax.set_yscale("log"); ax.set_xlabel("Number of interaction layers"); ax.set_ylabel("Metabolites (log)")
    ax.set_title("Cross-type metabolite overlap", fontsize=7.5)
    for i, v in enumerate(vs):
        ax.text(i, v * 1.2, f"{v:,}", ha="center", fontsize=5.2)
    ax.set_ylim(top=max(vs) * 4)
    ax.annotate(f"{STATS['cross_type_overlap']['metabolites_in_all_types']} in all 7 layers",
                xy=(6, vs[6]), xytext=(3.2, max(vs) * 1.5), fontsize=5.6, color=INK,
                arrowprops=dict(arrowstyle="->", lw=0.5, color="#718096"))

    # D cross-reference coverage
    ax = axes[1, 1]; _panel(ax, "D")
    xr = STATS["cross_references"]
    labels = ["HMDB", "SMILES", "KEGG", "PubChem", "ChEBI"]
    pct = [xr["HMDB"]["pct"], xr["SMILES"]["pct"], xr["KEGG"]["pct"], xr["PubChem"]["pct"], xr["ChEBI"]["pct"]]
    ax.bar(labels, pct, color=["#1a365d", "#2c7a7b", "#4a5568", "#718096", "#a0aec0"], width=0.66)
    ax.set_ylabel("Metabolite coverage (%)"); ax.set_ylim(0, 105)
    ax.set_title("Identifier cross-reference coverage", fontsize=7.5)
    for i, p in enumerate(pct):
        ax.text(i, p + 2, f"{p:g}%", ha="center", fontsize=5.6)

    _save(fig, "Figure2_content_coverage")


# ───────────────────────── Figure 3: use cases ─────────────────────────────
def _hub(ax, title, anchor_color, anchor_label, n_anchor, layers):
    """layers = list of (display, count, entity_color)."""
    ax.set_xlim(-1.25, 1.25); ax.set_ylim(-1.25, 1.25); ax.axis("off")
    ax.set_title(title, fontsize=7.5, pad=2)
    ax.add_patch(Circle((0, 0), 0.30, fc=anchor_color, ec="white", lw=1.0, zorder=5))
    ax.text(0, 0, anchor_label, ha="center", va="center", fontsize=5.6,
            color="white", fontweight="bold", zorder=6)
    if n_anchor is not None:
        ax.text(0, -0.40, f"{n_anchor} metabolites", ha="center", fontsize=5, style="italic")
    n = len(layers)
    for i, (lab, cnt, col) in enumerate(layers):
        ang = math.radians(90 + i * (360 / n))
        x, y = 0.92 * math.cos(ang), 0.92 * math.sin(ang)
        lw = 0.5 + 2.2 * (math.log10(cnt + 1) / 5)
        ax.plot([0, x * 0.78], [0, y * 0.78], color=col, lw=lw, zorder=1, alpha=0.7)
        r = 0.12 + 0.12 * (math.log10(cnt + 1) / 5)
        ax.add_patch(Circle((x, y), r, fc=col, ec="white", lw=0.7, zorder=4))
        ax.text(x, y, lab, ha="center", va="center", fontsize=4.8, color="white", zorder=5)
        ax.text(x * 1.34, y * 1.34, f"{cnt:,}", ha="center", va="center", fontsize=5.4,
                color=INK, fontweight="bold")


def figure3():
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.0))
    fig.subplots_adjust(hspace=0.18, wspace=0.06, left=0.02, right=0.98, top=0.95, bottom=0.02)

    uc = USECASES
    _panel(axes[0, 0], "A")
    _hub(axes[0, 0], "Hepatocellular carcinoma", EC["disease"], "HCC", uc["hcc"]["metabolites"],
         [("Pro", uc["hcc"]["MPI"], EC["protein"]), ("Gene", uc["hcc"]["MGI"], EC["gene"]),
          ("Mic", uc["hcc"]["MMI"], EC["microbe"]), ("Drug", uc["hcc"]["MDrI"], EC["drug"]),
          ("SNP", uc["hcc"]["mGWAS"], EC["SNP"])])

    _panel(axes[0, 1], "B")
    _hub(axes[0, 1], "Metformin", EC["drug"], "MET", uc["metformin"]["metabolites"],
         [("Dis", uc["metformin"]["MDI"], EC["disease"]), ("Mic", uc["metformin"]["MMI"], EC["microbe"]),
          ("Gene", uc["metformin"]["MGI"], EC["gene"]), ("Pro", uc["metformin"]["MPI"], EC["protein"])])

    _panel(axes[1, 0], "C")
    _hub(axes[1, 0], "Butyrate (HMDB0000039)", EC["metabolite"], "BUT", None,
         [("Pro", uc["butyrate"]["MPI"], EC["protein"]), ("Enz", uc["butyrate"]["MEI"], EC["enzyme"]),
          ("Dis", uc["butyrate"]["MDI"], EC["disease"]), ("Mic", uc["butyrate"]["MMI"], EC["microbe"]),
          ("Gene", uc["butyrate"]["MGI"], EC["gene"]), ("Drug", uc["butyrate"]["MDrI"], EC["drug"])])

    _panel(axes[1, 1], "D")
    _hub(axes[1, 1], "GCKR locus (rs1260326)", EC["SNP"], "GCKR", uc["gckr"]["metabolites"],
         [("Pro", uc["gckr"]["MPI"], EC["protein"]), ("Dis", uc["gckr"]["MDI"], EC["disease"]),
          ("Gene", uc["gckr"]["MGI"], EC["gene"]), ("Mic", uc["gckr"]["MMI"], EC["microbe"])])

    _save(fig, "Figure3_use_cases")


# ─────────── Supplementary Figure S1: optional model benchmark (toned down) ──
def figureS1():
    """Per-type link-prediction AUC-ROC for the optional CoreMet-FM models.
    Companion-model results; shown only to support the optional prediction layer."""
    sup = OUT.parent / "Supplementary"
    sup.mkdir(parents=True, exist_ok=True)
    # (type, AUC-ROC mean, SD) — companion model 5-fold CV (see Methods / separate report)
    data = [("MPI", 0.955, 0.004), ("MDI", 0.970, 0.002), ("MMI", 0.991, 0.001),
            ("MDrI", 0.913, 0.007), ("MGI", 0.997, 0.001), ("mGWAS", 0.997, 0.001)]
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    names = [d[0] for d in data]
    means = [d[1] for d in data]; sds = [d[2] for d in data]
    ax.bar(names, means, yerr=sds, color=[TC[n] for n in names], width=0.66,
           error_kw=dict(lw=0.6, capsize=2))
    ax.set_ylim(0.85, 1.0); ax.set_ylabel("AUC-ROC (5-fold CV)")
    ax.set_title("Optional CoreMet-FM link prediction", fontsize=7.5)
    for i, (m, s) in enumerate(zip(means, sds)):
        ax.text(i, m + s + 0.003, f"{m:.3f}", ha="center", fontsize=5.2)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(sup / f"FigureS1_model_benchmark.{ext}", dpi=DPI,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", sup / "FigureS1_model_benchmark.png")


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    figureS1()
    print("Done. Figure 4 (web screenshots) is generated after the UI redesign.")
