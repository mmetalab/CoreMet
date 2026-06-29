#!/usr/bin/env python3
"""
generate_graphical_abstract.py — Mandatory NAR graphical abstract for CoreMet.

5:2 landscape, 300 dpi, Arial, color, text sparing, reads left->right:
  fragmented sources  ->  CoreMet metabolite hub  ->  cross-domain query.
Numbers read from data/coremetdb_stats.json. Output: Submission_NAR/Figures/Graphical_Abstract.{pdf,png,tif}
"""
from __future__ import annotations
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle

ROOT = Path(__file__).resolve().parent.parent
S = json.loads((ROOT / "data" / "coremetdb_stats.json").read_text())
OUT = ROOT.parent / "Submission_NAR" / "Figures"
OUT.mkdir(parents=True, exist_ok=True)

matplotlib.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
EC = {"metabolite": "#e27a3f", "protein": "#3182ce", "enzyme": "#00a3c4", "disease": "#e53e3e",
      "microbe": "#38a169", "drug": "#805ad5", "gene": "#d69e2e", "SNP": "#319795"}
INK = "#1a202c"
T = S["totals"]


def box(ax, x, y, w, h, text, fc, fs=11, tc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=0.4",
                                fc=fc, ec=fc, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", color=tc, fontsize=fs, zorder=3)


def arrow(ax, p1, p2, color="#a0aec0", lw=2.2):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=20,
                                 color=color, lw=lw, zorder=1, shrinkA=4, shrinkB=4))


fig = plt.figure(figsize=(10, 4))           # 5:2 landscape
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")

# ── Left: fragmented sources ──
ax.text(15, 37, "Fragmented across\nspecialized databases", ha="center", fontsize=12,
        fontweight="bold", color=INK)
srcs = ["KEGG / Rhea", "CTD", "HMDB", "gutMGene", "DrugBank", "GWAS Catalog"]
pos = [(3, 26), (17, 28), (5, 19), (18, 20), (4, 11), (16, 12)]
for s, (x, y) in zip(srcs, pos):
    box(ax, x, y, 12, 4.2, s, "#94a3b8", fs=9)

# ── Center: CoreMet metabolite hub ──
cx, cy = 50, 20
ax.text(cx, 37, "CoreMet", ha="center", fontsize=17, fontweight="bold", color=INK)
ax.add_patch(Circle((cx, cy), 3.4, fc=EC["metabolite"], ec="white", lw=1.5, zorder=5))
ax.text(cx, cy, "MET", ha="center", va="center", color="white", fontsize=10,
        fontweight="bold", zorder=6)
spokes = [("protein", "Protein"), ("enzyme", "Enzyme"), ("disease", "Disease"), ("microbe", "Microbe"),
          ("drug", "Drug"), ("gene", "Gene"), ("SNP", "SNP")]
for i, (k, lab) in enumerate(spokes):
    a = math.radians(90 + i * (360 / 7))
    nx, ny = cx + 11 * math.cos(a), cy + 9.5 * math.sin(a)
    ax.plot([cx, nx], [cy, ny], color=EC[k], lw=2.2, zorder=2)
    ax.add_patch(Circle((nx, ny), 2.7, fc=EC[k], ec="white", lw=1.0, zorder=4))
    ax.text(nx, ny, lab, ha="center", va="center", color="white", fontsize=7.5, zorder=5)
ax.text(cx, 4.2, f"{T['interactions']:,} interactions · 7 layers · {T['unique_metabolites_union']:,} metabolites",
        ha="center", fontsize=10.5, color=INK, fontweight="bold")

# ── Right: unified query ──
ax.text(86, 37, "One evidence-aware\ncross-domain query", ha="center", fontsize=12,
        fontweight="bold", color=INK)
outs = [("Cross-layer metabolite profile", "#1a365d"),
        ("Mechanism tracing across scales", "#2c7a7b"),
        ("ML-ready graph · REST API", "#4a5568")]
for i, (t, c) in enumerate(outs):
    box(ax, 74, 25 - i * 6.5, 24, 5, t, c, fs=9.5)

arrow(ax, (31, 20), (45, 20))
arrow(ax, (61, 20), (73, 20))

for ext in ("pdf", "png", "tif"):
    fig.savefig(OUT / f"Graphical_Abstract.{ext}", dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)
print("wrote", OUT / "Graphical_Abstract.png")
