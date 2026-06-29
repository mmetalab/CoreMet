#!/usr/bin/env python3
"""
capture_figure4.py — Auto-capture the CoreMet web-interface figure (Figure 4).

Requires a running app (default http://127.0.0.1:8095) and playwright+chromium.
Captures three panels and composites them into Submission_NAR/Figures/Figure4_web_interface.{png,pdf}.

    conda run -n mpi-vgae python scripts/capture_figure4.py [BASE_URL]
"""
from __future__ import annotations
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://127.0.0.1:8095"
OUT = Path(__file__).resolve().parent.parent.parent / "Submission_NAR" / "Figures"
SHOTS = Path(__file__).resolve().parent.parent / "data" / "figures" / "screenshots"
SHOTS.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

VW, VH = 1440, 960
PANELS = [
    ("A", "/database", "Database browser", SHOTS / "panelA_database.png"),
    ("B", "/metabolite?id=HMDB0000039", "Cross-layer metabolite profile", SHOTS / "panelB_profile.png"),
    ("C", "/network", "Network exploration", SHOTS / "panelC_network.png"),
]


def capture():
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": VW, "height": VH}, device_scale_factor=2)
        for label, route, _title, path in PANELS:
            url = BASE + route
            print("capturing", url)
            page.goto(url, wait_until="networkidle", timeout=60000)
            # give Dash callbacks time to populate tables / charts / cytoscape
            try:
                page.wait_for_selector(".cm-page-container, .dash-table-container, #cytoscape-network",
                                       timeout=15000)
            except Exception:
                pass
            time.sleep(5)
            # Network panel: run an example query so the graph is populated.
            if route == "/network":
                try:
                    page.fill("#net-query-input", "Glucose")
                    page.click("#net-search-btn")
                    time.sleep(9)  # let the 6-DB search + cytoscape layout render
                except Exception as e:
                    print("  network query failed:", e)
            page.screenshot(path=str(path))  # viewport-only for uniform panels
        browser.close()


def composite():
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 12.6))
    for ax, (label, _route, title, path) in zip(axes, PANELS):
        ax.axis("off")
        if path.exists():
            ax.imshow(mpimg.imread(str(path)))
        ax.text(-0.01, 1.01, label, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="bottom", ha="left")
        ax.text(0.5, 1.005, title, transform=ax.transAxes, fontsize=9,
                va="bottom", ha="center", color="#4a5568")
    fig.subplots_adjust(hspace=0.06, left=0.01, right=0.99, top=0.98, bottom=0.01)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"Figure4_web_interface.{ext}", dpi=200,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", OUT / "Figure4_web_interface.png")


if __name__ == "__main__":
    capture()
    composite()
