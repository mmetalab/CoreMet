"""
Main application module for CoreMet Web Server
"""

import logging

from dash import Dash, dcc, html, Input, Output, State, no_update, ctx
import dash_bootstrap_components as dbc
import os
from pathlib import Path

logger = logging.getLogger(__name__)

# Module-level storage for dynamically generated module layouts
_MODULE_LAYOUTS = {}

from app.config import config
from app.services.performance_service import PerformanceMonitor, monitor_performance
from components.navbar import make_navbar
from components.footer import make_footer


@monitor_performance("app_creation")
def create_app(config_name='default'):
    """Create and configure the Dash application"""

    app_config = config[config_name]
    performance_monitor = PerformanceMonitor()

    project_root = Path(__file__).parent.parent
    assets_path = project_root / 'assets'

    app = Dash(
        __name__,
        title=app_config.APP_NAME,
        external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME],
        assets_folder=str(assets_path),
    )

    # Add Open Graph meta tags in <head> for SEO / social sharing
    app.index_string = '''<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <meta property="og:title" content="CoreMet: Metabolite Interaction Database" />
        <meta property="og:description" content="Explore 1.95M+ curated metabolite interactions across seven layers (MPI, MEI, MDI, MMI, MDrI, MGI, mGWAS), proteins, enzymes, diseases, microbes, drugs, genes, and SNPs." />

        <meta property="og:type" content="website" />
        <meta property="og:site_name" content="CoreMet" />
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>'''

    app.config.suppress_callback_exceptions = True
    app.performance_monitor = performance_monitor

    # ── Eagerly import ALL page modules so their @callback decorators
    #    register with Dash BEFORE the server starts serving requests.
    #    Without this, page-level callbacks (filters, loaders, etc.) are
    #    never picked up by the Dash dependency graph.
    from pages import (home, predict, results, help,        # noqa: F401
                       enrichment, disease, database_v2,
                       blankPage, network, metabolite_detail,
                       batch_search, search, explore, downloads,
                       coremetai)

    # ── Register modular database pages ──────────────────────
    from pages.modules.module_factory import (
        create_browse_layout, create_predict_layout, register_module_callbacks)
    from pages.modules.mpi_config import MPI_CONFIG, load_mpi_db
    from pages.modules.mei_config import MEI_CONFIG, load_mei_db
    from pages.modules.mdi_config import MDI_CONFIG, load_mdi_db
    from pages.modules.mmi_config import MMI_CONFIG, load_mmi_db
    from pages.modules.mdri_config import MDRI_CONFIG, load_mdri_db
    from pages.modules.mgi_config import MGI_CONFIG, load_mgi_db
    from pages.modules.mgwas_config import MGWAS_CONFIG, load_mgwas_db

    _module_configs = [
        (MPI_CONFIG, load_mpi_db),
        (MEI_CONFIG, load_mei_db),
        (MDI_CONFIG, load_mdi_db),
        (MMI_CONFIG, load_mmi_db),
        (MDRI_CONFIG, load_mdri_db),
        (MGI_CONFIG, load_mgi_db),
        (MGWAS_CONFIG, load_mgwas_db),
    ]
    _module_layouts = {}
    for cfg, loader in _module_configs:
        key = cfg['type_key']
        _module_layouts[f'/{key}'] = create_browse_layout(cfg)
        _module_layouts[f'/{key}/browse'] = _module_layouts[f'/{key}']
        _module_layouts[f'/{key}/predict'] = create_predict_layout(cfg)
        register_module_callbacks(app, cfg, loader)
    # Store module layouts, update the existing dict (mutable, shared ref)
    _MODULE_LAYOUTS.update(_module_layouts)

    # ── CoreMet ID registry ─────────────────────────────────
    # Keep this lazy in production: loading the JSON registry as Python dicts adds
    # a large resident-memory block, and most routes do not need it at startup.
    if os.getenv("COREMET_EAGER_REGISTRY", "").lower() in {"1", "true", "yes"}:
        try:
            from app.services.entity_registry import build_registry
            build_registry()
        except Exception as e:
            logger.warning("Failed to build entity registry: %s", e)

    # Main layout, navbar and footer are OUTSIDE the page-swapped area
    app.layout = html.Div([
        dcc.Location(id='url', refresh=False),
        dcc.Location(id='global-search-redirect', refresh=True),
        dcc.Location(id='home-search-redirect', refresh=True),
        make_navbar(),
        dcc.Loading(
            html.Div(id='page-content'),
            type="circle",
            color="#3182ce",
            fullscreen=False,
        ),
        make_footer(),
        dcc.Interval(id='performance-interval', interval=30000, n_intervals=0),
    ])

    # Register API blueprint if available
    try:
        from api.routes import api_bp
        app.server.register_blueprint(api_bp, url_prefix='/api/v1')
    except Exception as e:
        logger.warning(f"API blueprint not registered: {e}")

    # SEO routes, sitemap.xml and robots.txt
    _register_seo_routes(app.server)

    register_callbacks(app)

    performance_monitor.log_performance(
        "app_creation_complete", 0,
        config_name=config_name, debug_mode=app_config.DEBUG,
    )

    return app


def _register_seo_routes(server):
    """Add /sitemap.xml and /robots.txt to the Flask server."""
    from flask import Response

    BASE_URL = "https://www.coremet.org"  # live custom domain

    STATIC_PAGES = [
        "/home", "/database", "/predict", "/disease", "/enrichment",
        "/network", "/profile", "/gallery", "/help",
        "/documentation", "/about", "/batch-search", "/api-docs",
    ]

    @server.route("/robots.txt")
    def robots_txt():
        body = (
            "User-agent: *\n"
            "Allow: /\n"
            f"Sitemap: {BASE_URL}/sitemap.xml\n"
        )
        return Response(body, mimetype="text/plain")

    @server.route("/sitemap.xml")
    def sitemap_xml():
        urls = []
        for page in STATIC_PAGES:
            urls.append(f"  <url><loc>{BASE_URL}{page}</loc><priority>0.8</priority></url>")

        # Metabolite detail pages, unique HMDB IDs from MPI
        try:
            from app.config import Config
            import pandas as pd
            cfg = Config()
            mpi = pd.read_csv(cfg.MPI_DB_PATH, usecols=["HMDB ID"])
            for hid in mpi["HMDB ID"].dropna().unique():
                urls.append(f"  <url><loc>{BASE_URL}/metabolite?id={hid}</loc></url>")
        except Exception:
            pass

        xml = (
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
            + "\n".join(urls) + "\n"
            '</urlset>'
        )
        return Response(xml, mimetype="application/xml")


def register_callbacks(app):
    """Register application callbacks, including shared navbar/cite callbacks"""

    # ── Page routing ──────────────────────────────────────────
    @app.callback(
        Output('page-content', 'children'),
        Input('url', 'pathname'),
        Input('url', 'search'),
    )
    def display_page(pathname, search):
        from pages import (home, predict, results, help,
                           enrichment, disease, database_v2,
                           blankPage, metabolite_profile, network,
                           documentation, about, gallery, metabolite_detail,
                           api_docs, batch_search, search as search_page,
                           explore, downloads, entity_detail, coremetai)

        if not pathname:
            pathname = '/home'

        # Dynamic route: metabolite detail page needs query params
        if pathname == '/metabolite':
            return metabolite_detail.build_page(search)

        # Dynamic routes: non-metabolite entity detail pages
        _ENTITY_ROUTES = {
            '/disease-detail': 'disease',
            '/gene': 'gene',
            '/protein': 'protein',
            '/drug': 'drug',
            '/microbe': 'microbe',
            '/snp': 'snp',
        }
        if pathname in _ENTITY_ROUTES:
            return entity_detail.build_page(_ENTITY_ROUTES[pathname], search)

        # Dynamic routes: browse entities by type (/explore/<type>)
        if pathname and pathname.startswith('/explore/'):
            entity_key = pathname.split('/explore/')[-1].strip('/')
            if entity_key in explore.ENTITY_BROWSE_CONFIG:
                return explore.build_entity_browse_page(entity_key)

        # Dynamic route: search page with optional query
        if pathname == '/search':
            return search_page.build_search_page(search)

        page_map = {
            '/': home.page_content,
            '/home': home.page_content,
            '/search': search_page.page_content,
            '/explore': explore.page_content,
            '/downloads': downloads.page_content,
            '/database': database_v2.page_content,
            '/interactions': predict.page_content,
            '/predict': predict.page_content,
            '/help': help.page_content,
            '/enrichment': enrichment.page_content,
            '/disease': disease.page_content,
            '/profile': metabolite_profile.page_content,
            '/network': network.page_content,
            '/documentation': documentation.page_content,
            '/about': about.page_content,
            '/gallery': gallery.page_content,
            '/api-docs': api_docs.page_content,
            '/batch-search': batch_search.page_content,
            '/coremetai': coremetai.page_content,
        }
        if pathname in page_map:
            return page_map[pathname]
        # Module routes (/{type}, /{type}/browse, /{type}/predict)
        from app.main import _MODULE_LAYOUTS
        if pathname in _MODULE_LAYOUTS:
            return _MODULE_LAYOUTS[pathname]
        if pathname and pathname.startswith('/results'):
            return results.page_content
        return blankPage.layout

    # ── Shared navbar toggle (ONCE for the whole app) ─────────
    @app.callback(
        Output("navbar-collapse", "is_open"),
        Input("navbar-toggler", "n_clicks"),
        State("navbar-collapse", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_navbar(n_clicks, is_open):
        if n_clicks:
            return not is_open
        return is_open

    # ── Home search callbacks are defined in pages/home.py via @callback ──

    # ── Re-trigger Cytoscape layout when Network tab is activated ──
    # Cytoscape inside a hidden dbc.Tab can't compute layout (zero-size
    # container).  When the user switches to the Network tab we must
    # re-run the layout so nodes spread correctly.
    app.clientside_callback(
        """
        function(activeTab) {
            if (activeTab === 'tab-network') {
                // Small delay to let the tab become visible first
                setTimeout(function() {
                    window.dispatchEvent(new Event('resize'));
                }, 150);
                return {
                    "name": "cose", "animate": true, "animationDuration": 400,
                    "nodeRepulsion": 10000, "idealEdgeLength": 100,
                    "nodeOverlap": 25, "gravity": 0.25, "nestingFactor": 1.2,
                    "padding": 40
                };
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output("entity-cytoscape", "layout", allow_duplicate=True),
        Input("entity-detail-tabs", "active_tab"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        """
        function(activeTab) {
            if (activeTab === 'tab-network') {
                setTimeout(function() {
                    window.dispatchEvent(new Event('resize'));
                }, 150);
                return {
                    "name": "cose", "animate": true, "animationDuration": 400,
                    "nodeRepulsion": 8000, "idealEdgeLength": 90,
                    "nodeOverlap": 20, "gravity": 0.3, "padding": 40
                };
            }
            return window.dash_clientside.no_update;
        }
        """,
        Output("metdetail-cytoscape", "layout", allow_duplicate=True),
        Input("entity-detail-tabs", "active_tab"),
        prevent_initial_call=True,
    )

    # ── Cytoscape network controls ───────────────────────────
    # Zoom in / out via the Cytoscape zoom prop
    @app.callback(
        Output("metdetail-cytoscape", "zoom"),
        Input("cyto-zoom-in", "n_clicks"),
        Input("cyto-zoom-out", "n_clicks"),
        State("metdetail-cytoscape", "zoom"),
        prevent_initial_call=True,
    )
    def cyto_zoom(zoom_in, zoom_out, cur_zoom):
        if cur_zoom is None:
            cur_zoom = 1.0
        trigger = ctx.triggered_id
        if trigger == "cyto-zoom-in":
            return min(cur_zoom * 1.3, 5.0)
        elif trigger == "cyto-zoom-out":
            return max(cur_zoom / 1.3, 0.2)
        return no_update

    # Fit-to-view: reset zoom + pan to centre
    @app.callback(
        Output("metdetail-cytoscape", "layout", allow_duplicate=True),
        Input("cyto-fit", "n_clicks"),
        Input("cyto-reset", "n_clicks"),
        prevent_initial_call=True,
    )
    def cyto_fit_reset(fit_clicks, reset_clicks):
        trigger = ctx.triggered_id
        if trigger == "cyto-fit":
            # preset layout keeps current positions
            return {"name": "preset", "fit": True, "animate": True}
        elif trigger == "cyto-reset":
            # re-run CoSE layout from scratch
            return {"name": "cose", "animate": True, "nodeRepulsion": 6000,
                    "idealEdgeLength": 80, "nodeOverlap": 20}
        return no_update

    # Expand / collapse height
    @app.callback(
        Output("metdetail-cytoscape", "style"),
        Input("cyto-expand", "n_clicks"),
        State("metdetail-cytoscape", "style"),
        prevent_initial_call=True,
    )
    def cyto_toggle_height(n_clicks, cur_style):
        if cur_style is None:
            cur_style = {}
        h = cur_style.get("height", "380px")
        new_h = "680px" if h == "380px" else "380px"
        return {**cur_style, "height": new_h}

    # ── Cytoscape node hover tooltip (content + visibility) ──
    # Positioning is handled by assets/cyto_tooltip.js (mouse-following)
    @app.callback(
        Output("cyto-tooltip", "children"),
        Output("cyto-tooltip", "style"),
        Input("metdetail-cytoscape", "mouseoverNodeData"),
        State("cyto-tooltip", "style"),
        prevent_initial_call=True,
    )
    def cyto_node_tooltip(hover_data, cur_style):
        if cur_style is None:
            cur_style = {}
        if not hover_data:
            return "", {**cur_style, "display": "none"}
        full = hover_data.get("full_label", hover_data.get("label", ""))
        return full, {**cur_style, "display": "block"}
