"""MGI module configuration, Metabolite–Gene Interactions."""


def _mgi_transform(df):
    """CTD interactions backed by PubMed = Experimental."""
    df['evidence_type'] = 'Experimental'
    if 'confidence' not in df.columns:
        df['confidence'] = 0.85
    return df


MGI_CONFIG = {
    'type_key': 'mgi',
    'type_name': 'Metabolite–Gene Interactions (MGI)',
    'short_name': 'MGI',
    'icon': 'fas fa-dna',
    'color': '#CC8400',
    'description': '1,658,745 metabolite–gene interactions from CTD (2025) covering 50,164 genes across 7 organisms',
    'data_transform': _mgi_transform,
    'columns': {
        'display': ['Metabolite_Name', 'HMDB_ID', 'Gene_Symbol', 'Gene_ID',
                     'Organism', 'Interaction_Type', 'Interaction_Actions',
                     'evidence_type', 'Source', 'PMID'],
        'search': ['Metabolite_Name', 'HMDB_ID', 'Gene_Symbol', 'Gene_ID'],
    },
    'filters': {
        'organism': True,
        'organism_column': 'Organism',
        'extra': [
            {'label': 'Interaction Type', 'column': 'Interaction_Type', 'type': 'dropdown'},
            {'label': 'Interaction Actions', 'column': 'Interaction_Actions', 'type': 'dropdown'},
        ],
    },
    'viz': [
        {'type': 'donut', 'column': 'evidence_type', 'title': 'Evidence Type'},
        {'type': 'donut', 'column': 'Organism', 'title': 'Organism Distribution'},
        {'type': 'bar_top', 'column': 'Interaction_Type', 'title': 'Top Interaction Types', 'n': 15},
        {'type': 'bar_top', 'column': 'Gene_Symbol', 'title': 'Top Genes by Metabolite Count', 'n': 15},
    ],
    'external_links': {
        'Gene_ID': 'https://www.ncbi.nlm.nih.gov/gene/{value}',
        'HMDB_ID': 'https://hmdb.ca/metabolites/{value}',
    },
    'source_links': {
        'Source': {
            'CTD': 'https://ctdbase.org/',
        },
    },
    'source_versions': {
        'CTD': '2025 release',
    },
    'pmid_column': 'PMID',
    'entity_links': {
        'Metabolite_Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB_ID'},
        'Gene_Symbol': {'route': '/gene', 'param': 'name'},
    },
    'references': [
        {'title': 'CTD: The Comparative Toxicogenomics Database: update 2023', 'authors': 'Davis AP et al.', 'journal': 'Nucleic Acids Res', 'year': 2023, 'pmid': '36169237'},
        {'title': 'HMDB 5.0: the Human Metabolome Database for 2022', 'authors': 'Wishart DS et al.', 'journal': 'Nucleic Acids Res', 'year': 2022, 'pmid': '34986597'},
    ],
    'stat_cards': [
        {'label': 'Total Interactions', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB_ID', 'icon': 'fas fa-atom'},
        {'label': 'Genes', 'column': 'Gene_Symbol', 'icon': 'fas fa-dna'},
        {'label': 'Organisms', 'column': 'Organism', 'icon': 'fas fa-globe'},
    ],
}


def load_mgi_db():
    from pathlib import Path
    import pandas as pd
    rel = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mgi.csv"
    v2 = rel if rel.exists() else Path(__file__).parent.parent.parent / "data" / "databases_v2" / "mgi_database_v2.csv"
    if v2.exists():
        from app.services.csv_loader import load_optimized
        return load_optimized(v2)
    from app.services.mgi_service import get_mgi_db
    return get_mgi_db()
