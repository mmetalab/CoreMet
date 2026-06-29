"""MDI module configuration, Metabolite–Disease Interactions."""


def _mdi_transform(df):
    """Classify evidence: CTD/HMDB = Experimental, case_study/inferred = Predicted."""
    import pandas as pd
    if 'Source' in df.columns:
        df['evidence_type'] = df['Source'].apply(lambda x: (
            'Predicted' if any(k in str(x).lower() for k in ['case_study', 'inferred', 'predicted'])
            else 'Experimental'
        ))
    else:
        df['evidence_type'] = 'Experimental'
    return df


MDI_CONFIG = {
    'type_key': 'mdi',
    'type_name': 'Metabolite–Disease Interactions (MDI)',
    'short_name': 'MDI',
    'icon': 'fas fa-disease',
    'color': '#D55E00',
    'description': '82,882 metabolite–disease associations from CTD (2025) and HMDB (v5.0, 2024) across 130 diseases and 22 categories',
    'data_transform': _mdi_transform,
    'columns': {
        'display': ['HMDB_ID', 'Metabolite_Name', 'Disease_Name', 'Category',
                     'Association_Type', 'evidence_type', 'confidence',
                     'Source', 'pmid'],
        'search': ['Metabolite_Name', 'HMDB_ID', 'Disease_Name', 'Disease_ID'],
    },
    'filters': {
        'category': True,
        'category_column': 'Category',
        'extra': [
            {'label': 'Disease', 'column': 'Disease_Name', 'type': 'dropdown'},
            {'label': 'Association Type', 'column': 'Association_Type', 'type': 'dropdown'},
            {'label': 'Evidence Level', 'column': 'Evidence_Level', 'type': 'dropdown'},
        ],
    },
    'viz': [
        {'type': 'donut', 'column': 'Category', 'title': 'Disease Categories'},
        {'type': 'donut', 'column': 'evidence_type', 'title': 'Evidence Type'},
        {'type': 'bar_top', 'column': 'Disease_Name', 'title': 'Top Diseases by Metabolite Count', 'n': 15},
        {'type': 'heatmap', 'row_col': 'Disease_Name', 'col_col': 'Metabolite_Name',
         'title': 'Disease–Metabolite Association Matrix (Top 10×10)', 'n': 10},
    ],
    'stat_cards': [
        {'label': 'Total Associations', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB_ID', 'icon': 'fas fa-atom'},
        {'label': 'Diseases', 'column': 'Disease_Name', 'icon': 'fas fa-heartbeat'},
        {'label': 'Categories', 'column': 'Category', 'icon': 'fas fa-tags'},
    ],
    'external_links': {
        'HMDB_ID': 'https://hmdb.ca/metabolites/{value}',
    },
    'source_links': {
        'Source': {
            'CTD': 'https://ctdbase.org/',
            'HMDB': 'https://hmdb.ca/',
            'CoreMet_curated': '#',
            'CoreMet_disease_mpi': '#',
            'CoreMet_case_study': '#',
        },
    },
    'source_versions': {
        'CTD': '2025 release',
        'HMDB': 'v5.0 (2024)',
    },
    'pmid_column': 'pmid',
    'entity_links': {
        'Metabolite_Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB_ID'},
        'Disease_Name': {'route': '/disease-detail', 'param': 'name'},
    },
    'references': [
        {'title': 'CTD: The Comparative Toxicogenomics Database: update 2023', 'authors': 'Davis AP et al.', 'journal': 'Nucleic Acids Res', 'year': 2023, 'pmid': '36169237'},
        {'title': 'HMDB 5.0: the Human Metabolome Database for 2022', 'authors': 'Wishart DS et al.', 'journal': 'Nucleic Acids Res', 'year': 2022, 'pmid': '34986597'},
        {'title': 'DisGeNET: a comprehensive platform for disease-associated genes and variants', 'authors': 'Pinero J et al.', 'journal': 'Nucleic Acids Res', 'year': 2020, 'pmid': '31680165'},
    ],
}


def load_mdi_db():
    from pathlib import Path
    import pandas as pd
    rel = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mdi.csv"
    if rel.exists():
        from app.services.csv_loader import load_optimized
        return load_optimized(rel)
    from app.services.mdi_service import get_mdi_db
    return get_mdi_db()
