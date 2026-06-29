"""MDrI module configuration, Metabolite–Drug Interactions."""


def _mdri_transform(df):
    """All DrugBank interactions are expert-curated = Experimental."""
    df['evidence_type'] = 'Experimental'
    if 'confidence' not in df.columns:
        df['confidence'] = 0.85
    return df


MDRI_CONFIG = {
    'type_key': 'mdri',
    'type_name': 'Metabolite–Drug Interactions (MDrI)',
    'short_name': 'MDrI',
    'icon': 'fas fa-pills',
    'color': '#7B2D8E',
    'description': '3,500 metabolite–drug interactions from DrugBank (v5.1, 2024) covering 2,162 drugs',
    'data_transform': _mdri_transform,
    'columns': {
        'display': ['Metabolite_Name', 'HMDB_ID', 'Drug_Name', 'DrugBank_ID',
                     'Interaction_Type', 'Tissue', 'Evidence_Level',
                     'evidence_type', 'Source', 'PMID'],
        'search': ['Metabolite_Name', 'HMDB_ID', 'Drug_Name', 'DrugBank_ID'],
    },
    'filters': {
        'extra': [
            {'label': 'Interaction Type', 'column': 'Interaction_Type', 'type': 'dropdown'},
            {'label': 'Tissue', 'column': 'Tissue', 'type': 'dropdown'},
            {'label': 'Evidence Level', 'column': 'Evidence_Level', 'type': 'dropdown'},
        ],
    },
    'viz': [
        {'type': 'donut', 'column': 'evidence_type', 'title': 'Evidence Type'},
        {'type': 'donut', 'column': 'Interaction_Type', 'title': 'Interaction Types'},
        {'type': 'bar_top', 'column': 'Drug_Name', 'title': 'Top Drugs', 'n': 15},
        {'type': 'donut', 'column': 'Source', 'title': 'Source Breakdown'},
    ],
    'external_links': {
        'DrugBank_ID': 'https://go.drugbank.com/drugs/{value}',
        'HMDB_ID': 'https://hmdb.ca/metabolites/{value}',
    },
    'source_links': {
        'Source': {
            'DrugBank_cross_ref': 'https://go.drugbank.com/releases/latest',
            'DrugBank_enzyme_bridge': 'https://go.drugbank.com/releases/latest',
            'CoreMet_curated': '#',
        },
    },
    'source_versions': {
        'DrugBank': 'v5.1 (2024)',
    },
    'pmid_column': 'PMID',
    'entity_links': {
        'Metabolite_Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB_ID'},
        'Drug_Name': {'route': '/drug', 'param': 'name'},
    },
    'references': [
        {'title': 'DrugBank 5.0: a major update to the DrugBank database for 2018', 'authors': 'Wishart DS et al.', 'journal': 'Nucleic Acids Res', 'year': 2018, 'pmid': '29126136'},
        {'title': 'PharmGKB: The Pharmacogenomics Knowledge Base', 'authors': 'Whirl-Carrillo M et al.', 'journal': 'Clin Pharmacol Ther', 'year': 2012, 'pmid': '22992668'},
    ],
    'stat_cards': [
        {'label': 'Total Interactions', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB_ID', 'icon': 'fas fa-atom'},
        {'label': 'Drugs', 'column': 'Drug_Name', 'icon': 'fas fa-pills'},
        {'label': 'DrugBank IDs', 'column': 'DrugBank_ID', 'icon': 'fas fa-barcode'},
    ],
}


def load_mdri_db():
    from pathlib import Path
    import pandas as pd
    rel = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mdri.csv"
    v2 = rel if rel.exists() else Path(__file__).parent.parent.parent / "data" / "databases_v2" / "mdri_database_v2.csv"
    if v2.exists():
        return pd.read_csv(v2, low_memory=False)
    from app.services.mdri_service import get_mdri_db
    return get_mdri_db()
