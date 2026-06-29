"""MEI module configuration, Metabolite–Enzyme Interactions."""


def _mei_transform(df):
    """All MEI entries are curated KEGG/Rhea catalytic relationships."""
    if 'evidence_type' not in df.columns:
        df['evidence_type'] = 'Curated'
    if 'confidence' not in df.columns:
        df['confidence'] = 0.9
    return df


MEI_CONFIG = {
    'type_key': 'mei',
    'type_name': 'Metabolite–Enzyme Interactions (MEI)',
    'short_name': 'MEI',
    'icon': 'fas fa-vials',
    'color': '#00a3c4',
    'description': '47,551 curated metabolite–enzyme (EC) catalytic relationships from KEGG and Rhea, spanning 2,406 enzyme classes across 37 organisms',
    'data_transform': _mei_transform,
    'columns': {
        'display': ['Metabolite_Name', 'HMDB_ID', 'Enzyme_Name', 'EC_Number',
                    'Gene_Name', 'Species', 'Pathway_Name', 'Evidence_Source'],
        'search': ['Metabolite_Name', 'HMDB_ID', 'Enzyme_Name', 'EC_Number', 'Gene_Name'],
    },
    'filters': {
        'organism': True,
        'organism_column': 'Species',
        'extra': [
            {'label': 'Pathway', 'column': 'Pathway_Name', 'type': 'dropdown'},
            {'label': 'Evidence Source', 'column': 'Evidence_Source', 'type': 'dropdown'},
        ],
    },
    'viz': [
        {'type': 'donut', 'column': 'Species', 'title': 'Organism Distribution'},
        {'type': 'donut', 'column': 'Evidence_Source', 'title': 'Evidence Source'},
        {'type': 'bar_top', 'column': 'Metabolite_Name', 'title': 'Top Metabolites by Enzyme Count', 'n': 15},
        {'type': 'coverage_bar', 'column': 'Pathway_Name', 'title': 'KEGG Pathway Annotation Coverage'},
    ],
    'stat_cards': [
        {'label': 'Total Interactions', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB_ID', 'icon': 'fas fa-atom'},
        {'label': 'EC Numbers', 'column': 'EC_Number', 'icon': 'fas fa-vials'},
        {'label': 'Organisms', 'column': 'Species', 'icon': 'fas fa-globe'},
    ],
    'external_links': {
        'HMDB_ID': 'https://hmdb.ca/metabolites/{value}',
        'Uniprot_ID': 'https://www.uniprot.org/uniprot/{value}',
        'EC_Number': 'https://enzyme.expasy.org/EC/{value}',
    },
    'source_links': {
        'Evidence_Source': {
            'KEGG': 'https://www.genome.jp/kegg/',
            'Rhea': 'https://www.rhea-db.org/',
        },
    },
    'source_versions': {
        'KEGG': 'Release 2025',
        'Rhea': 'v139 (2024)',
    },
    'entity_links': {
        'Metabolite_Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB_ID'},
        'Gene_Name': {'route': '/gene', 'param': 'name'},
    },
    'references': [
        {'title': 'KEGG: Kyoto Encyclopedia of Genes and Genomes', 'authors': 'Kanehisa M, Goto S', 'journal': 'Nucleic Acids Res', 'year': 2000, 'pmid': '10592173'},
        {'title': 'Rhea, the reaction knowledgebase in 2022', 'authors': 'Bansal P et al.', 'journal': 'Nucleic Acids Res', 'year': 2022, 'pmid': '34755880'},
    ],
}


def load_mei_db():
    """Load the deduplicated MEI release table (counts match the manuscript)."""
    from pathlib import Path
    import pandas as pd
    base = Path(__file__).parent.parent.parent / "data" / "databases"
    release = base / "release" / "coremetdb_mei.csv"
    path = release if release.exists() else base / "mei_database_v2_enriched.csv"
    from app.services.csv_loader import load_optimized
    return load_optimized(path)
