"""MPI module configuration, Metabolite–Protein Interactions."""


def _mpi_transform(df):
    """Add evidence_type: all MPI entries are from curated databases (KEGG/Rhea) = Experimental."""
    import pandas as pd
    df['evidence_type'] = 'Experimental'
    if 'confidence' not in df.columns:
        df['confidence'] = 0.85
    return df


MPI_CONFIG = {
    'type_key': 'mpi',
    'type_name': 'Metabolite–Protein Interactions (MPI)',
    'short_name': 'MPI',
    'icon': 'fas fa-dna',
    'color': '#0072B2',
    'description': '38,061 curated metabolite–protein interactions from KEGG (Release 2025), Rhea (v139, 2024), and MPIDB across 10 organisms',
    'data_transform': _mpi_transform,
    'columns': {
        'display': ['Metabolite Name', 'HMDB ID', 'Protein Name', 'Uniprot ID',
                     'Gene Name', 'Species', 'Pathway_Name', 'evidence_type',
                     'Evidence_Source', 'confidence'],
        'search': ['Metabolite Name', 'HMDB ID', 'Protein Name', 'Gene Name', 'Uniprot ID'],
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
        {'type': 'donut', 'column': 'evidence_type', 'title': 'Evidence Type'},
        {'type': 'bar_top', 'column': 'Metabolite Name', 'title': 'Top Metabolites by Interaction Count', 'n': 15},
        {'type': 'coverage_bar', 'column': 'Pathway_Name', 'title': 'KEGG Pathway Annotation Coverage'},
    ],
    'stat_cards': [
        {'label': 'Total Interactions', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB ID', 'icon': 'fas fa-atom'},
        {'label': 'Proteins', 'column': 'Uniprot ID', 'icon': 'fas fa-microscope'},
        {'label': 'Organisms', 'column': 'Species', 'icon': 'fas fa-globe'},
    ],
    'external_links': {
        'HMDB ID': 'https://hmdb.ca/metabolites/{value}',
        'Uniprot ID': 'https://www.uniprot.org/uniprot/{value}',
    },
    'source_links': {
        'Evidence_Source': {
            'KEGG': 'https://www.genome.jp/kegg/reaction/',
            'Rhea': 'https://www.rhea-db.org/',
        },
    },
    'source_versions': {
        'KEGG': 'Release 2025',
        'Rhea': 'v139 (2024)',
    },
    'entity_links': {
        'Metabolite Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB ID'},
        'Protein Name': {'route': '/protein', 'param': 'name'},
        'Gene Name': {'route': '/gene', 'param': 'name'},
    },
    'references': [
        {'title': 'KEGG: Kyoto Encyclopedia of Genes and Genomes', 'authors': 'Kanehisa M, Goto S', 'journal': 'Nucleic Acids Res', 'year': 2000, 'pmid': '10592173'},
        {'title': 'Rhea, the reaction knowledgebase in 2022', 'authors': 'Bansal P et al.', 'journal': 'Nucleic Acids Res', 'year': 2022, 'pmid': '34755880'},
        {'title': 'HMDB 5.0: the Human Metabolome Database for 2022', 'authors': 'Wishart DS et al.', 'journal': 'Nucleic Acids Res', 'year': 2022, 'pmid': '34986597'},
    ],
}


def load_mpi_db():
    """Load the MPI v2 database with standardized evidence."""
    from pathlib import Path
    import pandas as pd
    rel = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mpi.csv"
    v2_path = rel if rel.exists() else Path(__file__).parent.parent.parent / "data" / "databases_v2" / "mpi_database_v2.csv"
    if v2_path.exists():
        return pd.read_csv(v2_path, low_memory=False)
    from app.services.data_service import get_db
    return get_db()
