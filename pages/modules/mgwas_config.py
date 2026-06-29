"""mGWAS module configuration, Metabolite–SNP Associations."""


def _mgwas_transform(df):
    """All mGWAS are genome-wide significant (p < 5e-8) experimental associations."""
    df['evidence_type'] = 'Experimental'
    if 'confidence' not in df.columns:
        import pandas as pd
        # Higher confidence for lower p-values
        pv = pd.to_numeric(df.get('P_Value', 1), errors='coerce').fillna(1)
        df['confidence'] = pv.apply(lambda p: 0.95 if p < 1e-20 else (0.85 if p < 1e-10 else 0.75))
    return df


MGWAS_CONFIG = {
    'type_key': 'mgwas',
    'type_name': 'Metabolite–SNP Associations (mGWAS)',
    'short_name': 'mGWAS',
    'icon': 'fas fa-map-marker-alt',
    'color': '#319795',
    'description': (
        '44,344 genome-wide significant metabolite–SNP associations curated from '
        'NHGRI-EBI GWAS Catalog (2025; 43,628), Yin et al. 2024 (621), and Shin et al. 2014 (95). '
        'All associations reach p < 5×10⁻⁸ with PubMed references.'
    ),
    'data_transform': _mgwas_transform,
    'columns': {
        'display': ['Metabolite_Name', 'HMDB_ID', 'rsID', 'Chromosome',
                     'Mapped_Gene', 'P_Value', 'Beta', 'Trait',
                     'evidence_type', 'Source', 'PMID'],
        'search': ['Metabolite_Name', 'HMDB_ID', 'rsID', 'Mapped_Gene', 'Source', 'Trait'],
    },
    'filters': {
        'extra': [
            {'label': 'Chromosome', 'column': 'Chromosome', 'type': 'dropdown'},
            {'label': 'Source', 'column': 'Source', 'type': 'dropdown'},
            {'label': 'Trait', 'column': 'Trait', 'type': 'dropdown'},
            {'label': 'Trait Category', 'column': 'trait_category', 'type': 'dropdown'},
            {'label': 'Significance', 'column': 'significance_tier', 'type': 'dropdown'},
        ],
    },
    'viz': [
        {'type': 'donut', 'column': 'Source', 'title': 'Data Source Distribution'},
        {'type': 'donut', 'column': 'evidence_type', 'title': 'Evidence Type'},
        {'type': 'donut', 'column': 'Chromosome', 'title': 'Chromosome Distribution'},
        {'type': 'bar_top', 'column': 'Mapped_Gene', 'title': 'Top Mapped Genes', 'n': 15},
        {'type': 'log_histogram', 'column': 'P_Value', 'title': 'P-Value Distribution (-log10)'},
    ],
    'stat_cards': [
        {'label': 'Total Associations', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB_ID', 'icon': 'fas fa-atom'},
        {'label': 'SNPs', 'column': 'rsID', 'icon': 'fas fa-map-marker-alt'},
        {'label': 'Mapped Genes', 'column': 'Mapped_Gene', 'icon': 'fas fa-dna'},
        {'label': 'With PMID', 'column': 'PMID', 'icon': 'fas fa-book'},
    ],
    'external_links': {
        'rsID': 'https://www.ncbi.nlm.nih.gov/snp/{value}',
        'HMDB_ID': 'https://hmdb.ca/metabolites/{value}',
    },
    'source_links': {
        'Source': {
            'GWAS_Catalog': 'https://www.ebi.ac.uk/gwas/',
            'Shin2014': 'https://pubmed.ncbi.nlm.nih.gov/24816252/',
            'Yin2024': 'https://pubmed.ncbi.nlm.nih.gov/38120091/',
            'GWAS_Catalog; Yin2024': 'https://www.ebi.ac.uk/gwas/',
            'GWAS_Catalog; Shin2014': 'https://www.ebi.ac.uk/gwas/',
        },
    },
    'source_versions': {
        'GWAS_Catalog': 'Release 2025-03-17, GRCh38.p14, dbSNP 156',
        'Shin2014': 'Nat Genet 46:543–550 (PMID: 24816252)',
        'Yin2024': 'Nat Genet (PMID: 38120091)',
    },
    'pmid_column': 'PMID',
    'entity_links': {
        'Metabolite_Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB_ID'},
        'rsID': {'route': '/snp', 'param': 'name'},
        'Mapped_Gene': {'route': '/gene', 'param': 'name'},
    },
    'references': [
        {'title': 'The NHGRI-EBI GWAS Catalog of published genome-wide association studies', 'authors': 'Sollis E et al.', 'journal': 'Nucleic Acids Res', 'year': 2023, 'pmid': '36350656'},
        {'title': 'An atlas of genetic influences on human blood metabolites', 'authors': 'Shin SY et al.', 'journal': 'Nat Genet', 'year': 2014, 'pmid': '24816252'},
        {'title': 'Genome-wide association study of metabolite levels', 'authors': 'Yin X et al.', 'journal': 'Nat Genet', 'year': 2024, 'pmid': ''},
    ],
}


def load_mgwas_db():
    from pathlib import Path
    import pandas as pd
    rel = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mgwas.csv"
    v2 = rel if rel.exists() else Path(__file__).parent.parent.parent / "data" / "databases_v2" / "mgwas_database_v2.csv"
    if v2.exists():
        return pd.read_csv(v2, low_memory=False)
    from app.services.mgwas_service import get_mgwas_db
    return get_mgwas_db()
