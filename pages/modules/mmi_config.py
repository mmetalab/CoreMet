"""MMI module configuration, Metabolite–Microbe Interactions."""


def _mmi_transform(df):
    """Classify evidence: gutMGene (wet-lab validated) = Experimental, AGORA2 (flux modeling) = Predicted."""
    df = df.copy()
    if 'Source' in df.columns:
        sources = df['Source'].astype('string').fillna('').astype(str)
        df['evidence_type'] = sources.map({
            'gutMGene': 'Experimental',
            'AGORA2': 'Predicted',
        }).fillna('Predicted').astype(str)
    else:
        df['evidence_type'] = 'Predicted'
    if 'confidence' not in df.columns:
        df['confidence'] = df['evidence_type'].astype(str).map({
            'Experimental': 0.9, 'Predicted': 0.6,
        }).fillna(0.5).astype(float)
    return df


MMI_CONFIG = {
    'type_key': 'mmi',
    'type_name': 'Metabolite–Microbe Interactions (MMI)',
    'short_name': 'MMI',
    'icon': 'fas fa-bacterium',
    'color': '#009E73',
    'description': '77,605 release metabolite–microbe interactions from gutMGene (v2.0, 2023) and AGORA2 (VMH, 2023) covering 1,262 microbes',
    'data_transform': _mmi_transform,
    'precomputed_summary': 'data/module_summaries/mmi_summary.json',
    'columns': {
        'display': ['Metabolite_Name', 'HMDB_ID', 'Microbe_Name', 'Rank',
                     'Relationship_Type', 'Tissue', 'Organism',
                     'Evidence_Level', 'evidence_type', 'Source', 'PMID'],
        'search': ['Metabolite_Name', 'HMDB_ID', 'Microbe_Name', 'Taxonomy_ID'],
    },
    'filters': {
        'organism': True,
        'organism_column': 'Organism',
        'extra': [
            {'label': 'Microbe Rank', 'column': 'Rank', 'type': 'dropdown'},
            {'label': 'Tissue', 'column': 'Tissue', 'type': 'dropdown'},
            {'label': 'Relationship', 'column': 'Relationship_Type', 'type': 'dropdown'},
            {'label': 'Experimental Method', 'column': 'Experimental_Method', 'type': 'dropdown'},
        ],
    },
    'viz': [
        {'type': 'donut', 'column': 'evidence_type', 'title': 'Evidence Type'},
        {'type': 'donut', 'column': 'Source', 'title': 'Data Source'},
        {'type': 'bar_top', 'column': 'Microbe_Name', 'title': 'Top Microbes', 'n': 15},
        {'type': 'donut', 'column': 'Relationship_Type', 'title': 'Relationship Types'},
    ],
    'stat_cards': [
        {'label': 'Total Interactions', 'key': 'total', 'icon': 'fas fa-link'},
        {'label': 'Metabolites', 'column': 'HMDB_ID', 'icon': 'fas fa-atom'},
        {'label': 'Microbes', 'column': 'Microbe_Name', 'icon': 'fas fa-bacterium'},
        {'label': 'Organisms', 'column': 'Organism', 'icon': 'fas fa-globe'},
    ],
    'source_links': {
        'Source': {
            'gutMGene': 'http://bio-annotation.cn/gutmgene/',
            'AGORA2': 'https://www.vmh.life/#microbe',
        },
    },
    'source_versions': {
        'gutMGene': 'v2.0 (2023)',
        'AGORA2': 'VMH 2023',
    },
    'pmid_column': 'PMID',
    'entity_links': {
        'Metabolite_Name': {'route': '/metabolite', 'param': 'id', 'id_col': 'HMDB_ID'},
        'Microbe_Name': {'route': '/microbe', 'param': 'name'},
    },
    'references': [
        {'title': 'gutMGene: a comprehensive database for target genes of gut microbes and microbial metabolites', 'authors': 'Cheng L et al.', 'journal': 'Nucleic Acids Res', 'year': 2022, 'pmid': '34500458'},
        {'title': 'AGORA2: Large scale reconstruction of the microbiome highlights wide-spread drug-metabolising capacities', 'authors': 'Heinken A et al.', 'journal': 'bioRxiv', 'year': 2023, 'pmid': ''},
        {'title': 'VMH: Virtual Metabolic Human database', 'authors': 'Noronha A et al.', 'journal': 'Nucleic Acids Res', 'year': 2019, 'pmid': '30371894'},
    ],
}


def load_mmi_db():
    from pathlib import Path
    import pandas as pd
    rel = Path(__file__).parent.parent.parent / "data" / "databases" / "release" / "coremetdb_mmi.csv"
    v2 = rel if rel.exists() else Path(__file__).parent.parent.parent / "data" / "databases_v2" / "mmi_database_v2.csv"
    if v2.exists():
        from app.services.csv_loader import load_optimized
        return load_optimized(v2)
    from app.services.mmi_service import get_mmi_db
    return get_mmi_db()
