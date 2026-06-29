# CoreMet, Dataset (v1, 2026)

A metabolite-centered cross-domain interaction knowledge graph. This deposit contains
the curated, deduplicated interaction tables underlying the CoreMet web resource.
Record counts match exactly those reported in the manuscript and on the website.

**License:** CC BY 4.0 · **Web resource:** https://www.coremet.org · **Code:** [GitHub URL]

## Contents

| File | Interactions | Description |
|------|-------------:|-------------|
| `coremetdb_mpi.csv` | 38,061 | Metabolite–protein (KEGG/Rhea; curated, directly-evidenced only) |
| `coremetdb_mei.csv` | 47,551 | Metabolite–enzyme (EC) catalytic relationships (KEGG/Rhea) |
| `coremetdb_mdi.csv` | 82,882 | Metabolite–disease associations (CTD, HMDB) |
| `coremetdb_mmi.csv` | 77,605 | Metabolite–microbe relationships (gutMGene, AGORA2) |
| `coremetdb_mdri.csv` | 3,500 | Metabolite–drug interactions (DrugBank) |
| `coremetdb_mgi.csv` | 1,658,745 | Metabolite–gene interactions (CTD) |
| `coremetdb_mgwas.csv` | 44,344 | Metabolite–SNP associations (GWAS Catalog) |
| `coremetdb_stats.json` |, | Machine-readable summary statistics (single source of truth) |
| **Total** | **1,952,688** | 30,674 unique metabolites; 8 node types; 7 edge types |

## Common columns

All tables are HMDB-anchored and share a provenance model:

- **HMDB_ID / `HMDB ID`**, metabolite primary key (HMDB accession).
- **Metabolite_Name / SMILES**, name and structure (SMILES where available).
- **Source / Evidence_Source**, originating database(s).
- **Evidence_Level / evidence_type**, evidence tier (experimental, curated, inferred, computational).
- **confidence**, normalized [0,1] score (compare within a layer, not across).
- **PMID / pmid**, PubMed reference where available.

## Per-table key columns

- **MPI** (`HMDB ID`, `Uniprot ID`, `Species`, `Protein Name`, `Gene Name`, `Pathway_Name`, `Evidence_Source`, `interaction_subtype`). Curated layer only; the larger unsourced "pathway_participant" co-membership set is intentionally excluded.
- **MEI** (`HMDB_ID`, `EC_Number`, `Enzyme_Name`, `Uniprot_ID`, `Gene_Name`, `Species`, `Pathway_Name`, `enzyme_role`, `reaction_id`).
- **MDI** (`HMDB_ID`, `Disease_Name`, `Disease_ID`, `MeSH_ID`, `evidence_type`, `confidence`, `pmid`, `Source`).
- **MMI** (`HMDB_ID`, `Microbe_Name`, `Taxonomy_ID`, `Organism`, `Relationship_Type`, `Evidence_Level`, `PMID`, plus KEGG/PubChem/ChEBI cross-refs).
- **MDrI** (`HMDB_ID`, `Drug_Name`, `DrugBank_ID`, `Interaction_Type`, `mechanism_category`, `Evidence_Level`, `PMID`).
- **MGI** (`HMDB_ID`, `Gene_Symbol`, `Gene_ID`, `Organism`, `Interaction_Type`, `Interaction_Actions`, `direction`, `PMID`).
- **mGWAS** (`HMDB_ID`, `rsID`, `Chromosome`, `Position`, `Mapped_Gene`, `P_Value`, `Beta`, `Trait`, `genome_build`, `PMID`).

## Reproducibility

These files are produced from the source databases by the scripts in the code repository:
`compute_db_stats.py` (statistics), `build_release_csvs.py` (these deduplicated tables),
`compute_use_cases.py`, `compute_xrefs.py`. Deduplication is by metabolite–target(–species)
edge key within each layer; MPI is restricted to directly-evidenced interactions.

## How to cite

A formal citation will be provided once the accompanying manuscript is published
(in preparation, 2026). Please cite this dataset by its DOI and link to https://www.coremet.org.
