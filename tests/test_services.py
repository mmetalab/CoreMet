"""Integration tests for CoreMet services.

Covers: data_service, export_service, enrichment_service, molecular_service,
        job_service, mmi_service, mdi_service, mdri_service, mei_service,
        mgi_service, mgwas_service.
"""

import json
import pandas as pd
import pytest
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# ═══════════════════════════════════════════════════════════════════════════
# DataService
# ═══════════════════════════════════════════════════════════════════════════

class TestDataService:

    @pytest.fixture(scope="class")
    def ds(self):
        from app.services.data_service import DataService
        return DataService()

    def test_mpi_db_loaded(self, ds):
        """MPI database should load with >30k interactions."""
        assert ds.mpi_db is not None
        assert len(ds.mpi_db) > 30_000

    def test_mpi_db_columns(self, ds):
        """Verify expected v2 columns exist."""
        required = {"Species", "Metabolite Name", "HMDB ID", "SMILES",
                     "Uniprot ID", "Protein Name", "Gene Name"}
        assert required.issubset(set(ds.mpi_db.columns))

    def test_get_mpi_database_info(self, ds):
        info = ds.get_mpi_database_info()
        assert info["total_interactions"] > 30_000
        assert info["unique_metabolites"] > 1000
        assert info["unique_proteins"] > 10_000
        assert len(info["organisms"]) >= 5

    def test_validate_metabolite_data_valid(self, ds):
        df = pd.DataFrame({
            "Metabolite Name": ["Glucose"],
            "HMDB ID": ["HMDB0000122"],
            "SMILES": ["OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"],
        })
        valid, errors = ds.validate_metabolite_data(df)
        assert valid
        assert errors == []

    def test_validate_metabolite_data_missing_cols(self, ds):
        df = pd.DataFrame({"Name": ["X"]})
        valid, errors = ds.validate_metabolite_data(df)
        assert not valid
        assert any("Missing" in e for e in errors)

    def test_validate_metabolite_data_nulls(self, ds):
        df = pd.DataFrame({
            "Metabolite Name": [None],
            "HMDB ID": ["HMDB0000122"],
            "SMILES": ["CCO"],
        })
        valid, errors = ds.validate_metabolite_data(df)
        assert not valid

    def test_validate_protein_data_valid(self, ds):
        df = pd.DataFrame({
            "UniprotID": ["P12345"],
            "Protein Name": ["Test Protein"],
            "Gene Name": ["TP1"],
            "Organism": ["Homo sapiens"],
            "Sequence": ["ACDEFGHIKLMNPQRSTVWY"],
        })
        valid, errors = ds.validate_protein_data(df)
        assert valid

    def test_validate_protein_data_missing_cols(self, ds):
        df = pd.DataFrame({"UniprotID": ["P12345"]})
        valid, errors = ds.validate_protein_data(df)
        assert not valid

    def test_parse_text_input_metabolite(self, ds):
        text = "Glucose,HMDB0000122,OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        df = ds.parse_text_input(text, "metabolite")
        assert len(df) == 1
        assert df.iloc[0]["HMDB ID"] == "HMDB0000122"

    def test_parse_text_input_protein(self, ds):
        text = "P12345,Albumin,ALB,Homo sapiens,ACDEFGHIKLMNPQRSTVWY"
        df = ds.parse_text_input(text, "protein")
        assert len(df) == 1
        assert df.iloc[0]["UniprotID"] == "P12345"

    def test_parse_text_input_empty(self, ds):
        df = ds.parse_text_input("", "metabolite")
        assert len(df) == 0

    def test_parse_text_input_incomplete_lines(self, ds):
        """Lines with fewer fields than required should be skipped."""
        text = "Glucose,HMDB0000122"  # Only 2 parts, metabolite needs 3
        df = ds.parse_text_input(text, "metabolite")
        assert len(df) == 0

    def test_get_example_data_metabolite(self, ds):
        df = ds.get_example_data("metabolite")
        assert len(df) == 1
        assert "SMILES" in df.columns

    def test_get_example_data_protein(self, ds):
        df = ds.get_example_data("protein")
        assert len(df) == 1
        assert "UniprotID" in df.columns

    def test_get_example_data_unknown_type(self, ds):
        df = ds.get_example_data("unknown")
        assert df.empty

    def test_parse_uploaded_csv(self, ds):
        """Simulate a base64-encoded CSV upload."""
        import base64
        csv_content = "Metabolite Name,HMDB ID,SMILES\nGlucose,HMDB0000122,CCO\n"
        encoded = base64.b64encode(csv_content.encode()).decode()
        contents = f"data:text/csv;base64,{encoded}"
        df = ds.parse_uploaded_file(contents, "test.csv")
        assert len(df) == 1
        assert df.iloc[0]["HMDB ID"] == "HMDB0000122"

    def test_parse_uploaded_unsupported(self, ds):
        """Unsupported file types should raise ValueError."""
        import base64
        encoded = base64.b64encode(b"binary data").decode()
        contents = f"data:application/octet-stream;base64,{encoded}"
        with pytest.raises(ValueError):
            ds.parse_uploaded_file(contents, "test.xyz")


# ═══════════════════════════════════════════════════════════════════════════
# ExportService
# ═══════════════════════════════════════════════════════════════════════════

class TestExportService:

    @pytest.fixture
    def predictions_df(self):
        return pd.DataFrame({
            "Metabolite": ["HMDB0000122", "HMDB0000191"],
            "Protein": ["P12345", "P67890"],
            "Protein Name": ["Albumin", "Hemoglobin"],
            "Prediction Score": [0.95, 0.80],
            "Existing": ["Yes", "No"],
        })

    def test_export_csv(self, predictions_df):
        from app.services.export_service import export_csv
        result = export_csv(predictions_df)
        assert "Metabolite" in result
        assert "HMDB0000122" in result
        lines = result.strip().split("\n")
        assert len(lines) == 3  # header + 2 rows

    def test_export_json(self, predictions_df):
        from app.services.export_service import export_json
        result = export_json(predictions_df)
        data = json.loads(result)
        assert isinstance(data, list)
        assert len(data) == 2
        assert data[0]["Metabolite"] == "HMDB0000122"

    def test_export_graphml(self, predictions_df):
        from app.services.export_service import export_graphml
        result = export_graphml(predictions_df)
        assert "<graphml" in result
        assert "HMDB0000122" in result

    def test_export_sif(self, predictions_df):
        from app.services.export_service import export_sif
        result = export_sif(predictions_df)
        lines = result.strip().split("\n")
        assert len(lines) == 2
        assert "interacts_with" in lines[0]

    def test_export_csv_empty_df(self):
        from app.services.export_service import export_csv
        result = export_csv(pd.DataFrame())
        assert result.strip() == ""

    def test_export_graphml_empty_df(self):
        from app.services.export_service import export_graphml
        result = export_graphml(pd.DataFrame(columns=["Metabolite", "Protein", "Prediction Score"]))
        assert "<graphml" in result


# ═══════════════════════════════════════════════════════════════════════════
# JobService
# ═══════════════════════════════════════════════════════════════════════════

class TestJobService:

    def test_create_and_get_job(self):
        from app.services.job_service import create_job, get_job
        job_id = create_job(
            metabolites_json='[{"name": "Glucose"}]',
            proteins_json='[{"uniprot_id": "P12345"}]',
            organism="Homo sapiens",
        )
        assert len(job_id) == 12

        job = get_job(job_id)
        assert job is not None
        assert job["status"] == "running"
        assert job["organism"] == "Homo sapiens"

    def test_update_job(self):
        from app.services.job_service import create_job, update_job, get_job
        job_id = create_job(
            metabolites_json="[]",
            proteins_json="[]",
            organism="All",
        )
        update_job(job_id, "completed", '{"results": []}')
        job = get_job(job_id)
        assert job["status"] == "completed"

    def test_get_nonexistent_job(self):
        from app.services.job_service import get_job
        assert get_job("nonexistent") is None

    def test_cleanup_expired(self):
        from app.services.job_service import cleanup_expired
        # Should not raise
        cleanup_expired()

    def test_job_has_expiry(self):
        from app.services.job_service import create_job, get_job
        job_id = create_job("[]", "[]", "All")
        job = get_job(job_id)
        assert job["expires_at"] is not None
        assert job["created_at"] < job["expires_at"]


# ═══════════════════════════════════════════════════════════════════════════
# MolecularService
# ═══════════════════════════════════════════════════════════════════════════

class TestMolecularService:

    @pytest.fixture(scope="class")
    def mol_svc(self):
        from app.services.molecular_service import MolecularService
        return MolecularService()

    def test_validate_smiles_valid(self, mol_svc):
        assert mol_svc.validate_smiles("CCO") is True
        assert mol_svc.validate_smiles("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O") is True

    def test_validate_smiles_invalid(self, mol_svc):
        # Note: RDKit considers empty string a valid (empty) molecule
        assert mol_svc.validate_smiles("NOT_A_SMILES_XYZ123") is False

    def test_validate_protein_sequence_valid(self, mol_svc):
        assert mol_svc.validate_protein_sequence("ACDEFGHIKLMNPQRSTVWY") is True

    def test_validate_protein_sequence_too_short(self, mol_svc):
        assert mol_svc.validate_protein_sequence("ACDE") is False

    def test_validate_protein_sequence_invalid_chars(self, mol_svc):
        assert mol_svc.validate_protein_sequence("ACDEFGHIKLX1234NOPQ") is False

    def test_validate_protein_sequence_empty(self, mol_svc):
        assert mol_svc.validate_protein_sequence("") is False
        assert mol_svc.validate_protein_sequence(None) is False

    def test_extract_metabolite_features(self, mol_svc):
        features = mol_svc.extract_metabolite_features("CCO")
        assert features is not None
        assert len(features) > 0

    def test_extract_metabolite_features_invalid(self, mol_svc):
        # Empty string produces all-zero fingerprint → valid PCA output (RDKit quirk)
        assert mol_svc.extract_metabolite_features("NOT_SMILES") is None


# ═══════════════════════════════════════════════════════════════════════════
# MMI / MDI / MDrI / MEI Services (database loading and stats)
# ═══════════════════════════════════════════════════════════════════════════

class TestMMIService:

    def test_get_mmi_db_returns_df(self):
        from app.services.mmi_service import get_mmi_db
        df = get_mmi_db()
        assert isinstance(df, pd.DataFrame)

    def test_get_mmi_stats(self):
        from app.services.mmi_service import get_mmi_stats
        stats = get_mmi_stats()
        assert "total" in stats
        assert isinstance(stats["total"], int)

    def test_mmi_db_columns(self):
        from app.services.mmi_service import get_mmi_db
        df = get_mmi_db()
        if not df.empty:
            assert "Metabolite_Name" in df.columns
            assert "Microbe_Name" in df.columns

    def test_search_mmi_by_query(self):
        from app.services.mmi_service import search_mmi
        result = search_mmi(query="Glucose")
        assert isinstance(result, pd.DataFrame)

    def test_search_mmi_empty(self):
        from app.services.mmi_service import search_mmi
        result = search_mmi()
        assert isinstance(result, pd.DataFrame)

    def test_search_mmi_no_results(self):
        from app.services.mmi_service import search_mmi
        result = search_mmi(query="ZZZZZ_NONEXISTENT_99999")
        assert len(result) == 0

    def test_get_microbes_for_metabolite(self):
        from app.services.mmi_service import get_microbes_for_metabolite
        result = get_microbes_for_metabolite("Glucose")
        assert isinstance(result, pd.DataFrame)

    def test_get_metabolites_for_microbe(self):
        from app.services.mmi_service import get_metabolites_for_microbe
        result = get_metabolites_for_microbe("Lactobacillus")
        assert isinstance(result, pd.DataFrame)

    def test_get_unique_microbes(self):
        from app.services.mmi_service import get_unique_microbes
        result = get_unique_microbes()
        assert isinstance(result, list)

    def test_get_unique_metabolites(self):
        from app.services.mmi_service import get_unique_metabolites
        result = get_unique_metabolites()
        assert isinstance(result, list)

    def test_annotate_metabolites_with_microbes(self):
        from app.services.mmi_service import annotate_metabolites_with_microbes
        result = annotate_metabolites_with_microbes(["Glucose"])
        assert isinstance(result, pd.DataFrame)


class TestMDIService:

    def test_get_mdi_db_returns_df(self):
        from app.services.mdi_service import get_mdi_db
        df = get_mdi_db()
        assert isinstance(df, pd.DataFrame)

    def test_get_mdi_stats(self):
        from app.services.mdi_service import get_mdi_stats
        stats = get_mdi_stats()
        assert "total" in stats

    def test_search_mdi_by_metabolite(self):
        from app.services.mdi_service import search_mdi
        result = search_mdi(metabolite="Glucose")
        assert isinstance(result, pd.DataFrame)

    def test_search_mdi_by_disease(self):
        from app.services.mdi_service import search_mdi
        result = search_mdi(disease="cancer")
        assert isinstance(result, pd.DataFrame)

    def test_search_mdi_no_results(self):
        from app.services.mdi_service import search_mdi
        result = search_mdi(metabolite="ZZZZZ_NONEXISTENT")
        assert len(result) == 0

    def test_get_metabolites_for_disease(self):
        from app.services.mdi_service import get_metabolites_for_disease, get_mdi_db
        df = get_mdi_db()
        if df.empty:
            pytest.skip("MDI database not available")
        disease = df["Disease_Name"].iloc[0]
        result = get_metabolites_for_disease(disease)
        assert len(result) > 0

    def test_get_diseases_for_metabolite(self):
        from app.services.mdi_service import get_diseases_for_metabolite, get_mdi_db
        df = get_mdi_db()
        if df.empty:
            pytest.skip("MDI database not available")
        hmdb = df["HMDB_ID"].iloc[0]
        result = get_diseases_for_metabolite(hmdb_id=hmdb)
        assert len(result) > 0

    def test_annotate_metabolites_with_diseases(self):
        from app.services.mdi_service import annotate_metabolites_with_diseases, get_mdi_db
        df = get_mdi_db()
        if df.empty:
            pytest.skip("MDI database not available")
        hmdb_ids = df["HMDB_ID"].dropna().unique()[:3].tolist()
        result = annotate_metabolites_with_diseases(hmdb_ids)
        assert isinstance(result, dict)
        assert len(result) > 0
        first = list(result.values())[0]
        assert "diseases" in first
        assert "count" in first


class TestMDrIService:

    def test_get_mdri_db_returns_df(self):
        from app.services.mdri_service import get_mdri_db
        df = get_mdri_db()
        assert isinstance(df, pd.DataFrame)

    def test_get_mdri_stats(self):
        from app.services.mdri_service import get_mdri_stats
        stats = get_mdri_stats()
        assert "total" in stats

    def test_search_mdri_by_query(self):
        from app.services.mdri_service import search_mdri
        result = search_mdri(query="aspirin")
        assert isinstance(result, pd.DataFrame)

    def test_search_mdri_no_results(self):
        from app.services.mdri_service import search_mdri
        result = search_mdri(query="ZZZZZ_NONEXISTENT")
        assert len(result) == 0

    def test_search_mdri_with_filters(self):
        from app.services.mdri_service import search_mdri, get_mdri_db
        df = get_mdri_db()
        if df.empty:
            pytest.skip("MDrI database not available")
        result = search_mdri(filters={"Interaction_Type": "PK"})
        assert isinstance(result, pd.DataFrame)


class TestMEIService:

    def test_get_mei_db_returns_df(self):
        from app.services.mei_service import get_mei_db
        df = get_mei_db()
        assert isinstance(df, pd.DataFrame)

    def test_get_mei_stats(self):
        from app.services.mei_service import get_mei_stats
        stats = get_mei_stats()
        assert "total" in stats

    def test_search_mei_by_metabolite(self):
        from app.services.mei_service import search_mei
        result = search_mei(metabolite="Glucose")
        assert isinstance(result, pd.DataFrame)

    def test_search_mei_by_enzyme(self):
        from app.services.mei_service import search_mei
        result = search_mei(enzyme="kinase")
        assert isinstance(result, pd.DataFrame)

    def test_search_mei_no_results(self):
        from app.services.mei_service import search_mei
        result = search_mei(metabolite="ZZZZZ_NONEXISTENT")
        assert len(result) == 0


# ═══════════════════════════════════════════════════════════════════════════
# EnrichmentService
# ═══════════════════════════════════════════════════════════════════════════

class TestEnrichmentService:

    def test_load_pathway_annotations(self):
        from app.services.enrichment_service import load_pathway_annotations
        ec_pathways, pw_names = load_pathway_annotations()
        assert isinstance(ec_pathways, dict)
        assert isinstance(pw_names, dict)

    def test_load_uniprot_to_ec(self):
        from app.services.enrichment_service import load_uniprot_to_ec
        mapping = load_uniprot_to_ec()
        assert isinstance(mapping, dict)

    def test_run_enrichment_empty_df(self):
        from app.services.enrichment_service import run_enrichment
        empty = pd.DataFrame(columns=["Protein"])
        result = run_enrichment(empty)
        assert isinstance(result, pd.DataFrame)

    def test_run_enrichment_no_proteins(self):
        from app.services.enrichment_service import run_enrichment
        df = pd.DataFrame({"Protein": [], "Prediction Score": []})
        result = run_enrichment(df)
        assert isinstance(result, pd.DataFrame)

    def test_run_disease_enrichment_empty(self):
        from app.services.enrichment_service import run_disease_enrichment
        result = run_disease_enrichment([])
        assert isinstance(result, pd.DataFrame)

    def test_run_disease_enrichment_nonexistent(self):
        from app.services.enrichment_service import run_disease_enrichment
        result = run_disease_enrichment(["HMDB_FAKE_999999"])
        assert isinstance(result, pd.DataFrame)

    def test_run_microbe_enrichment_empty(self):
        from app.services.enrichment_service import run_microbe_enrichment
        result = run_microbe_enrichment([])
        assert isinstance(result, pd.DataFrame)

    def test_run_drug_enrichment_empty(self):
        from app.services.enrichment_service import run_drug_enrichment
        result = run_drug_enrichment([])
        assert isinstance(result, pd.DataFrame)

    def test_run_disease_enrichment_with_real_metabolite(self):
        """Use a known metabolite to verify enrichment produces results."""
        from app.services.enrichment_service import run_disease_enrichment
        from app.services.mdi_service import get_mdi_db
        mdi = get_mdi_db()
        if mdi.empty:
            pytest.skip("MDI database not available")
        # Pick a metabolite that appears frequently in MDI
        top_met = mdi["HMDB_ID"].value_counts().index[0]
        result = run_disease_enrichment([top_met])
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "Disease_Name" in result.columns
        assert "FDR" in result.columns


# ═══════════════════════════════════════════════════════════════════════════
# MGI Service (metabolite–gene interactions)
# ═══════════════════════════════════════════════════════════════════════════

class TestMGIService:

    def test_get_mgi_db_returns_df(self):
        from app.services.mgi_service import get_mgi_db
        df = get_mgi_db()
        assert isinstance(df, pd.DataFrame)

    def test_mgi_db_not_empty(self):
        from app.services.mgi_service import get_mgi_db
        df = get_mgi_db()
        assert len(df) > 100_000, "MGI database should have >100K interactions"

    def test_mgi_db_columns(self):
        from app.services.mgi_service import get_mgi_db
        df = get_mgi_db()
        if not df.empty:
            required = {"HMDB_ID", "Metabolite_Name", "Gene_Symbol", "Gene_ID",
                        "Organism", "Interaction_Type"}
            assert required.issubset(set(df.columns))

    def test_get_mgi_stats(self):
        from app.services.mgi_service import get_mgi_stats
        stats = get_mgi_stats()
        assert "total" in stats
        assert isinstance(stats["total"], int)
        assert stats["total"] > 100_000
        assert stats["available"] is True

    def test_get_mgi_stats_keys(self):
        from app.services.mgi_service import get_mgi_stats
        stats = get_mgi_stats()
        for key in ("total", "metabolites", "genes", "organisms",
                     "interaction_types", "available"):
            assert key in stats

    def test_search_mgi_by_metabolite(self):
        from app.services.mgi_service import search_mgi
        result = search_mgi(metabolite="Glucose")
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_search_mgi_by_gene(self):
        from app.services.mgi_service import search_mgi
        result = search_mgi(gene="CYP1A2")
        assert isinstance(result, pd.DataFrame)

    def test_search_mgi_by_organism(self):
        from app.services.mgi_service import search_mgi
        result = search_mgi(organism="Homo sapiens")
        assert isinstance(result, pd.DataFrame)
        if len(result) > 0:
            assert all(result["Organism"].str.lower() == "homo sapiens")

    def test_search_mgi_no_results(self):
        from app.services.mgi_service import search_mgi
        result = search_mgi(metabolite="ZZZZZ_NONEXISTENT_99999")
        assert len(result) == 0

    def test_search_mgi_with_limit(self):
        from app.services.mgi_service import search_mgi
        result = search_mgi(metabolite="Glucose", limit=5)
        assert len(result) <= 5

    def test_get_genes_for_metabolite_by_hmdb(self):
        from app.services.mgi_service import get_genes_for_metabolite, get_mgi_db
        df = get_mgi_db()
        if df.empty:
            pytest.skip("MGI database not available")
        hmdb = df["HMDB_ID"].iloc[0]
        result = get_genes_for_metabolite(hmdb_id=hmdb)
        assert len(result) > 0

    def test_get_genes_for_metabolite_by_name(self):
        from app.services.mgi_service import get_genes_for_metabolite, get_mgi_db
        df = get_mgi_db()
        if df.empty:
            pytest.skip("MGI database not available")
        name = df["Metabolite_Name"].iloc[0]
        result = get_genes_for_metabolite(metabolite_name=name)
        assert len(result) > 0

    def test_get_metabolites_for_gene(self):
        from app.services.mgi_service import get_metabolites_for_gene, get_mgi_db
        df = get_mgi_db()
        if df.empty:
            pytest.skip("MGI database not available")
        gene = df["Gene_Symbol"].iloc[0]
        result = get_metabolites_for_gene(gene)
        assert len(result) > 0

    def test_annotate_metabolites_with_genes(self):
        from app.services.mgi_service import annotate_metabolites_with_genes, get_mgi_db
        df = get_mgi_db()
        if df.empty:
            pytest.skip("MGI database not available")
        hmdb_ids = df["HMDB_ID"].dropna().unique()[:3].tolist()
        result = annotate_metabolites_with_genes(hmdb_ids)
        assert isinstance(result, dict)
        assert len(result) > 0
        first = list(result.values())[0]
        assert "genes" in first
        assert "count" in first
        assert "organisms" in first

    def test_annotate_metabolites_with_genes_empty(self):
        from app.services.mgi_service import annotate_metabolites_with_genes
        result = annotate_metabolites_with_genes([])
        assert result == {}

    def test_annotate_metabolites_with_genes_nonexistent(self):
        from app.services.mgi_service import annotate_metabolites_with_genes
        result = annotate_metabolites_with_genes(["HMDB_FAKE_999999"])
        assert result == {}


# ═══════════════════════════════════════════════════════════════════════════
# mGWAS Service (metabolite–SNP associations)
# ═══════════════════════════════════════════════════════════════════════════

class TestMGWASService:

    def test_get_mgwas_db_returns_df(self):
        from app.services.mgwas_service import get_mgwas_db
        df = get_mgwas_db()
        assert isinstance(df, pd.DataFrame)

    def test_mgwas_db_not_empty(self):
        from app.services.mgwas_service import get_mgwas_db
        df = get_mgwas_db()
        assert len(df) > 10_000, "mGWAS database should have >10K associations"

    def test_mgwas_db_columns(self):
        from app.services.mgwas_service import get_mgwas_db
        df = get_mgwas_db()
        if not df.empty:
            required = {"HMDB_ID", "Metabolite_Name", "rsID", "Chromosome",
                        "Mapped_Gene"}
            assert required.issubset(set(df.columns))

    def test_get_mgwas_stats(self):
        from app.services.mgwas_service import get_mgwas_stats
        stats = get_mgwas_stats()
        assert "total" in stats
        assert isinstance(stats["total"], int)
        assert stats["total"] > 10_000
        assert stats["available"] is True

    def test_get_mgwas_stats_keys(self):
        from app.services.mgwas_service import get_mgwas_stats
        stats = get_mgwas_stats()
        for key in ("total", "metabolites", "snps", "genes",
                     "chromosomes", "available"):
            assert key in stats

    def test_search_mgwas_by_metabolite(self):
        from app.services.mgwas_service import search_mgwas
        result = search_mgwas(metabolite="Glucose")
        assert isinstance(result, pd.DataFrame)

    def test_search_mgwas_by_snp(self):
        from app.services.mgwas_service import search_mgwas, get_mgwas_db
        df = get_mgwas_db()
        if df.empty:
            pytest.skip("mGWAS database not available")
        rsid = df["rsID"].iloc[0]
        result = search_mgwas(snp=rsid)
        assert len(result) > 0

    def test_search_mgwas_by_gene(self):
        from app.services.mgwas_service import search_mgwas, get_mgwas_db
        df = get_mgwas_db()
        if df.empty:
            pytest.skip("mGWAS database not available")
        gene = df["Mapped_Gene"].dropna().iloc[0]
        result = search_mgwas(gene=gene)
        assert len(result) > 0

    def test_search_mgwas_by_chromosome(self):
        from app.services.mgwas_service import search_mgwas
        result = search_mgwas(chromosome="1")
        assert isinstance(result, pd.DataFrame)

    def test_search_mgwas_no_results(self):
        from app.services.mgwas_service import search_mgwas
        result = search_mgwas(metabolite="ZZZZZ_NONEXISTENT_99999")
        assert len(result) == 0

    def test_search_mgwas_with_limit(self):
        from app.services.mgwas_service import search_mgwas
        result = search_mgwas(limit=5)
        assert len(result) <= 5

    def test_get_snps_for_metabolite_by_hmdb(self):
        from app.services.mgwas_service import get_snps_for_metabolite, get_mgwas_db
        df = get_mgwas_db()
        if df.empty:
            pytest.skip("mGWAS database not available")
        hmdb = df["HMDB_ID"].iloc[0]
        result = get_snps_for_metabolite(hmdb_id=hmdb)
        assert len(result) > 0

    def test_get_snps_for_metabolite_by_name(self):
        from app.services.mgwas_service import get_snps_for_metabolite, get_mgwas_db
        df = get_mgwas_db()
        if df.empty:
            pytest.skip("mGWAS database not available")
        name = df["Metabolite_Name"].iloc[0]
        result = get_snps_for_metabolite(metabolite_name=name)
        assert len(result) > 0

    def test_get_metabolites_for_snp(self):
        from app.services.mgwas_service import get_metabolites_for_snp, get_mgwas_db
        df = get_mgwas_db()
        if df.empty:
            pytest.skip("mGWAS database not available")
        rsid = df["rsID"].iloc[0]
        result = get_metabolites_for_snp(rsid)
        assert len(result) > 0

    def test_annotate_metabolites_with_snps(self):
        from app.services.mgwas_service import annotate_metabolites_with_snps, get_mgwas_db
        df = get_mgwas_db()
        if df.empty:
            pytest.skip("mGWAS database not available")
        hmdb_ids = df["HMDB_ID"].dropna().unique()[:3].tolist()
        result = annotate_metabolites_with_snps(hmdb_ids)
        assert isinstance(result, dict)
        assert len(result) > 0
        first = list(result.values())[0]
        assert "snps" in first
        assert "count" in first
        assert "genes" in first

    def test_annotate_metabolites_with_snps_empty(self):
        from app.services.mgwas_service import annotate_metabolites_with_snps
        result = annotate_metabolites_with_snps([])
        assert result == {}

    def test_annotate_metabolites_with_snps_nonexistent(self):
        from app.services.mgwas_service import annotate_metabolites_with_snps
        result = annotate_metabolites_with_snps(["HMDB_FAKE_999999"])
        assert result == {}


# ═══════════════════════════════════════════════════════════════════════════
# Gene Enrichment (cross-database enrichment using MGI)
# ═══════════════════════════════════════════════════════════════════════════

class TestGeneEnrichment:

    def test_run_gene_enrichment_empty(self):
        from app.services.enrichment_service import run_gene_enrichment
        result = run_gene_enrichment([])
        assert isinstance(result, pd.DataFrame)

    def test_run_gene_enrichment_nonexistent(self):
        from app.services.enrichment_service import run_gene_enrichment
        result = run_gene_enrichment(["HMDB_FAKE_999999"])
        assert isinstance(result, pd.DataFrame)

    def test_run_gene_enrichment_with_real_metabolite(self):
        from app.services.enrichment_service import run_gene_enrichment
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if mgi.empty:
            pytest.skip("MGI database not available")
        # Pick a metabolite that appears frequently in MGI
        top_met = mgi["HMDB_ID"].value_counts().index[0]
        result = run_gene_enrichment([top_met])
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "Gene_Symbol" in result.columns
        assert "FDR" in result.columns

    def test_run_gene_enrichment_returns_expected_columns(self):
        from app.services.enrichment_service import run_gene_enrichment
        from app.services.mgi_service import get_mgi_db
        mgi = get_mgi_db()
        if mgi.empty:
            pytest.skip("MGI database not available")
        top_met = mgi["HMDB_ID"].value_counts().index[0]
        result = run_gene_enrichment([top_met])
        if not result.empty:
            expected_cols = {"Gene_Symbol", "Fold_Enrichment", "P_value",
                             "Metabolite_Count", "Background_Count", "FDR"}
            assert expected_cols.issubset(set(result.columns))


# ═══════════════════════════════════════════════════════════════════════════
# MultiPredictionService
# ═══════════════════════════════════════════════════════════════════════════

class TestMultiPredictionService:

    @pytest.fixture(scope="class")
    def svc(self):
        from app.services.multi_prediction_service import MultiPredictionService
        return MultiPredictionService()

    def test_get_available_models(self, svc):
        models = svc.get_available_models()
        assert isinstance(models, list)
        assert len(models) >= 4  # at least mpi, mdi, mmi, mdri

    def test_available_models_have_keys(self, svc):
        models = svc.get_available_models()
        for m in models:
            assert "key" in m
            assert "label" in m
            assert "src_type" in m
            assert "dst_type" in m

    def test_mpi_model_available(self, svc):
        keys = [m["key"] for m in svc.get_available_models()]
        assert "mpi" in keys

    def test_mgi_model_available(self, svc):
        keys = [m["key"] for m in svc.get_available_models()]
        assert "mgi" in keys

    def test_mgwas_model_available(self, svc):
        keys = [m["key"] for m in svc.get_available_models()]
        assert "mgwas" in keys

    def test_predict_mpi_empty_sources(self, svc):
        result = svc.predict("mpi", [])
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_predict_mpi_invalid_source(self, svc):
        result = svc.predict("mpi", ["NONEXISTENT_ID_12345"])
        assert isinstance(result, pd.DataFrame)
        assert result.empty

    def test_predict_mpi_valid(self, svc):
        # Use a metabolite known to be in the MPI graph
        feat = pd.read_pickle(
            PROJECT_ROOT / "data" / "processed" / "features" / "pca_feature_df_All.pkl"
        )
        if "node" in feat.columns:
            mets = [n for n in feat["node"] if str(n).startswith("HMDB")]
            if mets:
                result = svc.predict("mpi", [mets[0]])
                assert isinstance(result, pd.DataFrame)
                assert len(result) > 0
                assert "Source" in result.columns
                assert "Target" in result.columns
                assert "Score" in result.columns

    def test_predict_mgi_valid(self, svc):
        feat = pd.read_pickle(
            PROJECT_ROOT / "data" / "processed" / "features" / "mgi_feature_df.pkl"
        )
        mets = [n for n in feat["node"] if str(n).startswith("HMDB")]
        if not mets:
            pytest.skip("No metabolites in MGI graph")
        result = svc.predict("mgi", [mets[0]], ["CYP1A2"])
        assert isinstance(result, pd.DataFrame)
        # CYP1A2 should resolve to composite keys like CYP1A2|Homo sapiens
        assert len(result) > 0

    def test_predict_mgwas_valid(self, svc):
        feat = pd.read_pickle(
            PROJECT_ROOT / "data" / "processed" / "features" / "mgwas_feature_df.pkl"
        )
        mets = [n for n in feat["node"] if str(n).startswith("HMDB")]
        if not mets:
            pytest.skip("No metabolites in mGWAS graph")
        result = svc.predict("mgwas", [mets[0]], ["rs1260326"])
        assert isinstance(result, pd.DataFrame)

    def test_predict_with_metadata_mpi(self, svc):
        feat = pd.read_pickle(
            PROJECT_ROOT / "data" / "processed" / "features" / "pca_feature_df_All.pkl"
        )
        if "node" in feat.columns:
            mets = [n for n in feat["node"] if str(n).startswith("HMDB")]
            if mets:
                result = svc.predict_with_metadata("mpi", [mets[0]])
                assert isinstance(result, pd.DataFrame)
                if not result.empty:
                    assert "Metabolite" in result.columns or "Source" in result.columns

    def test_predict_with_metadata_mgi(self, svc):
        result = svc.predict_with_metadata("mgi", ["HMDB0001877"], ["CYP1A2"])
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert "Gene_Symbol" in result.columns
            assert "Organism" in result.columns

    def test_predict_with_metadata_mgwas(self, svc):
        result = svc.predict_with_metadata("mgwas", ["HMDB0000064"], ["rs1260326"])
        assert isinstance(result, pd.DataFrame)
        if not result.empty:
            assert "rsID" in result.columns or "Target" in result.columns

    def test_predict_returns_known_flag(self, svc):
        result = svc.predict("mgi", ["HMDB0001877"], ["CYP1A2"])
        if not result.empty:
            assert "Known" in result.columns
            assert result["Known"].isin(["Yes", "No"]).all()

    def test_predict_scores_between_0_and_1(self, svc):
        result = svc.predict("mgi", ["HMDB0001877"], ["CYP1A2"])
        if not result.empty:
            assert (result["Score"] >= 0).all()
            assert (result["Score"] <= 1).all()

    def test_predict_sorted_by_score(self, svc):
        result = svc.predict("mgi", ["HMDB0001877"], ["CYP1A2", "TP53", "TNF"])
        if len(result) > 1:
            scores = result["Score"].values
            assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))
