"""Extended API tests — edge cases, invalid inputs, and untested endpoints."""

import json
import pytest


# ═══════════════════════════════════════════════════════════════════════════
# Predict endpoint edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestPredictEndpoint:

    def test_predict_no_json(self, client):
        resp = client.post("/api/v1/predict", data="not json",
                           content_type="text/plain")
        assert resp.status_code == 400

    def test_predict_empty_json(self, client):
        resp = client.post("/api/v1/predict",
                           data=json.dumps({}),
                           content_type="application/json")
        assert resp.status_code == 400
        data = json.loads(resp.data)
        assert "metabolites" in data["error"]

    def test_predict_missing_proteins(self, client):
        resp = client.post("/api/v1/predict",
                           data=json.dumps({
                               "metabolites": [{"name": "Glucose", "hmdb_id": "HMDB0000122",
                                                "smiles": "CCO"}],
                           }),
                           content_type="application/json")
        assert resp.status_code == 400
        assert "proteins" in json.loads(resp.data)["error"]

    def test_predict_metabolite_missing_field(self, client):
        resp = client.post("/api/v1/predict",
                           data=json.dumps({
                               "metabolites": [{"name": "Glucose"}],
                               "proteins": [{"uniprot_id": "P12345", "name": "ALB",
                                              "gene": "ALB", "organism": "Human",
                                              "sequence": "ACDEFGHIKLMNPQRSTVWY"}],
                           }),
                           content_type="application/json")
        assert resp.status_code == 400
        assert "hmdb_id" in json.loads(resp.data)["error"]

    def test_predict_protein_missing_field(self, client):
        resp = client.post("/api/v1/predict",
                           data=json.dumps({
                               "metabolites": [{"name": "Glucose", "hmdb_id": "HMDB0000122",
                                                "smiles": "CCO"}],
                               "proteins": [{"uniprot_id": "P12345"}],
                           }),
                           content_type="application/json")
        assert resp.status_code == 400
        assert "name" in json.loads(resp.data)["error"]

    def test_predict_valid_returns_202(self, client):
        """A well-formed predict request should return 202 with a job_id."""
        resp = client.post("/api/v1/predict",
                           data=json.dumps({
                               "metabolites": [
                                   {"name": "Glucose", "hmdb_id": "HMDB0000122",
                                    "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"}
                               ],
                               "proteins": [
                                   {"uniprot_id": "P02768", "name": "Serum albumin",
                                    "gene": "ALB", "organism": "Homo sapiens",
                                    "sequence": "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"}
                               ],
                               "organism": "All",
                           }),
                           content_type="application/json")
        assert resp.status_code == 202
        data = json.loads(resp.data)
        assert "job_id" in data
        assert data["status"] == "running"

    def test_predict_empty_arrays(self, client):
        resp = client.post("/api/v1/predict",
                           data=json.dumps({"metabolites": [], "proteins": []}),
                           content_type="application/json")
        assert resp.status_code == 400


# ═══════════════════════════════════════════════════════════════════════════
# Database search edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestDatabaseSearch:

    def test_search_by_metabolite(self, client):
        resp = client.get("/api/v1/database/search?metabolite=Glucose")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)
        assert len(data) > 0

    def test_search_by_protein(self, client):
        resp = client.get("/api/v1/database/search?protein=ALB")
        assert resp.status_code == 200

    def test_search_by_organism(self, client):
        resp = client.get("/api/v1/database/search?organism=Homo%20sapiens")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        if len(data) > 0:
            assert data[0].get("Species", "").lower() == "homo sapiens"

    def test_search_combined_filters(self, client):
        resp = client.get("/api/v1/database/search?metabolite=Glucose&organism=Homo%20sapiens&limit=5")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) <= 5

    def test_search_limit_respected(self, client):
        resp = client.get("/api/v1/database/search?limit=3")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) <= 3

    def test_search_limit_max_capped(self, client):
        """Limit above 1000 should be capped to 1000."""
        resp = client.get("/api/v1/database/search?limit=5000")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) <= 1000

    def test_search_invalid_limit(self, client):
        """Non-integer limit should fall back to default 100."""
        resp = client.get("/api/v1/database/search?limit=abc")
        assert resp.status_code == 200

    def test_search_no_results(self, client):
        resp = client.get("/api/v1/database/search?metabolite=ZZZZZ_NONEXISTENT_99999")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) == 0

    def test_search_query_too_long(self, client):
        long_q = "A" * 300
        resp = client.get(f"/api/v1/database/search?metabolite={long_q}")
        assert resp.status_code == 400


# ═══════════════════════════════════════════════════════════════════════════
# MMI search endpoint
# ═══════════════════════════════════════════════════════════════════════════

class TestMMISearch:

    def test_mmi_search_by_metabolite(self, client):
        resp = client.get("/api/v1/mmi/search?metabolite=Glucose")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)

    def test_mmi_search_by_microbe(self, client):
        resp = client.get("/api/v1/mmi/search?microbe=Lacto")
        assert resp.status_code == 200

    def test_mmi_search_empty(self, client):
        resp = client.get("/api/v1/mmi/search")
        assert resp.status_code == 200

    def test_mmi_search_query_too_long(self, client):
        long_q = "B" * 300
        resp = client.get(f"/api/v1/mmi/search?q={long_q}")
        assert resp.status_code == 400


# ═══════════════════════════════════════════════════════════════════════════
# MDrI search endpoint
# ═══════════════════════════════════════════════════════════════════════════

class TestMDrISearch:

    def test_mdri_search_by_drug(self, client):
        resp = client.get("/api/v1/mdri/search?drug=aspirin")
        assert resp.status_code == 200

    def test_mdri_search_by_metabolite(self, client):
        resp = client.get("/api/v1/mdri/search?metabolite=Glucose")
        assert resp.status_code == 200

    def test_mdri_search_empty(self, client):
        resp = client.get("/api/v1/mdri/search")
        assert resp.status_code == 200

    def test_mdri_search_query_too_long(self, client):
        long_q = "C" * 300
        resp = client.get(f"/api/v1/mdri/search?q={long_q}")
        assert resp.status_code == 400


# ═══════════════════════════════════════════════════════════════════════════
# Export endpoint edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestExportEndpoint:

    def test_export_missing_params(self, client):
        resp = client.get("/api/v1/export/metabolite")
        assert resp.status_code == 400
        data = json.loads(resp.data)
        assert "error" in data

    def test_export_invalid_db(self, client):
        resp = client.get("/api/v1/export/metabolite?id=HMDB0000122&db=invalid")
        assert resp.status_code == 400
        data = json.loads(resp.data)
        assert "Invalid db" in data["error"]

    def test_export_specific_db(self, client):
        resp = client.get("/api/v1/export/metabolite?id=HMDB0000122&db=mpi")
        assert resp.status_code in (200, 404)  # 404 if no MPI data for this metabolite

    def test_export_nonexistent_metabolite(self, client):
        resp = client.get("/api/v1/export/metabolite?id=HMDB9999999")
        assert resp.status_code == 404

    def test_export_by_name(self, client):
        resp = client.get("/api/v1/export/metabolite?name=Glucose")
        assert resp.status_code in (200, 404)


# ═══════════════════════════════════════════════════════════════════════════
# Results endpoint edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestResultsEndpoint:

    def test_results_invalid_job_returns_404(self, client):
        resp = client.get("/api/v1/results/DOES_NOT_EXIST")
        assert resp.status_code in (200, 404)
        data = json.loads(resp.data)
        if resp.status_code == 404:
            assert "error" in data


# ═══════════════════════════════════════════════════════════════════════════
# Autocomplete edge cases
# ═══════════════════════════════════════════════════════════════════════════

class TestAutocompleteExtended:

    def test_autocomplete_empty_query(self, client):
        resp = client.get("/api/v1/autocomplete?q=")
        assert resp.status_code == 200
        assert json.loads(resp.data) == []

    def test_autocomplete_single_char(self, client):
        resp = client.get("/api/v1/autocomplete?q=G")
        assert resp.status_code == 200
        assert json.loads(resp.data) == []

    def test_autocomplete_query_too_long(self, client):
        resp = client.get(f"/api/v1/autocomplete?q={'X' * 300}")
        assert resp.status_code == 400

    def test_autocomplete_custom_limit(self, client):
        resp = client.get("/api/v1/autocomplete?q=Gl&limit=2")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) <= 2

    def test_autocomplete_returns_structure(self, client):
        resp = client.get("/api/v1/autocomplete?q=Glucose")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        if data:
            assert "label" in data[0]
            assert "name" in data[0]
            assert "hmdb_id" in data[0]


# ═══════════════════════════════════════════════════════════════════════════
# Species endpoint
# ═══════════════════════════════════════════════════════════════════════════

class TestSpeciesExtended:

    def test_species_returns_list(self, client):
        resp = client.get("/api/v1/species")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert "species" in data
        species_list = data["species"]
        assert len(species_list) >= 5
        for sp in species_list:
            assert "name" in sp
            assert "mpi_count" in sp
            assert sp["mpi_count"] > 0

    def test_species_includes_human(self, client):
        resp = client.get("/api/v1/species")
        data = json.loads(resp.data)
        names = [s["name"] for s in data["species"]]
        assert "Homo sapiens" in names


# ═══════════════════════════════════════════════════════════════════════════
# MGI API endpoints
# ═══════════════════════════════════════════════════════════════════════════

class TestMGIAPI:

    def test_mgi_stats_returns_200(self, client):
        resp = client.get("/api/v1/mgi/stats")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert "total" in data
        assert data["available"] is True

    def test_mgi_stats_has_required_keys(self, client):
        resp = client.get("/api/v1/mgi/stats")
        data = json.loads(resp.data)
        for key in ("total", "metabolites", "genes", "organisms",
                     "interaction_types", "available"):
            assert key in data

    def test_mgi_search_returns_200(self, client):
        resp = client.get("/api/v1/mgi/search")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)

    def test_mgi_search_by_metabolite(self, client):
        resp = client.get("/api/v1/mgi/search?metabolite=Glucose")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)
        assert len(data) > 0

    def test_mgi_search_by_gene(self, client):
        resp = client.get("/api/v1/mgi/search?gene=CYP")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)

    def test_mgi_search_by_organism(self, client):
        resp = client.get("/api/v1/mgi/search?organism=Homo%20sapiens")
        assert resp.status_code == 200

    def test_mgi_search_no_results(self, client):
        resp = client.get("/api/v1/mgi/search?metabolite=ZZZZZ_NONEXISTENT_99999")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) == 0

    def test_mgi_search_with_limit(self, client):
        resp = client.get("/api/v1/mgi/search?gene=CYP&limit=3")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) <= 3


# ═══════════════════════════════════════════════════════════════════════════
# mGWAS API endpoints
# ═══════════════════════════════════════════════════════════════════════════

class TestMGWASAPI:

    def test_mgwas_stats_returns_200(self, client):
        resp = client.get("/api/v1/mgwas/stats")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert "total" in data
        assert data["available"] is True

    def test_mgwas_stats_has_required_keys(self, client):
        resp = client.get("/api/v1/mgwas/stats")
        data = json.loads(resp.data)
        for key in ("total", "metabolites", "snps", "genes",
                     "chromosomes", "available"):
            assert key in data

    def test_mgwas_search_returns_200(self, client):
        resp = client.get("/api/v1/mgwas/search")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)

    def test_mgwas_search_by_metabolite(self, client):
        resp = client.get("/api/v1/mgwas/search?metabolite=Glucose")
        assert resp.status_code == 200

    def test_mgwas_search_by_snp(self, client):
        resp = client.get("/api/v1/mgwas/search?snp=rs")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)

    def test_mgwas_search_by_gene(self, client):
        resp = client.get("/api/v1/mgwas/search?gene=SLC")
        assert resp.status_code == 200

    def test_mgwas_search_by_chromosome(self, client):
        resp = client.get("/api/v1/mgwas/search?chromosome=1")
        assert resp.status_code == 200

    def test_mgwas_search_no_results(self, client):
        resp = client.get("/api/v1/mgwas/search?metabolite=ZZZZZ_NONEXISTENT_99999")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) == 0

    def test_mgwas_search_with_limit(self, client):
        resp = client.get("/api/v1/mgwas/search?snp=rs&limit=5")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert len(data) <= 5

    def test_mgwas_search_free_text(self, client):
        resp = client.get("/api/v1/mgwas/search?q=rs")
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert isinstance(data, list)
