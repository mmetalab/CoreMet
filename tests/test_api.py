"""API endpoint tests."""

import json
import pytest


def test_health_endpoint(client):
    resp = client.get("/api/v1/health")
    assert resp.status_code == 200
    data = json.loads(resp.data)
    assert data.get("status") in ("ok", "healthy")


def test_species_endpoint(client):
    resp = client.get("/api/v1/species")
    assert resp.status_code == 200
    data = json.loads(resp.data)
    assert isinstance(data, (list, dict))


def test_database_search_no_params(client):
    resp = client.get("/api/v1/database/search")
    assert resp.status_code in (200, 400)


def test_mmi_stats(client):
    resp = client.get("/api/v1/mmi/stats")
    assert resp.status_code == 200


def test_mdri_stats(client):
    resp = client.get("/api/v1/mdri/stats")
    assert resp.status_code == 200


def test_export_metabolite_by_id(client):
    resp = client.get("/api/v1/export/metabolite?id=HMDB0000122")
    assert resp.status_code == 200
    assert "text/csv" in resp.content_type or "application/octet-stream" in resp.content_type


def test_results_invalid_job(client):
    resp = client.get("/api/v1/results/nonexistent-job-id")
    assert resp.status_code in (200, 404)


def test_sitemap(client):
    resp = client.get("/sitemap.xml")
    assert resp.status_code == 200
    assert b"<urlset" in resp.data


def test_robots_txt(client):
    resp = client.get("/robots.txt")
    assert resp.status_code == 200
    assert b"Sitemap" in resp.data


def test_autocomplete_short_query(client):
    """Queries shorter than 2 chars should return empty list."""
    resp = client.get("/api/v1/autocomplete?q=A")
    assert resp.status_code == 200
    data = json.loads(resp.data)
    assert data == []


def test_autocomplete_valid_query(client):
    """A reasonable prefix should return results with label/name/hmdb_id."""
    resp = client.get("/api/v1/autocomplete?q=Glut")
    assert resp.status_code == 200
    data = json.loads(resp.data)
    assert isinstance(data, list)
    if len(data) > 0:
        assert "label" in data[0]
        assert "hmdb_id" in data[0]


def test_autocomplete_hmdb_prefix(client):
    """Search by HMDB ID prefix."""
    resp = client.get("/api/v1/autocomplete?q=HMDB00001")
    assert resp.status_code == 200
    data = json.loads(resp.data)
    assert isinstance(data, list)
