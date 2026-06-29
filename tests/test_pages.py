"""Smoke tests — every page must return HTTP 200."""

import pytest

STATIC_PAGES = [
    "/",
    "/home",
    "/database",
    "/predict",
    "/disease",
    "/enrichment",
    "/network",
    "/profile",
    "/gallery",
    "/help",
    "/documentation",
    "/about",
    "/batch-search",
    "/api-docs",
]


@pytest.mark.parametrize("path", STATIC_PAGES)
def test_static_page_200(client, path):
    resp = client.get(path)
    assert resp.status_code == 200, f"{path} returned {resp.status_code}"


def test_metabolite_detail_by_id(client):
    resp = client.get("/metabolite?id=HMDB0000122")
    assert resp.status_code == 200


def test_metabolite_detail_by_name(client):
    resp = client.get("/metabolite?name=Glucose")
    assert resp.status_code == 200


def test_metabolite_detail_empty(client):
    resp = client.get("/metabolite")
    assert resp.status_code == 200


def test_unknown_page_returns_200(client):
    """Unknown paths still return the SPA shell (200) with a 404 message."""
    resp = client.get("/nonexistent-page")
    assert resp.status_code == 200
