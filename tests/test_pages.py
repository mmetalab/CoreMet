"""Smoke tests — every page must return HTTP 200."""

import pytest
import pandas as pd

from pages.metabolite_detail import _linkify_entities, _linkify_pmids

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


def test_metabolite_detail_linkifiers_accept_categorical_columns():
    df = pd.DataFrame({
        "Uniprot_ID": pd.Categorical(["P12345"]),
        "PMID": pd.Categorical(["9742976"]),
    })

    linked = _linkify_pmids(_linkify_entities(df))

    assert linked.loc[0, "Uniprot_ID"] == "[P12345](https://www.uniprot.org/uniprot/P12345)"
    assert linked.loc[0, "PMID"] == "[9742976](https://pubmed.ncbi.nlm.nih.gov/9742976)"


def test_unknown_page_returns_200(client):
    """Unknown paths still return the SPA shell (200) with a 404 message."""
    resp = client.get("/nonexistent-page")
    assert resp.status_code == 200
