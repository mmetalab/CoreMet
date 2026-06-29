"""Unit tests for metabolite lookup logic."""

import pytest
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def test_lookup_by_hmdb_id():
    from pages.metabolite_detail import _lookup_metabolite
    data = _lookup_metabolite(hmdb_id="HMDB0000122")
    assert data["hmdb_id"] == "HMDB0000122"
    assert data["name"]  # should have a name
    assert data["smiles"]  # should have SMILES
    total = sum(len(data[k]) for k in ["mpi", "mei", "mdi", "mmi", "mdri"])
    assert total > 0, "Glucose should have interactions"


def test_lookup_by_name():
    from pages.metabolite_detail import _lookup_metabolite
    data = _lookup_metabolite(name="Glucose")
    assert data["name"].lower() == "glucose" or "glucose" in data["name"].lower()
    assert data["hmdb_id"]


def test_lookup_not_found():
    from pages.metabolite_detail import _lookup_metabolite
    data = _lookup_metabolite(hmdb_id="HMDB9999999")
    total = sum(len(data[k]) for k in ["mpi", "mei", "mdi", "mmi", "mdri"])
    assert total == 0


def test_lookup_empty_query():
    from pages.metabolite_detail import _lookup_metabolite
    data = _lookup_metabolite()
    total = sum(len(data[k]) for k in ["mpi", "mei", "mdi", "mmi", "mdri"])
    assert total == 0


def test_chem_properties():
    from pages.metabolite_detail import _compute_chem_properties
    props = _compute_chem_properties("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
    assert props["molecular_formula"] == "C6H12O6"
    assert props["mol_weight"] > 170
    assert "logp" in props


def test_chem_properties_invalid():
    from pages.metabolite_detail import _compute_chem_properties
    assert _compute_chem_properties("") == {}
    assert _compute_chem_properties("INVALID_SMILES_XYZ") == {}
