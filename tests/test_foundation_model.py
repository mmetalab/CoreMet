"""Tests for CoreMet Foundation Model.

Covers: model architecture, graph construction, forward pass, decoder outputs.
"""

import sys
from pathlib import Path

import pytest
import torch

PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


class TestFoundationModelArchitecture:
    """Test MetaboliteFoundationModel class."""

    @pytest.fixture(scope="class")
    def model(self):
        from app.models.foundation_model import MetaboliteFoundationModel
        return MetaboliteFoundationModel()

    @pytest.fixture(scope="class")
    def test_graph(self):
        import dgl
        return dgl.heterograph({
            ('metabolite', 'mpi', 'protein'): ([0, 1, 2], [0, 1, 2]),
            ('protein', 'rev_mpi', 'metabolite'): ([0, 1, 2], [0, 1, 2]),
            ('metabolite', 'mdi', 'disease'): ([0, 1], [0, 1]),
            ('disease', 'rev_mdi', 'metabolite'): ([0, 1], [0, 1]),
            ('metabolite', 'mmi', 'microbe'): ([0, 1], [0, 1]),
            ('microbe', 'rev_mmi', 'metabolite'): ([0, 1], [0, 1]),
            ('metabolite', 'mdri', 'drug'): ([0, 1], [0, 1]),
            ('drug', 'rev_mdri', 'metabolite'): ([0, 1], [0, 1]),
            ('metabolite', 'mgi', 'gene'): ([0, 1, 2], [0, 1, 2]),
            ('gene', 'rev_mgi', 'metabolite'): ([0, 1, 2], [0, 1, 2]),
            ('metabolite', 'mgwas', 'snp'): ([0, 1], [0, 1]),
            ('snp', 'rev_mgwas', 'metabolite'): ([0, 1], [0, 1]),
        })

    def test_model_has_7_input_projectors(self, model):
        assert len(model.input_proj) == 7

    def test_model_has_6_decoders(self, model):
        expected = {'mpi', 'mdi', 'mmi', 'mdri', 'mgi', 'mgwas'}
        assert set(model.decoders.keys()) == expected

    def test_parameter_count_reasonable(self, model):
        total = sum(p.numel() for p in model.parameters())
        assert 1_000_000 < total < 5_000_000

    def test_encode_output_shapes(self, model, test_graph):
        feat_dict = {nt: torch.randn(test_graph.num_nodes(nt), 128)
                     for nt in test_graph.ntypes}
        model.eval()
        with torch.no_grad():
            embs = model.encode(test_graph, feat_dict)
        for ntype in test_graph.ntypes:
            assert embs[ntype].shape == (test_graph.num_nodes(ntype), 128)

    def test_forward_returns_scores(self, model, test_graph):
        feat_dict = {nt: torch.randn(test_graph.num_nodes(nt), 128)
                     for nt in test_graph.ntypes}
        pos = {'mpi': (torch.tensor([0, 1]), torch.tensor([0, 1]))}
        neg = {'mpi': (torch.tensor([0, 1]), torch.tensor([2, 0]))}
        model.eval()
        with torch.no_grad():
            pos_s, neg_s, _ = model(test_graph, feat_dict, pos, neg)
        assert pos_s['mpi'].shape == (2,)
        assert neg_s['mpi'].shape == (2,)

    def test_scores_are_finite(self, model, test_graph):
        feat_dict = {nt: torch.randn(test_graph.num_nodes(nt), 128)
                     for nt in test_graph.ntypes}
        pos = {'mdi': (torch.tensor([0]), torch.tensor([0]))}
        neg = {'mdi': (torch.tensor([0]), torch.tensor([1]))}
        model.eval()
        with torch.no_grad():
            pos_s, neg_s, _ = model(test_graph, feat_dict, pos, neg)
        assert torch.isfinite(pos_s['mdi']).all()

    def test_sigmoid_scores_in_range(self, model, test_graph):
        feat_dict = {nt: torch.randn(test_graph.num_nodes(nt), 128)
                     for nt in test_graph.ntypes}
        pos = {'mmi': (torch.tensor([0, 1]), torch.tensor([0, 1]))}
        neg = {'mmi': (torch.tensor([0, 1]), torch.tensor([1, 0]))}
        model.eval()
        with torch.no_grad():
            pos_s, _, _ = model(test_graph, feat_dict, pos, neg)
        probs = torch.sigmoid(pos_s['mmi'])
        assert (probs >= 0).all() and (probs <= 1).all()

    def test_multi_type_forward(self, model, test_graph):
        """Test forward pass with multiple edge types simultaneously."""
        feat_dict = {nt: torch.randn(test_graph.num_nodes(nt), 128)
                     for nt in test_graph.ntypes}
        pos = {
            'mpi': (torch.tensor([0]), torch.tensor([0])),
            'mdi': (torch.tensor([0]), torch.tensor([0])),
            'mgi': (torch.tensor([0]), torch.tensor([0])),
        }
        neg = {
            'mpi': (torch.tensor([0]), torch.tensor([1])),
            'mdi': (torch.tensor([0]), torch.tensor([1])),
            'mgi': (torch.tensor([0]), torch.tensor([1])),
        }
        model.eval()
        with torch.no_grad():
            pos_s, neg_s, embs = model(test_graph, feat_dict, pos, neg)
        assert len(pos_s) == 3
        assert 'metabolite' in embs


class TestGraphConstruction:
    """Test unified heterograph construction from databases."""

    @pytest.fixture(scope="class")
    def graph_data(self):
        from scripts.train_foundation_model import build_unified_heterograph
        return build_unified_heterograph()

    def test_metabolite_nodes_exist(self, graph_data):
        g, _, _, _ = graph_data
        assert 'metabolite' in g.ntypes
        assert g.num_nodes('metabolite') > 1000

    def test_edges_are_bidirectional(self, graph_data):
        g, _, _, _ = graph_data
        etypes = [et[1] for et in g.canonical_etypes]
        for etype in ['mdi', 'mmi', 'mdri', 'mgi', 'mgwas']:
            if etype in etypes:
                assert f'rev_{etype}' in etypes

    def test_node_maps_match_graph(self, graph_data):
        g, node_maps, _, _ = graph_data
        for ntype, nmap in node_maps.items():
            if nmap and ntype in g.ntypes:
                assert len(nmap) == g.num_nodes(ntype)

    def test_smiles_map_not_empty(self, graph_data):
        _, _, _, smiles_map = graph_data
        assert len(smiles_map) > 1000

    def test_total_edges_reasonable(self, graph_data):
        g, _, _, _ = graph_data
        total = sum(g.num_edges(et) for et in g.canonical_etypes) // 2
        assert total > 100_000
