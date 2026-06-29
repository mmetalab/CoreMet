"""
CoreMet Foundation Model — Unified Heterogeneous Graph Neural Network
=====================================================================

Unifies 6 separate VGAE models (MPI, MDI, MMI, MDrI, MGI, mGWAS) into
a single HeteroGraphConv encoder with type-specific decoders.

Architecture:
    7 node types × 12 edge types (6 forward + 6 reverse)
    → Per-type input projector (Linear + LayerNorm + ReLU)
    → HeteroGraphConv Layer 1 (SAGEConv 128→256 per edge type)
    → HeteroGraphConv Layer 2 (SAGEConv 256→128 per edge type)
    → 6 type-specific MLP decoders (reuse MLPPredictor)

~1.94M parameters, fits in ~500 MB RAM.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

import dgl
import dgl.nn as dglnn

from app.models.vgae_model import MLPPredictor


# Node and edge type constants
NODE_TYPES = ['metabolite', 'protein', 'disease', 'microbe', 'drug', 'gene', 'snp']

EDGE_TYPES = ['mpi', 'mdi', 'mmi', 'mdri', 'mgi', 'mgwas']

# Maps edge type → (src_ntype, etype, dst_ntype)
CANONICAL_ETYPES = {
    'mpi':   ('metabolite', 'mpi',   'protein'),
    'mdi':   ('metabolite', 'mdi',   'disease'),
    'mmi':   ('metabolite', 'mmi',   'microbe'),
    'mdri':  ('metabolite', 'mdri',  'drug'),
    'mgi':   ('metabolite', 'mgi',   'gene'),
    'mgwas': ('metabolite', 'mgwas', 'snp'),
}

# Reverse edge types for bidirectional message passing
REVERSE_ETYPES = {
    'rev_mpi':   ('protein',  'rev_mpi',   'metabolite'),
    'rev_mdi':   ('disease',  'rev_mdi',   'metabolite'),
    'rev_mmi':   ('microbe',  'rev_mmi',   'metabolite'),
    'rev_mdri':  ('drug',     'rev_mdri',  'metabolite'),
    'rev_mgi':   ('gene',     'rev_mgi',   'metabolite'),
    'rev_mgwas': ('snp',      'rev_mgwas', 'metabolite'),
}

ALL_CANONICAL_ETYPES = list(CANONICAL_ETYPES.values()) + list(REVERSE_ETYPES.values())


class MetaboliteFoundationModel(nn.Module):
    """Unified heterogeneous GNN for multi-type metabolite interaction prediction.

    Args:
        in_feats: Input feature dimension (default 128).
        hidden_feats: Hidden layer dimension (default 256).
        out_feats: Output embedding dimension (default 128).
        dropout: Dropout rate (default 0.2).
        node_types: List of node type names.
        canonical_etypes: List of (src, etype, dst) tuples.
    """

    def __init__(
        self,
        in_feats: int = 128,
        hidden_feats: int = 256,
        out_feats: int = 128,
        dropout: float = 0.2,
        node_types: list = None,
        canonical_etypes: list = None,
        use_input_proj: bool = True,
        shared_decoder: bool = False,
        encoder_type: str = 'sage',
    ):
        super().__init__()
        self.in_feats = in_feats
        self.hidden_feats = hidden_feats
        self.out_feats = out_feats
        self.use_input_proj = use_input_proj
        self.shared_decoder = shared_decoder
        self.encoder_type = encoder_type

        if node_types is None:
            node_types = NODE_TYPES
        if canonical_etypes is None:
            canonical_etypes = ALL_CANONICAL_ETYPES

        # Per-type input projector: raw features → shared space
        if use_input_proj:
            self.input_proj = nn.ModuleDict({
                ntype: nn.Sequential(
                    nn.Linear(in_feats, in_feats),
                    nn.LayerNorm(in_feats),
                    nn.ReLU(),
                )
                for ntype in node_types
            })
        else:
            self.input_proj = nn.ModuleDict()  # empty, skip projection

        # Build encoder layers based on encoder_type
        def _make_conv(d_in, d_out):
            if encoder_type == 'sage':
                return dglnn.HeteroGraphConv(
                    {etype: dglnn.SAGEConv(d_in, d_out, 'mean')
                     for _, etype, _ in canonical_etypes},
                    aggregate='sum')
            elif encoder_type == 'gat':
                num_heads = 4
                d_head = d_out // num_heads
                return dglnn.HeteroGraphConv(
                    {etype: dglnn.GATConv(d_in, d_head, num_heads,
                                           allow_zero_in_degree=True)
                     for _, etype, _ in canonical_etypes},
                    aggregate='sum')
            elif encoder_type == 'gcn':
                return dglnn.HeteroGraphConv(
                    {etype: dglnn.GraphConv(d_in, d_out,
                                             allow_zero_in_degree=True)
                     for _, etype, _ in canonical_etypes},
                    aggregate='sum')
            else:
                raise ValueError(f"Unknown encoder_type: {encoder_type}")

        self.conv1 = _make_conv(in_feats, hidden_feats)
        self.bn1 = nn.ModuleDict({
            ntype: nn.BatchNorm1d(hidden_feats) for ntype in node_types
        })

        self.conv2 = _make_conv(hidden_feats, out_feats)
        self.bn2 = nn.ModuleDict({
            ntype: nn.BatchNorm1d(out_feats) for ntype in node_types
        })

        self.dropout = nn.Dropout(dropout)

        # Decoders: type-specific (default) or shared single decoder
        if shared_decoder:
            _shared = MLPPredictor(out_feats)
            self.decoders = nn.ModuleDict({
                etype: _shared for etype in EDGE_TYPES
            })
        else:
            self.decoders = nn.ModuleDict({
                etype: MLPPredictor(out_feats) for etype in EDGE_TYPES
            })

    def encode(self, g, feat_dict):
        """Encode all nodes into unified embedding space.

        Args:
            g: DGL heterogeneous graph.
            feat_dict: {ntype: Tensor[N, in_feats]} raw features.

        Returns:
            {ntype: Tensor[N, out_feats]} embeddings.
        """
        # Project inputs per type (skip if ablation disables projectors)
        h = {}
        for ntype, feat in feat_dict.items():
            if self.use_input_proj and ntype in self.input_proj:
                h[ntype] = self.input_proj[ntype](feat)
            else:
                h[ntype] = feat

        # Keep original features for nodes that may not receive messages
        h_original = {k: v.clone() for k, v in h.items()}

        # Layer 1 — fill missing node types with projected input (for unidirectional graphs)
        h1 = self.conv1(g, h)
        for ntype in h_original:
            if ntype not in h1:
                h1[ntype] = torch.zeros(h_original[ntype].shape[0], self.hidden_feats,
                                         device=h_original[ntype].device)
        h1_out = {}
        for ntype in h1:
            x = h1[ntype]
            # GAT returns [N, num_heads, d_head] — flatten to [N, hidden_feats]
            if x.dim() == 3:
                x = x.flatten(1)
            if x.shape[0] > 1:
                x = self.bn1[ntype](x)
            h1_out[ntype] = self.dropout(F.relu(x))
        h1 = h1_out

        # Layer 2
        h2 = self.conv2(g, h1)
        for ntype in h_original:
            if ntype not in h2:
                h2[ntype] = torch.zeros(h_original[ntype].shape[0], self.out_feats,
                                         device=h_original[ntype].device)
        h2_out = {}
        for ntype in h2:
            x = h2[ntype]
            # GAT returns [N, num_heads, d_head] — flatten
            if x.dim() == 3:
                x = x.flatten(1)
            if x.shape[0] > 1:
                x = self.bn2[ntype](x)
            h2_out[ntype] = x
        h2 = h2_out

        return h2

    def decode(self, etype, h_src, h_dst, edge_graph=None):
        """Score edges of a given type.

        Args:
            etype: Edge type key (e.g., 'mpi', 'mdi').
            h_src: Tensor[E, out_feats] source node embeddings.
            h_dst: Tensor[E, out_feats] destination node embeddings.
            edge_graph: Optional DGL graph for edge-level scoring.
                        If None, builds one from h_src/h_dst.

        Returns:
            Tensor[E] edge scores (logits, not sigmoid).
        """
        n_edges = h_src.shape[0]
        if edge_graph is None:
            edge_graph = dgl.graph(
                (torch.arange(n_edges), torch.arange(n_edges, 2 * n_edges)),
                num_nodes=2 * n_edges,
            )
        h = torch.cat([h_src, h_dst], dim=0)
        return self.decoders[etype](edge_graph, h)

    def forward(self, g, feat_dict, pos_edges, neg_edges):
        """Full forward pass: encode + decode for all edge types.

        Args:
            g: DGL heterogeneous graph.
            feat_dict: {ntype: Tensor[N, in_feats]}.
            pos_edges: {etype: (src_idx, dst_idx)} positive edges.
            neg_edges: {etype: (src_idx, dst_idx)} negative edges.

        Returns:
            pos_scores: {etype: Tensor[E]} positive edge scores.
            neg_scores: {etype: Tensor[E]} negative edge scores.
        """
        embeddings = self.encode(g, feat_dict)

        pos_scores = {}
        neg_scores = {}

        for etype in pos_edges:
            src_type = CANONICAL_ETYPES[etype][0]
            dst_type = CANONICAL_ETYPES[etype][2]

            # Positive edges
            p_src, p_dst = pos_edges[etype]
            h_src_pos = embeddings[src_type][p_src]
            h_dst_pos = embeddings[dst_type][p_dst]
            n_pos = len(p_src)
            pos_g = dgl.graph(
                (torch.arange(n_pos), torch.arange(n_pos, 2 * n_pos)),
                num_nodes=2 * n_pos,
            )
            pos_scores[etype] = self.decoders[etype](pos_g, torch.cat([h_src_pos, h_dst_pos]))

            # Negative edges
            n_src, n_dst = neg_edges[etype]
            h_src_neg = embeddings[src_type][n_src]
            h_dst_neg = embeddings[dst_type][n_dst]
            n_neg = len(n_src)
            neg_g = dgl.graph(
                (torch.arange(n_neg), torch.arange(n_neg, 2 * n_neg)),
                num_nodes=2 * n_neg,
            )
            neg_scores[etype] = self.decoders[etype](neg_g, torch.cat([h_src_neg, h_dst_neg]))

        return pos_scores, neg_scores, embeddings
