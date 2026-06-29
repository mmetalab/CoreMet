"""
VGAE model definitions for CoreMet Application
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class MLPPredictor(nn.Module):
    """
    3-layer MLP predictor for edge prediction in VGAE.
    Architecture: 2h→h→h/2→1 with dropout=0.3
    """

    def __init__(self, h_feats: int):
        super().__init__()
        self.W1 = nn.Linear(h_feats * 2, h_feats)
        self.W2 = nn.Linear(h_feats, h_feats // 2)
        self.W3 = nn.Linear(h_feats // 2, 1)
        self.dropout = nn.Dropout(0.3)

    def apply_edges(self, edges):
        h = torch.cat([edges.src['h'], edges.dst['h']], 1)
        h = self.dropout(F.relu(self.W1(h)))
        h = self.dropout(F.relu(self.W2(h)))
        return {'score': self.W3(h).squeeze(1)}

    def forward(self, g, h):
        with g.local_scope():
            g.ndata['h'] = h
            g.apply_edges(self.apply_edges)
            return g.edata['score']


class GraphSAGEEncoder(nn.Module):
    """
    2-layer GraphSAGE encoder that aggregates neighbor information.
    Uses 'mean' aggregation: h_v = W * CONCAT(h_v, MEAN_{u in N(v)} h_u)
    """

    def __init__(self, in_feats: int, hidden_feats: int, out_feats: int, dropout: float = 0.2):
        super().__init__()
        from dgl.nn import SAGEConv
        self.conv1 = SAGEConv(in_feats, hidden_feats, 'mean')
        self.conv2 = SAGEConv(hidden_feats, out_feats, 'mean')
        self.bn1 = nn.BatchNorm1d(hidden_feats)
        self.dropout = nn.Dropout(dropout)

    def forward(self, g, x):
        x = F.relu(self.bn1(self.conv1(g, x)))
        x = self.dropout(x)
        x = self.conv2(g, x)
        return x


class VGAEEncoder(nn.Module):
    """
    Legacy MLP-based encoder (kept for backward compatibility).
    For new training use GraphSAGEEncoder.
    """

    def __init__(self, in_feats: int, hidden_feats: int, out_feats: int):
        super().__init__()
        self.conv1 = nn.Linear(in_feats, hidden_feats)
        self.conv2 = nn.Linear(hidden_feats, out_feats)
        self.dropout = nn.Dropout(0.1)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = self.dropout(x)
        x = self.conv2(x)
        return x


class VGAEDecoder(nn.Module):
    """
    Variational Graph Autoencoder Decoder
    """

    def __init__(self, in_feats: int, hidden_feats: int, out_feats: int):
        super().__init__()
        self.conv1 = nn.Linear(in_feats, hidden_feats)
        self.conv2 = nn.Linear(hidden_feats, out_feats)
        self.dropout = nn.Dropout(0.1)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = self.dropout(x)
        x = torch.sigmoid(self.conv2(x))
        return x


class VGAE(nn.Module):
    """
    Complete Variational Graph Autoencoder
    """

    def __init__(self, in_feats: int, hidden_feats: int, out_feats: int):
        super().__init__()
        self.encoder = GraphSAGEEncoder(in_feats, hidden_feats, out_feats)
        self.decoder = VGAEDecoder(out_feats, hidden_feats, in_feats)
        self.predictor = MLPPredictor(out_feats)

    def forward(self, g, features):
        encoded = self.encoder(g, features)
        decoded = self.decoder(encoded)
        edge_scores = self.predictor(g, encoded)
        return decoded, edge_scores, encoded
