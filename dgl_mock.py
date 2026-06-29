"""
Mock DGL module for when DGL is not available
"""

class MockGraph:
    def __init__(self, *args, **kwargs):
        self.num_nodes = 0
        self.ndata = {}
    
    def number_of_nodes(self):
        return self.num_nodes
    
    def to(self, device):
        return self

def from_networkx(graph):
    mock_graph = MockGraph()
    mock_graph.num_nodes = len(graph.nodes())
    return mock_graph

def graph(edges, num_nodes, device=None):
    mock_graph = MockGraph()
    mock_graph.num_nodes = num_nodes
    return mock_graph

# Create a mock dgl module
import sys
import types

mock_dgl = types.ModuleType('dgl')
mock_dgl.from_networkx = from_networkx
mock_dgl.graph = graph

# Add to sys.modules so imports work
sys.modules['dgl'] = mock_dgl
