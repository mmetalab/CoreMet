"""
Export service — CSV, JSON, GraphML, SIF formats.
"""

import json
import io

import pandas as pd
import networkx as nx


def export_csv(predictions_df: pd.DataFrame) -> str:
    """Export predictions as CSV string."""
    return predictions_df.to_csv(index=False)


def export_json(predictions_df: pd.DataFrame) -> str:
    """Export predictions as JSON string."""
    return predictions_df.to_json(orient='records', indent=2)


def export_graphml(predictions_df: pd.DataFrame) -> str:
    """Export predictions as GraphML for Cytoscape."""
    G = nx.Graph()

    for _, row in predictions_df.iterrows():
        met = str(row.get('Metabolite', ''))
        prot = str(row.get('Protein', row.get('Protein Name', '')))
        score = float(row.get('Prediction Score', 0))

        G.add_node(met, node_type='metabolite')
        G.add_node(prot, node_type='protein')
        G.add_edge(met, prot, score=score)

    buf = io.BytesIO()
    nx.write_graphml(G, buf)
    return buf.getvalue().decode('utf-8')


def export_sif(predictions_df: pd.DataFrame) -> str:
    """Export predictions as Simple Interaction Format."""
    lines = []
    for _, row in predictions_df.iterrows():
        met = str(row.get('Metabolite', ''))
        prot = str(row.get('Protein', row.get('Protein Name', '')))
        score = row.get('Prediction Score', 0)
        lines.append(f"{met}\tinteracts_with\t{prot}")
    return '\n'.join(lines)
