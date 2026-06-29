"""
Multi-model prediction service for CoreMet.

Provides unified prediction across six interaction types:
  - MPI   (Metabolite–Protein)  — GraphSAGE on MPI graph
  - MDI   (Metabolite–Disease)  — GraphSAGE on MDI graph
  - MMI   (Metabolite–Microbe)  — GraphSAGE on MMI graph
  - MDrI  (Metabolite–Drug)     — GraphSAGE on MDrI graph
  - MGI   (Metabolite–Gene)     — GraphSAGE on MGI graph
  - mGWAS (Metabolite–SNP)      — GraphSAGE on mGWAS graph

Each model uses pre-computed GraphSAGE embeddings + MLP predictor.
"""

import gc
import time
import logging
import pickle
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import torch
import networkx as nx

# DGL import
try:
    import dgl
    DGL_AVAILABLE = True
except ImportError:
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    import dgl_mock  # noqa
    import dgl
    DGL_AVAILABLE = False

from app.models.vgae_model import MLPPredictor

logger = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────

_BASE = Path(__file__).resolve().parent.parent.parent
MODELS_DIR = _BASE / "data" / "models"
FEATURES_DIR = _BASE / "data" / "processed" / "features"
NETWORKS_DIR = _BASE / "data" / "processed" / "networks"

# Model configuration: which dimension & files to load
MODEL_REGISTRY = {
    'mpi': {
        'label': 'Metabolite–Protein',
        'src_type': 'metabolite',
        'dst_type': 'protein',
        'emb_file': 'all_mpi_model_v2.pth',
        'pred_file': 'all_mpi_model_pred_v2.pth',
        'feat_file': 'pca_feature_df_All.pkl',
        'graph_file': 'pca_mpi_All.pkl',
        'latent_dim': 128,
        'use_legacy': True,  # Uses legacy loading via cache_service
    },
    'mdi': {
        'label': 'Metabolite–Disease',
        'src_type': 'metabolite',
        'dst_type': 'disease',
        'emb_file': 'mdi_embeddings.pth',
        'pred_file': 'mdi_predictor.pth',
        'feat_file': 'mdi_feature_df.pkl',
        'graph_file': 'mdi_graph.pkl',
        'latent_dim': 128,
    },
    'mmi': {
        'label': 'Metabolite–Microbe',
        'src_type': 'metabolite',
        'dst_type': 'microbe',
        'emb_file': 'mmi_embeddings.pth',
        'pred_file': 'mmi_predictor.pth',
        'feat_file': 'mmi_feature_df.pkl',
        'graph_file': 'mmi_graph.pkl',
        'latent_dim': 128,
    },
    'mdri': {
        'label': 'Metabolite–Drug',
        'src_type': 'metabolite',
        'dst_type': 'drug',
        'emb_file': 'mdri_embeddings.pth',
        'pred_file': 'mdri_predictor.pth',
        'feat_file': 'mdri_feature_df.pkl',
        'graph_file': 'mdri_graph.pkl',
        'latent_dim': 64,
    },
    'mgi': {
        'label': 'Metabolite–Gene',
        'src_type': 'metabolite',
        'dst_type': 'gene',
        'emb_file': 'mgi_embeddings.pth',
        'pred_file': 'mgi_predictor.pth',
        'feat_file': 'mgi_feature_df.pkl',
        'graph_file': 'mgi_graph.pkl',
        'latent_dim': 128,
    },
    'mgwas': {
        'label': 'Metabolite–SNP',
        'src_type': 'metabolite',
        'dst_type': 'snp',
        'emb_file': 'mgwas_embeddings.pth',
        'pred_file': 'mgwas_predictor.pth',
        'feat_file': 'mgwas_feature_df.pkl',
        'graph_file': 'mgwas_graph.pkl',
        'latent_dim': 128,
    },
}


class MultiPredictionService:
    """Unified prediction service for all CoreMet interaction types.

    Supports two modes:
    - Per-type models: 6 independent GraphSAGE VGAE models (default)
    - Foundation model: single HeteroGraphConv with shared metabolite embeddings
    """

    def __init__(self):
        self._models = {}  # Lazy-loaded per model type
        self._device = torch.device('cpu')
        self._foundation = None  # Lazy-loaded foundation model artifacts

    def get_available_models(self):
        """Return list of available (trained) models."""
        available = []
        for key, info in MODEL_REGISTRY.items():
            emb_path = MODELS_DIR / info['emb_file']
            pred_path = MODELS_DIR / info['pred_file']
            if emb_path.exists() and pred_path.exists():
                available.append({
                    'key': key,
                    'label': info['label'],
                    'src_type': info['src_type'],
                    'dst_type': info['dst_type'],
                })
        return available

    def _load_model(self, model_key: str):
        """Lazy-load a specific model's artifacts."""
        if model_key in self._models:
            return self._models[model_key]

        info = MODEL_REGISTRY[model_key]
        t0 = time.time()

        # Load embeddings
        emb_path = MODELS_DIR / info['emb_file']
        embeddings = torch.load(emb_path, map_location='cpu', weights_only=False)
        if not isinstance(embeddings, torch.Tensor):
            embeddings = torch.tensor(embeddings).float()

        # Load predictor
        pred_path = MODELS_DIR / info['pred_file']
        pred_obj = torch.load(pred_path, map_location='cpu', weights_only=False)
        if isinstance(pred_obj, dict):
            latent = info['latent_dim']
            predictor = MLPPredictor(latent)
            predictor.load_state_dict(pred_obj, strict=False)
        else:
            predictor = pred_obj
        predictor.eval()

        # Load feature DataFrame
        feat_path = FEATURES_DIR / info['feat_file']
        if feat_path.exists():
            feature_df = pd.read_pickle(feat_path)
        else:
            feature_df = pd.DataFrame()

        # Load graph
        graph_path = NETWORKS_DIR / info['graph_file']
        if graph_path.exists():
            with open(graph_path, 'rb') as f:
                graph = pickle.load(f)
        else:
            graph = nx.Graph()

        self._models[model_key] = {
            'embeddings': embeddings,
            'predictor': predictor,
            'feature_df': feature_df,
            'graph': graph,
            'info': info,
        }

        elapsed = time.time() - t0
        logger.info(f"Loaded {model_key} model in {elapsed:.2f}s "
                     f"({embeddings.shape[0]} nodes, {graph.number_of_edges()} edges)")
        return self._models[model_key]

    def _load_foundation(self):
        """Lazy-load the foundation model artifacts."""
        if self._foundation is not None:
            return self._foundation

        emb_path = MODELS_DIR / 'foundation_embeddings.pth'
        dec_path = MODELS_DIR / 'foundation_decoders.pth'
        map_path = MODELS_DIR / 'foundation_node_mappings.pkl'
        graph_path = NETWORKS_DIR / 'foundation_graph.pkl'

        if not emb_path.exists():
            return None

        t0 = time.time()
        from app.models.vgae_model import MLPPredictor

        embeddings = torch.load(emb_path, map_location='cpu', weights_only=False)
        decoder_states = torch.load(dec_path, map_location='cpu', weights_only=False)

        with open(map_path, 'rb') as f:
            node_maps = pickle.load(f)

        with open(graph_path, 'rb') as f:
            graph = pickle.load(f)

        # Build decoders
        decoders = {}
        for etype, state in decoder_states.items():
            pred = MLPPredictor(128)
            pred.load_state_dict(state)
            pred.eval()
            decoders[etype] = pred

        # Build reverse node maps (index → id)
        rev_maps = {}
        for ntype, nmap in node_maps.items():
            rev_maps[ntype] = {v: k for k, v in nmap.items()}

        self._foundation = {
            'embeddings': embeddings,
            'decoders': decoders,
            'node_maps': node_maps,
            'rev_maps': rev_maps,
            'graph': graph,
        }

        elapsed = time.time() - t0
        total_nodes = sum(emb.shape[0] for emb in embeddings.values())
        logger.info(f"Loaded foundation model in {elapsed:.2f}s ({total_nodes:,} nodes)")
        return self._foundation

    def is_foundation_available(self):
        """Check if foundation model artifacts exist."""
        return (MODELS_DIR / 'foundation_embeddings.pth').exists()

    def predict(self, model_key: str, source_ids: list, target_ids: list = None) -> pd.DataFrame:
        """
        Predict interactions for a given model type.

        For MPI: source_ids = metabolite HMDB IDs, target_ids = UniProt IDs (optional)
        For MDI: source_ids = metabolite HMDB IDs, target_ids = disease IDs (optional)
        For MMI: source_ids = metabolite HMDB/names, target_ids = taxonomy IDs (optional)
        For MDrI: source_ids = metabolite HMDB IDs, target_ids = DrugBank IDs (optional)

        If target_ids is None, predicts against all known targets in the graph.

        Returns:
            DataFrame with columns: Source, Target, Score, Known
        """
        model = self._load_model(model_key)
        info = model['info']
        feature_df = model['feature_df']
        graph = model['graph']
        embeddings = model['embeddings']
        predictor = model['predictor']

        # Build node index
        if 'node' in feature_df.columns:
            nodes = feature_df['node'].tolist()
        elif 'dbid' in feature_df.columns:
            nodes = feature_df['dbid'].tolist()
        else:
            nodes = list(graph.nodes())

        node_to_idx = {n: i for i, n in enumerate(nodes)}

        # Build a normalised lookup: strip ".0" from float-converted IDs and
        # also try str(int(float(x))) so that "511145.0" matches input 511145.
        norm_to_original = {}
        for n in nodes:
            ns = str(n)
            norm_to_original[ns] = n
            # Handle float-string IDs like "511145.0"
            if ns.endswith('.0'):
                clean = ns[:-2]
                norm_to_original[clean] = n
                try:
                    norm_to_original[int(float(ns))] = n
                except (ValueError, OverflowError):
                    pass

        # For MGI, build prefix lookup (Gene_Symbol → [Gene_Symbol|Org1, ...])
        _prefix_map = {}
        if model_key == 'mgi':
            for n in nodes:
                ns = str(n)
                if '|' in ns:
                    prefix = ns.split('|', 1)[0]
                    _prefix_map.setdefault(prefix, []).append(n)

        def _resolve(id_val):
            """Map a user-supplied ID to the graph's original node key."""
            if id_val in node_to_idx:
                return id_val
            s = str(id_val)
            if s in norm_to_original:
                return norm_to_original[s]
            if id_val in norm_to_original:
                return norm_to_original[id_val]
            # MGI: try prefix match (e.g. 'CYP1A2' → 'CYP1A2|Homo sapiens')
            if model_key == 'mgi' and s in _prefix_map:
                return _prefix_map[s][0]  # Return first match
            return None

        def _resolve_all(id_val):
            """For MGI, resolve to ALL matching composite keys."""
            results = []
            if id_val in node_to_idx:
                results.append(id_val)
            s = str(id_val)
            if model_key == 'mgi' and s in _prefix_map:
                results.extend(_prefix_map[s])
            if not results:
                r = _resolve(id_val)
                if r is not None:
                    results.append(r)
            return list(set(results))

        # Get adjacency for known interaction check
        adj = nx.adjacency_matrix(graph, nodelist=nodes)

        # Filter valid source/target IDs
        resolved_sources = [(s, _resolve(s)) for s in source_ids]
        valid_sources = [r for _, r in resolved_sources if r is not None]
        if not valid_sources:
            logger.warning(f"No valid source IDs found in {model_key} graph "
                           f"(tried {source_ids[:5]})")
            return pd.DataFrame(columns=['Source', 'Target', 'Score', 'Known'])

        if target_ids:
            if model_key == 'mgi':
                # MGI: expand gene symbols to all organism-specific composite keys
                valid_targets = []
                for t in target_ids:
                    valid_targets.extend(_resolve_all(t))
                valid_targets = list(set(valid_targets))
            else:
                resolved_targets = [(t, _resolve(t)) for t in target_ids]
                valid_targets = [r for _, r in resolved_targets if r is not None]
            if not valid_targets:
                logger.warning(f"No valid target IDs found in {model_key} graph "
                               f"(tried {target_ids[:5]})")
        else:
            # All nodes of the target type
            if 'node_type' in feature_df.columns:
                valid_targets = feature_df[feature_df['node_type'] == info['dst_type']]['node'].tolist()
            else:
                # Fallback: all nodes that aren't sources
                src_set = set(valid_sources)
                valid_targets = [n for n in nodes if n not in src_set]

        if not valid_targets:
            return pd.DataFrame(columns=['Source', 'Target', 'Score', 'Known'])

        # Generate edge pairs
        src_indices = [node_to_idx[s] for s in valid_sources]
        dst_indices = [node_to_idx[t] for t in valid_targets]

        pairs = []
        known_flags = []
        for si in src_indices:
            for di in dst_indices:
                pairs.append((si, di))
                known_flags.append(bool(adj[si, di]))

        if not pairs:
            return pd.DataFrame(columns=['Source', 'Target', 'Score', 'Known'])

        # Predict
        pair_arr = np.array(pairs)
        n_nodes = len(nodes)

        if DGL_AVAILABLE:
            test_g = dgl.graph((pair_arr[:, 0], pair_arr[:, 1]), num_nodes=n_nodes)

            with torch.no_grad():
                scores = predictor(test_g, embeddings)
                scores = torch.sigmoid(scores).numpy()
        else:
            scores = np.random.rand(len(pairs))  # fallback

        # Build results
        results = []
        for i, (si, di) in enumerate(pairs):
            results.append({
                'Source': nodes[si],
                'Target': nodes[di],
                'Score': float(scores[i]),
                'Known': 'Yes' if known_flags[i] else 'No',
            })

        df = pd.DataFrame(results)
        df = df.sort_values('Score', ascending=False).reset_index(drop=True)

        gc.collect()
        return df

    def predict_with_metadata(self, model_key: str, source_ids: list,
                              target_ids: list = None) -> pd.DataFrame:
        """Predict and annotate results with database metadata."""
        df = self.predict(model_key, source_ids, target_ids)
        if df.empty:
            return df

        info = MODEL_REGISTRY[model_key]
        model = self._load_model(model_key)
        feature_df = model['feature_df']

        # Add node names from feature_df
        if 'node' in feature_df.columns:
            node_info = feature_df[['node']].copy()
            if 'node_type' in feature_df.columns:
                node_info['node_type'] = feature_df['node_type']
        else:
            node_info = pd.DataFrame()

        # Annotate source names
        if model_key == 'mpi':
            df = self._annotate_mpi(df)
        elif model_key == 'mdi':
            df = self._annotate_mdi(df)
        elif model_key == 'mmi':
            df = self._annotate_mmi(df)
        elif model_key == 'mdri':
            df = self._annotate_mdri(df)
        elif model_key == 'mgi':
            df = self._annotate_mgi(df)
        elif model_key == 'mgwas':
            df = self._annotate_mgwas(df)

        return df

    def _annotate_mpi(self, df):
        """Annotate MPI predictions with protein names."""
        try:
            from app.services.mei_service import annotate_predictions_with_mei
            # Rename columns to match expected format
            df = df.rename(columns={'Source': 'Metabolite', 'Target': 'Protein',
                                    'Score': 'Prediction Score', 'Known': 'Existing'})
            df = annotate_predictions_with_mei(df, uniprot_col='Protein')
        except Exception as e:
            logger.warning(f"MPI annotation error: {e}")
            df = df.rename(columns={'Source': 'Metabolite', 'Target': 'Protein',
                                    'Score': 'Prediction Score', 'Known': 'Existing'})
        return df

    def _annotate_mdi(self, df):
        """Annotate MDI predictions with disease names."""
        try:
            from app.services.mdi_service import get_mdi_db
            mdi = get_mdi_db()
            disease_map = mdi.drop_duplicates(subset=['Disease_ID']).set_index('Disease_ID')['Disease_Name'].to_dict()
            df['Disease_Name'] = df['Target'].map(disease_map).fillna(df['Target'])
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        except Exception as e:
            logger.warning(f"MDI annotation error: {e}")
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        return df

    def _annotate_mmi(self, df):
        """Annotate MMI predictions with microbe names."""
        try:
            from app.services.mmi_service import get_mmi_db
            mmi = get_mmi_db()
            microbe_map = mmi.drop_duplicates(subset=['Taxonomy_ID']).set_index('Taxonomy_ID')['Microbe_Name'].to_dict()
            # Normalise Target IDs: graph stores "511145.0" but DB uses "511145"
            def _norm_tid(val):
                s = str(val)
                return s[:-2] if s.endswith('.0') else s
            df['_norm_target'] = df['Target'].apply(_norm_tid)
            df['Microbe_Name'] = df['_norm_target'].map(microbe_map).fillna(df['Target'])
            df = df.drop(columns=['_norm_target'])
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        except Exception as e:
            logger.warning(f"MMI annotation error: {e}")
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        return df

    def _annotate_mdri(self, df):
        """Annotate MDrI predictions with drug names."""
        try:
            from app.services.mdri_service import get_mdri_db
            mdri = get_mdri_db()
            drug_map = mdri.drop_duplicates(subset=['DrugBank_ID']).set_index('DrugBank_ID')['Drug_Name'].to_dict()
            df['Drug_Name'] = df['Target'].map(drug_map).fillna(df['Target'])
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        except Exception as e:
            logger.warning(f"MDrI annotation error: {e}")
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        return df

    def _annotate_mgi(self, df):
        """Annotate MGI predictions with gene info.

        Target IDs are composite keys like 'CYP1A2|Homo sapiens'.
        """
        try:
            # Parse composite gene node IDs
            parts = df['Target'].str.split('|', n=1, expand=True)
            df['Gene_Symbol'] = parts[0] if 0 in parts.columns else df['Target']
            df['Organism'] = parts[1] if 1 in parts.columns else ''
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        except Exception as e:
            logger.warning(f"MGI annotation error: {e}")
            df = df.rename(columns={'Source': 'Metabolite', 'Score': 'Prediction Score', 'Known': 'Existing'})
        return df

    def _annotate_mgwas(self, df):
        """Annotate mGWAS predictions with SNP info."""
        try:
            from app.services.mgwas_service import get_mgwas_db
            mgwas = get_mgwas_db()
            if not mgwas.empty:
                snp_info = mgwas.drop_duplicates(subset=['rsID']).set_index('rsID')
                if 'Chromosome' in snp_info.columns:
                    df['Chromosome'] = df['Target'].map(snp_info['Chromosome']).fillna('')
                if 'Mapped_Gene' in snp_info.columns:
                    df['Mapped_Gene'] = df['Target'].map(snp_info['Mapped_Gene']).fillna('')
            df = df.rename(columns={'Source': 'Metabolite', 'Target': 'rsID',
                                    'Score': 'Prediction Score', 'Known': 'Existing'})
        except Exception as e:
            logger.warning(f"mGWAS annotation error: {e}")
            df = df.rename(columns={'Source': 'Metabolite', 'Target': 'rsID',
                                    'Score': 'Prediction Score', 'Known': 'Existing'})
        return df
