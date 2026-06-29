"""
Optimized prediction service for CoreMet Application
"""

import logging
import numpy as np
import pandas as pd
import torch
import networkx as nx
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import gc
import time

logger = logging.getLogger(__name__)

# Handle DGL import gracefully
try:
    import dgl
    DGL_AVAILABLE = True
except ImportError as e:
    logger.warning(f"DGL not available, using mock: {e}")
    # Import mock DGL
    import sys
    import os
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    import dgl_mock
    import dgl
    DGL_AVAILABLE = False

from app.config import Config
from app.models.vgae_model import MLPPredictor
from app.utils.data_utils import cosine_similarity
from app.services.cache_service import OptimizedDataLoader


class PredictionService:
    """Optimized service for metabolite-protein interaction prediction"""
    
    def __init__(self):
        self.config = Config()
        self.data_loader = OptimizedDataLoader()
        self.model = None
        self.predictor = None
        self.node_features = None
        self.graph = None
        self.dgl_graph = None
        self._device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self._load_models()
    
    def _load_models(self):
        """Load pre-trained models and data with optimization"""
        try:
            start_time = time.time()
            
            # Load models with caching
            self.model = self.data_loader.load_mpi_model()
            self.predictor = self.data_loader.load_mpi_predictor()
            
            # Load network data with caching
            self.graph = self.data_loader.load_network_graph("All")
            self.node_features = self.data_loader.load_feature_dataframe("All")
            
            # Create DGL graph efficiently
            self._create_dgl_graph()
            
            load_time = time.time() - start_time
            logger.info(f"Models loaded successfully in {load_time:.2f} seconds")

        except Exception as e:
            logger.error(f"Error loading models: {e}")
            raise
    
    def _create_dgl_graph(self):
        """Create DGL graph efficiently"""
        if not DGL_AVAILABLE:
            logger.warning("DGL not available, skipping graph creation")
            self.dgl_graph = None
            return
            
        try:
            # Create DGL graph
            self.dgl_graph = dgl.from_networkx(self.graph)
            
            # Optimize feature loading
            features = np.stack(self.node_features['pca_128'].values, axis=0)
            features = torch.from_numpy(features).float()
            
            # Move to device if available
            if torch.cuda.is_available():
                features = features.to(self._device)
                self.dgl_graph = self.dgl_graph.to(self._device)
            
            self.dgl_graph.ndata['h'] = features
            
        except Exception as e:
            logger.error(f"Error creating DGL graph: {e}")
            self.dgl_graph = None
    
    def predict_interactions(self, 
                           metabolites_df: pd.DataFrame, 
                           proteins_df: pd.DataFrame,
                           organism: str = "All") -> pd.DataFrame:
        """
        Predict metabolite-protein interactions with optimization
        
        Args:
            metabolites_df: DataFrame with metabolite data
            proteins_df: DataFrame with protein data  
            organism: Target organism for prediction
            
        Returns:
            DataFrame with prediction results
        """
        try:
            start_time = time.time()
            
            # Load organism-specific data if needed
            if organism != "All":
                self._load_organism_data(organism)
            
            # Process input data efficiently
            feature_df = self._process_input_data_optimized(metabolites_df, proteins_df)
            
            # Generate predictions with memory optimization
            results = self._generate_predictions_optimized(feature_df, proteins_df)
            
            # Clean up memory
            gc.collect()
            
            prediction_time = time.time() - start_time
            logger.info(f"Prediction completed in {prediction_time:.2f} seconds")

            return results

        except Exception as e:
            logger.error(f"Error in prediction: {e}")
            # Clean up on error
            gc.collect()
            raise
    
    def _load_organism_data(self, organism: str):
        """Load organism-specific network and features with caching"""
        try:
            # Load organism-specific data with caching
            self.graph = self.data_loader.load_network_graph(organism)
            self.node_features = self.data_loader.load_feature_dataframe(organism)
            
            # Recreate DGL graph for new organism
            self._create_dgl_graph()
                
        except Exception as e:
            logger.warning(f"Could not load organism-specific data: {e}")
    
    def _process_input_data_optimized(self, 
                                    metabolites_df: pd.DataFrame, 
                                    proteins_df: pd.DataFrame) -> pd.DataFrame:
        """Process input data with optimization"""
        from app.services.molecular_service import MolecularService
        
        molecular_service = MolecularService()
        feature_data = []
        
        # Process metabolites in batches for better memory usage
        batch_size = 50
        for i in range(0, len(metabolites_df), batch_size):
            batch = metabolites_df.iloc[i:i+batch_size]
            for _, row in batch.iterrows():
                features = molecular_service.extract_metabolite_features(row['SMILES'])
                if features is not None:
                    feature_data.append({
                        'Name': row['Metabolite Name'],
                        'Type': 'Metabolite',
                        'HMDB ID': row['HMDB ID'],
                        'Feature': features
                    })
        
        # Process proteins in batches
        for i in range(0, len(proteins_df), batch_size):
            batch = proteins_df.iloc[i:i+batch_size]
            for _, row in batch.iterrows():
                features = molecular_service.extract_protein_features(row['UniprotID'])
                if features is not None:
                    feature_data.append({
                        'Name': row['UniprotID'],
                        'Type': 'Protein',
                        'Feature': features
                    })
        
        return pd.DataFrame(feature_data)
    
    def _generate_predictions_optimized(self, 
                                      feature_df: pd.DataFrame, 
                                      proteins_df: pd.DataFrame) -> pd.DataFrame:
        """Generate predictions with memory optimization"""
        
        # Get adjacency matrix efficiently
        adj_true = nx.adjacency_matrix(self.graph, nodelist=self.node_features['node'].tolist())
        
        # Extract metabolites and proteins
        metabolites = feature_df[feature_df['Type'] == 'Metabolite']['HMDB ID'].tolist()
        proteins_df_filtered = feature_df[feature_df['Type'] == 'Protein']
        proteins = proteins_df_filtered['Name'].tolist()
        
        # Map proteins to network nodes efficiently
        protein_mapping = self._map_proteins_to_nodes_optimized(proteins_df_filtered)
        mapped_proteins = [mapping[0] for mapping in protein_mapping.values()]
        
        # Generate edge pairs efficiently
        node_dict, test_pairs, true_pairs = self._generate_edge_pairs_optimized(
            metabolites, mapped_proteins, adj_true
        )
        
        if len(test_pairs) == 0:
            return pd.DataFrame(columns=['Metabolite', 'Protein', 'Prediction Score', 'Existing'])
        
        # Make predictions with GPU optimization
        results = self._make_predictions_gpu_optimized(test_pairs, true_pairs)
        
        # Merge with protein information efficiently
        results = self._merge_protein_info_optimized(results, protein_mapping, proteins_df)
        
        # Annotate with enzyme information from MEI database
        try:
            from app.services.mei_service import annotate_predictions_with_mei
            results = annotate_predictions_with_mei(results, uniprot_col="Protein")
        except Exception as e:
            logger.warning(f"MEI annotation skipped: {e}")
        
        return results
    
    def _map_proteins_to_nodes_optimized(self, proteins_df: pd.DataFrame) -> Dict:
        """Map proteins to closest network nodes using fully vectorized cosine similarity."""
        result_dict = dict(zip(proteins_df['Name'].tolist(), proteins_df['Feature'].tolist()))

        # Build node matrix once (N_nodes × dim)
        node_ids = list(self.node_features.dbid)
        node_matrix = np.array(list(self.node_features.pca_128))     # (N_nodes, dim)
        # Pre-compute node norms
        node_norms = np.linalg.norm(node_matrix, axis=1, keepdims=True)
        node_norms[node_norms == 0] = 1.0  # avoid div-by-zero
        node_matrix_normed = node_matrix / node_norms

        mapping_dict = {}
        for protein_id, protein_features in result_dict.items():
            protein_vec = np.array(protein_features).reshape(1, -1)
            p_norm = np.linalg.norm(protein_vec)
            if p_norm == 0:
                mapping_dict[protein_id] = [node_ids[0], 0.0]
                continue
            protein_vec_normed = protein_vec / p_norm
            # Vectorized cosine similarity against ALL nodes at once
            similarities = (protein_vec_normed @ node_matrix_normed.T).ravel()
            best_idx = int(np.argmax(similarities))
            mapping_dict[protein_id] = [node_ids[best_idx], float(similarities[best_idx])]

        return mapping_dict
    
    def _generate_edge_pairs_optimized(self, 
                                     metabolites: List[str], 
                                     proteins: List[str], 
                                     adj_true) -> Tuple[Dict, np.ndarray, List]:
        """Generate edge pairs with optimization"""
        node_dict = dict(zip(list(self.node_features.dbid), list(self.node_features.index)))
        
        # Use set operations for faster lookup
        metabolite_indices = [node_dict[met] for met in metabolites if met in node_dict]
        protein_indices = [node_dict[prot] for prot in proteins if prot in node_dict]
        
        # Generate test pairs efficiently
        test_pairs = []
        true_pairs = []
        
        for met_idx in metabolite_indices:
            for prot_idx in protein_indices:
                test_pairs.append([met_idx, prot_idx])
                true_pairs.append(adj_true[met_idx, prot_idx])
        
        return node_dict, np.array(test_pairs), true_pairs
    
    def _make_predictions_gpu_optimized(self, test_pairs: np.ndarray, true_pairs: List) -> pd.DataFrame:
        """Make predictions with GPU optimization and OOM fallback to CPU."""
        if not DGL_AVAILABLE or self.dgl_graph is None:
            logger.warning("DGL not available, returning empty results")
            return pd.DataFrame(columns=['Metabolite', 'Protein', 'Prediction Score', 'Existing'])

        # Build node feature tensor from the current (organism-specific) feature DataFrame
        node_feat_array = np.stack(self.node_features['pca_128'].values, axis=0)
        h = torch.from_numpy(node_feat_array).float()
        n_nodes = self.dgl_graph.number_of_nodes()

        # Try GPU first, fall back to CPU on OOM
        device = self._device
        try:
            real_apply_g = dgl.graph((test_pairs[:, 0], test_pairs[:, 1]),
                                    num_nodes=n_nodes, device=device)
            if device.type == "cuda":
                h = h.to(device)
                self.predictor = self.predictor.to(device)

            with torch.no_grad():
                pos_scores = self.predictor(real_apply_g, h)

        except (torch.cuda.OutOfMemoryError, RuntimeError) as e:
            if "out of memory" in str(e).lower() or isinstance(e, torch.cuda.OutOfMemoryError):
                logger.warning("GPU OOM during prediction — falling back to CPU")
                torch.cuda.empty_cache()
                gc.collect()
                device = torch.device("cpu")
                real_apply_g = dgl.graph((test_pairs[:, 0], test_pairs[:, 1]),
                                        num_nodes=n_nodes, device=device)
                h = h.to(device)
                self.predictor = self.predictor.to(device)
                with torch.no_grad():
                    pos_scores = self.predictor(real_apply_g, h)
            else:
                raise

        # Format results efficiently
        results = []
        for i in range(len(test_pairs)):
            met_idx, prot_idx = test_pairs[i][0], test_pairs[i][1]
            metabolite = self.node_features.loc[met_idx]['node']
            protein = self.node_features.loc[prot_idx]['node']
            score = f"{torch.sigmoid(pos_scores[i]).item():.5f}"
            existing = "Yes" if true_pairs[i] == 1 else "No"

            results.append({
                'Metabolite': metabolite,
                'Protein': protein,
                'Prediction Score': score,
                'Existing': existing
            })

        return pd.DataFrame(results)
    
    def _merge_protein_info_optimized(self, 
                                    results: pd.DataFrame, 
                                    protein_mapping: Dict, 
                                    proteins_df: pd.DataFrame) -> pd.DataFrame:
        """Merge protein information efficiently"""
        # Create mapping dataframe
        mapping_df = pd.DataFrame([
            {'source_protein': k, 'dbid': v[0], 'similarity_score': v[1]}
            for k, v in protein_mapping.items()
        ])
        
        # Merge results with mapping using efficient pandas operations
        merged_results = results.merge(self.node_features, 
                                     left_on='Protein', right_on='node', how='left')
        merged_results = merged_results.merge(mapping_df, 
                                            left_on='dbid', right_on='dbid', how='left')
        merged_results = merged_results.merge(proteins_df, 
                                            left_on='source_protein', right_on='UniprotID', how='left')
        
        # Select final columns
        final_columns = ['Metabolite', 'Protein Name', 'Prediction Score', 
                        'Existing', 'similarity_score', 'Protein']
        return merged_results[final_columns]
    
