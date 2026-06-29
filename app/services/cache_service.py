"""
Cache service for optimized data loading and management
"""

import logging
import pickle
import pandas as pd
import numpy as np
import torch
from typing import Any, Dict, Optional, Union
from functools import lru_cache
from pathlib import Path
import threading
import time

from app.config import Config

logger = logging.getLogger(__name__)


class CacheService:
    """Service for caching frequently accessed data and models"""
    
    def __init__(self):
        self.config = Config()
        self._cache: Dict[str, Any] = {}
        self._cache_timestamps: Dict[str, float] = {}
        self._cache_lock = threading.RLock()
        self._cache_ttl = 3600  # 1 hour TTL for cached data
    
    def get_or_load(self, key: str, loader_func, *args, **kwargs) -> Any:
        """
        Get data from cache or load it using the provided function
        
        Args:
            key: Cache key
            loader_func: Function to load data if not in cache
            *args: Arguments for loader function
            **kwargs: Keyword arguments for loader function
            
        Returns:
            Cached or newly loaded data
        """
        with self._cache_lock:
            # Check if data is in cache and not expired
            if self._is_cached(key):
                return self._cache[key]
            
            # Load data using provided function
            data = loader_func(*args, **kwargs)
            
            # Cache the data
            self._cache[key] = data
            self._cache_timestamps[key] = time.time()
            
            return data
    
    def _is_cached(self, key: str) -> bool:
        """Check if data is cached and not expired"""
        if key not in self._cache:
            return False
        
        # Check TTL
        if time.time() - self._cache_timestamps[key] > self._cache_ttl:
            self._remove_from_cache(key)
            return False
        
        return True
    
    def _remove_from_cache(self, key: str):
        """Remove data from cache"""
        if key in self._cache:
            del self._cache[key]
        if key in self._cache_timestamps:
            del self._cache_timestamps[key]
    
    def clear_cache(self):
        """Clear all cached data"""
        with self._cache_lock:
            self._cache.clear()
            self._cache_timestamps.clear()
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        with self._cache_lock:
            return {
                'cache_size': len(self._cache),
                'cached_keys': list(self._cache.keys()),
                'memory_usage': sum(self._estimate_size(v) for v in self._cache.values())
            }
    
    def _estimate_size(self, obj: Any) -> int:
        """Estimate memory size of an object"""
        try:
            if isinstance(obj, (pd.DataFrame, pd.Series)):
                return obj.memory_usage(deep=True).sum()
            elif isinstance(obj, np.ndarray):
                return obj.nbytes
            elif isinstance(obj, torch.Tensor):
                return obj.element_size() * obj.nelement()
            else:
                return len(str(obj))
        except Exception:
            logging.getLogger(__name__).debug("Could not estimate object size", exc_info=True)
            return 0


class OptimizedDataLoader:
    """Optimized data loader with caching and lazy loading"""
    
    def __init__(self):
        self.config = Config()
        self.cache = CacheService()
    
    @lru_cache(maxsize=1)
    def load_metabolite_pca(self) -> Any:
        """Load metabolite PCA model with caching"""
        return self.cache.get_or_load(
            'mets_pca',
            self._load_pickle_file,
            self.config.METS_PCA_PATH
        )
    
    @lru_cache(maxsize=1)
    def load_protein_pca(self) -> Any:
        """Load protein PCA model with caching"""
        return self.cache.get_or_load(
            'protein_pca',
            self._load_pickle_file,
            self.config.PROTEIN_PCA_PATH
        )
    
    @lru_cache(maxsize=1)
    def load_protein_vectors(self) -> Dict[str, np.ndarray]:
        """Load protein vectors with caching"""
        return self.cache.get_or_load(
            'protein_vectors',
            self._load_pickle_file,
            self.config.PROTEIN_VECTOR_PATH
        )
    
    @lru_cache(maxsize=1)
    def load_feature_dataframe(self, organism: str = "All") -> pd.DataFrame:
        """Load feature dataframe with caching"""
        key = f'feature_df_{organism}'
        path = self.config.get_feature_df_path(organism)
        return self.cache.get_or_load(
            key,
            self._load_pickle_file,
            path
        )
    
    @lru_cache(maxsize=1)
    def load_network_graph(self, organism: str = "All") -> Any:
        """Load network graph with caching"""
        key = f'network_graph_{organism}'
        path = self.config.get_network_path(organism)
        return self.cache.get_or_load(
            key,
            self._load_pickle_file,
            path
        )
    
    @lru_cache(maxsize=1)
    def load_mpi_model(self) -> torch.nn.Module:
        """Load MPI model with caching"""
        return self.cache.get_or_load(
            'mpi_model',
            self._load_torch_model,
            self.config.MPI_MODEL_PATH
        )
    
    @lru_cache(maxsize=1)
    def load_mpi_predictor(self) -> torch.nn.Module:
        """Load MPI predictor — handles both full model and state_dict formats"""
        return self.cache.get_or_load(
            'mpi_predictor',
            self._load_predictor_model,
            self.config.MPI_PREDICTOR_PATH
        )
    
    @lru_cache(maxsize=1)
    def load_docking_database(self) -> pd.DataFrame:
        """Load docking database with caching"""
        return self.cache.get_or_load(
            'docking_db',
            self._load_csv_file,
            self.config.DOCKING_DB_PATH
        )
    
    @lru_cache(maxsize=1)
    def load_mpi_database(self) -> pd.DataFrame:
        """Load MPI database with caching"""
        return self.cache.get_or_load(
            'mpi_db',
            self._load_csv_file,
            self.config.MPI_DB_PATH
        )
    
    def _load_pickle_file(self, path: Path) -> Any:
        """Load pickle file with error handling"""
        try:
            with open(path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            logger.error(f"Error loading pickle file {path}: {e}")
            raise
    
    def _load_torch_model(self, path: Path) -> torch.nn.Module:
        """Load PyTorch model with error handling"""
        try:
            return torch.load(path, map_location='cpu', weights_only=False)
        except Exception as e:
            logger.error(f"Error loading torch model {path}: {e}")
            raise

    def _load_predictor_model(self, path: Path) -> torch.nn.Module:
        """Load MLPPredictor — reconstructs from state_dict if needed.

        Handles both legacy 2-layer (W1, W2) and current 3-layer (W1, W2, W3) formats.
        W1.weight.shape[0] == h_feats in both formats.
        """
        from app.models.vgae_model import MLPPredictor
        try:
            obj = torch.load(path, map_location='cpu', weights_only=False)
            if isinstance(obj, torch.nn.Module):
                return obj
            # It's a state dict — infer h_feats from W1.weight shape
            if isinstance(obj, dict) and 'W1.weight' in obj:
                h_feats = obj['W1.weight'].shape[0]
                model = MLPPredictor(h_feats)
                try:
                    model.load_state_dict(obj)  # strict: new 3-layer format
                except RuntimeError:
                    model.load_state_dict(obj, strict=False)  # legacy 2-layer (missing W3)
                model.eval()
                return model
            raise ValueError(f"Unknown predictor format: {type(obj)}")
        except Exception as e:
            logger.error(f"Error loading predictor model {path}: {e}")
            raise
    
    def _load_csv_file(self, path: Path) -> pd.DataFrame:
        """Load CSV file with error handling"""
        try:
            return pd.read_csv(path)
        except Exception as e:
            logger.error(f"Error loading CSV file {path}: {e}")
            raise
    
    def preload_essential_data(self):
        """Preload essential data for faster access"""
        try:
            # Load most commonly used data
            self.load_metabolite_pca()
            self.load_protein_pca()
            self.load_protein_vectors()
            self.load_feature_dataframe("All")
            self.load_network_graph("All")
            logger.info("Essential data preloaded successfully")
        except Exception as e:
            logger.error(f"Error preloading essential data: {e}")
    
    def get_memory_usage(self) -> Dict[str, Any]:
        """Get memory usage statistics"""
        return self.cache.get_cache_stats()
