"""
Configuration management for CoreMet Web Application
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional
from functools import lru_cache

class Config:
    """Base configuration class with optimized path resolution"""
    
    # Application settings
    APP_NAME = "CoreMet"
    APP_VERSION = "3.0.0"
    DEBUG = False
    
    def __init__(self):
        """Initialize configuration with lazy path resolution"""
        self._base_dir = Path(__file__).parent.parent
        self._data_dir = self._base_dir / "data"
        self._path_cache: Dict[str, Path] = {}
    
    @property
    def BASE_DIR(self) -> Path:
        """Get base directory"""
        return self._base_dir
    
    @property
    def DATA_DIR(self) -> Path:
        """Get data directory"""
        return self._data_dir
    
    @property
    def MODELS_DIR(self) -> Path:
        """Get models directory"""
        return self._get_cached_path("models", self._data_dir / "models")
    
    @property
    def FEATURES_DIR(self) -> Path:
        """Get features directory"""
        return self._get_cached_path("features", self._data_dir / "processed" / "features")
    
    @property
    def NETWORKS_DIR(self) -> Path:
        """Get networks directory"""
        return self._get_cached_path("networks", self._data_dir / "processed" / "networks")
    
    @property
    def DATABASE_DIR(self) -> Path:
        """Get database directory, check both data/mpidatabase and data/raw/mpidatabase"""
        primary = self._data_dir / "mpidatabase"
        fallback = self._data_dir / "raw" / "mpidatabase"
        # Prefer the directory that has MPIDB files
        if (primary / "MPIDB_v2.csv").exists() or (primary / "MPIDB_May2024.csv").exists():
            return self._get_cached_path("database", primary)
        return self._get_cached_path("database", fallback)
    
    @property
    def DOCKING_DB_PATH(self) -> Path:
        """Get docking database path"""
        return self._get_cached_path("docking_db", self.DATABASE_DIR / "Docking_DB.csv")
    
    @property
    def MPI_DB_PATH(self) -> Path:
        """Get MPI database path. Prefer the curated, backfilled release file (the same
        38,061-edge table the browse modules and manuscript use) so the API, data_service,
        and sitemap all serve identical, fully-annotated MPI records; then v3 > v2 > original."""
        release_path = self._data_dir / "databases" / "release" / "coremetdb_mpi.csv"
        if release_path.exists():
            return self._get_cached_path("mpi_db_release", release_path)
        v3_path = self.DATABASE_DIR / "MPIDB_v3.csv"
        if v3_path.exists():
            return self._get_cached_path("mpi_db_v3", v3_path)
        v2_path = self.DATABASE_DIR / "MPIDB_v2.csv"
        if v2_path.exists():
            return self._get_cached_path("mpi_db_v2", v2_path)
        return self._get_cached_path("mpi_db", self.DATABASE_DIR / "MPIDB_May2024.csv")
    
    @property
    def METS_PCA_PATH(self) -> Path:
        """Get metabolite PCA path"""
        return self._get_cached_path("mets_pca", self.FEATURES_DIR / "mets_pca.pkl")
    
    @property
    def PROTEIN_PCA_PATH(self) -> Path:
        """Get protein PCA path"""
        return self._get_cached_path("protein_pca", self.FEATURES_DIR / "protein_pca.pkl")
    
    @property
    def PROTEIN_VECTOR_PATH(self) -> Path:
        """Get protein vector path"""
        return self._get_cached_path("protein_vector", self.FEATURES_DIR / "protein_vector.p")
    
    @property
    def FEATURE_DF_PATH(self) -> Path:
        """Get feature dataframe path"""
        return self._get_cached_path("feature_df", self.FEATURES_DIR / "pca_feature_df_All.pkl")
    
    @property
    def MPI_NETWORK_PATH(self) -> Path:
        """Get MPI network path"""
        return self._get_cached_path("mpi_network", self.NETWORKS_DIR / "pca_mpi_All.pkl")
    
    @property
    def MPI_MODEL_PATH(self) -> Path:
        """Get MPI model path, prefer v2 if available"""
        v2_path = self.MODELS_DIR / "all_mpi_model_v2.pth"
        if v2_path.exists():
            return self._get_cached_path("mpi_model_v2", v2_path)
        return self._get_cached_path("mpi_model", self.MODELS_DIR / "all_mpi_model.pth")

    @property
    def MPI_PREDICTOR_PATH(self) -> Path:
        """Get MPI predictor path, prefer v2 if available"""
        v2_path = self.MODELS_DIR / "all_mpi_model_pred_v2.pth"
        if v2_path.exists():
            return self._get_cached_path("mpi_predictor_v2", v2_path)
        return self._get_cached_path("mpi_predictor", self.MODELS_DIR / "all_mpi_model_pred.pth")
    
    def _get_cached_path(self, key: str, path: Path) -> Path:
        """Get cached path or create and cache it"""
        if key not in self._path_cache:
            self._path_cache[key] = path
        return self._path_cache[key]
    
    # Genome mapping
    GENOME_MAPPING = {
        'Homo sapiens': "Homo_sapiens",
        'Mus musculus': "Mus_musculus", 
        'Rattus norvegicus': "Rattus_norvegicus",
        'Escherichia coli': "Escherichia_coli",
        'Bos taurus': "Bos_taurus",
        'Pseudomonas aeruginosa': "Pseudomonas_aeruginosa",
        'Arabidopsis thaliana': "Arabidopsis_thaliana",
        'Saccharomyces cerevisiae': "Saccharomyces_cerevisiae",
        'Drosophila melanogaster': "Drosophila_melanogaster",
        'Caenorhabditis elegans': "Caenorhabditis_elegans"
    }
    
    # File upload settings
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB
    ALLOWED_EXTENSIONS = {'csv', 'xlsx', 'xls'}
    
    def get_network_path(self, organism: str) -> Path:
        """Get network file path for specific organism"""
        if organism == "All":
            return self.MPI_NETWORK_PATH
        organism_key = self.GENOME_MAPPING.get(organism, organism)
        return self.NETWORKS_DIR / f"pca_mpi_{organism_key}.pkl"

    def get_feature_df_path(self, organism: str) -> Path:
        """Get feature dataframe path for specific organism"""
        if organism == "All":
            return self.FEATURE_DF_PATH
        organism_key = self.GENOME_MAPPING.get(organism, organism)
        return self.FEATURES_DIR / f"pca_feature_df_{organism_key}.pkl"


class DevelopmentConfig(Config):
    """Development configuration"""
    DEBUG = True


class ProductionConfig(Config):
    """Production configuration"""
    DEBUG = False


# Configuration mapping
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'default': DevelopmentConfig
}
