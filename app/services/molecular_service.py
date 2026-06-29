"""
Optimized molecular service for feature extraction and processing
"""

import logging
import numpy as np
import pandas as pd
from typing import Optional, Dict, Any
from pathlib import Path
from functools import lru_cache
import time

from app.config import Config
from app.services.cache_service import OptimizedDataLoader


logger = logging.getLogger(__name__)


class MolecularService:
    """Optimized service for molecular feature extraction and processing"""

    def __init__(self):
        self.config = Config()
        self.data_loader = OptimizedDataLoader()
        self.mets_pca = None
        self.protein_pca = None
        self.protein_vectors = None
        self._load_models()
    
    def _load_models(self):
        """Load pre-trained feature extraction models with caching"""
        try:
            start_time = time.time()
            
            # Load models with caching
            self.mets_pca = self.data_loader.load_metabolite_pca()
            self.protein_pca = self.data_loader.load_protein_pca()
            self.protein_vectors = self.data_loader.load_protein_vectors()
            
            load_time = time.time() - start_time
            logger.info(f"Molecular feature models loaded successfully in {load_time:.2f} seconds")

        except Exception as e:
            logger.error(f"Error loading molecular models: {e}")
            raise
    
    @lru_cache(maxsize=1000)
    def extract_metabolite_features(self, smiles: str) -> Optional[np.ndarray]:
        """
        Extract features from SMILES string with caching
        
        Args:
            smiles: SMILES string of the metabolite
            
        Returns:
            PCA-transformed features or None if extraction fails
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, DataStructs
            
            # Convert SMILES to molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # Generate Morgan fingerprint efficiently
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            array = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, array)
            
            # Apply PCA transformation
            if self.mets_pca is not None:
                pca_features = self.mets_pca.transform(array.reshape(1, -1))
                return pca_features[0]
            else:
                return array
                
        except Exception as e:
            logger.error(f"Error extracting metabolite features: {e}")
            return None
    
    def extract_protein_features(self, uniprot_id: str) -> Optional[np.ndarray]:
        """
        Extract features from protein Uniprot ID
        
        Args:
            uniprot_id: Uniprot ID of the protein
            
        Returns:
            PCA-transformed features or None if not found
        """
        try:
            if self.protein_vectors is None:
                logger.warning("Protein vectors not loaded")
                return None

            if uniprot_id not in self.protein_vectors:
                logger.warning(f"Protein {uniprot_id} not found in database")
                return None
            
            # Get protein vector
            protein_vector = self.protein_vectors[uniprot_id]
            
            # Apply PCA transformation
            if self.protein_pca is not None:
                pca_features = self.protein_pca.transform(protein_vector.reshape(1, -1))
                return pca_features[0]
            else:
                return protein_vector
                
        except Exception as e:
            logger.error(f"Error extracting protein features: {e}")
            return None
    
    def create_molecular_feature_dataframe(self, 
                                         metabolites_df: pd.DataFrame,
                                         proteins_df: pd.DataFrame) -> pd.DataFrame:
        """
        Create feature dataframe from metabolite and protein data
        
        Args:
            metabolites_df: DataFrame with metabolite data
            proteins_df: DataFrame with protein data
            
        Returns:
            DataFrame with molecular features
        """
        feature_data = []
        
        # Process metabolites
        for _, row in metabolites_df.iterrows():
            features = self.extract_metabolite_features(row['SMILES'])
            if features is not None:
                feature_data.append({
                    'Name': row['Metabolite Name'],
                    'Type': 'Metabolite',
                    'HMDB ID': row['HMDB ID'],
                    'Feature': features
                })
        
        # Process proteins
        for _, row in proteins_df.iterrows():
            features = self.extract_protein_features(row['UniprotID'])
            if features is not None:
                feature_data.append({
                    'Name': row['UniprotID'],
                    'Type': 'Protein',
                    'Feature': features
                })
        
        return pd.DataFrame(feature_data)
    
    def validate_smiles(self, smiles: str) -> bool:
        """
        Validate SMILES string
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            True if valid, False otherwise
        """
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except Exception:
            logging.getLogger(__name__).debug("SMILES validation failed", exc_info=True)
            return False
    
    def validate_protein_sequence(self, sequence: str) -> bool:
        """
        Validate protein sequence
        
        Args:
            sequence: Protein sequence to validate
            
        Returns:
            True if valid, False otherwise
        """
        if not sequence or len(sequence) < 10:
            return False
        
        # Check for valid amino acid characters
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        return all(aa in valid_aa for aa in sequence.upper())
    
    def get_protein_info(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """
        Get protein information from database
        
        Args:
            uniprot_id: Uniprot ID of the protein
            
        Returns:
            Dictionary with protein information or None
        """
        try:
            if self.protein_vectors is None or uniprot_id not in self.protein_vectors:
                return None
            
            return {
                'uniprot_id': uniprot_id,
                'has_features': True,
                'feature_dimension': len(self.protein_vectors[uniprot_id])
            }
            
        except Exception as e:
            logger.error(f"Error getting protein info: {e}")
            return None
