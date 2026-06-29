"""
Feature extraction models for CoreMet Application
"""

import numpy as np
from typing import Optional, Dict, Any
from functools import lru_cache

class FeatureExtractor:
    """Feature extraction utility class"""
    
    def __init__(self):
        self.cache = {}
    
    @lru_cache(maxsize=1000)
    def extract_smiles_features(self, smiles: str) -> Optional[np.ndarray]:
        """Extract features from SMILES string"""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, DataStructs
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            array = np.zeros((0,), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, array)
            return array
        except Exception as e:
            print(f"Error extracting SMILES features: {e}")
            return None
    
    def extract_protein_features(self, sequence: str) -> Optional[np.ndarray]:
        """Extract features from protein sequence"""
        # Simple feature extraction - can be enhanced
        if not sequence:
            return None
        
        # Basic features: length, amino acid composition
        features = np.array([
            len(sequence),
            sequence.count('A') / len(sequence),
            sequence.count('C') / len(sequence),
            sequence.count('D') / len(sequence),
            sequence.count('E') / len(sequence),
            sequence.count('F') / len(sequence),
            sequence.count('G') / len(sequence),
            sequence.count('H') / len(sequence),
            sequence.count('I') / len(sequence),
            sequence.count('K') / len(sequence),
            sequence.count('L') / len(sequence),
            sequence.count('M') / len(sequence),
            sequence.count('N') / len(sequence),
            sequence.count('P') / len(sequence),
            sequence.count('Q') / len(sequence),
            sequence.count('R') / len(sequence),
            sequence.count('S') / len(sequence),
            sequence.count('T') / len(sequence),
            sequence.count('V') / len(sequence),
            sequence.count('W') / len(sequence),
            sequence.count('Y') / len(sequence),
        ])
        return features
