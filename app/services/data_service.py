"""
Data service for loading and processing data files
"""

import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import base64
import io

from app.config import Config

logger = logging.getLogger(__name__)


class DataService:
    """Service for data loading and processing"""
    
    def __init__(self):
        self.config = Config()
        self._docking_db = None
        self._mpi_db = None

    @property
    def docking_db(self):
        """Load Docking DB only when a docking route needs it."""
        if self._docking_db is None:
            self._load_docking_db()
        return self._docking_db

    @property
    def mpi_db(self):
        """Load MPI DB only when an MPI route needs it."""
        if self._mpi_db is None:
            self._load_mpi_db()
        return self._mpi_db

    def _load_docking_db(self):
        try:
            if self.config.DOCKING_DB_PATH.exists():
                self._docking_db = pd.read_csv(self.config.DOCKING_DB_PATH)
            else:
                self._docking_db = pd.DataFrame()
        except Exception as e:
            logger.error(f"Error loading docking database: {e}")
            self._docking_db = pd.DataFrame()

    def _load_mpi_db(self):
        try:
            if self.config.MPI_DB_PATH.exists():
                from app.services.csv_loader import load_optimized
                self._mpi_db = load_optimized(self.config.MPI_DB_PATH)
            else:
                self._mpi_db = pd.DataFrame()
        except Exception as e:
            logger.error(f"Error loading MPI database: {e}")
            self._mpi_db = pd.DataFrame()
    
    def _load_databases(self):
        """Load database files"""
        try:
            self._load_docking_db()
            self._load_mpi_db()
            logger.info("Databases loaded successfully")

        except Exception as e:
            logger.error(f"Error loading databases: {e}")
    
    def parse_uploaded_file(self, contents: str, filename: str) -> pd.DataFrame:
        """
        Parse uploaded file content
        
        Args:
            contents: Base64 encoded file content
            filename: Name of the uploaded file
            
        Returns:
            DataFrame with parsed data
        """
        try:
            content_type, content_string = contents.split(',')
            decoded = base64.b64decode(content_string)
            
            if 'csv' in filename.lower():
                df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            elif 'xls' in filename.lower():
                df = pd.read_excel(io.BytesIO(decoded))
            else:
                raise ValueError(f"Unsupported file type: {filename}")
            
            return df
            
        except Exception as e:
            logger.error(f"Error parsing file {filename}: {e}")
            raise
    
    def parse_text_input(self, text: str, data_type: str) -> pd.DataFrame:
        """
        Parse manually entered text data
        
        Args:
            text: Text input from user
            data_type: Type of data ('metabolite' or 'protein')
            
        Returns:
            DataFrame with parsed data
        """
        try:
            lines = [line.strip() for line in text.split('\n') if line.strip()]
            data = []
            
            if data_type == 'metabolite':
                columns = ['Metabolite Name', 'HMDB ID', 'SMILES']
                for line in lines:
                    parts = line.split(',')
                    if len(parts) >= 3:
                        data.append({
                            'Metabolite Name': parts[0].strip(),
                            'HMDB ID': parts[1].strip(),
                            'SMILES': parts[2].strip()
                        })
            
            elif data_type == 'protein':
                columns = ['UniprotID', 'Protein Name', 'Gene Name', 'Organism', 'Sequence']
                for line in lines:
                    parts = line.split(',')
                    if len(parts) >= 5:
                        data.append({
                            'UniprotID': parts[0].strip(),
                            'Protein Name': parts[1].strip(),
                            'Gene Name': parts[2].strip(),
                            'Organism': parts[3].strip(),
                            'Sequence': parts[4].strip()
                        })
            
            return pd.DataFrame(data)
            
        except Exception as e:
            logger.error(f"Error parsing text input: {e}")
            raise
    
    def get_docking_options(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """
        Get available docking options
        
        Returns:
            Tuple of (metabolite_dict, protein_dict)
        """
        if self.docking_db is None:
            return {}, {}
        
        # Create metabolite dictionary
        metabolite_dict = {}
        for _, row in self.docking_db.iterrows():
            metabolite_dict[row['Metabolite Name']] = row['HMDB ID']
            metabolite_dict[row['HMDB ID']] = row['HMDB ID']
        
        # Create protein dictionary
        protein_dict = {}
        for _, row in self.docking_db.iterrows():
            protein_dict[row['Protein Name']] = row['Uniprot ID']
            protein_dict[row['Uniprot ID']] = row['Uniprot ID']
            protein_dict[row['PDB']] = row['Uniprot ID']
        
        return metabolite_dict, protein_dict
    
    def get_docking_file_paths(self, metabolite_id: str, protein_id: str) -> Tuple[Optional[Path], Optional[Path]]:
        """
        Get file paths for docking visualization
        
        Args:
            metabolite_id: HMDB ID of metabolite
            protein_id: Uniprot ID of protein
            
        Returns:
            Tuple of (pdb_path, docking_path)
        """
        if self.docking_db is None:
            return None, None
        
        # Find matching row
        row = self.docking_db[
            (self.docking_db['HMDB ID'] == metabolite_id) & 
            (self.docking_db['Uniprot ID'] == protein_id)
        ]
        
        if row.empty:
            return None, None
        
        pdb_id = row.iloc[0]['PDB']
        
        # Construct file paths
        pdb_path = self.config.DATABASE_DIR / "DockingDB" / "PDB" / f"{pdb_id}.pdb"
        docking_path = self.config.DATABASE_DIR / "DockingDB" / "Docking" / f"{metabolite_id}_{protein_id}_{pdb_id}.pdbqt"
        
        return pdb_path, docking_path
    
    def get_mpi_database_info(self) -> Dict[str, any]:
        """
        Get MPI database information

        Returns:
            Dictionary with database statistics
        """
        if self.mpi_db is None:
            return {}

        info = {
            'total_interactions': len(self.mpi_db),
            'unique_metabolites': self.mpi_db['HMDB ID'].nunique() if 'HMDB ID' in self.mpi_db.columns else 0,
            'unique_proteins': self.mpi_db['Uniprot ID'].nunique() if 'Uniprot ID' in self.mpi_db.columns else 0,
            'organisms': self.mpi_db['Species'].unique().tolist() if 'Species' in self.mpi_db.columns else [],
        }

        # New v2 columns
        if 'Evidence_Source' in self.mpi_db.columns:
            info['evidence_sources'] = self.mpi_db['Evidence_Source'].value_counts().to_dict()
        if 'Pathway_ID' in self.mpi_db.columns:
            info['rows_with_pathways'] = int(self.mpi_db['Pathway_ID'].notna().sum() & (self.mpi_db['Pathway_ID'] != '').sum())

        return info
    
    def validate_metabolite_data(self, df: pd.DataFrame) -> Tuple[bool, List[str]]:
        """
        Validate metabolite data format
        
        Args:
            df: DataFrame to validate
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        required_columns = ['Metabolite Name', 'HMDB ID', 'SMILES']
        
        # Check required columns
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            errors.append(f"Missing required columns: {missing_columns}")
        
        # Check for empty values
        for col in required_columns:
            if col in df.columns and df[col].isnull().any():
                errors.append(f"Column '{col}' contains empty values")
        
        return len(errors) == 0, errors
    
    def validate_protein_data(self, df: pd.DataFrame) -> Tuple[bool, List[str]]:
        """
        Validate protein data format
        
        Args:
            df: DataFrame to validate
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        required_columns = ['UniprotID', 'Protein Name', 'Gene Name', 'Organism', 'Sequence']
        
        # Check required columns
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            errors.append(f"Missing required columns: {missing_columns}")
        
        # Check for empty values
        for col in required_columns:
            if col in df.columns and df[col].isnull().any():
                errors.append(f"Column '{col}' contains empty values")
        
        return len(errors) == 0, errors
    
    def get_example_data(self, data_type: str) -> pd.DataFrame:
        """
        Get example data for templates
        
        Args:
            data_type: Type of example data ('metabolite' or 'protein')
            
        Returns:
            DataFrame with example data
        """
        if data_type == 'metabolite':
            return pd.DataFrame({
                'Metabolite Name': ['1-methyl-L-histidine'],
                'HMDB ID': ['HMDB0000001'],
                'SMILES': ['CN1C=NC(C[C@H](N)C(O)=O)=C1']
            })
        
        elif data_type == 'protein':
            return pd.DataFrame({
                'UniprotID': ['A0A075B6H7'],
                'Protein Name': ['Probable non-functional immunoglobulin kappa variable 3-7'],
                'Gene Name': ['IGKV3-7'],
                'Organism': ['Homo sapiens (Human)'],
                'Sequence': ['MEAPAQLLFLLLLWLPDTTREIVMTQSPPTLSLSPGERVTLSCRASQSVSSSYLTWYQQKPGQAPRLLIYGASTRATSIPARFSGSGSGTDFTLTISSLQPEDFAVYYCQQDYNLP']
            })
        
        return pd.DataFrame()
