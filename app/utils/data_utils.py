"""
Data utility functions for CoreMet Application
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Any, Tuple


def cosine_similarity(arr1: np.ndarray, arr2: np.ndarray) -> float:
    """
    Calculate cosine similarity between two arrays
    
    Args:
        arr1: First array
        arr2: Second array
        
    Returns:
        Cosine similarity score
    """
    dot_product = np.dot(arr1, arr2)
    norm_arr1 = np.linalg.norm(arr1)
    norm_arr2 = np.linalg.norm(arr2)
    
    if norm_arr1 == 0 or norm_arr2 == 0:
        return 0.0
    
    similarity = dot_product / (norm_arr1 * norm_arr2)
    return similarity


def validate_data_format(data: pd.DataFrame, data_type: str) -> Tuple[bool, List[str]]:
    """
    Validate data format for metabolites or proteins
    
    Args:
        data: DataFrame to validate
        data_type: Type of data ('metabolite' or 'protein')
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    
    if data_type == 'metabolite':
        required_columns = ['Metabolite Name', 'HMDB ID', 'SMILES']
    elif data_type == 'protein':
        required_columns = ['UniprotID', 'Protein Name', 'Gene Name', 'Organism', 'Sequence']
    else:
        return False, [f"Unknown data type: {data_type}"]
    
    # Check required columns
    missing_columns = [col for col in required_columns if col not in data.columns]
    if missing_columns:
        errors.append(f"Missing required columns: {missing_columns}")
    
    # Check for empty values
    for col in required_columns:
        if col in data.columns:
            empty_count = data[col].isnull().sum()
            if empty_count > 0:
                errors.append(f"Column '{col}' has {empty_count} empty values")
    
    # Check for duplicate entries
    if data_type == 'metabolite' and 'HMDB ID' in data.columns:
        duplicate_count = data['HMDB ID'].duplicated().sum()
        if duplicate_count > 0:
            errors.append(f"Found {duplicate_count} duplicate HMDB IDs")
    
    elif data_type == 'protein' and 'UniprotID' in data.columns:
        duplicate_count = data['UniprotID'].duplicated().sum()
        if duplicate_count > 0:
            errors.append(f"Found {duplicate_count} duplicate Uniprot IDs")
    
    return len(errors) == 0, errors


def clean_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Clean and standardize data
    
    Args:
        data: DataFrame to clean
        
    Returns:
        Cleaned DataFrame
    """
    # Remove leading/trailing whitespace from string columns
    string_columns = data.select_dtypes(include=['object']).columns
    for col in string_columns:
        data[col] = data[col].astype(str).str.strip()
    
    # Replace empty strings with NaN
    data = data.replace('', np.nan)
    
    # Remove rows where all values are NaN
    data = data.dropna(how='all')
    
    return data


def merge_metabolite_protein_data(metabolites_df: pd.DataFrame, 
                                proteins_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge metabolite and protein data for processing
    
    Args:
        metabolites_df: Metabolite DataFrame
        proteins_df: Protein DataFrame
        
    Returns:
        Merged DataFrame
    """
    # Add type column to distinguish metabolites and proteins
    metabolites_df = metabolites_df.copy()
    proteins_df = proteins_df.copy()
    
    metabolites_df['Type'] = 'Metabolite'
    proteins_df['Type'] = 'Protein'
    
    # Standardize column names for merging
    metabolite_columns = {
        'Metabolite Name': 'Name',
        'HMDB ID': 'ID'
    }
    protein_columns = {
        'UniprotID': 'ID',
        'Protein Name': 'Name'
    }
    
    metabolites_df = metabolites_df.rename(columns=metabolite_columns)
    proteins_df = proteins_df.rename(columns=protein_columns)
    
    # Select common columns
    common_columns = ['Name', 'ID', 'Type']
    metabolites_subset = metabolites_df[common_columns]
    proteins_subset = proteins_df[common_columns]
    
    # Merge dataframes
    merged_df = pd.concat([metabolites_subset, proteins_subset], ignore_index=True)
    
    return merged_df


def calculate_statistics(data: pd.DataFrame) -> Dict[str, Any]:
    """
    Calculate basic statistics for data
    
    Args:
        data: DataFrame to analyze
        
    Returns:
        Dictionary with statistics
    """
    stats = {
        'total_rows': len(data),
        'total_columns': len(data.columns),
        'missing_values': data.isnull().sum().to_dict(),
        'data_types': data.dtypes.to_dict()
    }
    
    # Add type-specific statistics
    if 'Type' in data.columns:
        stats['type_counts'] = data['Type'].value_counts().to_dict()
    
    return stats


def filter_data_by_organism(data: pd.DataFrame, organism: str) -> pd.DataFrame:
    """
    Filter data by organism
    
    Args:
        data: DataFrame to filter
        organism: Target organism
        
    Returns:
        Filtered DataFrame
    """
    if 'Organism' not in data.columns:
        return data
    
    if organism == 'All':
        return data
    
    return data[data['Organism'].str.contains(organism, case=False, na=False)]


def sample_data(data: pd.DataFrame, n_samples: int = 10) -> pd.DataFrame:
    """
    Sample data for preview
    
    Args:
        data: DataFrame to sample
        n_samples: Number of samples to return
        
    Returns:
        Sampled DataFrame
    """
    if len(data) <= n_samples:
        return data
    
    return data.sample(n=n_samples, random_state=42)
