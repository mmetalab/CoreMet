"""
Optimized validation utility functions for CoreMet Application
"""

import re
from typing import List, Tuple, Optional
from functools import lru_cache

# Pre-compiled regex patterns for better performance
SMILES_INVALID_CHARS_PATTERN = re.compile(r'[^A-Za-z0-9@+\-\[\]()=#\\/]')
SMILES_ATOM_PATTERN = re.compile(r'[A-Za-z]')
HMDB_ID_PATTERN = re.compile(r'^HMDB\d{7}$')
UNIPROT_ID_PATTERN = re.compile(r'^[A-Z0-9]{1,5}_[A-Z0-9]{1,5}$')
ORGANISM_NAME_PATTERN = re.compile(r'^[A-Za-z\s()\-]+$')


@lru_cache(maxsize=1000)
def validate_smiles(smiles: str) -> Tuple[bool, Optional[str]]:
    """
    Validate SMILES string format
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not smiles or not isinstance(smiles, str):
        return False, "SMILES must be a non-empty string"
    
    # Basic SMILES validation patterns using pre-compiled regex
    # Check for invalid characters
    invalid_chars = SMILES_INVALID_CHARS_PATTERN.findall(smiles)
    if invalid_chars:
        return False, f"Invalid characters in SMILES: {set(invalid_chars)}"
    
    # Check for balanced parentheses and brackets
    if not _balanced_brackets(smiles):
        return False, "Unbalanced parentheses or brackets in SMILES"
    
    # Check for basic structure (at least one atom)
    if not SMILES_ATOM_PATTERN.search(smiles):
        return False, "SMILES must contain at least one atom"
    
    return True, None


def validate_protein_sequence(sequence: str) -> Tuple[bool, Optional[str]]:
    """
    Validate protein sequence format
    
    Args:
        sequence: Protein sequence to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not sequence or not isinstance(sequence, str):
        return False, "Sequence must be a non-empty string"
    
    # Check sequence length
    if len(sequence) < 10:
        return False, "Protein sequence must be at least 10 amino acids long"
    
    if len(sequence) > 10000:
        return False, "Protein sequence too long (max 10,000 amino acids)"
    
    # Check for valid amino acid characters
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    sequence_upper = sequence.upper()
    invalid_aa = set(sequence_upper) - valid_aa
    
    if invalid_aa:
        return False, f"Invalid amino acid characters: {invalid_aa}"
    
    return True, None


@lru_cache(maxsize=1000)
def validate_hmdb_id(hmdb_id: str) -> Tuple[bool, Optional[str]]:
    """
    Validate HMDB ID format with caching
    
    Args:
        hmdb_id: HMDB ID to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not hmdb_id or not isinstance(hmdb_id, str):
        return False, "HMDB ID must be a non-empty string"
    
    # HMDB ID format: HMDB followed by 7 digits using pre-compiled pattern
    if not HMDB_ID_PATTERN.match(hmdb_id.upper()):
        return False, "HMDB ID must be in format HMDB followed by 7 digits (e.g., HMDB0000001)"
    
    return True, None


@lru_cache(maxsize=1000)
def validate_uniprot_id(uniprot_id: str) -> Tuple[bool, Optional[str]]:
    """
    Validate Uniprot ID format with caching
    
    Args:
        uniprot_id: Uniprot ID to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not uniprot_id or not isinstance(uniprot_id, str):
        return False, "Uniprot ID must be a non-empty string"
    
    # Uniprot ID format using pre-compiled pattern
    if not UNIPROT_ID_PATTERN.match(uniprot_id.upper()):
        return False, "Uniprot ID must be in format XXXX_XXXX (e.g., P12345_1)"
    
    return True, None


@lru_cache(maxsize=1000)
def validate_organism_name(organism: str) -> Tuple[bool, Optional[str]]:
    """
    Validate organism name format with caching
    
    Args:
        organism: Organism name to validate
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not organism or not isinstance(organism, str):
        return False, "Organism name must be a non-empty string"
    
    # Check for valid organism name format (genus species)
    if len(organism.split()) < 2:
        return False, "Organism name should include genus and species"
    
    # Check for valid characters using pre-compiled pattern
    if not ORGANISM_NAME_PATTERN.match(organism):
        return False, "Organism name contains invalid characters"
    
    return True, None


def validate_file_extension(filename: str, allowed_extensions: List[str]) -> Tuple[bool, Optional[str]]:
    """
    Validate file extension
    
    Args:
        filename: Name of the file
        allowed_extensions: List of allowed extensions
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if not filename:
        return False, "Filename cannot be empty"
    
    # Get file extension
    if '.' not in filename:
        return False, "File must have an extension"
    
    extension = filename.split('.')[-1].lower()
    
    if extension not in allowed_extensions:
        return False, f"File extension '{extension}' not allowed. Allowed extensions: {allowed_extensions}"
    
    return True, None


def validate_file_size(file_size: int, max_size: int) -> Tuple[bool, Optional[str]]:
    """
    Validate file size
    
    Args:
        file_size: Size of the file in bytes
        max_size: Maximum allowed size in bytes
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if file_size <= 0:
        return False, "File size must be greater than 0"
    
    if file_size > max_size:
        max_size_mb = max_size / (1024 * 1024)
        file_size_mb = file_size / (1024 * 1024)
        return False, f"File size ({file_size_mb:.1f} MB) exceeds maximum allowed size ({max_size_mb:.1f} MB)"
    
    return True, None


def _balanced_brackets(text: str) -> bool:
    """
    Check if brackets and parentheses are balanced in text
    
    Args:
        text: Text to check
        
    Returns:
        True if balanced, False otherwise
    """
    stack = []
    brackets = {'(': ')', '[': ']', '{': '}'}
    
    for char in text:
        if char in brackets:
            stack.append(char)
        elif char in brackets.values():
            if not stack:
                return False
            if brackets[stack.pop()] != char:
                return False
    
    return len(stack) == 0


def validate_upload_data(data_type: str, data: dict) -> Tuple[bool, List[str]]:
    """
    Validate uploaded data based on type
    
    Args:
        data_type: Type of data ('metabolite' or 'protein')
        data: Dictionary with uploaded data
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    
    if data_type == 'metabolite':
        errors.extend(_validate_metabolite_data(data))
    elif data_type == 'protein':
        errors.extend(_validate_protein_data(data))
    else:
        errors.append(f"Unknown data type: {data_type}")
    
    return len(errors) == 0, errors


def _validate_metabolite_data(data: dict) -> List[str]:
    """Validate metabolite data"""
    errors = []
    
    required_fields = ['Metabolite Name', 'HMDB ID', 'SMILES']
    for field in required_fields:
        if field not in data or not data[field]:
            errors.append(f"Missing required field: {field}")
    
    # Validate HMDB ID if present
    if 'HMDB ID' in data and data['HMDB ID']:
        is_valid, error_msg = validate_hmdb_id(data['HMDB ID'])
        if not is_valid:
            errors.append(f"Invalid HMDB ID: {error_msg}")
    
    # Validate SMILES if present
    if 'SMILES' in data and data['SMILES']:
        is_valid, error_msg = validate_smiles(data['SMILES'])
        if not is_valid:
            errors.append(f"Invalid SMILES: {error_msg}")
    
    return errors


def _validate_protein_data(data: dict) -> List[str]:
    """Validate protein data"""
    errors = []
    
    required_fields = ['UniprotID', 'Protein Name', 'Gene Name', 'Organism', 'Sequence']
    for field in required_fields:
        if field not in data or not data[field]:
            errors.append(f"Missing required field: {field}")
    
    # Validate Uniprot ID if present
    if 'UniprotID' in data and data['UniprotID']:
        is_valid, error_msg = validate_uniprot_id(data['UniprotID'])
        if not is_valid:
            errors.append(f"Invalid Uniprot ID: {error_msg}")
    
    # Validate sequence if present
    if 'Sequence' in data and data['Sequence']:
        is_valid, error_msg = validate_protein_sequence(data['Sequence'])
        if not is_valid:
            errors.append(f"Invalid protein sequence: {error_msg}")
    
    # Validate organism if present
    if 'Organism' in data and data['Organism']:
        is_valid, error_msg = validate_organism_name(data['Organism'])
        if not is_valid:
            errors.append(f"Invalid organism name: {error_msg}")
    
    return errors
