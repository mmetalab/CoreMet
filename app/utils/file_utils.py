"""
File utility functions for CoreMet Application
"""

import os
import pandas as pd
from pathlib import Path
from typing import Optional, Union
import tempfile

def save_uploaded_file(uploaded_file, upload_dir: str = "uploads") -> Optional[str]:
    """Save uploaded file to specified directory"""
    try:
        # Create upload directory if it doesn't exist
        upload_path = Path(upload_dir)
        upload_path.mkdir(exist_ok=True)
        
        # Generate unique filename
        file_path = upload_path / uploaded_file.filename
        
        # Save file
        with open(file_path, 'wb') as f:
            f.write(uploaded_file.read())
        
        return str(file_path)
    except Exception as e:
        print(f"Error saving uploaded file: {e}")
        return None

def load_data_file(file_path: Union[str, Path], file_type: str = "csv") -> Optional[pd.DataFrame]:
    """Load data file based on file type"""
    try:
        file_path = Path(file_path)
        
        if not file_path.exists():
            print(f"File not found: {file_path}")
            return None
        
        if file_type.lower() == "csv":
            return pd.read_csv(file_path)
        elif file_type.lower() in ["xlsx", "xls"]:
            return pd.read_excel(file_path)
        elif file_type.lower() == "json":
            return pd.read_json(file_path)
        else:
            print(f"Unsupported file type: {file_type}")
            return None
            
    except Exception as e:
        print(f"Error loading data file: {e}")
        return None

def validate_file_extension(filename: str, allowed_extensions: list = None) -> bool:
    """Validate file extension"""
    if allowed_extensions is None:
        allowed_extensions = ['.csv', '.xlsx', '.xls', '.json']
    
    file_ext = Path(filename).suffix.lower()
    return file_ext in allowed_extensions

def get_temp_file_path(extension: str = ".csv") -> str:
    """Get temporary file path"""
    temp_file = tempfile.NamedTemporaryFile(suffix=extension, delete=False)
    temp_file.close()
    return temp_file.name
