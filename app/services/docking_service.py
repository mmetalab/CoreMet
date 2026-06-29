"""
Docking service for molecular visualization
"""

# Handle py3Dmol import gracefully
try:
    import py3Dmol
    PY3DMOL_AVAILABLE = True
except ImportError as e:
    print(f"Warning: py3Dmol not available: {e}")
    PY3DMOL_AVAILABLE = False
    py3Dmol = None

from typing import List, Optional, Tuple
from pathlib import Path

from app.config import Config


class DockingService:
    """Service for molecular docking visualization"""
    
    def __init__(self):
        self.config = Config()
    
    def parse_vina_output(self, output_file_path: Path) -> List[str]:
        """
        Parse AutoDock Vina output file to extract docked poses
        
        Args:
            output_file_path: Path to the Vina output file
            
        Returns:
            List of pose strings
        """
        poses = []
        current_pose = ""
        
        try:
            with open(output_file_path, 'r') as f:
                for line in f:
                    if line.startswith('MODEL'):
                        current_pose = line
                    elif line.startswith('ENDMDL'):
                        if current_pose:
                            poses.append(current_pose)
                            current_pose = ""
                    else:
                        current_pose += line
            
            # Add the last pose if file doesn't end with ENDMDL
            if current_pose and not current_pose.strip().endswith('ENDMDL'):
                poses.append(current_pose)
                
        except Exception as e:
            print(f"Error parsing Vina output: {e}")
            return []
        
        return poses
    
    def create_molecular_viewer(self, 
                              pdb_path: Path, 
                              poses: List[str], 
                              width: int = 800, 
                              height: int = 600) -> str:
        """
        Create 3D molecular viewer HTML
        
        Args:
            pdb_path: Path to PDB file
            poses: List of docked poses
            width: Viewer width
            height: Viewer height
            
        Returns:
            HTML string for the viewer
        """
        try:
            # Create py3Dmol view
            view = py3Dmol.view(width=width, height=height)
            
            # Add protein structure
            if pdb_path.exists():
                with open(pdb_path, 'r') as f:
                    pdb_content = f.read()
                view.addModel(pdb_content, 'pdb')
                view.setStyle({'cartoon': {'color': 'spectrum'}})
            
            # Add ligand poses
            if poses:
                for i, pose in enumerate(poses):
                    view.addModel(pose, 'pdb')
                    view.setStyle({'model': i + 1}, {'stick': {'colorscheme': 'greenCarbon'}})
            
            # Center and zoom
            view.zoomTo()
            
            # Generate HTML
            return view._make_html()
            
        except Exception as e:
            print(f"Error creating molecular viewer: {e}")
            return self._create_error_viewer(width, height)
    
    def _create_error_viewer(self, width: int, height: int) -> str:
        """Create error viewer when molecular visualization fails"""
        return f"""
        <div style="width: {width}px; height: {height}px; 
                    display: flex; align-items: center; justify-content: center;
                    border: 1px solid #ccc; background-color: #f8f9fa;">
            <p style="color: #6c757d; text-align: center;">
                Unable to load molecular visualization.<br>
                Please check that the required files are available.
            </p>
        </div>
        """
    
    def get_available_docking_pairs(self) -> List[Tuple[str, str, str]]:
        """
        Get list of available docking pairs
        
        Returns:
            List of tuples (metabolite_name, protein_name, pdb_id)
        """
        try:
            from app.services.data_service import DataService
            
            data_service = DataService()
            if data_service.docking_db is None:
                return []
            
            pairs = []
            for _, row in data_service.docking_db.iterrows():
                pairs.append((
                    row['Metabolite Name'],
                    row['Protein Name'], 
                    row['PDB']
                ))
            
            return pairs
            
        except Exception as e:
            print(f"Error getting docking pairs: {e}")
            return []
    
    def validate_docking_files(self, metabolite_id: str, protein_id: str) -> Tuple[bool, str]:
        """
        Validate that docking files exist
        
        Args:
            metabolite_id: HMDB ID of metabolite
            protein_id: Uniprot ID of protein
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            from app.services.data_service import DataService
            
            data_service = DataService()
            pdb_path, docking_path = data_service.get_docking_file_paths(metabolite_id, protein_id)
            
            if pdb_path is None or docking_path is None:
                return False, "Docking pair not found in database"
            
            if not pdb_path.exists():
                return False, f"PDB file not found: {pdb_path}"
            
            if not docking_path.exists():
                return False, f"Docking file not found: {docking_path}"
            
            return True, "Files validated successfully"
            
        except Exception as e:
            return False, f"Validation error: {e}"
    
    def get_docking_statistics(self) -> dict:
        """
        Get docking database statistics
        
        Returns:
            Dictionary with statistics
        """
        try:
            from app.services.data_service import DataService
            
            data_service = DataService()
            if data_service.docking_db is None:
                return {}
            
            return {
                'total_pairs': len(data_service.docking_db),
                'unique_metabolites': data_service.docking_db['Metabolite Name'].nunique(),
                'unique_proteins': data_service.docking_db['Protein Name'].nunique(),
                'unique_pdb_structures': data_service.docking_db['PDB'].nunique()
            }
            
        except Exception as e:
            print(f"Error getting docking statistics: {e}")
            return {}
