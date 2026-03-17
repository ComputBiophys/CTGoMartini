"""
Author: Song Yang
Date: 2025-05-24
Version: 0.1.0
"""
import numpy as np
from openmmtools.multistate import MultiStateReporter
from openmm import unit
import warnings
from typing import List, Tuple, Optional

warnings.filterwarnings("ignore")

class PositionExtractor:
    def __init__(self, main_reporter_nc_path: str, checkpoint_nc_path: Optional[str] = None):
        """
        Initialize the PositionExtractor with structure and NetCDF files.
        
        Args:
            main_reporter_nc_path: Path to main NetCDF reporter file
            checkpoint_nc_path: Path to checkpoint NetCDF file (optional)
        """
        self.main_reporter_nc_path = main_reporter_nc_path
        self.checkpoint_nc_path = checkpoint_nc_path
        
        # Initialize reporter and keep it open
        self._reporter = MultiStateReporter(self.main_reporter_nc_path, open_mode='r', 
                                          checkpoint_storage=self.checkpoint_nc_path)
        self._initialize_metadata()
        
    def _initialize_metadata(self):
        """Initialize metadata from the reporter file."""
        trajectory_storage = self._reporter._storage_checkpoint
        
        # Store basic dimensions
        self._n_frames = trajectory_storage.variables['positions'].shape[0]
        self._n_replicas = trajectory_storage.variables['positions'].shape[1]
        self._n_atoms = trajectory_storage.variables['positions'].shape[2]

        try:
            mcmove = self._reporter.read_mcmc_moves()[0] 
            time_interval_per_iteration = mcmove.n_steps * mcmove.timestep
            self._dt = time_interval_per_iteration.value_in_unit(unit.picosecond)
        except Exception as e:
            print(f"Warning: Could not determine dt from reporter MCMC moves ({e}).")
            self._dt = 1.0

    def close(self):
        """Close the reporter explicitly."""
        if self._reporter is not None:
            self._reporter.close()
            self._reporter = None
    
    def __del__(self):
        """Destructor to ensure reporter is closed when object is garbage collected."""
        self.close()
    
    @property
    def n_frames(self) -> int:
        """Number of frames in the trajectory."""
        return self._n_frames
    
    @property
    def n_replicas(self) -> int:
        """Number of replicas in the simulation."""
        return self._n_replicas
    
    @property
    def n_atoms(self) -> int:
        """Number of atoms in the system."""
        return self._n_atoms
    
    @property
    def dt(self) -> float:
        """Time interval between frames in picoseconds"""
        return self._dt

    def get_positions(self, frame_indices: np.ndarray, replica_indices: np.ndarray, 
                     atom_indices: np.ndarray) -> np.ndarray:
        """
        Read positions for specified frames, replicas, and atoms from NetCDF storage.
        
        Args:
            frame_indices: Array of frame indices to read
            replica_indices: Array of replica indices to read
            atom_indices: Array of atom indices to read
            
        Returns:
            Array of positions in Angstroms (shape: n_frames x n_replicas x n_atoms x 3)
        """
        try:
            if self._reporter is None:
                raise RuntimeError("Reporter is closed")
                
            trajectory_storage = self._reporter._storage_checkpoint
            positions = trajectory_storage.variables['positions'][
                frame_indices, replica_indices, atom_indices, :
            ] * 10.0  # Convert nm to Angstrom

            return positions
            
        except Exception as e:
            print(f"Error processing frames {frame_indices[0]}-{frame_indices[-1]}: {e}")
            return np.empty([0, 4])

    def get_positions_for_atom_groups(self, atom_groups: List[np.ndarray], 
                                    frame_indices: np.ndarray, 
                                    replica_indices: np.ndarray) -> List[np.ndarray]:
        """
        Get positions for multiple atom groups with minimal I/O.
        
        Args:
            atom_groups: List of arrays containing atom indices for each group
            frame_indices: Array of frame indices to read
            replica_indices: Array of replica indices to read
            
        Returns:
            List of position arrays for each atom group
        """
        unique_atom_indices = np.unique(np.concatenate(atom_groups))
        
        if len(unique_atom_indices) == 0:
            print("Error: No atoms selected. Please check atom selection strings.")
            return [np.empty([0, 4]) for _ in atom_groups]
        
        # Create mapping from original indices to compact indices
        map_orig_to_compact = {orig_idx: compact_idx for compact_idx, orig_idx in enumerate(unique_atom_indices)}
        group_indices_list = [
            np.array([map_orig_to_compact[idx] for idx in atom_group]) for atom_group in atom_groups
        ]
        positions = self.get_positions(frame_indices, replica_indices, unique_atom_indices)
        
        # Extract positions for each group using compact indices
        # print(positions.shape)
        positions_groups = [positions[:, :, group_indices, :] for group_indices in group_indices_list]
        return positions_groups

