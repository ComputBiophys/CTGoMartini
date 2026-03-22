#!/usr/bin/env python3
"""
REMD dRMS trajectory analysis and plotting.

Analyzes distance root-mean-square (dRMS) trajectories from replica exchange
molecular dynamics simulations and generates publication-quality plots.

This module provides both calculation (from NetCDF) and plotting functionality.

"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import sys
import time
import warnings
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np

# Suppress warnings
warnings.filterwarnings("ignore")

# Lazy imports for openmmtools
_MultiStateReporter = None


def _import_openmmtools():
    """Lazy import openmmtools to avoid import-time warnings."""
    global _MultiStateReporter
    if _MultiStateReporter is None:
        from openmmtools.multistate import MultiStateReporter
        _MultiStateReporter = MultiStateReporter
    return _MultiStateReporter


def _estimate_chunk_size(
    n_process_frames: int,
    n_replicas: int,
    n_pairs: int,
    num_workers: int,
    target_chunk_time: float = 0.1,
    min_chunks_per_worker: int = 2,
) -> int:
    """
    Automatically estimate optimal chunk_size based on workload characteristics.
    
    The optimization considers:
    1. Load balancing: Ensure enough chunks for all workers
    2. Task overhead: Each chunk should take ~target_chunk_time seconds
    3. Memory efficiency: Avoid excessive memory usage
    
    Args:
        n_process_frames: Number of frames to process
        n_replicas: Number of replicas
        n_pairs: Number of atom pairs (for dRMS calculation)
        num_workers: Number of worker processes
        target_chunk_time: Target computation time per chunk (seconds)
        min_chunks_per_worker: Minimum chunks per worker for load balancing
        
    Returns:
        Optimal chunk_size
    """
    if n_process_frames <= 0:
        return 1
    
    # Constraint 1: Minimum chunks for load balancing
    # Each worker should get at least min_chunks_per_worker chunks
    min_chunks = num_workers * min_chunks_per_worker
    max_chunk_size_load = max(1, n_process_frames // min_chunks)
    
    # Constraint 2: Target computation time per chunk
    # Empirical estimate: ~1e-6 seconds per (frame * replica * pair) operation
    # This is a rough estimate based on typical dRMS calculation performance
    operations_per_frame = n_replicas * n_pairs
    if operations_per_frame > 0:
        # Estimate: 1e-7 seconds per operation (empirical)
        time_per_frame = operations_per_frame * 1e-7
        target_frames = max(1, int(target_chunk_time / time_per_frame))
    else:
        target_frames = 100
    
    # Constraint 3: Memory efficiency
    # Each frame reads: n_replicas * n_pairs * 3 coords * 4 bytes
    # Target: < 100MB per chunk
    bytes_per_frame = n_replicas * n_pairs * 3 * 4  # float32
    target_memory = 100 * 1024 * 1024  # 100MB
    max_chunk_size_memory = max(1, int(target_memory / max(bytes_per_frame, 1)))
    
    # Constraint 4: Absolute bounds
    # Minimum: 10 frames (avoid too many small tasks)
    # Maximum: 1000 frames or 10% of total (avoid huge chunks)
    min_chunk_size = 10
    max_chunk_size_absolute = min(1000, max(1, n_process_frames // 10))
    
    # Combine constraints: take the most restrictive
    chunk_size = min(target_frames, max_chunk_size_load, max_chunk_size_memory)
    chunk_size = max(min_chunk_size, min(chunk_size, max_chunk_size_absolute))
    
    # Final adjustment: ensure we don't exceed total frames
    chunk_size = min(chunk_size, n_process_frames)
    
    return chunk_size


def _format_bytes(size_bytes: int) -> str:
    """Format bytes to human readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


class PositionExtractor:
    """Extract positions from NetCDF trajectory files."""
    
    def __init__(self, main_reporter_nc_path: str, checkpoint_nc_path: Optional[str] = None):
        """
        Initialize the PositionExtractor.
        
        Args:
            main_reporter_nc_path: Path to main NetCDF reporter file
            checkpoint_nc_path: Path to checkpoint NetCDF file (optional)
        """
        MultiStateReporter = _import_openmmtools()
        
        self.main_reporter_nc_path = main_reporter_nc_path
        self.checkpoint_nc_path = checkpoint_nc_path
        
        # Initialize reporter
        self._reporter = MultiStateReporter(
            self.main_reporter_nc_path, 
            open_mode='r',
            checkpoint_storage=self.checkpoint_nc_path
        )
        self._initialize_metadata()
        
    def _initialize_metadata(self):
        """Initialize metadata from the reporter file."""
        try:
            trajectory_storage = self._reporter._storage_checkpoint
            
            # Store basic dimensions
            self._n_frames = trajectory_storage.variables['positions'].shape[0]
            self._n_replicas = trajectory_storage.variables['positions'].shape[1]
            self._n_atoms = trajectory_storage.variables['positions'].shape[2]
            
            # Try to get dt from MCMC moves
            try:
                mcmove = self._reporter.read_mcmc_moves()[0]
                time_interval = mcmove.n_steps * mcmove.timestep
                # Convert to ps
                import openmm
                self._dt = time_interval.value_in_unit(openmm.unit.picosecond)
            except Exception:
                self._dt = 1.0  # Default: 1 ps per frame
        except Exception as e:
            print(f"Error initializing metadata: {e}")
            self._n_frames = self._n_replicas = self._n_atoms = 0
            self._dt = 1.0
    
    def close(self):
        """Close the reporter."""
        if hasattr(self, '_reporter') and self._reporter is not None:
            self._reporter.close()
            self._reporter = None
    
    def __del__(self):
        """Destructor."""
        self.close()
    
    @property
    def n_frames(self) -> int:
        return self._n_frames
    
    @property
    def n_replicas(self) -> int:
        return self._n_replicas
    
    @property
    def n_atoms(self) -> int:
        return self._n_atoms
    
    @property
    def dt(self) -> float:
        return self._dt
    
    def get_positions(self, frame_idx: int, replica_idx: int) -> np.ndarray:
        """
        Get positions for a specific frame and replica.
        
        Args:
            frame_idx: Frame index
            replica_idx: Replica index
            
        Returns:
            Positions array (n_atoms x 3) in Angstroms
        """
        trajectory_storage = self._reporter._storage_checkpoint
        positions = trajectory_storage.variables['positions'][frame_idx, replica_idx, :, :]
        return positions * 10.0  # Convert nm to Angstrom
    
    def get_positions_batch(
        self,
        frame_indices: np.ndarray,
        replica_indices: np.ndarray,
        atom_indices: np.ndarray
    ) -> np.ndarray:
        """
        Get positions for multiple frames, replicas, and atoms.
        
        Uses frame-by-frame reading which is efficient for NetCDF4 and allows
        reading arbitrary combinations of frames, replicas, and atoms.
        
        Args:
            frame_indices: Array of frame indices
            replica_indices: Array of replica indices  
            atom_indices: Array of atom indices
            
        Returns:
            Positions array (n_frames x n_replicas x n_atoms x 3) in Angstroms
        """
        trajectory_storage = self._reporter._storage_checkpoint
        var = trajectory_storage.variables['positions']
        
        n_frames = len(frame_indices)
        n_replicas = len(replica_indices)
        n_atoms = len(atom_indices)
        
        # Pre-allocate result array
        positions = np.empty((n_frames, n_replicas, n_atoms, 3), dtype=np.float32)
        
        # Read frame by frame (netCDF4 supports list indexing for replicas/atoms)
        for i, f_idx in enumerate(frame_indices):
            # For each frame, read all replicas and selected atoms
            # netCDF4 supports: scalar, list, list, slice
            frame_data = var[f_idx, replica_indices.tolist(), atom_indices.tolist(), :]
            positions[i] = frame_data
        
        return positions * 10.0  # Convert nm to Angstrom
    
    def get_positions_for_atom_groups(
        self,
        atom_groups: List[np.ndarray],
        frame_indices: np.ndarray,
        replica_indices: np.ndarray
    ) -> List[np.ndarray]:
        """
        Get positions for multiple atom groups with minimal I/O.
        
        Args:
            atom_groups: List of arrays containing atom indices for each group
            frame_indices: Array of frame indices to read
            replica_indices: Array of replica indices to read
            
        Returns:
            List of position arrays for each atom group, each with shape
            (n_frames x n_replicas x n_atoms_in_group x 3)
        """
        unique_atom_indices = np.unique(np.concatenate(atom_groups))
        
        if len(unique_atom_indices) == 0:
            print("Error: No atoms selected.")
            return [np.empty((0, 4)) for _ in atom_groups]
        
        # Create mapping from original indices to compact indices
        map_orig_to_compact = {
            orig_idx: compact_idx 
            for compact_idx, orig_idx in enumerate(unique_atom_indices)
        }
        group_indices_list = [
            np.array([map_orig_to_compact[idx] for idx in atom_group])
            for atom_group in atom_groups
        ]
        
        # Single read for all unique atoms
        positions = self.get_positions_batch(
            frame_indices, replica_indices, unique_atom_indices
        )
        
        # Extract positions for each group
        return [
            positions[:, :, group_indices, :]
            for group_indices in group_indices_list
        ]


def calculate_reference_distances(
    ref_files: list[str],
    selected_atom: str = "name BB",
    min_dist: float = 6.0,
    max_dist: float = 50.0,
    min_diff: float = 5.0,
    excl_res: int = 4
) -> list[np.ndarray]:
    """
    Calculate reference distances from multiple structure files.
    
    This function mirrors the behavior of drms_analysis.DRMSAnalyzer.
    It filters atom pairs based on distance differences across multiple reference states,
    only keeping pairs where the minimum difference between any two states is >= min_diff.
    Each state's distances are filtered independently for the distance range.
    
    Args:
        ref_files: List of reference structure files (pdb/gro)
        selected_atom: Atom selection string
        min_dist: Minimum distance in Angstrom
        max_dist: Maximum distance in Angstrom
        min_diff: Minimum difference between states (pairs with smaller differences are excluded)
        excl_res: Exclude residues within this separation
        
    Returns:
        List of reference distance arrays, one per state (n_pairs x 3): [atom_i, atom_j, distance]
    """
    import MDAnalysis as mda
    
    # Load all reference structures
    ref_universes = [mda.Universe(f) for f in ref_files]
    ref_atoms = [u.select_atoms(selected_atom) for u in ref_universes]
    
    # Validate that all references have the same atoms
    base_count = len(ref_atoms[0])
    for i, atoms in enumerate(ref_atoms[1:], 1):
        if len(atoms) != base_count:
            raise ValueError(f"Reference structure {i+1} has different atom count!")
        # Check segid and resid consistency
        for j in range(base_count):
            if (atoms[j].segid != ref_atoms[0][j].segid or 
                atoms[j].resid != ref_atoms[0][j].resid):
                raise ValueError(f"Atom {j} segid/resid mismatch in references")
    
    n_atoms = base_count
    n_states = len(ref_atoms)
    
    # Results for each state
    results = [[] for _ in range(n_states)]
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            atom_i, atom_j = ref_atoms[0][i], ref_atoms[0][j]
            
            # Check residue exclusion (same logic as drms_analysis)
            if atom_i.segid == atom_j.segid and atom_i.resid + excl_res > atom_j.resid:
                continue
            
            # Calculate distances for all states
            dists = []
            for atoms in ref_atoms:
                pos_i = atoms[i].position
                pos_j = atoms[j].position
                dist = np.linalg.norm(pos_i - pos_j)
                dists.append(dist)
            
            # Calculate minimum difference between any two states
            min_difference = min(abs(dists[k] - dists[l]) 
                                for k in range(n_states) 
                                for l in range(k + 1, n_states))
            
            # Only process pairs with sufficient difference between states
            if min_difference >= min_diff:
                # For each state, independently check distance range
                for state_idx in range(n_states):
                    dist = dists[state_idx]
                    if min_dist <= dist <= max_dist:
                        results[state_idx].append([atoms[i].ix, atoms[j].ix, dist])
    
    # Convert to numpy arrays
    return [np.array(state_data) if state_data else np.array([]) for state_data in results]


def _worker_process_chunk(
    netcdf_path: str,
    checkpoint_path: str,
    ref_distances_list: List[np.ndarray],
    frame_chunk_indices: List[int],
    replica_indices: np.ndarray,
    dt: float
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Worker function to process a chunk of frames and calculate dRMS.
    
    Creates its own PositionExtractor instance in the worker process to avoid
    pickling issues with file handles.
    
    Args:
        netcdf_path: Path to the NetCDF trajectory file
        checkpoint_path: Path to the checkpoint file
        ref_distances_list: List of reference distance arrays
        frame_chunk_indices: List of frame indices to process
        replica_indices: Array of replica indices to analyze
        dt: Time step between frames (ps)
        
    Returns:
        List of (times, drms) tuples for each reference state
    """
    chunk_start = time.time()
    extractor = PositionExtractor(netcdf_path, checkpoint_path)
    n_states = len(ref_distances_list)
    results = []
    
    try:
        # Prepare atom groups for batch reading
        # Each state needs two atom groups: atom_i and atom_j
        atom_groups = []
        for ref_distances in ref_distances_list:
            atom_groups.append(ref_distances[:, 0].astype(int))
            atom_groups.append(ref_distances[:, 1].astype(int))
        
        frame_indices = np.array(frame_chunk_indices)
        
        # Batch read all positions
        positions_groups = extractor.get_positions_for_atom_groups(
            atom_groups, frame_indices, replica_indices
        )
        
        # Process each state
        for state_idx in range(n_states):
            pos_i = positions_groups[2 * state_idx]      # (n_frames, n_replicas, n_pairs, 3)
            pos_j = positions_groups[2 * state_idx + 1]
            ref_dist = ref_distances_list[state_idx][:, 2]
            
            # Calculate distances: (n_frames, n_replicas, n_pairs)
            distances = np.linalg.norm(pos_i - pos_j, axis=-1)
            
            # Calculate dRMS: (n_frames, n_replicas)
            drms = np.sqrt(np.mean((distances - ref_dist) ** 2, axis=-1))
            
            # Convert frame indices to time
            times = frame_indices * dt
            
            results.append((times, drms))
    finally:
        extractor.close()
    
    chunk_time = time.time() - chunk_start
    print(f"  Frames {frame_chunk_indices[0]}-{frame_chunk_indices[-1]}: {chunk_time:.2f}s")
    return results


def calculate_trajectory_drms(
    netcdf_file: str,
    checkpoint_file: str,
    ref_distances: np.ndarray,
    replica_indices: Optional[np.ndarray] = None,
    skip: int = 1,
    num_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    auto_chunk: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate dRMS for entire trajectory using parallel processing.
    
    Args:
        netcdf_file: Path to NetCDF file
        checkpoint_file: Path to checkpoint file
        ref_distances: Reference distances
        replica_indices: Replica indices to analyze (None for all)
        skip: Process every N-th frame
        num_workers: Number of parallel workers (default: CPU count)
        chunk_size: Number of frames per worker task (auto-optimized if None)
        auto_chunk: Whether to automatically optimize chunk_size
        
    Returns:
        Tuple of (times, drms_array)
    """
    # Get metadata first (single extractor for metadata)
    extractor = PositionExtractor(netcdf_file, checkpoint_file)
    try:
        n_frames = extractor.n_frames
        dt = extractor.dt
        if replica_indices is None:
            replica_indices = np.arange(extractor.n_replicas)
    finally:
        extractor.close()
    
    # Frame indices to process
    frame_indices = np.arange(0, n_frames, skip)
    n_process_frames = len(frame_indices)
    
    if n_process_frames == 0:
        return np.array([]), np.array([])
    
    print(f"Processing {n_process_frames} frames (every {skip} of {n_frames})...")
    print(f"Analyzing {len(replica_indices)} replicas")
    
    # Determine number of workers
    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    num_workers = min(num_workers, n_process_frames)
    
    # Determine chunk_size
    n_pairs = len(ref_distances)
    if chunk_size is None and auto_chunk:
        chunk_size = _estimate_chunk_size(
            n_process_frames=n_process_frames,
            n_replicas=len(replica_indices),
            n_pairs=n_pairs,
            num_workers=num_workers,
        )
        mem_per_chunk = len(replica_indices) * n_pairs * 3 * 4 * chunk_size
        print(f"Auto-optimized chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    elif chunk_size is None:
        chunk_size = 100  # Default
        print(f"Using default chunk_size={chunk_size}")
    else:
        mem_per_chunk = len(replica_indices) * n_pairs * 3 * 4 * chunk_size
        print(f"Using specified chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    
    # Prepare tasks (chunk frames for parallel processing)
    tasks = []
    for i in range(0, n_process_frames, chunk_size):
        chunk_idx_end = min(i + chunk_size, n_process_frames)
        chunk_frame_indices = frame_indices[i:chunk_idx_end].tolist()
        tasks.append((
            netcdf_file,
            checkpoint_file,
            [ref_distances],  # Wrap in list for single state
            chunk_frame_indices,
            replica_indices,
            dt
        ))
    
    print(f"Using {num_workers} workers for {len(tasks)} chunks")
    
    start_time = time.time()
    
    # Parallel processing
    with multiprocessing.Pool(num_workers) as pool:
        chunk_results = pool.starmap(_worker_process_chunk, tasks)
    
    # Combine results
    all_times = np.concatenate([r[0][0] for r in chunk_results])
    all_drms = np.vstack([r[0][1] for r in chunk_results])
    
    elapsed = time.time() - start_time
    print(f"Total time: {elapsed:.1f}s ({n_process_frames/elapsed:.1f} fps)")
    
    return all_times, all_drms


def save_drms_results(
    times: np.ndarray,
    drms_data: np.ndarray,
    output_file: str,
    netcdf_file: str = "",
    checkpoint_file: str = "",
    replica_indices: Optional[np.ndarray] = None
):
    """Save dRMS results to file."""
    header_lines = [
        f"# time(ps), dRMS values for each replica",
        f"# NetCDF: {netcdf_file}",
        f"# Checkpoint: {checkpoint_file}",
    ]
    if replica_indices is not None:
        header_lines.append(f"# Replicas: {','.join(map(str, replica_indices))}")
    header = "\n".join(header_lines)
    
    data = np.column_stack([times, drms_data])
    
    fmt = ['%.1f'] + ['%.4f'] * drms_data.shape[1]
    np.savetxt(output_file, data, header=header, fmt=fmt)
    print(f"Saved dRMS results to: {output_file}")


def main():
    """Command-line interface for REMD dRMS analysis."""
    parser = argparse.ArgumentParser(
        description="Calculate REMD dRMS trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    parser.add_argument("-nc", "--netcdf", required=True,
                        help="NetCDF trajectory file")
    parser.add_argument("-c", "--checkpoint", default="output_checkpoint.nc",
                        help="Checkpoint NetCDF file")
    parser.add_argument("-ref", "--ref-files", nargs="+", required=True,
                        help="Reference structure files (pdb/gro)")
    
    # Output
    parser.add_argument("-o", "--output", default="dRMSTraj_nc",
                        help="Output prefix")
    
    # Analysis parameters
    parser.add_argument("--sel", "--selected-atom", default="name BB",
                        help="Atom selection string")
    parser.add_argument("--min-dist", type=float, default=6.0,
                        help="Minimum distance in Angstrom")
    parser.add_argument("--max-dist", type=float, default=50.0,
                        help="Maximum distance in Angstrom")
    parser.add_argument("--min-diff", type=float, default=5.0,
                        help="Minimum difference between references")
    parser.add_argument("--excl-res", type=int, default=4,
                        help="Residue exclusion distance")
    parser.add_argument("--skip", type=int, default=1,
                        help="Process every N-th frame")
    parser.add_argument("--num-workers", type=int, default=None,
                        help="Number of parallel workers (default: CPU count)")
    parser.add_argument("--chunk-size", type=int, default=None,
                        help="Frames per chunk (auto-optimized if not specified)")
    parser.add_argument("--auto-chunk", action="store_true", default=True,
                        help="Automatically optimize chunk size (default: enabled)")
    parser.add_argument("--no-auto-chunk", action="store_false", dest="auto_chunk",
                        help="Disable automatic chunk size optimization")
    parser.add_argument("--replicas", type=str, default="all",
                        help="Comma-separated replica indices or 'all'")
    
    args = parser.parse_args()
    
    # Calculate reference distances for each state
    print("\nCalculating reference distances...")
    states_str = 'ABCDEFGHIJKLMN'
    
    # Calculate all reference distances with multi-state filtering
    all_ref_distances = calculate_reference_distances(
        args.ref_files,
        selected_atom=args.sel,
        min_dist=args.min_dist,
        max_dist=args.max_dist,
        min_diff=args.min_diff,
        excl_res=args.excl_res
    )
    
    for i, ref_dist in enumerate(all_ref_distances):
        print(f"  State {states_str[i]}: {len(ref_dist)} distance pairs")
    
    # Parse replica indices
    if args.replicas == "all":
        replica_indices = None
    else:
        replica_indices = np.array([int(x) for x in args.replicas.split(",")])
    
    # Get metadata once
    extractor = PositionExtractor(args.netcdf, args.checkpoint)
    try:
        n_frames = extractor.n_frames
        dt = extractor.dt
        if replica_indices is None:
            replica_indices = np.arange(extractor.n_replicas)
    finally:
        extractor.close()
    
    frame_indices = np.arange(0, n_frames, args.skip)
    n_process_frames = len(frame_indices)
    
    # Determine workers and chunk size
    num_workers = args.num_workers if args.num_workers else multiprocessing.cpu_count()
    chunk_size = args.chunk_size
    
    print(f"\nProcessing {n_process_frames} frames (every {args.skip} of {n_frames})...")
    
    # Determine chunk_size
    if args.chunk_size is None and args.auto_chunk:
        # Use the maximum pairs across all states for estimation
        max_pairs = max(len(ref_dist) for ref_dist in all_ref_distances)
        chunk_size = _estimate_chunk_size(
            n_process_frames=n_process_frames,
            n_replicas=len(replica_indices),
            n_pairs=max_pairs,
            num_workers=num_workers,
        )
        mem_per_chunk = len(replica_indices) * max_pairs * 3 * 4 * chunk_size
        print(f"Auto-optimized chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    elif args.chunk_size is not None:
        chunk_size = args.chunk_size
        max_pairs = max(len(ref_dist) for ref_dist in all_ref_distances)
        mem_per_chunk = len(replica_indices) * max_pairs * 3 * 4 * chunk_size
        print(f"Using specified chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    else:
        chunk_size = 100  # Default fallback
        print(f"Using default chunk_size={chunk_size}")
    
    print(f"Using {num_workers} workers")
    
    # Prepare tasks for all states
    tasks = []
    for ref_dist in all_ref_distances:
        for i in range(0, n_process_frames, chunk_size):
            chunk_end = min(i + chunk_size, n_process_frames)
            chunk_frame_indices = frame_indices[i:chunk_end].tolist()
            tasks.append((
                args.netcdf,
                args.checkpoint,
                [ref_dist],
                chunk_frame_indices,
                replica_indices,
                dt
            ))
    
    # Parallel process all chunks
    start_time = time.time()
    with multiprocessing.Pool(num_workers) as pool:
        all_chunk_results = pool.starmap(_worker_process_chunk, tasks)
    
    # Organize results by state
    n_states = len(all_ref_distances)
    chunks_per_state = len(tasks) // n_states
    
    for state_idx in range(n_states):
        print(f"\nSaving State {states_str[state_idx]}...")
        state_chunks = all_chunk_results[
            state_idx * chunks_per_state : (state_idx + 1) * chunks_per_state
        ]
        
        times = np.concatenate([r[0][0] for r in state_chunks])
        drms_data = np.vstack([r[0][1] for r in state_chunks])
        
        output_file = f"{args.output}_State{states_str[state_idx]}.dat"
        save_drms_results(
            times, drms_data, output_file,
            args.netcdf, args.checkpoint, replica_indices
        )


if __name__ == "__main__":
    main()
