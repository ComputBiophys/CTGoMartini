#!/usr/bin/env python3
"""
REMD dRMS trajectory analysis and plotting.

Analyzes distance root-mean-square (dRMS) trajectories from replica exchange
molecular dynamics simulations and generates publication-quality plots.

This module provides both calculation (from NetCDF) and plotting functionality.

Author: Song Yang
Date: 2025
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
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

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
                self._dt = time_interval.value_in_unit(
                    __import__('openmm').unit.picosecond
                )
            except Exception:
                self._dt = 5.0  # Default: 5 ps per frame
        except Exception as e:
            print(f"Error initializing metadata: {e}")
            self._n_frames = self._n_replicas = self._n_atoms = 0
            self._dt = 5.0
    
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
    ref_file: str,
    selected_atom: str = "name BB",
    min_dist: float = 6.0,
    max_dist: float = 50.0,
    min_diff: float = 5.0,
    excl_res: int = 4
) -> np.ndarray:
    """
    Calculate reference distances from a structure file.
    
    Args:
        ref_file: Reference structure file (pdb/gro)
        selected_atom: Atom selection string
        min_dist: Minimum distance in Angstrom
        max_dist: Maximum distance in Angstrom
        min_diff: Minimum difference between states
        excl_res: Exclude residues within this separation
        
    Returns:
        Array of reference distances (n_pairs x 3): [atom_i, atom_j, distance]
    """
    import MDAnalysis as mda
    
    u = mda.Universe(ref_file)
    atoms = u.select_atoms(selected_atom)
    n_atoms = len(atoms)
    
    positions = atoms.positions
    resids = atoms.resids
    
    results = []
    
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            # Check residue exclusion
            if abs(resids[j] - resids[i]) < excl_res:
                continue
            
            # Calculate distance
            dist = np.linalg.norm(positions[i] - positions[j])
            
            # Check distance range
            if dist < min_dist or dist > max_dist:
                continue
            
            results.append([atoms[i].ix, atoms[j].ix, dist])
    
    return np.array(results)


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
    chunk_size: int = 100
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
        chunk_size: Number of frames per worker task
        
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
    
    print(f"Using {num_workers} workers for {len(tasks)} chunks (chunk_size={chunk_size})")
    
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


def plot_drms_trajectory(
    data_file: str | Path,
    output_file: str | Path = "drms_trajectory.pdf",
    dt: float = 0.005,
    skip: int = 10,
    colormap: list | None = None,
    figsize: tuple | None = None,
) -> None:
    """
    Plot dRMS trajectory for all replicas.
    
    Args:
        data_file: Path to dRMS trajectory data file.
        output_file: Path for output plot file.
        dt: Time step in microseconds.
        skip: Frame skipping interval for plotting.
        colormap: Custom colormap for replicas.
        figsize: Figure size tuple (width, height).
    """
    data = np.loadtxt(data_file)
    
    # Default colormap
    if colormap is None:
        colormap = [
            np.array([255, 81, 81]) / 255,
            np.array([0, 162, 210]) / 255,
        ]

    # Get number of replicas
    n_replica = data.shape[1] - 1
    
    # Calculate time in microseconds
    time = data[::skip, 0] * dt

    # Create figure
    if figsize is None:
        figsize = (18, 2 * n_replica)
    
    fig, axes = plt.subplots(n_replica, 1, figsize=figsize, sharex=True)
    fig.subplots_adjust(hspace=0)

    if n_replica == 1:
        axes = [axes]

    labels = [f"Replica {i}" for i in range(n_replica)]

    for i, ax in enumerate(axes):
        color = colormap[i % len(colormap)]
        ax.fill_between(time, data[::skip, i + 1], color=color, alpha=0.5)
        ax.plot(time, data[::skip, i + 1], color=color, linewidth=0.5)
        ax.set_ylabel(labels[i], fontsize=12)
        ax.set_ylim(0, max(data[:, i + 1]) * 1.1)
        ax.set_xlim(0, time[-1])

    axes[-1].set_xlabel("Time (μs)", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_file}")


def main():
    """Command-line interface for REMD dRMS analysis."""
    parser = argparse.ArgumentParser(
        description="Calculate and plot REMD dRMS trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Mode selection
    parser.add_argument(
        "--mode",
        choices=["calculate", "plot", "both"],
        default="both",
        help="Analysis mode: calculate dRMS, plot only, or both"
    )
    
    # Input files for calculation
    parser.add_argument("-nc", "--netcdf", help="NetCDF trajectory file")
    parser.add_argument("-c", "--checkpoint", default="output_checkpoint.nc",
                        help="Checkpoint NetCDF file")
    parser.add_argument("-ref", "--ref-files", nargs="+",
                        help="Reference structure files (pdb/gro)")
    
    # Input for plotting
    parser.add_argument("-f", "--file", default="dRMStraj_nc_StateA.dat",
                        help="Input dRMS data file (for plot mode)")
    
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
    parser.add_argument("--chunk-size", type=int, default=100,
                        help="Frames per chunk for parallel processing")
    parser.add_argument("--replicas", type=str, default="all",
                        help="Comma-separated replica indices or 'all'")
    
    # Plotting parameters
    parser.add_argument("--dt", type=float, default=0.005,
                        help="Time step in microseconds (for plotting)")
    parser.add_argument("--plot-skip", type=int, default=10,
                        help="Frame skipping for plotting")
    
    args = parser.parse_args()
    
    if args.mode in ["calculate", "both"]:
        if not args.netcdf or not args.ref_files:
            print("Error: --netcdf and --ref-files required for calculation mode")
            sys.exit(1)
        
        # Calculate reference distances for each state
        print("\nCalculating reference distances...")
        all_ref_distances = []
        states_str = 'ABCDEFGHIJKLMN'
        
        for i, ref_file in enumerate(args.ref_files):
            ref_dist = calculate_reference_distances(
                ref_file,
                selected_atom=args.sel,
                min_dist=args.min_dist,
                max_dist=args.max_dist,
                min_diff=args.min_diff,
                excl_res=args.excl_res
            )
            all_ref_distances.append(ref_dist)
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
        print(f"Using {num_workers} workers, chunk_size={chunk_size}")
        
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
    
    if args.mode in ["plot", "both"]:
        # Plot for each state file
        states_str = 'ABCDEFGHIJKLMN'
        # Determine number of states to plot
        if args.mode == "plot" and args.ref_files is None:
            # Try to find existing state files
            n_states = 0
            for i in range(len(states_str)):
                if os.path.exists(f"{args.output}_State{states_str[i]}.dat"):
                    n_states += 1
                else:
                    break
            if n_states == 0:
                # Try default state files
                for i in range(len(states_str)):
                    if os.path.exists(f"dRMSTraj_nc_State{states_str[i]}.dat"):
                        n_states += 1
                    else:
                        break
        else:
            n_states = len(args.ref_files) if args.ref_files else 1
        
        for i in range(n_states):
            data_file = f"{args.output}_State{states_str[i]}.dat"
            if os.path.exists(data_file):
                plot_file = f"{args.output}_State{states_str[i]}.pdf"
                print(f"\nPlotting {data_file}...")
                plot_drms_trajectory(
                    data_file, plot_file,
                    dt=args.dt, skip=args.plot_skip
                )


if __name__ == "__main__":
    main()
