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


def calculate_drms_for_frame(
    frame_idx: int,
    replica_indices: np.ndarray,
    extractor: PositionExtractor,
    ref_distances: np.ndarray
) -> Tuple[int, np.ndarray]:
    """
    Calculate dRMS for a single frame across specified replicas.
    
    Args:
        frame_idx: Frame index
        replica_indices: Array of replica indices
        extractor: PositionExtractor instance
        ref_distances: Reference distances array
        
    Returns:
        Tuple of (frame_idx, drms_values)
    """
    drms_values = np.zeros(len(replica_indices))
    
    atom_i = ref_distances[:, 0].astype(int)
    atom_j = ref_distances[:, 1].astype(int)
    ref_dist = ref_distances[:, 2]
    
    for r_idx, replica_idx in enumerate(replica_indices):
        positions = extractor.get_positions(frame_idx, replica_idx)
        
        # Calculate distances
        distances = np.linalg.norm(positions[atom_i] - positions[atom_j], axis=1)
        
        # Calculate dRMS
        drms = np.sqrt(np.mean((distances - ref_dist) ** 2))
        drms_values[r_idx] = drms
    
    return frame_idx, drms_values


def calculate_trajectory_drms(
    netcdf_file: str,
    checkpoint_file: str,
    ref_distances: np.ndarray,
    replica_indices: Optional[np.ndarray] = None,
    skip: int = 1,
    num_workers: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate dRMS for entire trajectory.
    
    Args:
        netcdf_file: Path to NetCDF file
        checkpoint_file: Path to checkpoint file
        ref_distances: Reference distances
        replica_indices: Replica indices to analyze (None for all)
        skip: Process every N-th frame
        num_workers: Number of parallel workers
        
    Returns:
        Tuple of (times, drms_array)
    """
    extractor = PositionExtractor(netcdf_file, checkpoint_file)
    
    try:
        n_frames = extractor.n_frames
        
        if replica_indices is None:
            replica_indices = np.arange(extractor.n_replicas)
        
        # Frame indices to process
        frame_indices = np.arange(0, n_frames, skip)
        n_process_frames = len(frame_indices)
        
        print(f"Processing {n_process_frames} frames (every {skip} of {n_frames})...")
        print(f"Analyzing {len(replica_indices)} replicas")
        
        # Sequential processing (parallel requires pickling)
        drms_results = np.zeros((n_process_frames, len(replica_indices)))
        
        start_time = time.time()
        
        for i, frame_idx in enumerate(frame_indices):
            _, drms = calculate_drms_for_frame(
                frame_idx, replica_indices, extractor, ref_distances
            )
            drms_results[i] = drms
            
            if (i + 1) % 100 == 0 or i == n_process_frames - 1:
                elapsed = time.time() - start_time
                fps = (i + 1) / elapsed
                remaining = (n_process_frames - i - 1) / fps
                print(f"  Processed {i+1}/{n_process_frames} frames "
                      f"({fps:.1f} fps, {remaining:.0f}s remaining)")
        
        # Convert frame indices to time
        times = frame_indices * extractor.dt
        
        return times, drms_results
        
    finally:
        extractor.close()


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
                        help="Number of parallel workers")
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
        
        # Calculate dRMS for each state
        for i, ref_dist in enumerate(all_ref_distances):
            print(f"\nCalculating dRMS for State {states_str[i]}...")
            times, drms_data = calculate_trajectory_drms(
                args.netcdf,
                args.checkpoint,
                ref_dist,
                replica_indices=replica_indices,
                skip=args.skip
            )
            
            output_file = f"{args.output}_State{states_str[i]}.dat"
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
