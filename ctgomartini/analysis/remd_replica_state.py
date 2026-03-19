#!/usr/bin/env python3
"""
REMD replica state trajectory analysis and plotting.

Analyzes and visualizes replica state indices over time from replica exchange
molecular dynamics simulations.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from openmm import unit


kB = (unit.MOLAR_GAS_CONSTANT_R).in_units_of(
    unit.kilojoule / (unit.kelvin * unit.mole)
).value_in_unit(unit.kilojoule / (unit.kelvin * unit.mole))


def _import_openmmtools_state():
    """Lazy import openmmtools to avoid import-time warnings."""
    from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer
    return MultiStateReporter, ReplicaExchangeAnalyzer


def load_replica_states(output_file: str | Path) -> np.ndarray:
    """Load replica state indices from REMD output.

    Args:
        output_file: Path to REMD output NetCDF file.

    Returns:
        Array of replica state indices with shape (n_replicas, n_timesteps).
    """
    MultiStateReporter, ReplicaExchangeAnalyzer = _import_openmmtools_state()
    reporter = MultiStateReporter(str(output_file), open_mode="r")
    analyzer = ReplicaExchangeAnalyzer(reporter)
    
    # Read energies and states
    _, _, _, replica_state_indices = analyzer.read_energies()
    reporter.close()
    
    return replica_state_indices


def plot_replica_states(
    output_file: str | Path,
    output_plot: str | Path = "replica_states.pdf",
    dt: float = 0.005,
    skip: int = 100,
    figsize: tuple | None = None,
) -> None:
    """Plot replica state trajectories.

    Args:
        output_file: Path to REMD output NetCDF file.
        output_plot: Path for output plot file.
        dt: Time step in microseconds.
        skip: Frame skipping interval for plotting.
        figsize: Figure size tuple.
    """
    replica_state_indices = load_replica_states(output_file)
    n_replica = replica_state_indices.shape[0]
    
    # Calculate time
    n_frames = replica_state_indices.shape[1]
    time = np.arange(0, n_frames * dt, dt * skip)[: n_frames // skip]

    # Create figure
    if figsize is None:
        figsize = (18, 1 * n_replica)
    
    fig, axes = plt.subplots(n_replica, 1, figsize=figsize, sharex=True)
    fig.subplots_adjust(hspace=0)

    if n_replica == 1:
        axes = [axes]

    colormap = plt.cm.viridis(np.linspace(0, 1, n_replica))

    for i, ax in enumerate(axes):
        states = replica_state_indices[i, ::skip]
        ax.scatter(time, states, c=colormap[states.astype(int) % n_replica], 
                   s=1, alpha=0.5)
        ax.set_ylabel(f"R{i}", fontsize=10)
        ax.set_ylim(-0.5, n_replica - 0.5)
        ax.set_yticks(range(n_replica))

    axes[-1].set_xlabel("Time (μs)", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_plot}")


def compute_state_occupancies(output_file: str | Path) -> dict:
    """Compute state occupancies for each replica.

    Args:
        output_file: Path to REMD output NetCDF file.

    Returns:
        Dictionary with occupancy statistics.
    """
    replica_state_indices = load_replica_states(output_file)
    n_replicas, n_frames = replica_state_indices.shape
    
    occupancies = {}
    for i in range(n_replicas):
        states, counts = np.unique(replica_state_indices[i], return_counts=True)
        occupancies[f"replica_{i}"] = {
            "states": states.tolist(),
            "counts": counts.tolist(),
            "fractions": (counts / n_frames).tolist(),
        }
    
    return occupancies


def main():
    """Command-line interface for REMD replica state analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze and plot REMD replica state trajectories"
    )
    parser.add_argument(
        "-f", "--file",
        default="output.nc",
        help="REMD output NetCDF file (default: output.nc)"
    )
    parser.add_argument(
        "-o", "--output",
        default="replica_states.pdf",
        help="Output plot file (default: replica_states.pdf)"
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=0.005,
        help="Time step in microseconds (default: 0.005)"
    )
    parser.add_argument(
        "--skip",
        type=int,
        default=100,
        help="Frame skipping interval (default: 100)"
    )
    parser.add_argument(
        "--occupancies",
        action="store_true",
        help="Print state occupancy statistics"
    )
    
    args = parser.parse_args()
    
    if args.occupancies:
        occ = compute_state_occupancies(args.file)
        for replica, stats in occ.items():
            print(f"\n{replica}:")
            for state, frac in zip(stats["states"], stats["fractions"]):
                print(f"  State {int(state)}: {frac:.3f}")
    else:
        plot_replica_states(
            output_file=args.file,
            output_plot=args.output,
            dt=args.dt,
            skip=args.skip,
        )


if __name__ == "__main__":
    main()
