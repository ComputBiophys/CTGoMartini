#!/usr/bin/env python3
"""
REMD replica state trajectory analysis.

Analyzes replica state indices over time from replica exchange
molecular dynamics simulations and saves results to files.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from openmm import unit


kB = (unit.MOLAR_GAS_CONSTANT_R).in_units_of(
    unit.kilojoule / (unit.kelvin * unit.mole)
).value_in_unit(unit.kilojoule / (unit.kelvin * unit.mole))


def _import_openmmtools_state():
    """Lazy import openmmtools to avoid import-time warnings."""
    from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer
    return MultiStateReporter, ReplicaExchangeAnalyzer


def _get_dt_from_netcdf(output_file: str | Path) -> float | None:
    """Try to get time step (dt) from NetCDF file.
    
    Args:
        output_file: Path to REMD output NetCDF file.
        
    Returns:
        Time step in ps, or None if cannot be determined.
    """
    try:
        from openmmtools.multistate import MultiStateReporter
        import openmm
        
        reporter = MultiStateReporter(str(output_file), open_mode="r")
        try:
            mcmove = reporter.read_mcmc_moves()[0]
            time_interval = mcmove.n_steps * mcmove.timestep
            dt_per_move = time_interval.value_in_unit(openmm.unit.picosecond)
            # Multiply by checkpoint_interval to get actual time between frames
            dt = dt_per_move * reporter.checkpoint_interval
            reporter.close()
            return dt
        except Exception:
            reporter.close()
            return None
    except Exception:
        return None


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


def save_replica_states(
    output_file: str | Path,
    output_data: str | Path = "replica_states.dat",
    dt: float | None = None,
) -> None:
    """Save replica state trajectories to file.

    Args:
        output_file: Path to REMD output NetCDF file.
        output_data: Path for output data file.
        dt: Time step in picoseconds (ps). If None, will try to auto-detect
            from NetCDF, otherwise defaults to 1.0 ps.
    """
    # Try to auto-detect dt from NetCDF
    dt_detected = _get_dt_from_netcdf(output_file)
    
    if dt is not None:
        # User specified dt, use it but warn if different from detected
        if dt_detected is not None and abs(dt - dt_detected) > 0.01:
            print(f"Warning: Specified dt={dt} ps differs from detected dt={dt_detected:.2f} ps")
        dt_use = dt
    else:
        # Use detected dt or default
        if dt_detected is not None:
            dt_use = dt_detected
            print(f"Auto-detected time step: {dt_use:.2f} ps")
        else:
            dt_use = 1.0
            print(f"Warning: Could not auto-detect time step from NetCDF. Using default: {dt_use:.2f} ps")
    
    replica_state_indices = load_replica_states(output_file)
    n_replica = replica_state_indices.shape[0]
    
    # Calculate time
    n_frames = replica_state_indices.shape[1]
    time = np.arange(0, n_frames * dt_use, dt_use)

    # Save to file
    header = "Time(ps)" + "".join([f"\tReplica_{i}" for i in range(n_replica)])
    data = np.column_stack([time, replica_state_indices.T])
    np.savetxt(output_data, data, header=header, comments="", fmt="%.6f")
    print(f"Saved replica state data to: {output_data}")


def compute_state_occupancies(
    output_file: str | Path, 
    output_file_path: str | Path | None = None
) -> dict:
    """Compute state occupancies for each replica.

    Args:
        output_file: Path to REMD output NetCDF file.
        output_file_path: Optional path to save occupancy results.

    Returns:
        Dictionary with occupancy statistics.
    """
    replica_state_indices = load_replica_states(output_file)
    n_replicas, n_frames = replica_state_indices.shape
    
    occupancies = {}
    lines = ["# State Occupancy Statistics\n"]
    lines.append(f"# Total frames: {n_frames}\n")
    lines.append("#\n")
    
    for i in range(n_replicas):
        states, counts = np.unique(replica_state_indices[i], return_counts=True)
        occupancies[f"replica_{i}"] = {
            "states": states.tolist(),
            "counts": counts.tolist(),
            "fractions": (counts / n_frames).tolist(),
        }
        
        lines.append(f"# Replica {i}:\n")
        for state, count, frac in zip(states, counts, counts / n_frames):
            lines.append(f"{i}\t{int(state)}\t{count}\t{frac:.6f}\n")
    
    if output_file_path is not None:
        with open(output_file_path, "w") as f:
            f.writelines(lines)
        print(f"Saved occupancy statistics to: {output_file_path}")
    
    return occupancies


def main():
    """Command-line interface for REMD replica state analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze REMD replica state trajectories and save results"
    )
    parser.add_argument(
        "-f", "--file",
        default="output.nc",
        help="REMD output NetCDF file (default: output.nc)"
    )
    parser.add_argument(
        "-o", "--output",
        default="replica_states.dat",
        help="Output data file (default: replica_states.dat)"
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=None,
        help="Time step in picoseconds (default: auto-detect, fallback to 1.0)"
    )
    parser.add_argument(
        "--occupancies",
        action="store_true",
        help="Compute and save state occupancy statistics"
    )
    parser.add_argument(
        "--occupancy-output",
        default="state_occupancies.dat",
        help="Output file for occupancy statistics (default: state_occupancies.dat)"
    )
    
    args = parser.parse_args()
    
    if args.occupancies:
        compute_state_occupancies(args.file, output_file_path=args.occupancy_output)
    else:
        save_replica_states(
            output_file=args.file,
            output_data=args.output,
            dt=args.dt,
        )


if __name__ == "__main__":
    main()
