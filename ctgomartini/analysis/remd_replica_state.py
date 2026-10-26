#!/usr/bin/env python3
"""REMD replica state trajectory analysis.

Analyzes replica state indices over time from replica exchange
molecular dynamics simulations (OpenMM >= 8.6 built-in backend output
directory) and saves results to files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _get_dt_from_metadata(output_dir: str | Path) -> float | None:
    """Try to get the time step between logged frames from metadata.json.

    The built-in REMD reporter logs one row every ``reportInterval``
    iterations of ``stepsPerIteration`` MD steps each, so the time between
    frames is ``dt * stepsPerIteration * reportInterval``.

    Args:
        output_dir: Path to the built-in REMD output directory.

    Returns:
        Time step between frames in ps, or None if it cannot be determined.
    """
    try:
        metadata = json.loads((Path(output_dir) / "metadata.json").read_text())
        dt = float(metadata["dt"])
        steps_per_iteration = float(
            metadata.get("stepsPerIteration", metadata["exc_freq"])
        )
        report_interval = float(metadata["reportInterval"])
        return dt * steps_per_iteration * report_interval
    except Exception:
        return None


def load_replica_states(output_dir: str | Path) -> np.ndarray:
    """Load replica state indices from a built-in REMD output directory.

    Reads ``log.csv`` (columns ``Iteration,Step,Replica_0_State,...``).

    Args:
        output_dir: Path to the built-in REMD output directory.

    Returns:
        Array of replica state indices with shape (n_replicas, n_frames).

    Raises:
        ValueError: If output_dir is not a directory (openmmtools NetCDF
            output is no longer supported).
        FileNotFoundError: If log.csv is missing.
    """
    directory = Path(output_dir)
    if not directory.is_dir():
        raise ValueError(
            f"load_replica_states requires a built-in REMD output directory "
            f"(containing log.csv); got {output_dir}. openmmtools NetCDF "
            f"(.nc) output is no longer supported."
        )
    log_csv = directory / "log.csv"
    if not log_csv.is_file():
        raise FileNotFoundError(
            f"Cannot load replica states: {log_csv} does not exist."
        )
    log = np.loadtxt(log_csv, delimiter=",", skiprows=1, ndmin=2)
    # Columns 2: hold Replica_i_State, one row per reported iteration.
    return log[:, 2:].astype(int).T


def save_replica_states(
    output_dir: str | Path,
    output_data: str | Path = "replica_states.dat",
    dt: float | None = None,
) -> None:
    """Save replica state trajectories to file.

    Args:
        output_dir: Path to the built-in REMD output directory.
        output_data: Path for output data file.
        dt: Time step in picoseconds (ps) between logged frames. If None,
            auto-detect from metadata.json, otherwise default to 1.0 ps.
    """
    dt_detected = _get_dt_from_metadata(output_dir)

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
            print(
                "Warning: Could not auto-detect time step from metadata.json. "
                f"Using default: {dt_use:.2f} ps"
            )

    replica_state_indices = load_replica_states(output_dir)
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
    output_dir: str | Path,
    output_file_path: str | Path | None = None
) -> dict:
    """Compute state occupancies for each replica.

    Args:
        output_dir: Path to the built-in REMD output directory.
        output_file_path: Optional path to save occupancy results.

    Returns:
        Dictionary with occupancy statistics.
    """
    replica_state_indices = load_replica_states(output_dir)
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
        for state, count, frac in zip(states, counts, counts / n_frames, strict=True):
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
        default="output_remd",
        help="REMD output directory (default: output_remd)"
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
            output_dir=args.file,
            output_data=args.output,
            dt=args.dt,
        )


if __name__ == "__main__":
    main()
