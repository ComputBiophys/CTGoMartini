#!/usr/bin/env python3
"""
REMD dRMS trajectory analysis and plotting.

Analyzes distance root-mean-square (dRMS) trajectories from replica exchange
molecular dynamics simulations and generates publication-quality plots.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def plot_drms_trajectory(
    data_file: str | Path,
    output_file: str | Path = "drms_trajectory.pdf",
    dt: float = 0.005,
    skip: int = 10,
    colormap: list | None = None,
    figsize: tuple | None = None,
) -> None:
    """Plot dRMS trajectory for all replicas.

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

    # Create figure with vertically stacked subplots
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
        description="Analyze and plot REMD dRMS trajectories"
    )
    parser.add_argument(
        "-f", "--file",
        default="dRMStraj_nc_StateA.dat",
        help="Input dRMS trajectory file (default: dRMStraj_nc_StateA.dat)"
    )
    parser.add_argument(
        "-o", "--output",
        default="drms_trajectory.pdf",
        help="Output plot file (default: drms_trajectory.pdf)"
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
        default=10,
        help="Frame skipping interval (default: 10)"
    )
    
    args = parser.parse_args()
    
    plot_drms_trajectory(
        data_file=args.file,
        output_file=args.output,
        dt=args.dt,
        skip=args.skip,
    )


if __name__ == "__main__":
    main()
