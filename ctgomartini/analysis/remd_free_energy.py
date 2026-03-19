#!/usr/bin/env python3
"""
REMD free energy analysis and comparison plotting.

Computes and visualizes free energy profiles from replica exchange simulations,
including comparisons with unbiased sampling data.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import matplotlib.pyplot as plt

if TYPE_CHECKING:
    from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer


def _import_openmmtools_analysis():
    """Lazy import openmmtools to avoid import-time warnings."""
    from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer
    return MultiStateReporter, ReplicaExchangeAnalyzer


def load_free_energy_data(
    output_file: str | Path,
    unbiased_file: str | Path | None = None,
) -> dict:
    """Load free energy data from REMD output.

    Args:
        output_file: Path to REMD output NetCDF file.
        unbiased_file: Optional path to unbiased sampling data.

    Returns:
        Dictionary containing free energy data.
    """
    MultiStateReporter, ReplicaExchangeAnalyzer = _import_openmmtools_analysis()
    reporter = MultiStateReporter(str(output_file), open_mode="r")
    analyzer = ReplicaExchangeAnalyzer(reporter)
    
    # Get free energy differences
    f_ij, df_ij = analyzer.get_free_energy_differences()
    reporter.close()
    
    data = {
        "f_ij": f_ij,
        "df_ij": df_ij,
        "n_states": f_ij.shape[0],
    }
    
    if unbiased_file and Path(unbiased_file).exists():
        data["unbiased"] = np.loadtxt(unbiased_file)
    
    return data


def plot_free_energy_comparison(
    data_files: list[str | Path],
    labels: list[str] | None = None,
    unbiased_file: str | Path | None = None,
    output_file: str | Path = "free_energy_comparison.pdf",
    colormap: np.ndarray | None = None,
) -> None:
    """Plot free energy comparison across multiple REMD runs.

    Args:
        data_files: List of REMD output NetCDF files.
        labels: List of labels for each dataset.
        unbiased_file: Optional unbiased sampling data file.
        output_file: Output plot file path.
        colormap: Custom colormap array.
    """
    if colormap is None:
        colormap = np.array([
            [128, 128, 128],
            [144, 191, 249],
            [249, 100, 149],
            [109, 219, 0],
        ]) / 255

    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i, data_file in enumerate(data_files):
        data = load_free_energy_data(data_file)
        f_ij = data["f_ij"]
        df_ij = data["df_ij"]
        
        # Calculate free energy relative to first state
        fe = f_ij[0, :] - f_ij[0, 0]
        err = df_ij[0, :]
        
        x = np.arange(len(fe))
        label = labels[i] if labels else f"Run {i+1}"
        color = colormap[i % len(colormap)]
        
        ax.errorbar(x, fe, yerr=err, fmt='o-', color=color, 
                   label=label, capsize=3)
    
    # Add unbiased sampling if available
    if unbiased_file and Path(unbiased_file).exists():
        unbiased = np.loadtxt(unbiased_file)
        ax.plot(unbiased[:, 0], unbiased[:, 1], 'k--', 
               label='Unbiased Sampling', linewidth=2)
    
    ax.set_xlabel("State Index", fontsize=14)
    ax.set_ylabel("Free Energy (kJ/mol)", fontsize=14)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_file}")


def plot_free_energy_profile(
    output_file: str | Path,
    output_plot: str | Path = "free_energy_profile.pdf",
    xlabel: str = "Reaction Coordinate",
) -> None:
    """Plot free energy profile from single REMD run.

    Args:
        output_file: Path to REMD output NetCDF file.
        output_plot: Output plot file path.
        xlabel: Label for x-axis.
    """
    data = load_free_energy_data(output_file)
    f_ij = data["f_ij"]
    df_ij = data["df_ij"]
    n_states = data["n_states"]
    
    # Free energy relative to first state
    fe = f_ij[0, :] - f_ij[0, 0]
    err = df_ij[0, :]
    
    x = np.linspace(-500, -100, n_states) if n_states == 11 else np.arange(n_states)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(x, fe, yerr=err, fmt='o-', capsize=5, linewidth=2, markersize=8)
    ax.fill_between(x, fe - err, fe + err, alpha=0.3)
    
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Free Energy (kJ/mol)", fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_plot}")


def main():
    """Command-line interface for REMD free energy analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze and plot REMD free energy profiles"
    )
    parser.add_argument(
        "-f", "--files",
        nargs="+",
        required=True,
        help="REMD output NetCDF file(s)"
    )
    parser.add_argument(
        "-o", "--output",
        default="free_energy_comparison.pdf",
        help="Output plot file"
    )
    parser.add_argument(
        "-l", "--labels",
        nargs="+",
        help="Labels for each dataset"
    )
    parser.add_argument(
        "-u", "--unbiased",
        help="Unbiased sampling data file for comparison"
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="Plot single free energy profile instead of comparison"
    )
    
    args = parser.parse_args()
    
    if args.single:
        plot_free_energy_profile(
            output_file=args.files[0],
            output_plot=args.output,
        )
    else:
        plot_free_energy_comparison(
            data_files=args.files,
            labels=args.labels,
            unbiased_file=args.unbiased,
            output_file=args.output,
        )


if __name__ == "__main__":
    main()
