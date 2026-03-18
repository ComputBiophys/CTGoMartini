#!/usr/bin/env python3
"""
MBAR parameter optimization analysis for REMD simulations.

Analyzes MBAR (Multistate Bennett Acceptance Ratio) parameter optimization
results including selected states, start ratios, and G-values.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def plot_selected_state_analysis(
    results_file: str | Path,
    output_file: str | Path = "selected_state_analysis.pdf",
) -> None:
    """Plot MBAR analysis for different numbers of selected states.

    Args:
        results_file: Path to pickle file with selected state results.
        output_file: Output plot file path.
    """
    with open(results_file, 'rb') as f:
        results = pickle.load(f)
    
    n_selected_states = list(results.keys())
    fe_values = [results[k]["free_energy"] for k in n_selected_states]
    errors = [results[k].get("error", 0) for k in n_selected_states]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(n_selected_states, fe_values, yerr=errors, fmt='o-', capsize=5)
    ax.set_xlabel("Number of Selected States", fontsize=14)
    ax.set_ylabel("Free Energy (kJ/mol)", fontsize=14)
    ax.set_title("MBAR Convergence with Selected States", fontsize=16)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_file}")


def plot_start_ratio_analysis(
    results_file: str | Path,
    output_file: str | Path = "start_ratio_analysis.pdf",
) -> None:
    """Plot MBAR analysis for different starting ratios.

    Args:
        results_file: Path to pickle file with start ratio results.
        output_file: Output plot file path.
    """
    with open(results_file, 'rb') as f:
        results = pickle.load(f)
    
    start_ratios = np.array(list(results.keys()))
    fe_values = np.array([results[k]["free_energy"] for k in start_ratios])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(start_ratios, fe_values, 'o-', linewidth=2, markersize=8)
    ax.set_xlabel("Start Ratio", fontsize=14)
    ax.set_ylabel("Free Energy (kJ/mol)", fontsize=14)
    ax.set_title("MBAR Convergence with Start Ratio", fontsize=16)
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_file}")


def plot_g_values_analysis(
    results_file: str | Path,
    output_file: str | Path = "g_values_analysis.pdf",
) -> None:
    """Plot G-values analysis from MBAR optimization.

    Args:
        results_file: Path to pickle file with G-values results.
        output_file: Output plot file path.
    """
    with open(results_file, 'rb') as f:
        results = pickle.load(f)
    
    g_values = np.array(list(results.keys()))
    fe_values = np.array([results[k]["free_energy"] for k in g_values])
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g_values, fe_values, 'o-', linewidth=2, markersize=8, color='green')
    ax.set_xlabel("G Value", fontsize=14)
    ax.set_ylabel("Free Energy (kJ/mol)", fontsize=14)
    ax.set_title("MBAR Convergence with G Values", fontsize=16)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved plot to: {output_file}")


def compute_keq_from_free_energy(
    delta_g: float,
    temperature: float = 310.0,
) -> float:
    """Compute equilibrium constant from free energy difference.

    Args:
        delta_g: Free energy difference in kJ/mol.
        temperature: Temperature in Kelvin.

    Returns:
        Equilibrium constant K_eq.
    """
    R = 8.314  # J/(mol*K)
    delta_g_j = delta_g * 1000  # Convert to J/mol
    keq = np.exp(-delta_g_j / (R * temperature))
    return keq


def main():
    """Command-line interface for MBAR parameter analysis."""
    parser = argparse.ArgumentParser(
        description="Analyze MBAR parameter optimization results"
    )
    parser.add_argument(
        "-t", "--type",
        choices=["selected_states", "start_ratio", "g_values"],
        required=True,
        help="Type of analysis to perform"
    )
    parser.add_argument(
        "-f", "--file",
        required=True,
        help="Path to pickle results file"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output plot file (default: auto-generated based on type)"
    )
    parser.add_argument(
        "--keq",
        action="store_true",
        help="Also compute and print equilibrium constant"
    )
    parser.add_argument(
        "--temp",
        type=float,
        default=310.0,
        help="Temperature for K_eq calculation (default: 310 K)"
    )
    
    args = parser.parse_args()
    
    # Set default output filename
    if args.output is None:
        args.output = f"{args.type}_analysis.pdf"
    
    if args.type == "selected_states":
        plot_selected_state_analysis(args.file, args.output)
    elif args.type == "start_ratio":
        plot_start_ratio_analysis(args.file, args.output)
    elif args.type == "g_values":
        plot_g_values_analysis(args.file, args.output)
    
    # Compute K_eq if requested
    if args.keq:
        with open(args.file, 'rb') as f:
            results = pickle.load(f)
        # Use first key's free energy
        first_key = list(results.keys())[0]
        delta_g = results[first_key]["free_energy"]
        keq = compute_keq_from_free_energy(delta_g, args.temp)
        print(f"\nFree Energy: {delta_g:.2f} kJ/mol")
        print(f"Equilibrium Constant (K_eq) at {args.temp}K: {keq:.4e}")
        print(f"ln(K_eq): {np.log(keq):.4f}")


if __name__ == "__main__":
    main()
