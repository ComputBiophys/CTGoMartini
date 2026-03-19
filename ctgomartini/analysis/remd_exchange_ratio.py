#!/usr/bin/env python3
"""
REMD exchange ratio analysis.

Calculates and reports replica exchange ratios from REMD simulation output.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
from pathlib import Path

from openmm import unit

kB = (unit.MOLAR_GAS_CONSTANT_R).in_units_of(
    unit.kilojoule / (unit.kelvin * unit.mole)
)


def _import_openmmtools_reporter():
    """Lazy import openmmtools to avoid import-time warnings."""
    from openmmtools.multistate import MultiStateReporter
    return MultiStateReporter


def ReportExchangeRatio(output_file: str | Path = 'output.nc') -> dict:
    """Calculate and report replica exchange ratios.

    Args:
        output_file: Path to REMD output NetCDF file.

    Returns:
        Dictionary with exchange ratios for each replica pair.
    """
    MultiStateReporter = _import_openmmtools_reporter()
    
    reporter = MultiStateReporter(str(output_file), open_mode="r")
    n_accepted_matrix, n_proposed_matrix = reporter.read_mixing_statistics()
    
    # Because use the swap-neighbors method to exchange states, 
    # only half of the exchanges happened for state i and i+1
    ratio = n_accepted_matrix.sum(axis=0) / (n_accepted_matrix.shape[0] / 2)
    
    results = {}
    for i in range(ratio.shape[0] - 1):
        j = i + 1
        results[f"{i}-{j}"] = float(ratio[i, j])
        print(f'Exchange ratio between replica {i} and {j}: {ratio[i, j]:.4f}')
    
    reporter.close()
    return results


def main():
    """Command-line interface for exchange ratio analysis."""
    parser = argparse.ArgumentParser(
        description="Calculate REMD replica exchange ratios"
    )
    parser.add_argument(
        "-f", "--file",
        default="output.nc",
        help="REMD output NetCDF file (default: output.nc)"
    )
    
    args = parser.parse_args()
    ReportExchangeRatio(args.file)


if __name__ == "__main__":
    main()
