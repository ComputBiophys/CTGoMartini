#!/usr/bin/env python3
"""
REMD exchange ratio analysis.

Calculates and reports replica exchange ratios from REMD simulation output.

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


def ReportExchangeRatio(
    output_file: str | Path = 'output.nc',
    output_txt: str | Path | None = None
) -> dict:
    """Calculate and report replica exchange ratios.

    Args:
        output_file: Path to REMD output: an openmmtools NetCDF file, or an
            OpenMM >= 8.6 built-in REMD output directory (containing
            exchange.csv).
        output_txt: Optional path to save exchange ratios to text file.

    Returns:
        Dictionary with exchange ratios for each replica pair.
    """
    output_path = Path(output_file)
    if output_path.is_dir():
        return _report_exchange_ratio_openmm(output_path, output_txt)

    MultiStateReporter = _import_openmmtools_reporter()

    reporter = MultiStateReporter(str(output_file), open_mode="r")
    n_accepted_matrix, n_proposed_matrix = reporter.read_mixing_statistics()

    # Because use the swap-neighbors method to exchange states,
    # only half of the exchanges happened for state i and i+1
    ratio = n_accepted_matrix.sum(axis=0) / (n_accepted_matrix.shape[0] / 2)

    results = {}
    lines = ["# REMD Exchange Ratios\n"]
    lines.append("# Replica Pair\tExchange Ratio\n")

    for i in range(ratio.shape[0] - 1):
        j = i + 1
        results[f"{i}-{j}"] = float(ratio[i, j])
        line = f"{i}-{j}\t{ratio[i, j]:.6f}\n"
        lines.append(line)
        print(f'Exchange ratio between replica {i} and {j}: {ratio[i, j]:.4f}')

    reporter.close()

    # Save to file if output path is provided
    if output_txt is not None:
        with open(output_txt, "w") as f:
            f.writelines(lines)
        print(f"Saved exchange ratios to: {output_txt}")

    return results


def _report_exchange_ratio_openmm(
    directory: Path,
    output_txt: str | Path | None = None
) -> dict:
    """Calculate exchange ratios from a built-in REMD output directory.

    Reads the cumulative counters in exchange.csv (last row per pair), so the
    ratio is simply accepted/proposed for the whole run.
    """
    exchange_csv = directory / "exchange.csv"
    if not exchange_csv.is_file():
        raise FileNotFoundError(
            f"Cannot compute exchange ratios: {exchange_csv} does not exist."
        )

    # exchange.csv is cumulative: the last row of each pair holds the totals
    accepted: dict[str, int] = {}
    proposed: dict[str, int] = {}
    with open(exchange_csv) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("Iteration"):
                continue
            _, pair, acc, prop = line.split(",")[:4]
            accepted[pair] = int(acc)
            proposed[pair] = int(prop)

    results = {}
    lines = ["# REMD Exchange Ratios\n"]
    lines.append("# Replica Pair\tExchange Ratio\n")

    for pair in sorted(accepted, key=lambda p: tuple(int(x) for x in p.split("-"))):
        ratio = accepted[pair] / proposed[pair] if proposed[pair] > 0 else float("nan")
        results[pair] = float(ratio)
        lines.append(f"{pair}\t{ratio:.6f}\n")
        print(f'Exchange ratio between replica {pair.replace("-", " and ")}: {ratio:.4f}')

    if output_txt is not None:
        with open(output_txt, "w") as f:
            f.writelines(lines)
        print(f"Saved exchange ratios to: {output_txt}")

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
    parser.add_argument(
        "-o", "--output",
        default="exchange_ratios.dat",
        help="Output file for exchange ratios (default: exchange_ratios.dat)"
    )
    
    args = parser.parse_args()
    ReportExchangeRatio(args.file, output_txt=args.output)


if __name__ == "__main__":
    main()
