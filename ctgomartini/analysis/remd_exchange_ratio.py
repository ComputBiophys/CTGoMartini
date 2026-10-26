#!/usr/bin/env python3
"""REMD exchange ratio analysis.

Calculates and reports replica exchange ratios from the built-in REMD
output directory (exchange.csv).

"""

from __future__ import annotations

import argparse
from pathlib import Path


def ReportExchangeRatio(  # noqa: N802
    output_file: str | Path = 'output_remd',
    output_txt: str | Path | None = None
) -> dict:
    """Calculate and report replica exchange ratios.

    Reads the cumulative counters in exchange.csv (last row per pair), so the
    ratio is simply accepted/proposed for the whole run.

    Args:
        output_file: Path to the OpenMM >= 8.6 built-in REMD output
            directory (containing exchange.csv).
        output_txt: Optional path to save exchange ratios to text file.

    Returns:
        Dictionary with exchange ratios for each replica pair.

    Raises:
        ValueError: If output_file is not a directory (NetCDF output from the
            removed openmmtools backend is no longer supported).
        FileNotFoundError: If exchange.csv is missing.
    """
    directory = Path(output_file)
    if not directory.is_dir():
        raise ValueError(
            f"ReportExchangeRatio requires a built-in REMD output directory "
            f"(containing exchange.csv); got {output_file}. openmmtools "
            f"NetCDF (.nc) output is no longer supported."
        )

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
        default="output_remd",
        help="REMD output directory (default: output_remd)"
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
