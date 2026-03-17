#!/usr/bin/env python3
"""
Unified simulation runner for CTGoMartini.

Provides a single entry point for running both standard MD and REMD simulations.
Automatically detects the simulation mode based on the REMD parameter in the
input file. Supports the --append flag for continuing existing simulations.

Authors: Song Yang
"""

from __future__ import annotations

import argparse
import sys
from typing import Any

from ctgomartini.core.simulation import MDRunner, REMDRunner
from ctgomartini.utils import read_inputs


def main() -> int:
    """Main entry point for CTGoMartini simulation runner.

    Parses command line arguments, detects simulation mode (MD/REMD),
    and dispatches to the appropriate runner.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    parser = argparse.ArgumentParser(
        description="Run CTGoMartini molecular dynamics simulation (MD or REMD)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run a new standard MD simulation
  run_ctgomartini -i md.inp

  # Continue an existing MD simulation (append mode)
  run_ctgomartini -i md.inp --append

  # Run a new REMD simulation (REMD = yes in input file)
  run_ctgomartini -i remd.inp

  # Continue an existing REMD simulation
  run_ctgomartini -i remd.inp --append
        """,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="inpfile",
        required=True,
        help="Input parameter file containing simulation settings",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to existing simulation (auto-detected if not specified)",
    )
    parser.add_argument(
        "--replica-params",
        type=str,
        default=None,
        help="JSON file with replica exchange parameters (REMD only)",
    )
    parser.add_argument(
        "--output-data",
        type=str,
        default="output.nc",
        help="Output NetCDF file for REMD data (default: output.nc)",
    )

    args = parser.parse_args()

    # Read input file to detect simulation mode
    try:
        inputs = read_inputs(args.inpfile)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Determine if this is REMD or standard MD
    is_remd = inputs.remd == "yes"

    if is_remd:
        # Run REMD simulation
        return _run_remd(args, inputs)
    else:
        # Run standard MD simulation
        return _run_md(args, inputs)


def _run_md(args: argparse.Namespace, inputs: Any) -> int:
    """Run standard MD simulation.

    Args:
        args: Parsed command line arguments.
        inputs: Input parameters object.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    print("=" * 60)
    print("CTGoMartini Standard MD Simulation")
    print("=" * 60)

    try:
        runner = MDRunner(args.inpfile)

        # Determine whether to append or start fresh
        if args.append or runner.detect_existing_output():
            if args.append:
                print("\nAppend mode: Continuing existing simulation")
            else:
                print("\nAuto-detected existing output, continuing simulation")
            runner.extend()
        else:
            print("\nStarting new simulation")
            runner.run()

        return 0

    except Exception as e:
        print(f"\nError during MD simulation: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


def _run_remd(args: argparse.Namespace, inputs: Any) -> int:
    """Run REMD simulation.

    Args:
        args: Parsed command line arguments.
        inputs: Input parameters object.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    print("=" * 60)
    print("CTGoMartini Replica Exchange MD Simulation")
    print("=" * 60)

    # Load replica parameters if provided
    replica_params: dict[str, Any] | None = None
    if args.replica_params:
        import json

        try:
            with open(args.replica_params, "r") as f:
                replica_params = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Error loading replica parameters: {e}", file=sys.stderr)
            return 1
    else:
        # Use default replica parameters for testing
        # In production, these should come from the input file or external config
        print("\nWarning: No replica parameters provided, using defaults")
        import numpy as np

        replica_params = {
            "molname": "MOL",
            "beta": [1 / 500] * 21,
            "C1": np.linspace(-2000, 2000, 21).tolist(),
            "C2": [0] * 21,
        }

    try:
        runner = REMDRunner(
            args.inpfile,
            replica_params=replica_params,
            output_data=args.output_data,
        )

        # Determine whether to append or start fresh
        if args.append or runner.detect_existing_output():
            if args.append:
                print("\nAppend mode: Continuing existing REMD simulation")
            else:
                print("\nAuto-detected existing output, continuing REMD simulation")
            runner.extend()
        else:
            print("\nStarting new REMD simulation")
            runner.run()

        return 0

    except ImportError as e:
        print(f"\nError: {e}", file=sys.stderr)
        print("REMD requires openmmtools. Install with: pip install openmmtools", file=sys.stderr)
        return 1

    except Exception as e:
        print(f"\nError during REMD simulation: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
