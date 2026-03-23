#!/usr/bin/env python3
"""
Unified simulation runner for CTGoMartini.

Provides a single entry point for running both standard MD and REMD simulations.
Automatically detects the simulation mode based on the REMD parameter in the
input file. Supports the --append flag for continuing existing simulations.

"""

from __future__ import annotations

import argparse
import sys
from typing import Any

from ctgomartini.simulation import MDRunner, REMDRunner, load_config


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
        help="Append to existing simulation",
    )

    args = parser.parse_args()

    # Read input file to detect simulation mode
    try:
        config = load_config(args.inpfile)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Determine if this is REMD or standard MD
    is_remd = config.remd == "yes"

    if is_remd:
        # Run REMD simulation
        return _run_remd(args, config)
    else:
        # Run standard MD simulation
        return _run_md(args, config)


def _run_md(args: argparse.Namespace, config: Any) -> int:
    """Run standard MD simulation.

    Args:
        args: Parsed command line arguments.
        config: Simulation configuration object.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    print("=" * 60)
    print("CTGoMartini Standard MD")
    print("=" * 60)

    try:
        runner = MDRunner(args.inpfile)

        # Determine whether to append or start fresh
        if args.append:
            print("\n[Mode] Continue (append)")
            runner.extend()
        else:
            print("\n[Mode] New simulation")
            runner.run()

        return 0

    except Exception as e:
        print(f"\nError during MD simulation: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        return 1


def _run_remd(args: argparse.Namespace, config: Any) -> int:
    """Run REMD simulation.

    Args:
        args: Parsed command line arguments.
        config: Simulation configuration object.

    Returns:
        Exit code (0 for success, non-zero for failure).
    """
    print("=" * 60)
    print("CTGoMartini Replica Exchange MD")
    print("=" * 60)

    # Build replica parameters from config
    # NOTE: For H-REMD, 'beta' is the coupling constant, not 1/(kB*T)
    # Original code: beta = [1/300] * n_replicas
    replica_params: dict[str, Any] = {
        "molname": config.replica_molname,
        "beta": config.replica_coupling,  # Use coupling constant from input
        "C1": config.replica_c1,
        "C2": config.replica_c2,
    }

    try:
        runner = REMDRunner(
            args.inpfile,
            replica_params=replica_params,
            output_data=config.remd_output,
        )

        # Determine whether to append or start fresh
        if args.append:
            print("\nAppend mode: Continuing existing REMD simulation")
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
