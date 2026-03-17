#!/usr/bin/env python3
"""
Main simulation runner for CTGoMartini.

.. deprecated::
    This module is deprecated. Use ctgomartini.run_ctgomartini instead.
    The new unified runner supports both MD and REMD modes with automatic detection.

Provides functionality for running molecular dynamics simulations using
the OpenMM engine with Martini coarse-grained force fields.

Authors: Song Yang
"""

from __future__ import annotations

import warnings
import sys

# Emit deprecation warning
warnings.warn(
    "data/run_ctgomartini.py is deprecated. "
    "Use the unified runner: run_ctgomartini -i <input_file> "
    "or from ctgomartini.run_ctgomartini import main",
    DeprecationWarning,
    stacklevel=2,
)

# Import and use the new unified runner
from ctgomartini.run_ctgomartini import main as _new_main


def mdrun(inpfile: str) -> None:
    """Execute CTGoMartini molecular dynamics simulation.

    .. deprecated::
        Use run_ctgomartini.main() instead.

    Args:
        inpfile: Path to input parameter file.
    """
    import sys
    from unittest.mock import patch

    # Set up arguments as if called from command line
    test_args = ["run_ctgomartini", "-i", inpfile]
    with patch.object(sys, "argv", test_args):
        _new_main()


if __name__ == "__main__":
    sys.exit(_new_main())
