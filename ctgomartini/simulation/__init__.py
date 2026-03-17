"""
Simulation module for CTGoMartini.

This module provides classes and functions for running molecular dynamics simulations,
including standard MD and replica exchange MD (REMD).

Example:
    >>> from ctgomartini.simulation import MDRunner, load_config
    >>> config = load_config("md.inp")
    >>> runner = MDRunner(config)
    >>> runner.run()
"""

from ctgomartini.simulation.config import SimulationConfig, load_config
from ctgomartini.simulation.base import SimulationRunner, add_restraints, generate_restraints
from ctgomartini.simulation.md import MDRunner
from ctgomartini.simulation.remd import REMDRunner

__all__ = [
    # Configuration
    "SimulationConfig",
    "load_config",
    # Runners
    "SimulationRunner",
    "MDRunner",
    "REMDRunner",
    # Utility functions
    "add_restraints",
    "generate_restraints",
]
