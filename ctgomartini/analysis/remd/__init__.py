"""REMD analysis tools for CTGoMartini.

This module provides analysis scripts for replica exchange molecular dynamics
simulations, including free energy calculations, replica state analysis, and
dRMS trajectory analysis.

Note: Some functions require openmmtools and pymbar packages.
"""

from .remd_drms_analysis import plot_drms_trajectory
from .remd_exchange_ratio import ReportExchangeRatio

# Optional imports (require openmmtools)
try:
    from .remd_replica_state import load_replica_states, plot_replica_states
    from .remd_free_energy import load_free_energy_data, plot_free_energy_comparison
    from .remd_mbar_analysis import (
        plot_selected_state_analysis,
        plot_start_ratio_analysis,
        plot_g_values_analysis,
    )
    _has_openmmtools = True
except ImportError:
    _has_openmmtools = False
    load_replica_states = None  # type: ignore
    plot_replica_states = None  # type: ignore
    load_free_energy_data = None  # type: ignore
    plot_free_energy_comparison = None  # type: ignore
    plot_selected_state_analysis = None  # type: ignore
    plot_start_ratio_analysis = None  # type: ignore
    plot_g_values_analysis = None  # type: ignore

__all__ = [
    "plot_drms_trajectory",
    "ReportExchangeRatio",
]

if _has_openmmtools:
    __all__.extend([
        "load_replica_states",
        "plot_replica_states",
        "load_free_energy_data",
        "plot_free_energy_comparison",
        "plot_selected_state_analysis",
        "plot_start_ratio_analysis",
        "plot_g_values_analysis",
    ])
