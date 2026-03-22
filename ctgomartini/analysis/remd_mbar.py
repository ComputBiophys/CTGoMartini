"""
REMD MBAR Free Energy Surface (FES) Analysis.

Analyzes free energy surfaces from replica exchange molecular dynamics simulations
using the Multistate Bennett Acceptance Ratio (MBAR) method.

Supports the mixed-state analysis (EXP and HAM methods).

"""

from __future__ import annotations

import os
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

import numpy as np
from openmm import unit

# Constants
kB = (unit.MOLAR_GAS_CONSTANT_R).value_in_unit(unit.kilojoule / (unit.kelvin * unit.mole))


def _import_pymbar():
    """Lazy import pymbar to avoid import-time overhead."""
    import pymbar
    from pymbar import timeseries
    return pymbar, timeseries


def _import_openmmtools():
    """Lazy import openmmtools to avoid import-time warnings."""
    from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer
    return MultiStateReporter, ReplicaExchangeAnalyzer


class FESAnalyzer:
    """
    Free Energy Surface (FES) analyzer for REMD simulations using MBAR.
    
    This class provides tools to analyze free energy surfaces from REMD
    simulations, supporting different mixed-state analysis
    (EXP and HAM methods).
    
    Attributes:
        output_file: Path to REMD output NetCDF file.
        cv_file: Path to collective variable (CV) data file.
        interval: Interval between CV frames and energy frames.
        n_states: Number of thermodynamic states/replicas.
        temperatures_k: Array of temperatures for each state (K).
        beta_k: Inverse temperatures (1/kT) for each state.
    
    Example:
        >>> analyzer = FESAnalyzer('output.nc', 'distance.dat', interval=5)
        >>> analyzer.initialize_fes(g=1.0, start_ratio=0.2)
        >>> results = analyzer.analyze_onestate(selected_state=5)
        >>> print(results['metrics']['barrier'])
    """
    
    def __init__(self, output_file: Union[str, Path], cv_file: Union[str, Path], 
                 interval: int = 1):
        """
        Initialize the FES analyzer with simulation data files.
        
        Args:
            output_file: Path to simulation output NetCDF file.
            cv_file: Path to collective variable (CV) data file.
            interval: Interval for getting CV data compared to energies.
                      For example, if energy is saved every frame but CV
                      is saved every 5 frames, interval=5.
        
        Raises:
            FileNotFoundError: If input files do not exist.
            RuntimeError: If data loading fails.
        """
        self.output_file = Path(output_file)
        self.cv_file = Path(cv_file)
        self.interval = interval
        self.fes = None
        
        # Analysis results storage
        self.subsampled_indices: Optional[List[int]] = None
        self.u_kn: Optional[np.ndarray] = None
        self.concatenated_cv: Optional[np.ndarray] = None
        self.subsampled_unsampled_energies: Optional[np.ndarray] = None
        
        self._validate_files()
        self._load_simulation_data()
    
    def _validate_files(self) -> None:
        """Validate that input files exist."""
        if not self.output_file.exists():
            raise FileNotFoundError(f"Output file not found: {self.output_file}")
        if not self.cv_file.exists():
            raise FileNotFoundError(f"CV file not found: {self.cv_file}")
    
    def _load_simulation_data(self) -> None:
        """Load basic simulation data from files."""
        MultiStateReporter, ReplicaExchangeAnalyzer = _import_openmmtools()
        
        try:
            # Load reporter and analyzer
            reporter = MultiStateReporter(str(self.output_file), open_mode="r")
            analyzer = ReplicaExchangeAnalyzer(reporter)
            
            # Read simulation data
            states, _ = reporter.read_thermodynamic_states()
            temperature_list = [s.temperature for s in states]
            
            # Read energies and states
            replica_energies, unsampled_energies, _, replica_state_indices = analyzer.read_energies()
            
            # Downsample energies to match CV interval
            self.replica_energies = replica_energies[:, :, ::self.interval]
            self.unsampled_energies = unsampled_energies[:, :, ::self.interval]
            
            reporter.close()
            
            # Verify number of replicas matches number of states
            n_replicas = self.replica_energies.shape[0]
            n_states = len(states)
            if n_replicas != n_states:
                raise ValueError(
                    f"Number of replicas ({n_replicas}) does not match "
                    f"number of states ({n_states})"
                )
            
            # Process temperature data
            self.n_states = n_states
            self.temperatures_k = np.array([
                temp.value_in_unit(unit.kelvin) for temp in temperature_list
            ])
            self.beta_k = 1 / (kB * self.temperatures_k)
            
            # Load CV data
            cv_data = np.loadtxt(self.cv_file)
            if cv_data.shape[1] != n_states + 1:
                raise ValueError(
                    f"CV data has wrong shape. Expected {n_states + 1} columns, "
                    f"got {cv_data.shape[1]}"
                )
            # Transpose to get (n_states, n_frames)
            self.cv_values_replica = cv_data[:, 1:].T
            
        except Exception as e:
            raise RuntimeError(f"Error loading simulation data: {e}")
    
    def _get_subsampled_indices(self, g: float = 1.0, length_ratio: float = 1.0,
                               start_point: Optional[int] = None,
                               start_ratio: Optional[float] = None) -> np.ndarray:
        """
        Determine subsampled indices based on analysis parameters.
        
        Uses replica 0's CV data for subsampling. This assumes that the
        correlation time is similar across all replicas.
        
        Args:
            g: Statistical inefficiency value (g=1 means no subsampling).
            length_ratio: Fraction of trajectory length to use.
            start_point: Starting frame index (overrides start_ratio).
            start_ratio: Fraction of trajectory to skip at beginning.
        
        Returns:
            Array of subsampled indices.
        """
        _, timeseries = _import_pymbar()
        
        total_frames = self.replica_energies.shape[2]
        
        # Determine start index
        if start_point is not None:
            start_idx = start_point
        elif start_ratio is not None:
            start_idx = int(total_frames * start_ratio)
        else:
            start_idx = 0
        
        # Determine end index
        n_frames = int((total_frames - start_idx) * length_ratio)
        end_idx = start_idx + n_frames
        
        if end_idx <= start_idx:
            raise ValueError("Invalid frame selection: end index <= start index")
        
        # Subsample to remove correlation (using replica 0)
        subsampled = timeseries.subsample_correlated_data(
            self.cv_values_replica[0, start_idx:end_idx], g=g
        )
        self.subsampled_indices = [x + start_idx for x in subsampled]
        
        return self.subsampled_indices
    
    def initialize_fes(self, g: float = 1.0, length_ratio: float = 1.0,
                      start_point: Optional[int] = None,
                      start_ratio: Optional[float] = None) -> None:
        """
        Initialize free energy surface (FES) calculation using MBAR.
        
        Prepares the data by subsampling and creates the MBAR object.
        Must be called before any analyze* methods.
        
        Args:
            g: Statistical inefficiency value for subsampling.
            length_ratio: Fraction of trajectory to analyze.
            start_point: Starting frame index.
            start_ratio: Fraction of trajectory to skip at start.
        
        Raises:
            RuntimeError: If FES initialization fails.
        """
        pymbar, _ = _import_pymbar()
        
        try:
            # Get subsampled indices
            self._get_subsampled_indices(
                g=g, length_ratio=length_ratio,
                start_point=start_point, start_ratio=start_ratio
            )
            
            # Prepare energy and CV arrays
            u_kln = self.replica_energies[:, :, self.subsampled_indices]
            cv_kn = self.cv_values_replica[:, self.subsampled_indices]
            
            # Convert from kln to kn format
            # k: state, l: replica, n: frame
            u_kn = pymbar.utils.kln_to_kn(u_kln)
            
            # Calculate FES
            n_states_array = len(self.subsampled_indices) * np.ones(
                self.n_states, dtype=int
            )
            self.fes = pymbar.FES(
                u_kn, n_states_array, 
                mbar_options=dict(solver_protocol="robust")
            )
            
            # Store for later analysis
            self.u_kn = u_kn
            self.concatenated_cv = np.concatenate(cv_kn)
            self.subsampled_unsampled_energies = self.unsampled_energies[:, :, self.subsampled_indices]
            
        except Exception as e:
            raise RuntimeError(f"Failed to initialize FES: {e}")
    
    def analyze(self, u_n: np.ndarray, cv_n: np.ndarray, temperature: float,
               n_bins: int = 100, ranges: Optional[List[float]] = None,
               left_bound: float = 3, right_bound: float = 6) -> Dict:
        """
        Analyze free energy surface.
        
        Calculates the potential of mean force (PMF) as a function of CV.
        
        Args:
            u_n: Reduced potential energy array.
            cv_n: Collective variable array.
            temperature: Temperature for analysis (K).
            n_bins: Number of bins for histogram.
            ranges: [min, max] range for CV values.
            left_bound: Left boundary for barrier search.
            right_bound: Right boundary for barrier search.
        
        Returns:
            Dictionary containing:
            - cv_values: Bin centers
            - pmf: Potential of mean force (kJ/mol)
            - pmf_uncertainty: PMF uncertainty
            - metrics: Dict with barrier_pos, barrier, basin positions, keq
        
        Raises:
            RuntimeError: If FES not initialized.
        """
        if self.fes is None:
            raise RuntimeError("FES not initialized. Call initialize_fes() first.")
        
        # Determine bin ranges
        if ranges is None:
            # padding = 0.05 * (cv_n.max() - cv_n.min())  # Old padding
            padding = 3 / n_bins * (cv_n.max() - cv_n.min())  # v2: scale with bins
            ranges = [cv_n.min() + padding, cv_n.max() - padding]
        
        # Bin the CV
        counts, bin_edges = np.histogram(cv_n, bins=n_bins, range=ranges)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Calculate FES
        self.fes.generate_fes(
            u_n, cv_n, fes_type="histogram",
            histogram_parameters={"bin_edges": [bin_edges]}
        )
        results = self.fes.get_fes(
            bin_centers, reference_point="from-lowest",
            uncertainty_method="analytical"
        )
        
        # Convert to kJ/mol
        beta = 1 / (temperature * kB)
        pmf = results["f_i"] / beta
        pmf_uncertainty = results["df_i"] / beta
        
        # Calculate metrics
        metrics = self.calculate_metrics(
            bin_centers, pmf, temperature, left_bound, right_bound
        )
        
        return {
            "cv_values": bin_centers,
            "pmf": pmf,
            "pmf_uncertainty": pmf_uncertainty,
            "metrics": metrics
        }
    
    def analyze_onestate(self, selected_state: int = 0, **kwargs) -> Dict:
        """
        Perform FES analysis for a single thermodynamic state.
        
        Args:
            selected_state: Index of thermodynamic state to analyze.
            **kwargs: Additional arguments passed to analyze().
        
        Returns:
            Dictionary containing analysis results.
        
        Raises:
            ValueError: If selected_state is invalid.
        """
        if not (0 <= selected_state < self.n_states):
            raise ValueError(
                f"Invalid state index {selected_state}. "
                f"Must be between 0 and {self.n_states - 1}"
            )
        
        return self.analyze(
            self.u_kn[selected_state, :],
            self.concatenated_cv,
            self.temperatures_k[selected_state],
            **kwargs
        )
    
    def analyze_one_mixing_parameters(self, mixing_parameters: Dict,
                                     temperature: float, method: str = "EXP",
                                     **kwargs) -> Dict:
        """
        Perform FES analysis for a single mixing parameter set.
        
        Supports both EXP (exponential) and HAM (Hamiltonian) mixing methods
        for multiple-basin Go-Martini simulations.
        
        Args:
            mixing_parameters: Dictionary of mixing parameters.
                For EXP: {'beta': float, 'C1': float, 'C2': float}
                For HAM: {'delta': float, 'C1': float, 'C2': float}
            temperature: Temperature for analysis (K).
            method: Mixing method ('EXP' or 'HAM').
            **kwargs: Additional arguments passed to analyze().
        
        Returns:
            Dictionary containing analysis results.
        
        Raises:
            ValueError: If method is not 'EXP' or 'HAM'.
        """
        if method == "EXP":
            beta = mixing_parameters["beta"]
            C1 = mixing_parameters["C1"]
            C2 = mixing_parameters["C2"]
            u_n = self._exp_mixing(beta, C1, C2, temperature)
        elif method == "HAM":
            delta = mixing_parameters["delta"]
            C1 = mixing_parameters["C1"]
            C2 = mixing_parameters["C2"]
            u_n = self._ham_mixing(delta, C1, C2, temperature)
        else:
            raise ValueError(f"Invalid mixing method {method}. Must be 'EXP' or 'HAM'.")
        
        return self.analyze(u_n, self.concatenated_cv, temperature, **kwargs)
    
    def _exp_mixing(self, beta: float, C1: float, C2: float,
                   temperature: float) -> np.ndarray:
        """
        EXP mixing method for multiple-basin potential.
        
        Energy = -(1/beta) * ln[ exp(-beta*(E1+C1)) + exp(-beta*(E2+C2)) ]
        """
        E1 = np.concatenate(self.subsampled_unsampled_energies[:, 0, :]) * (kB * temperature)
        E2 = np.concatenate(self.subsampled_unsampled_energies[:, 1, :]) * (kB * temperature)
        
        term1 = np.exp(-beta * (E1 + C1))
        term2 = np.exp(-beta * (E2 + C2))
        
        return -(1 / beta) * np.log(term1 + term2) / (kB * temperature)
    
    def _ham_mixing(self, delta: float, C1: float, C2: float,
                   temperature: float) -> np.ndarray:
        """
        HAM (Hamiltonian) mixing method for multiple-basin potential.
        
        Energy = (E1+E2+dV)/2 - sqrt[((E1-E2-dV)/2)^2 + delta^2]
        where dV = C2 - C1
        """
        E1 = np.concatenate(self.subsampled_unsampled_energies[:, 0, :]) * (kB * temperature)
        E2 = np.concatenate(self.subsampled_unsampled_energies[:, 1, :]) * (kB * temperature)
        
        dV = C2 - C1
        part1 = (E1 + E2 + dV) / 2
        part2 = (E1 - E2 - dV) / 2
        
        energy = part1 - np.sqrt(part2**2 + delta**2)
        return energy / (kB * temperature)
    
    def parameter_sweep(self, param_name: str, values: List,
                       default_params: Optional[Dict] = None,
                       **analyze_kwargs) -> Dict:
        """
        Perform a parameter sweep analysis.
        
        Useful for checking convergence with respect to various parameters
        like start_ratio, g, length_ratio, or selected_state.
        
        Args:
            param_name: Name of parameter to vary.
                Options: 'g', 'length_ratio', 'start_ratio', 'start_point', 
                         'selected_state'
            values: List of values to test for the parameter.
            default_params: Dictionary of default parameter values.
            **analyze_kwargs: Additional arguments passed to analyze().
        
        Returns:
            Dictionary mapping parameter values to analysis results.
        """
        # Set defaults and merge with user-provided params
        full_defaults = {
            "g": 1,
            "length_ratio": 1,
            "start_point": None,
            "start_ratio": 0,
            "selected_state": 0
        }
        if default_params is not None:
            full_defaults.update(default_params)
        
        results = {}
        for value in values:
            params = full_defaults.copy()
            params[param_name] = value
            
            self.initialize_fes(
                g=params["g"],
                length_ratio=params["length_ratio"],
                start_point=params["start_point"],
                start_ratio=params["start_ratio"]
            )
            
            result = self.analyze_onestate(
                selected_state=params["selected_state"],
                **analyze_kwargs
            )
            results[value] = result
            print(f"{param_name}={value}: {result['metrics']}")
        
        return results
    
    # -------------------------------------------------------------------------
    # Static utility methods for calculation
    # -------------------------------------------------------------------------
    
    @staticmethod
    def save_results(results: Dict, filename: Union[str, Path]) -> None:
        """Save analysis results to a pickle file."""
        with open(filename, "wb") as f:
            pickle.dump(results, f)
    
    @staticmethod
    def load_results(filename: Union[str, Path]) -> Dict:
        """Load analysis results from a pickle file."""
        with open(filename, "rb") as f:
            return pickle.load(f)
    
    @staticmethod
    def calculate_barrier(cv_values: np.ndarray, pmf: np.ndarray,
                         left_bound: float = 3, right_bound: float = 6) -> Tuple[float, float]:
        """
        Calculate barrier position and height.
        
        Args:
            cv_values: Array of CV values.
            pmf: Potential of mean force values.
            left_bound: Left boundary for barrier search.
            right_bound: Right boundary for barrier search.
        
        Returns:
            Tuple of (barrier_position, barrier_height).
        """
        mask = (cv_values > left_bound) & (cv_values < right_bound)
        if not np.any(mask):
            raise ValueError("No data points found within specified bounds")
        
        barrier_idx = np.argmax(pmf[mask])
        barrier_pos = cv_values[mask][barrier_idx]
        barrier_height = pmf[mask][barrier_idx]
        
        return float(barrier_pos), float(barrier_height)
    
    @staticmethod
    def calculate_basins(cv_values: np.ndarray, pmf: np.ndarray,
                        barrier_position: float) -> Tuple[float, float, float, float]:
        """
        Calculate basin positions and free energies relative to barrier.
        
        Args:
            cv_values: Array of CV values.
            pmf: Potential of mean force values.
            barrier_position: Position dividing the two basins.
        
        Returns:
            Tuple of (basin1_pos, basin2_pos, basin1_pmf, basin2_pmf).
        """
        # Left basin (CV < barrier)
        left_basin = cv_values < barrier_position
        if not np.any(left_basin):
            raise ValueError("No data points found in left basin")
        basin1_idx = np.argmin(pmf[left_basin])
        basin1_pos = cv_values[left_basin][basin1_idx]
        basin1_pmf = pmf[left_basin][basin1_idx]
        
        # Right basin (CV > barrier)
        right_basin = cv_values > barrier_position
        if not np.any(right_basin):
            raise ValueError("No data points found in right basin")
        basin2_idx = np.argmin(pmf[right_basin])
        basin2_pos = cv_values[right_basin][basin2_idx]
        basin2_pmf = pmf[right_basin][basin2_idx]
        
        return float(basin1_pos), float(basin2_pos), float(basin1_pmf), float(basin2_pmf)
    
    @staticmethod
    def calculate_equilibrium_constant(cv_values: np.ndarray, pmf: np.ndarray,
                                      barrier_position: float,
                                      temperature: float = 310) -> float:
        """
        Calculate equilibrium constant between basins.
        
        K_eq = population1 / population2
        
        Args:
            cv_values: Array of CV values.
            pmf: Potential of mean force values.
            barrier_position: Position dividing basins.
            temperature: Temperature for Boltzmann weighting (K).
        
        Returns:
            Equilibrium constant (population1/population2).
        """
        beta = 1 / (kB * temperature)
        
        # Calculate populations in each basin
        left_basin = cv_values < barrier_position
        right_basin = cv_values > barrier_position
        
        pop1 = np.exp(-beta * pmf[left_basin]).sum()
        pop2 = np.exp(-beta * pmf[right_basin]).sum()
        
        if pop2 == 0:
            return float("inf")
        return float(pop1 / pop2)
    
    def calculate_metrics(self, cv_values: np.ndarray, pmf: np.ndarray,
                         temperature: float, left_bound: float = 3,
                         right_bound: float = 6) -> Dict:
        """
        Calculate key metrics from PMF and CV values.
        
        Args:
            cv_values: Array of CV values.
            pmf: Potential of mean force values.
            temperature: Temperature for equilibrium constant (K).
            left_bound: Left boundary for barrier search.
            right_bound: Right boundary for barrier search.
        
        Returns:
            Dictionary of calculated metrics including barrier position,
            barrier height, basin positions, and equilibrium constant.
        """
        try:
            barrier_pos, barrier = self.calculate_barrier(
                cv_values, pmf, left_bound, right_bound
            )
            basin1_pos, basin2_pos, basin1_pmf, basin2_pmf = self.calculate_basins(
                cv_values, pmf, barrier_pos
            )
            keq = self.calculate_equilibrium_constant(
                cv_values, pmf, barrier_pos, temperature
            )
            
            return {
                "barrier_pos": barrier_pos,
                "barrier": barrier,
                "basin1_pos": basin1_pos,
                "basin2_pos": basin2_pos,
                "basin1_pmf": basin1_pmf,
                "basin2_pmf": basin2_pmf,
                "keq": keq
            }
        except Exception:
            return {
                "barrier_pos": float("nan"),
                "barrier": float("nan"),
                "basin1_pos": float("nan"),
                "basin2_pos": float("nan"),
                "basin1_pmf": float("nan"),
                "basin2_pmf": float("nan"),
                "keq": float("nan")
            }
