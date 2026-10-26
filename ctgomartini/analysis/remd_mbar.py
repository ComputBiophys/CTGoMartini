"""
REMD MBAR Free Energy Surface (FES) Analysis.

Analyzes free energy surfaces from replica exchange molecular dynamics simulations
using the Multistate Bennett Acceptance Ratio (MBAR) method.

Supports the mixed-state analysis (EXP and HAM methods).

"""

from __future__ import annotations

import csv
import json
import os
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union

import numpy as np
from openmm import unit


# Constants
kB = (unit.MOLAR_GAS_CONSTANT_R).value_in_unit(unit.kilojoule / (unit.kelvin * unit.mole))

# 1 bar * nm^3 per particle, in kJ/mol. Kept in sync with
# ctgomartini.simulation.remd_openmm.BAR_NM3_TO_KJMOL (duplicated so this
# module does not depend on the simulation package).
BAR_NM3_TO_KJMOL = 0.0602214076


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
                 interval: Optional[int] = None):
        """
        Initialize the FES analyzer with simulation data files.
        
        Args:
            output_file: Path to simulation output: an openmmtools NetCDF file,
                or an OpenMM >= 8.6 built-in REMD output directory (containing
                metadata.json and energy.csv).
            cv_file: Path to collective variable (CV) data file.
            interval: Interval for getting CV data compared to energies.
                      For example, if energy is saved every frame but CV
                      is saved every 5 frames, interval=5.
                      If None (default), will auto-detect from checkpoint_interval
                      in the NetCDF file; for a built-in REMD directory the
                      data is already spaced by the report interval, so the
                      default is 1.
        
        Raises:
            FileNotFoundError: If input files do not exist.
            RuntimeError: If data loading fails.
        """
        self.output_file = Path(output_file)
        self.cv_file = Path(cv_file)
        self._user_interval = interval  # Store user-specified interval (None means auto-detect)
        self.fes = None
        
        # Analysis results storage
        self.subsampled_indices: Optional[List[int]] = None
        self.u_kn: Optional[np.ndarray] = None
        self.concatenated_cv: Optional[np.ndarray] = None
        self.subsampled_unsampled_energies: Optional[np.ndarray] = None
        # Raw basin energies in kJ/mol (built-in REMD directories only), in
        # the same [replica][basin][frame] layout as unsampled_energies. The
        # mixing analysis uses them directly so that no temperature rescaling
        # is involved (correct for mixed-temperature REMD).
        self._basin_energies_kj: np.ndarray | None = None
        self._subsampled_basin_energies_kj: np.ndarray | None = None
        # NPT work term W = p*V - gamma*A in kJ/mol (built-in REMD directories
        # with metadata npt.enabled only), [replica][frame]. It enters both
        # the per-state reduced energies (replica_energies) and the mixing
        # analysis target potential.
        self._work_energies_kj: np.ndarray | None = None
        self._subsampled_work_kj: np.ndarray | None = None
        
        self._validate_files()
        self._load_simulation_data()  # This will also set self.interval
    
    def _validate_files(self) -> None:
        """Validate that input files exist."""
        if not self.output_file.exists():
            raise FileNotFoundError(f"Output file not found: {self.output_file}")
        if not self.cv_file.exists():
            raise FileNotFoundError(f"CV file not found: {self.cv_file}")
    
    def _load_simulation_data(self) -> None:
        """Load basic simulation data from files.

        Dispatches on the output layout: an OpenMM >= 8.6 built-in REMD output
        directory (metadata.json + energy.csv) or an openmmtools NetCDF file.
        Also auto-detects interval if not specified by user.
        """
        if self.output_file.is_dir():
            self._load_openmm_data()
            return

        MultiStateReporter, ReplicaExchangeAnalyzer = _import_openmmtools()

        try:
            # Load reporter and analyzer
            reporter = MultiStateReporter(str(self.output_file), open_mode="r")
            analyzer = ReplicaExchangeAnalyzer(reporter)

            # Auto-detect interval if not specified by user
            if self._user_interval is None:
                try:
                    self.interval = reporter.checkpoint_interval
                except Exception:
                    self.interval = 1
            else:
                self.interval = self._user_interval

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
            self._load_cv_data()

        except Exception as e:
            raise RuntimeError(f"Error loading simulation data: {e}")

    def _load_cv_data(self) -> None:
        """Load the external CV data file and transpose to (n_states, n_frames).

        Expects n_states + 1 columns (first column is time/index).
        """
        cv_data = np.loadtxt(self.cv_file, ndmin=2)
        if cv_data.shape[1] != self.n_states + 1:
            raise ValueError(
                f"CV data has wrong shape. Expected {self.n_states + 1} columns, "
                f"got {cv_data.shape[1]}"
            )
        # Transpose to get (n_states, n_frames)
        self.cv_values_replica = cv_data[:, 1:].T

    def _load_openmm_data(self) -> None:
        """Load simulation data from an OpenMM >= 8.6 built-in REMD directory.

        Reads metadata.json, energy.csv, log.csv and cv_energies.csv as written
        by BuiltinREMDRunner. Does not require openmmtools. Frames are already
        spaced by the run's report interval, so the default interval is 1.

        NPT runs (metadata "npt" block with enabled=true) additionally read
        volume.csv (and area.csv for p_type=membrane): the reduced potential
        of state s is beta_s*(U + p*V - gamma*A), so the per-state reduced
        energies gain the work term W/(kB*T_s) with W = p*V - gamma*A.

        Raises:
            RuntimeError: If data loading fails.
        """
        try:
            metadata = json.loads((self.output_file / "metadata.json").read_text())
            n_basins = metadata.get("n_basins", 2)
            if n_basins != 2:
                raise ValueError(
                    f"Mixing analysis supports only two-basin setups "
                    f"(n_basins=2); metadata.json reports n_basins={n_basins}."
                )
            temperatures = [float(s["temperature"]) for s in metadata["states"]]
            n_states = len(temperatures)
            # Data is already written every reportInterval iterations.
            self.interval = self._user_interval if self._user_interval is not None else 1

            # energy.csv: one row per reported iteration holding every
            # replica x state reduced energy, [replica][state] flattened.
            energy = np.loadtxt(
                self.output_file / "energy.csv", delimiter=",", ndmin=2
            )
            n_frames = energy.shape[0]
            if energy.shape[1] != n_states * n_states:
                raise ValueError(
                    f"energy.csv has {energy.shape[1]} columns, expected "
                    f"{n_states * n_states} (replica x state)"
                )
            u_rsn = energy.reshape(n_frames, n_states, n_states)
            # u_kln: [state, replica, frame]
            self.replica_energies = u_rsn.transpose(2, 1, 0)[:, :, ::self.interval]

            # log.csv: per-iteration state assignment (Iteration,Step,Replica_i_State...)
            log = np.loadtxt(
                self.output_file / "log.csv", delimiter=",", skiprows=1, ndmin=2
            )
            state_assignment = log[:, 2 : 2 + n_states].astype(int)

            # cv_energies.csv: long table Iteration,Replica,Molecule,E1,E2[,...];
            # sum over molecules to get the total basin energies per replica.
            basin: dict[tuple[int, int], list[float]] = {}
            with open(self.output_file / "cv_energies.csv") as f:
                for row in csv.DictReader(f):
                    key = (int(row["Iteration"]), int(row["Replica"]))
                    e1, e2 = float(row["E1"]), float(row["E2"])
                    if key in basin:
                        basin[key][0] += e1
                        basin[key][1] += e2
                    else:
                        basin[key] = [e1, e2]
            cv_iterations = sorted({iteration for iteration, _ in basin})

            # volume.csv (and area.csv for membrane coupling) are written for
            # NPT runs only, one row per reported iteration, one value per
            # replica, row-aligned with energy.csv.
            npt = metadata.get("npt") or {}
            npt_enabled = bool(npt.get("enabled"))
            volume = area = None
            if npt_enabled:
                volume = np.loadtxt(
                    self.output_file / "volume.csv", delimiter=",", ndmin=2
                )
                if npt.get("p_type") == "membrane":
                    area = np.loadtxt(
                        self.output_file / "area.csv", delimiter=",", ndmin=2
                    )

            frame_counts = {
                "energy.csv": n_frames,
                "log.csv": log.shape[0],
                "cv_energies.csv": len(cv_iterations),
            }
            if volume is not None:
                frame_counts["volume.csv"] = volume.shape[0]
            if area is not None:
                frame_counts["area.csv"] = area.shape[0]
            if len(set(frame_counts.values())) != 1:
                details = ", ".join(
                    f"{name} has {count}" for name, count in frame_counts.items()
                )
                raise ValueError(f"Frame count mismatch: {details}")
            for name, table in (("volume.csv", volume), ("area.csv", area)):
                if table is not None and table.shape[1] != n_states:
                    raise ValueError(
                        f"{name} has {table.shape[1]} columns, expected "
                        f"{n_states} (one per replica)"
                    )

            # Reduce basin energies by the temperature of the state each
            # replica occupies at that frame, and keep the raw kJ/mol values
            # alongside: the mixing analysis combines the raw energies, which
            # is exact even when states differ in temperature.
            # Both arrays follow the openmmtools layout [replica][basin][frame].
            basin_energies_kj = np.zeros((n_states, 2, n_frames))
            unsampled_energies = np.zeros((n_states, 2, n_frames))
            for frame, iteration in enumerate(cv_iterations):
                for replica in range(n_states):
                    e1, e2 = basin[(iteration, replica)]
                    kt_replica = kB * temperatures[state_assignment[frame][replica]]
                    basin_energies_kj[replica, 0, frame] = e1
                    basin_energies_kj[replica, 1, frame] = e2
                    unsampled_energies[replica, 0, frame] = e1 / kt_replica
                    unsampled_energies[replica, 1, frame] = e2 / kt_replica
            self.unsampled_energies = unsampled_energies[:, :, ::self.interval]
            self._basin_energies_kj = basin_energies_kj[:, :, ::self.interval]

            # NPT work term W = p*V - gamma*A per replica frame (kJ/mol). The
            # reduced potential of state s gains W/(kB*T_s); W itself is
            # state-independent (all states share the reference pressure).
            if npt_enabled:
                pressure = float(npt["pressure_bar"])
                gamma = float(npt.get("surface_tension_bar_nm") or 0.0)
                work_fr = pressure * volume  # [frame][replica]
                if area is not None:
                    work_fr = work_fr - gamma * area
                # [replica][frame], downsampled like replica_energies
                self._work_energies_kj = (work_fr * BAR_NM3_TO_KJMOL).T[
                    :, ::self.interval
                ]
                self.replica_energies = self.replica_energies + (
                    self._work_energies_kj[None, :, :]
                    / (kB * np.asarray(temperatures)[:, None, None])
                )

            self.n_states = n_states
            self.temperatures_k = np.array(temperatures)
            self.beta_k = 1 / (kB * self.temperatures_k)

            # Load CV data
            self._load_cv_data()

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
            if self._basin_energies_kj is not None:
                self._subsampled_basin_energies_kj = self._basin_energies_kj[
                    :, :, self.subsampled_indices
                ]
            if self._work_energies_kj is not None:
                self._subsampled_work_kj = self._work_energies_kj[
                    :, self.subsampled_indices
                ]
            
        except Exception as e:
            raise RuntimeError(f"Failed to initialize FES: {e}")
    
    def analyze(self, u_n: np.ndarray, cv_n: np.ndarray, temperature: float,
               n_bins: int = 100, ranges: Optional[List[float]] = None,
               left_bound: float = 3, right_bound: float = 6) -> Dict:
        """
        Analyze free energy surface (FES).
        
        Calculates the FES as a function of CV, normalized to basin1 = 0.
        
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
            - fes: Free energy surface (kJ/mol), normalized to basin1 = 0
            - fes_uncertainty: FES uncertainty
            - metrics: Dict with barrier_pos, barrier, basin positions, keq
        
        Raises:
            RuntimeError: If FES not initialized.
        """
        if self.fes is None:
            raise RuntimeError("FES not initialized. Call initialize_fes() first.")
        
        # Determine bin ranges
        if ranges is None:
            padding = 0.05 * (cv_n.max() - cv_n.min())  # 5% padding on each side
            ranges = [cv_n.min() + padding, cv_n.max() - padding]
        
        # Validate ranges
        if ranges[0] >= ranges[1]:
            raise ValueError(f"Invalid ranges: {ranges[0]} >= {ranges[1]}")
        
        # Bin the CV
        counts, bin_edges = np.histogram(cv_n, bins=n_bins, range=ranges)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Calculate FES
        self.fes.generate_fes(
            u_n, cv_n, fes_type="histogram",
            histogram_parameters={"bin_edges": [bin_edges]}
        )
        
        # Get FES results with error handling for edge cases
        try:
            fes_results = self.fes.get_fes(
                bin_centers, reference_point="from-lowest",
                uncertainty_method="analytical"
            )
            results = fes_results
        except KeyError as e:
            # Handle PyMBAR KeyError when bins have no samples
            # This can happen with certain parameter combinations
            raise RuntimeError(
                f"FES calculation failed: {e}. "
                f"This often occurs when the data range is too narrow or "
                f"some bins have no samples. Try adjusting ranges "
                f"or using a different start_ratio."
            ) from e
        
        # Convert to kJ/mol
        beta = 1 / (temperature * kB)
        fes = results["f_i"] / beta
        fes_uncertainty = results["df_i"] / beta
        
        # Normalize to left basin (basin1) as zero
        metrics = self.calculate_metrics(
            bin_centers, fes, temperature, left_bound, right_bound
        )
        
        # Apply normalization: shift FES so basin1 is at 0
        if not np.isnan(metrics["basin1_fes"]):
            fes = fes - metrics["basin1_fes"]
            # Update metrics after normalization
            metrics = self.calculate_metrics(
                bin_centers, fes, temperature, left_bound, right_bound
            )
        
        return {
            "cv_values": bin_centers,
            "fes": fes,
            "fes_uncertainty": fes_uncertainty,
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

        # NPT: the target reduced potential includes W/(kB*T) with
        # W = p*V - gamma*A per sample (replica-major order, like the basin
        # energies). analyze_onestate needs no addition: u_kn already carries
        # the work term from the per-state reduced energies.
        if self._subsampled_work_kj is not None:
            u_n = u_n + np.concatenate(self._subsampled_work_kj) / (kB * temperature)

        return self.analyze(u_n, self.concatenated_cv, temperature, **kwargs)
    
    def _get_basin_energies(self, temperature: float) -> tuple[np.ndarray, np.ndarray]:
        """Return (E1, E2) basin energies in kJ/mol for all subsampled samples.

        Concatenated replica-major (replica 0 frames, then replica 1 frames,
        ...), matching the sample order of concatenated_cv. Built-in REMD
        directories store the raw basin energies in kJ/mol, so they are used
        directly (no temperature rescaling); openmmtools NetCDF files store
        reduced basin energies, which are converted back to kJ/mol at the
        analysis temperature.

        Args:
            temperature: Analysis temperature (K), used only to un-reduce
                openmmtools reduced energies.
        """
        if self._subsampled_basin_energies_kj is not None:
            E1 = np.concatenate(self._subsampled_basin_energies_kj[:, 0, :])
            E2 = np.concatenate(self._subsampled_basin_energies_kj[:, 1, :])
        else:
            E1 = np.concatenate(self.subsampled_unsampled_energies[:, 0, :]) * (
                kB * temperature
            )
            E2 = np.concatenate(self.subsampled_unsampled_energies[:, 1, :]) * (
                kB * temperature
            )
        return E1, E2

    def _exp_mixing(self, beta: float, C1: float, C2: float,
                   temperature: float) -> np.ndarray:
        """
        EXP mixing method for multiple-basin potential.
        
        Energy = -(1/beta) * ln[ exp(-beta*(E1+C1)) + exp(-beta*(E2+C2)) ]
        
        Numerical stability is improved by subtracting the mean energy
        to avoid large exponential values.
        """
        E1, E2 = self._get_basin_energies(temperature)
        
        # Avoid large exp values by subtracting mean energy (numerical stability)
        E1_mean = E1.mean()
        E1 = E1 - E1_mean
        E2 = E2 - E1_mean
        
        term1 = np.exp(-beta * (E1 + C1))
        term2 = np.exp(-beta * (E2 + C2))
        
        E = -(1 / beta) * np.log(term1 + term2)
        E = E + E1_mean  # Add back the mean
        
        return E / (kB * temperature)
    
    def _ham_mixing(self, delta: float, C1: float, C2: float,
                   temperature: float) -> np.ndarray:
        """
        HAM (Hamiltonian) mixing method for multiple-basin potential.
        
        Energy = (E1+E2+dV)/2 - sqrt[((E1-E2-dV)/2)^2 + delta^2]
        where dV = C2 - C1
        """
        E1, E2 = self._get_basin_energies(temperature)
        
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
            
            try:
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
            except Exception as e:
                import traceback
                print(f"\nWARNING: {param_name}={value} failed with error: {e}")
                print(f"  Full parameters used: g={params['g']}, "
                      f"length_ratio={params['length_ratio']}, "
                      f"start_ratio={params['start_ratio']}, "
                      f"selected_state={params['selected_state']}")
                print(f"  analyze_kwargs: {analyze_kwargs}")
                print(f"  Original error: {type(e).__name__}: {e}")
                print(f"  Traceback:\n{traceback.format_exc()}")
                
                # Store NaN results for failed parameters
                results[value] = {
                    "cv_values": np.array([]),
                    "fes": np.array([]),
                    "fes_uncertainty": np.array([]),
                    "metrics": {
                        "barrier_pos": float("nan"),
                        "barrier": float("nan"),
                        "basin1_pos": float("nan"),
                        "basin2_pos": float("nan"),
                        "basin1_fes": float("nan"),
                        "basin2_fes": float("nan"),
                        "keq": float("nan")
                    },
                    "error": str(e),
                    "traceback": traceback.format_exc()
                }
        
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
    def calculate_barrier(cv_values: np.ndarray, fes: np.ndarray,
                         left_bound: float = 3, right_bound: float = 6) -> Tuple[float, float]:
        """
        Calculate barrier position and height.
        
        Args:
            cv_values: Array of CV values.
            fes: Free energy surface values.
            left_bound: Left boundary for barrier search.
            right_bound: Right boundary for barrier search.
        
        Returns:
            Tuple of (barrier_position, barrier_height).
        """
        mask = (cv_values > left_bound) & (cv_values < right_bound)
        if not np.any(mask):
            raise ValueError("No data points found within specified bounds")
        
        # Handle NaN values in FES
        fes_masked = fes[mask]
        cv_masked = cv_values[mask]
        
        # Filter out NaN values
        valid_mask = ~np.isnan(fes_masked)
        if not np.any(valid_mask):
            raise ValueError("No valid (non-NaN) FES values within specified bounds")
        
        barrier_idx = np.argmax(fes_masked[valid_mask])
        barrier_pos = cv_masked[valid_mask][barrier_idx]
        barrier_height = fes_masked[valid_mask][barrier_idx]
        
        return float(barrier_pos), float(barrier_height)
    
    @staticmethod
    def calculate_basins(cv_values: np.ndarray, fes: np.ndarray,
                        barrier_position: float) -> Tuple[float, float, float, float]:
        """
        Calculate basin positions and free energies relative to barrier.
        
        Args:
            cv_values: Array of CV values.
            fes: Free energy surface values.
            barrier_position: Position dividing the two basins.
        
        Returns:
            Tuple of (basin1_pos, basin2_pos, basin1_fes, basin2_fes).
        """
        # Left basin (CV < barrier)
        left_basin = cv_values < barrier_position
        if not np.any(left_basin):
            raise ValueError("No data points found in left basin")
        
        fes_left = fes[left_basin]
        cv_left = cv_values[left_basin]
        valid_left = ~np.isnan(fes_left)
        if not np.any(valid_left):
            raise ValueError("No valid FES values in left basin")
        basin1_idx = np.argmin(fes_left[valid_left])
        basin1_pos = cv_left[valid_left][basin1_idx]
        basin1_fes = fes_left[valid_left][basin1_idx]
        
        # Right basin (CV > barrier)
        right_basin = cv_values > barrier_position
        if not np.any(right_basin):
            raise ValueError("No data points found in right basin")
        
        fes_right = fes[right_basin]
        cv_right = cv_values[right_basin]
        valid_right = ~np.isnan(fes_right)
        if not np.any(valid_right):
            raise ValueError("No valid FES values in right basin")
        basin2_idx = np.argmin(fes_right[valid_right])
        basin2_pos = cv_right[valid_right][basin2_idx]
        basin2_fes = fes_right[valid_right][basin2_idx]
        
        return float(basin1_pos), float(basin2_pos), float(basin1_fes), float(basin2_fes)
    
    @staticmethod
    def calculate_equilibrium_constant(cv_values: np.ndarray, fes: np.ndarray,
                                      barrier_position: float,
                                      temperature: float = 310) -> float:
        """
        Calculate equilibrium constant between basins.
        
        K_eq = population1 / population2
        
        Args:
            cv_values: Array of CV values.
            fes: Free energy surface values.
            barrier_position: Position dividing basins.
            temperature: Temperature for Boltzmann weighting (K).
        
        Returns:
            Equilibrium constant (population1/population2).
        """
        beta = 1 / (kB * temperature)
        
        # Calculate populations in each basin
        left_basin = cv_values < barrier_position
        right_basin = cv_values > barrier_position
        
        # Filter out NaN values
        fes_left = fes[left_basin]
        fes_right = fes[right_basin]
        
        pop1 = np.nansum(np.exp(-beta * fes_left))
        pop2 = np.nansum(np.exp(-beta * fes_right))
        
        if pop2 == 0:
            return float("inf")
        return float(pop1 / pop2)
    
    def calculate_metrics(self, cv_values: np.ndarray, fes: np.ndarray,
                         temperature: float, left_bound: float = 3,
                         right_bound: float = 6) -> Dict:
        """
        Calculate key metrics from FES and CV values.
        
        Args:
            cv_values: Array of CV values.
            fes: Free energy surface values.
            temperature: Temperature for equilibrium constant (K).
            left_bound: Left boundary for barrier search.
            right_bound: Right boundary for barrier search.
        
        Returns:
            Dictionary of calculated metrics including barrier position,
            barrier height, basin positions, and equilibrium constant.
        """
        try:
            barrier_pos, barrier = self.calculate_barrier(
                cv_values, fes, left_bound, right_bound
            )
            basin1_pos, basin2_pos, basin1_fes, basin2_fes = self.calculate_basins(
                cv_values, fes, barrier_pos
            )
            keq = self.calculate_equilibrium_constant(
                cv_values, fes, barrier_pos, temperature
            )
            
            return {
                "barrier_pos": barrier_pos,
                "barrier": barrier,
                "basin1_pos": basin1_pos,
                "basin2_pos": basin2_pos,
                "basin1_fes": basin1_fes,
                "basin2_fes": basin2_fes,
                "keq": keq
            }
        except Exception:
            return {
                "barrier_pos": float("nan"),
                "barrier": float("nan"),
                "basin1_pos": float("nan"),
                "basin2_pos": float("nan"),
                "basin1_fes": float("nan"),
                "basin2_fes": float("nan"),
                "keq": float("nan")
            }
