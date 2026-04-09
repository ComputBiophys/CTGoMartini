"""
Replica Exchange MD simulation runner module for CTGoMartini.

Provides the REMDRunner class for running replica exchange molecular dynamics
simulations with support for new simulations and extending existing simulations.
"""

from __future__ import annotations

import datetime
import os
import sys
from typing import Any

import numpy as np
import openmm as mm
import openmm.unit as u
from openmm.app import GromacsGroFile, PDBFile

from ctgomartini.simulation.base import (
    SimulationRunner,
    report_time,
    load_structure,
    generate_restraints,
    add_restraints,
)

# Import MartiniTopFile here to avoid circular imports
from ctgomartini.topology import MartiniTopFile

# Lazy imports for openmmtools to avoid warnings on module load
_global_context_cache = None
_MultiStateReporter = None
_ReplicaExchangeSampler = None
_XTCMultiStateReporter = None


def _import_remd_dependencies():
    """
    Lazy import all REMD dependencies to avoid import-time warnings.
    
    Imports openmmtools (which triggers pymbar/jax warnings) and 
    XTCMultiStateReporter only when actually needed.
    
    Returns:
        Tuple of (global_context_cache, MultiStateReporter, 
                  ReplicaExchangeSampler, XTCMultiStateReporter)
    """
    global _global_context_cache, _MultiStateReporter, _ReplicaExchangeSampler, _XTCMultiStateReporter
    
    if _MultiStateReporter is None:
        try:
            from openmmtools.cache import global_context_cache
            from openmmtools.multistate import (
                MultiStateReporter,
                ReplicaExchangeSampler,
            )
            _global_context_cache = global_context_cache
            _MultiStateReporter = MultiStateReporter
            _ReplicaExchangeSampler = ReplicaExchangeSampler
        except ImportError as e:
            raise ImportError(
                "openmmtools is required for REMD simulations. "
                "Install it with: conda install openmmtools"
            ) from e
    
    if _XTCMultiStateReporter is None:
        try:
            from ctgomartini.simulation.xtc_reporter import XTCMultiStateReporter
            _XTCMultiStateReporter = XTCMultiStateReporter
        except ImportError:
            _XTCMultiStateReporter = None
    
    return _global_context_cache, _MultiStateReporter, _ReplicaExchangeSampler, _XTCMultiStateReporter


class REMDRunner(SimulationRunner):
    """Replica exchange molecular dynamics simulation runner.

    This class handles both new REMD simulations and extending existing
    REMD simulations using openmmtools. It uses the 'append' parameter
    from the input file or command line to determine whether to start
    fresh or continue from existing output.

    Attributes:
        config: Parsed simulation configuration object.
        platform: Configured OpenMM Platform.
        platform_properties: Dictionary of platform properties.
        replica_params: Dictionary containing replica exchange parameters.
        sampler: ReplicaExchangeSampler object (set after initialization).
    """

    def __init__(
        self,
        inpfile: str,
        replica_params: dict[str, Any] | None = None,
        output_data: str = "output.nc",
    ) -> None:
        """Initialize the REMD runner with input file and parameters.

        Args:
            inpfile: Path to input parameter file.
            replica_params: Dictionary containing replica exchange parameters
                with keys 'molname', 'beta', 'C1', and 'C2'. If None, must be
                set before running.
            output_data: Path to output NetCDF file for simulation data.

        Raises:
            ImportError: If openmmtools is not installed.
        """
        # Lazy import openmmtools to avoid import-time warnings
        global_context_cache, MultiStateReporter, ReplicaExchangeSampler, _ = _import_remd_dependencies()

        super().__init__(inpfile)
        self.output_data = output_data
        self.unsampled_topfiles = self.config.remd_unsampled_topfiles
        self.replica_params = replica_params
        self.sampler: Any | None = None

        # Set global context cache
        global_context_cache.platform = self.platform
        global_context_cache.platform_properties = self.platform_properties

    def detect_existing_output(self) -> bool:
        """Detect if existing REMD output file is present.

        Checks for the existence of the NetCDF output file.

        Returns:
            True if output file exists and can be extended, False otherwise.
        """
        return os.path.isfile(self.output_data)

    def _set_replica_params(
        self,
        replica_params: dict[str, Any],
    ) -> None:
        """Set replica exchange parameters.

        Args:
            replica_params: Dictionary containing 'molname', 'beta', 'C1', 'C2'.
        """
        self.replica_params = replica_params

    def _build_systems(
        self,
    ) -> tuple[list[mm.System], list[mm.System] | None]:
        """Build OpenMM systems for all replicas.

        Creates a system for each replica with different parameters.

        Returns:
            Tuple of (system_list, unsampled_system_list).
        """
        if self.replica_params is None:
            raise ValueError("Replica parameters not set. Call _set_replica_params first.")

        # Load positions and box vectors
        conf, box_vectors = load_structure(self.config.input)
        positions = conf.getPositions()
        velocities = None

        if self.config.ichk:
            print(f"\n[Checkpoint] {self.config.ichk}")
            with open(self.config.ichk, "r") as f:
                states = mm.XmlSerializer.deserialize(f.read())
            positions = states.getPositions()
            velocities = states.getVelocities()
            box_vectors = states.getPeriodicBoxVectors()

        # Build systems for each replica
        system_list: list[mm.System] = []
        target_mol: str = self.replica_params["molname"]

        for i, (beta, C1, C2) in enumerate(zip(
            self.replica_params["beta"],
            self.replica_params["C1"],
            self.replica_params["C2"],
        )):
            print(f"\n  [Replica {i}] beta={beta}, C1={C1}, C2={C2}")
            # Load topology
            top = MartiniTopFile(
                self.config.topol,
                periodicBoxVectors=box_vectors,
                defines=self.config.defines,
            )

            # Set replica exchange parameters
            top.moleculeTypes[target_mol].multiple_basin = [
                ["True", "exp", "2", str(beta), str(C1), str(C2)]
            ]

            # Create system
            system = top.create_system(
                nonbonded_cutoff=self.config.nonbonded_cutoff * u.nanometer,
                epsilon_r=self.config.epsilon_r,
            )

            # Check system charges
            if top.charges != 0:
                print(f"[Warning] System charge: {top.charges} (expected 0)")

            # Check number of atoms
            assert len(positions) == top.topology.getNumAtoms(), (
                f"Error: Number of atoms in {self.config.input} is not the same as "
                f"that from {self.config.topol}!"
            )

            # Add restraints
            if self.config.gen_rest == "yes":
                generate_restraints(
                    self.config.input,
                    self.config.atomname,
                    self.config.fc,
                    self.config.gen_rest_file,
                )
            if self.config.rest == "yes":
                system = add_restraints(system, self.config)

            # Add plumed (not supported)
            if self.config.plumed == "yes":
                raise ValueError("Error: Plumed is not supported!")

            # Add barostat
            if self.config.pcouple == "yes":
                if self.config.p_type == "isotropic":
                    barostat = mm.MonteCarloBarostat(
                        self.config.p_ref * u.bar,
                        self.config.temp * u.kelvin,
                        self.config.p_freq,
                    )
                elif self.config.p_type == "membrane":
                    barostat = mm.MonteCarloMembraneBarostat(
                        self.config.p_ref * u.bar,
                        self.config.p_tens * u.bar * u.nanometers,
                        self.config.temp * u.kelvin,
                        self.config.p_XYMode,
                        self.config.p_ZMode,
                        self.config.p_freq,
                    )
                else:
                    raise Exception(
                        "Unsupported barostat type: ", self.config.p_type
                    )
                system.addForce(barostat)

            system_list.append(system)

        # Build unsampled systems if provided
        unsampled_system_list: list[mm.System] | None = None
        if self.unsampled_topfiles is not None:
            unsampled_system_list = []
            for topfile in self.unsampled_topfiles:
                top = MartiniTopFile(
                    topfile,
                    periodicBoxVectors=box_vectors,
                    defines=self.config.defines,
                )

                system = top.create_system(
                    nonbonded_cutoff=self.config.nonbonded_cutoff * u.nanometer,
                    epsilon_r=self.config.epsilon_r,
                )

                if top.charges != 0:
                    print(
                        f"Warning: The charges of the system are {top.charges} instead of 0."
                    )

                assert len(positions) == top.topology.getNumAtoms(), (
                    f"Error: Number of atoms in {self.config.input} is not the same as "
                    f"that from {self.config.topol}!"
                )

                if self.config.gen_rest == "yes":
                    generate_restraints(
                        self.config.input,
                        self.config.atomname,
                        self.config.fc,
                        self.config.gen_rest_file,
                    )
                if self.config.rest == "yes":
                    system = add_restraints(system, self.config)

                if self.config.plumed == "yes":
                    raise ValueError("Error: Plumed is not supported!")

                if self.config.pcouple == "yes":
                    if self.config.p_type == "isotropic":
                        barostat = mm.MonteCarloBarostat(
                            self.config.p_ref * u.bar,
                            self.config.temp * u.kelvin,
                            self.config.p_freq,
                        )
                    elif self.config.p_type == "membrane":
                        barostat = mm.MonteCarloMembraneBarostat(
                            self.config.p_ref * u.bar,
                            self.config.p_tens * u.bar * u.nanometers,
                            self.config.temp * u.kelvin,
                            self.config.p_XYMode,
                            self.config.p_ZMode,
                            self.config.p_freq,
                        )
                    else:
                        raise Exception(
                            "Unsupported barostat type: ", self.config.p_type
                        )
                    system.addForce(barostat)

                unsampled_system_list.append(system)

        return system_list, unsampled_system_list

    def run(self) -> None:
        """Execute a new REMD simulation.

        This method sets up the replica exchange simulation from scratch
        and runs for the specified number of iterations.
        """
        # Lazy import openmmtools to avoid import-time warnings
        global_context_cache, MultiStateReporter, ReplicaExchangeSampler, _ = _import_remd_dependencies()
        import openmmtools.states
        import openmmtools.mcmc

        start_time = datetime.datetime.now()

        # Load parameters
        print("Loading parameters")

        # Build systems
        system_list, unsampled_system_list = self._build_systems()

        report_time("REMD Setup", start_time)

        # Production MD
        if self.config.nstep > 0:
            total_simulation_time = self.config.nstep * self.config.dt * u.picosecond
            simulation_time_step = self.config.dt * u.picosecond

            simulation_steps = int(np.floor(total_simulation_time / simulation_time_step))
            exchange_frequency = self.config.exc_freq
            exchange_attempts = int(np.floor(simulation_steps / exchange_frequency))

            sampler_states: list[openmmtools.states.SamplerState] = []
            thermodynamic_states: list[openmmtools.states.ThermodynamicState] = []

            # Get positions and box vectors
            conf, box_vectors = load_structure(self.config.input)
            positions = conf.getPositions()
            velocities = None

            if self.config.ichk:
                with open(self.config.ichk, "r") as f:
                    states = mm.XmlSerializer.deserialize(f.read())
                positions = states.getPositions()
                velocities = states.getVelocities()
                box_vectors = states.getPeriodicBoxVectors()

            # Define thermodynamic states with replica-specific temperatures
            for i, system in enumerate(system_list):
                thermodynamic_state = openmmtools.states.ThermodynamicState(
                    system=system, temperature=self.config.replica_temp[i] * u.kelvin
                )
                thermodynamic_states.append(thermodynamic_state)
                sampler_states.append(
                    openmmtools.states.SamplerState(
                        positions, velocities, box_vectors=box_vectors
                    )
                )

            unsampled_thermodynamic_states: list[
                openmmtools.states.ThermodynamicState
            ] | None = None
            if unsampled_system_list is not None:
                unsampled_thermodynamic_states = []
                for system in unsampled_system_list:
                    thermodynamic_state = openmmtools.states.ThermodynamicState(
                        system=system, temperature=self.config.replica_temp[0] * u.kelvin
                    )
                    unsampled_thermodynamic_states.append(thermodynamic_state)

            # Create and configure simulation object
            if not self.config.const_tol:
                self.config.const_tol = 1e-5
            move = openmmtools.mcmc.LangevinDynamicsMove(
                timestep=simulation_time_step,
                collision_rate=self.config.fric_coeff / u.picosecond,
                n_steps=exchange_frequency,
                reassign_velocities=self.config.remd_reassign_velocities,
                constraint_tolerance=self.config.const_tol,
            )

            self.sampler = ReplicaExchangeSampler(
                mcmc_moves=move,
                number_of_iterations=exchange_attempts,
                replica_mixing_scheme=self.config.remd_mixing_scheme,
                online_analysis_interval=self.config.remd_online_analysis_interval,
            )

            # Remove existing output file if present
            if os.path.exists(self.output_data):
                os.remove(self.output_data)

            # Choose reporter based on configuration
            _, _, _, XTCMultiStateReporter = _import_remd_dependencies()
            if self.config.remd_xtc_output == 'yes' and XTCMultiStateReporter is not None:
                n_replicas = len(system_list)
                reporter = XTCMultiStateReporter(
                    self.output_data,
                    checkpoint_interval=self.config.remd_checkpoint_interval,
                    xtc_dir=self.config.remd_xtc_dir,
                    total_iterations=exchange_attempts,
                    n_replicas=n_replicas,
                    exc_freq=self.config.exc_freq,
                    dt=self.config.dt,
                )
            else:
                reporter = MultiStateReporter(
                    self.output_data,
                    checkpoint_interval=self.config.remd_checkpoint_interval
                )
            if unsampled_thermodynamic_states is not None:
                self.sampler.create(
                    thermodynamic_states,
                    sampler_states,
                    reporter,
                    unsampled_thermodynamic_states=unsampled_thermodynamic_states,
                )
            else:
                self.sampler.create(
                    thermodynamic_states, sampler_states, reporter
                )

            total_steps = exchange_attempts * exchange_frequency
            print(f"\n[Production REMD]")
            print(f"  Time step: {simulation_time_step}")
            print(f"  Iterations: {exchange_attempts} ({total_steps:,} steps)", flush=True)

            self.sampler.run()

    def extend(self, n_iterations: int | None = None) -> None:
        """Continue an existing REMD simulation.

        This method loads the existing REMD simulation from the NetCDF
        file and continues running.

        Args:
            n_iterations: Number of additional iterations to run. If None,
                runs until the original target is reached or indefinitely.
        """
        # Lazy import openmmtools to avoid import-time warnings
        _, _, ReplicaExchangeSampler, XTCMultiStateReporter = _import_remd_dependencies()

        if not os.path.exists(self.output_data):
            raise FileNotFoundError(
                f"Cannot extend REMD: output file {self.output_data} does not exist. "
                f"Use run() to start a new simulation."
            )

        print(f"\n[Extend REMD] {self.output_data}")

        # Calculate total simulation parameters first (needed for reporter)
        total_simulation_time = self.config.nstep * self.config.dt * u.picosecond
        simulation_time_step = self.config.dt * u.picosecond
        simulation_steps = int(np.floor(total_simulation_time / simulation_time_step))
        exchange_frequency = self.config.exc_freq
        exchange_attempts = int(np.floor(simulation_steps / exchange_frequency))

        # Load sampler from storage
        if self.config.remd_xtc_output == 'yes' and XTCMultiStateReporter is not None:
            print(f"  [XTC Mode] Resuming with XTC separate storage")
            n_replicas = self.config.replica_count
            xtc_reporter = XTCMultiStateReporter(
                self.output_data,
                checkpoint_interval=self.config.remd_checkpoint_interval,
                xtc_dir=self.config.remd_xtc_dir,
                total_iterations=exchange_attempts,
                n_replicas=n_replicas,
                exc_freq=self.config.exc_freq,
                dt=self.config.dt,
            )
            self.sampler = ReplicaExchangeSampler.from_storage(xtc_reporter)
        else:
            self.sampler = ReplicaExchangeSampler.from_storage(self.output_data)

        # Set online analysis interval
        self.sampler._online_analysis_interval = self.config.remd_online_analysis_interval

        # Determine number of iterations to run
        current_iteration = self.sampler.iteration
        current_steps = current_iteration * exchange_frequency
        
        if n_iterations is not None:
            target_iteration = current_iteration + n_iterations
            target_steps = target_iteration * exchange_frequency
            print(f"  Progress: iteration {current_iteration} → {target_iteration} ({current_steps:,} → {target_steps:,} steps)")
            self.sampler.extend(n_iterations=n_iterations)
        else:
            # Run until original target is reached
            n_iter_remain = exchange_attempts - current_iteration
            if n_iter_remain <= 0:
                total_steps = exchange_attempts * exchange_frequency
                print(f"  Simulation already complete: {current_iteration}/{exchange_attempts} iterations ({current_steps:,}/{total_steps:,} steps)")
                return
            target_iteration = current_iteration + n_iter_remain
            target_steps = target_iteration * exchange_frequency
            print(f"  Progress: iteration {current_iteration} → {target_iteration} ({current_steps:,} → {target_steps:,} steps)")
            self.sampler.extend(n_iterations=n_iter_remain)
