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

try:
    from openmmtools.cache import global_context_cache
    from openmmtools.multistate import (
        MultiStateReporter,
        ReplicaExchangeSampler,
    )
    HAS_OPENMMTOOLS = True
except ImportError:
    HAS_OPENMMTOOLS = False
    global_context_cache = None
    MultiStateReporter = None
    ReplicaExchangeSampler = None

from ctgomartini.simulation.base import (
    SimulationRunner,
    report_time,
    load_structure,
    generate_restraints,
    add_restraints,
)

# Import MartiniTopFile here to avoid circular imports
from ctgomartini.topology import MartiniTopFile


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
        if not HAS_OPENMMTOOLS:
            raise ImportError(
                "openmmtools is required for REMD simulations. "
                "Install it with: pip install openmmtools"
            )

        super().__init__(inpfile)
        self.output_data = output_data
        self.unsampled_topfiles = self.config.remd_unsampled_topfiles
        self.replica_params = replica_params
        self.sampler: ReplicaExchangeSampler | None = None

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
            print(f"\nLoad checkpoint file: {self.config.ichk}")
            with open(self.config.ichk, "r") as f:
                states = mm.XmlSerializer.deserialize(f.read())
            positions = states.getPositions()
            velocities = states.getVelocities()
            box_vectors = states.getPeriodicBoxVectors()

        # Build systems for each replica
        system_list: list[mm.System] = []
        target_mol: str = self.replica_params["molname"]

        for beta, C1, C2 in zip(
            self.replica_params["beta"],
            self.replica_params["C1"],
            self.replica_params["C2"],
        ):
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
                print(
                    f"Warning: The charges of the system are {top.charges} instead of 0."
                )

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
        import openmmtools.states
        import openmmtools.mcmc

        start_time = datetime.datetime.now()

        # Load parameters
        print("Loading parameters")

        # Build systems
        system_list, unsampled_system_list = self._build_systems()

        # System Loading Finishes!
        print("\nLoading system finishes!")
        report_time(start_time)

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

            # Define thermodynamic states
            for system in system_list:
                thermodynamic_state = openmmtools.states.ThermodynamicState(
                    system=system, temperature=self.config.temp * u.kelvin
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
                        system=system, temperature=self.config.temp * u.kelvin
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

            print("Running OpenMM replica exchange simulation...")
            print(f"Time step: {simulation_time_step}")
            print(f"Iterations: {exchange_attempts}", flush=True)

            self.sampler.run()

    def extend(self, n_iterations: int | None = None) -> None:
        """Continue an existing REMD simulation.

        This method loads the existing REMD simulation from the NetCDF
        file and continues running.

        Args:
            n_iterations: Number of additional iterations to run. If None,
                runs until the original target is reached or indefinitely.
        """
        from openmmtools.multistate import ReplicaExchangeSampler

        if not os.path.exists(self.output_data):
            raise FileNotFoundError(
                f"Cannot extend REMD: output file {self.output_data} does not exist. "
                f"Use run() to start a new simulation."
            )

        print(f"Extending REMD simulation from {self.output_data}")

        # Create a new sampler and resume from storage
        self.sampler = ReplicaExchangeSampler.from_storage(self.output_data)

        # Determine number of iterations to run
        if n_iterations is not None:
            current_iteration = self.sampler.iteration
            target_iteration = current_iteration + n_iterations
            print(
                f"Current iteration: {current_iteration}, target: {target_iteration}"
            )
            self.sampler.run(n_iterations=n_iterations)
        else:
            # Run until completion (original target)
            print("Continuing REMD simulation...")
            self.sampler.run()
