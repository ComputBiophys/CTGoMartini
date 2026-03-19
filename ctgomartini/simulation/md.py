"""
Standard MD simulation runner module for CTGoMartini.

Provides the MDRunner class for running standard molecular dynamics simulations
with support for new simulations and appending to existing simulations.
"""

from __future__ import annotations

import datetime
import os
import signal
import sys

import openmm as mm
import openmm.unit as u
from openmm.app import (
    DCDReporter,
    Simulation,
    StateDataReporter,
    XTCReporter,
    CheckpointReporter,
)

from ctgomartini.simulation.base import (
    SimulationRunner,
    report_time,
    backup_file,
    write_output,
    write_checkpoint,
    _cleanup_handler,
    load_structure,
)


class MDRunner(SimulationRunner):
    """Standard molecular dynamics simulation runner.

    This class handles both new MD simulations and appending to existing
    simulations. It uses the 'append' parameter from the input file to
    determine whether to start fresh or continue from existing output.

    Attributes:
        config: Parsed simulation configuration object.
        platform: Configured OpenMM Platform.
        platform_properties: Dictionary of platform properties.
        simulation: OpenMM Simulation object.
    """

    def __init__(self, inpfile: str) -> None:
        """Initialize the MD runner with input file.

        Args:
            inpfile: Path to input parameter file.
        """
        super().__init__(inpfile)
        self.simulation: Simulation | None = None

    def detect_existing_output(self) -> bool:
        """Detect if existing MD output files are present.

        Checks for the existence of trajectory files (DCD or XTC) and
        checkpoint files.

        Returns:
            True if output files exist and can be appended to, False otherwise.
        """
        # Check for trajectory file
        has_traj = False
        if self.config.odcd and os.path.isfile(self.config.odcd):
            has_traj = True
        if self.config.oxtc and os.path.isfile(self.config.oxtc):
            has_traj = True

        # Check for checkpoint file
        has_chk = self.config.ochk and os.path.isfile(self.config.ochk)

        return has_traj or has_chk

    def _setup_simulation(self) -> Simulation:
        """Set up the OpenMM simulation.

        Creates the system, integrator, and simulation object.

        Returns:
            Configured OpenMM Simulation object.
        """
        import datetime

        start_time = datetime.datetime.now()

        # Generate system and topology
        system, topology, box_vectors = self.generate_topology()

        # Create integrator
        integrator = mm.LangevinIntegrator(
            self.config.temp * u.kelvin,
            self.config.fric_coeff / u.picosecond,
            self.config.dt * u.picosecond,
        )

        if self.config.const_tol:
            integrator.setConstraintTolerance(self.config.const_tol)
            print("Set Constraint Tolerance:", integrator.getConstraintTolerance())

        # Create simulation
        from openmm.app import Simulation

        simulation = Simulation(
            topology,
            system,
            integrator,
            self.platform,
            self.platform_properties,
        )

        # Load structure
        conf, _ = load_structure(self.config.input)
        simulation.context.setPositions(conf.getPositions())

        # Load checkpoint if provided
        if self.config.ichk:
            try:
                with open(self.config.ichk, "r") as f:
                    simulation.context.setState(
                        mm.XmlSerializer.deserialize(f.read())
                    )
            except Exception:
                with open(self.config.ichk, "rb") as f:
                    simulation.context.loadCheckpoint(f.read())
            print(f"\nLoaded checkpoint file: {self.config.ichk}")

        report_time("Simulation Setup", start_time)

        return simulation

    def _setup_reporters(self, simulation: Simulation, append: bool) -> None:
        """Set up simulation reporters.

        Args:
            simulation: The OpenMM Simulation object.
            append: Whether to append to existing trajectory files.
        """
        # Determine trajectory format
        if self.config.odcd and not self.config.oxtc:
            traj_reporter = DCDReporter
            traj_file = self.config.odcd
        elif self.config.oxtc and not self.config.odcd:
            traj_reporter = XTCReporter
            traj_file = self.config.oxtc
        else:
            raise ValueError(
                "Error: Please specify either odcd or oxtc, not both!"
            )

        # Setup trajectory reporter
        if self.config.nstdcd > 0:
            if append:
                simulation.reporters.append(
                    traj_reporter(traj_file, self.config.nstdcd, append=True)
                )
            else:
                backup_file(traj_file)
                simulation.reporters.append(
                    traj_reporter(traj_file, self.config.nstdcd)
                )

            # Setup checkpoint reporter
            backup_file(self.config.ochk)
            simulation.reporters.append(
                CheckpointReporter(
                    self.config.ochk, self.config.nstdcd, writeState=True
                )
            )

        # Setup state data reporter
        simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                self.config.nstout,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=self.config.nstep,
                separator="\t",
            )
        )

    def run(self) -> None:
        """Execute a new MD simulation.

        This method sets up the simulation from scratch and runs for the
        specified number of steps.
        """
        start_time = datetime.datetime.now()

        # Setup simulation
        self.simulation = self._setup_simulation()

        # Calculate initial energy
        state = self.simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        print(f"\n[Initial Energy] {energy:.3f} kJ/mol")

        # Energy minimization
        if self.config.mini_nstep > 0:
            print(f"\n[Minimization] {self.config.mini_nstep} steps")
            self.simulation.minimizeEnergy(
                tolerance=self.config.mini_Tol,
                maxIterations=self.config.mini_nstep,
            )
            energy = (
                self.simulation.context.getState(getEnergy=True)
                .getPotentialEnergy()
                .value_in_unit(u.kilojoule_per_mole)
            )
            print(f"  Energy: {energy:.3f} kJ/mol")
            report_time("Minimization", start_time)

        # Generate initial velocities
        if self.config.gen_vel == "yes":
            print(f"\n[Initial Velocities] T = {self.config.gen_temp} K")
            if self.config.gen_seed:
                self.simulation.context.setVelocitiesToTemperature(
                    self.config.gen_temp * u.kelvin, self.config.gen_seed
                )
            else:
                self.simulation.context.setVelocitiesToTemperature(
                    self.config.gen_temp * u.kelvin
                )

        # Production MD
        if self.config.nstep > 0:
            start_time = datetime.datetime.now()

            # Set starting point
            begin_step = self.config.b_step
            self.simulation.context.setStepCount(begin_step)
            self.simulation.context.setTime(self.config.dt * begin_step)

            end_step = self.config.nstep
            remaining_steps = end_step - begin_step

            print(f"\n[Production MD] {remaining_steps} steps ({begin_step} → {end_step})")

            # Setup reporters (not appending for new simulation)
            self._setup_reporters(self.simulation, append=False)

            # Run simulation
            try:
                self.simulation.step(remaining_steps)
            except KeyboardInterrupt:
                print("\n[Interrupted]")
                report_time("Simulation", start_time)
                _cleanup_handler(signal.SIGINT, self.simulation, self.config)

        # Write final output
        if self.config.output:
            write_output(self.config.output, self.simulation, self.config.input)
        if self.config.output_pdb:
            write_output(
                self.config.output_pdb, self.simulation, self.config.input
            )

        if self.config.ochk:
            write_checkpoint(self.simulation, self.config.ochk)

        print(f"\n[Completed]")
        report_time("Total", start_time)

    def extend(self, n_iterations: int | None = None) -> None:
        """Continue an existing MD simulation (append mode).

        This method loads the existing checkpoint and continues the simulation,
        appending to the existing trajectory files.

        Args:
            n_iterations: Number of additional steps to run. If None, runs
                until the target nstep is reached.
        """
        start_time = datetime.datetime.now()

        # Setup simulation
        self.simulation = self._setup_simulation()

        # Load from checkpoint for continuation
        if self.config.ochk and os.path.isfile(self.config.ochk):
            try:
                with open(self.config.ochk, "r") as f:
                    self.simulation.context.setState(
                        mm.XmlSerializer.deserialize(f.read())
                    )
                print(f"\nLoaded checkpoint for continuation: {self.config.ochk}")
            except Exception:
                with open(self.config.ochk, "rb") as f:
                    self.simulation.context.loadCheckpoint(f.read())
                print(f"\nLoaded binary checkpoint for continuation: {self.config.ochk}")

        # Determine steps to run
        current_step = self.simulation.context.getStepCount()
        target_step = self.config.nstep

        if n_iterations is not None:
            remaining_steps = n_iterations
            end_step = current_step + n_iterations
        else:
            remaining_steps = target_step - current_step
            end_step = target_step

        if remaining_steps <= 0:
            print(f"\nSimulation already at step {current_step}, target is {target_step}")
            print("Nothing to extend.")
            return

        print(f"\nMD extend: current {current_step}, target {end_step}, remaining {remaining_steps}")

        # Setup reporters (appending mode)
        self._setup_reporters(self.simulation, append=True)

        # Run simulation
        try:
            self.simulation.step(remaining_steps)
        except KeyboardInterrupt:
            print("\n[Interrupted]")
            report_time("Extension", start_time)
            _cleanup_handler(signal.SIGINT, self.simulation, self.config)

        # Write final output
        if self.config.output:
            write_output(self.config.output, self.simulation, self.config.input)
        if self.config.output_pdb:
            write_output(
                self.config.output_pdb, self.simulation, self.config.input
            )

        if self.config.ochk:
            write_checkpoint(self.simulation, self.config.ochk)

        print("\n[Extension Complete]")
        report_time("Extension", start_time)
