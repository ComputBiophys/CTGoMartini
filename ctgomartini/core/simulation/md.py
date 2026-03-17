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
from typing import Any

import openmm as mm
import openmm.unit as u
from openmm.app import (
    DCDReporter,
    Simulation,
    StateDataReporter,
    XTCReporter,
    CheckpointReporter,
)

from .base import (
    BaseRunner,
    report_time,
    backup_file,
    write_output,
    write_checkpoint,
    cleanup,
)


class MDRunner(BaseRunner):
    """Standard molecular dynamics simulation runner.

    This class handles both new MD simulations and appending to existing
    simulations. It uses the 'append' parameter from the input file to
    determine whether to start fresh or continue from existing output.

    Attributes:
        inputs: Parsed input parameters object.
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
        if self.inputs.odcd and os.path.isfile(self.inputs.odcd):
            has_traj = True
        if self.inputs.oxtc and os.path.isfile(self.inputs.oxtc):
            has_traj = True

        # Check for checkpoint file
        has_chk = self.inputs.ochk and os.path.isfile(self.inputs.ochk)

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
            self.inputs.temp * u.kelvin,
            self.inputs.fric_coeff / u.picosecond,
            self.inputs.dt * u.picosecond,
        )

        if self.inputs.const_tol:
            integrator.setConstraintTolerance(self.inputs.const_tol)
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
        from .base import load_structure

        conf, _ = load_structure(self.inputs.input)
        simulation.context.setPositions(conf.getPositions())

        # Load checkpoint if provided
        if self.inputs.ichk:
            try:
                with open(self.inputs.ichk, "r") as f:
                    simulation.context.setState(
                        mm.XmlSerializer.deserialize(f.read())
                    )
            except Exception:
                with open(self.inputs.ichk, "rb") as f:
                    simulation.context.loadCheckpoint(f.read())
            print(f"\nLoaded checkpoint file: {self.inputs.ichk}")

        print("\nLoading system finishes!")
        report_time(start_time)

        return simulation

    def _setup_reporters(self, simulation: Simulation, append: bool) -> None:
        """Set up simulation reporters.

        Args:
            simulation: The OpenMM Simulation object.
            append: Whether to append to existing trajectory files.
        """
        # Determine trajectory format
        if self.inputs.odcd and not self.inputs.oxtc:
            traj_reporter = DCDReporter
            traj_file = self.inputs.odcd
        elif self.inputs.oxtc and not self.inputs.odcd:
            traj_reporter = XTCReporter
            traj_file = self.inputs.oxtc
        else:
            raise ValueError(
                "Error: Please specify either odcd or oxtc, not both!"
            )

        # Setup trajectory reporter
        if self.inputs.nstdcd > 0:
            if append:
                simulation.reporters.append(
                    traj_reporter(traj_file, self.inputs.nstdcd, append=True)
                )
            else:
                backup_file(traj_file)
                simulation.reporters.append(
                    traj_reporter(traj_file, self.inputs.nstdcd)
                )

            # Setup checkpoint reporter
            backup_file(self.inputs.ochk)
            simulation.reporters.append(
                CheckpointReporter(
                    self.inputs.ochk, self.inputs.nstdcd, writeState=True
                )
            )

        # Setup state data reporter
        simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                self.inputs.nstout,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=self.inputs.nstep,
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
        print(f"\nInitial system energy: {energy:.3f} kJ/mol")

        # Energy minimization
        if self.inputs.mini_nstep > 0:
            print(f"\nEnergy minimization: {self.inputs.mini_nstep} steps")
            self.simulation.minimizeEnergy(
                tolerance=self.inputs.mini_Tol,
                maxIterations=self.inputs.mini_nstep,
            )
            energy = (
                self.simulation.context.getState(getEnergy=True)
                .getPotentialEnergy()
                .value_in_unit(u.kilojoule_per_mole)
            )
            print(f"Minimized system energy: {energy:.3f} kJ/mol")
            report_time(start_time)

        # Generate initial velocities
        if self.inputs.gen_vel == "yes":
            print(f"\nGenerating initial velocities: {self.inputs.gen_temp} K")
            if self.inputs.gen_seed:
                self.simulation.context.setVelocitiesToTemperature(
                    self.inputs.gen_temp * u.kelvin, self.inputs.gen_seed
                )
            else:
                self.simulation.context.setVelocitiesToTemperature(
                    self.inputs.gen_temp * u.kelvin
                )

        # Production MD
        if self.inputs.nstep > 0:
            start_time = datetime.datetime.now()

            # Set starting point
            begin_step = self.inputs.b_step
            self.simulation.context.setStepCount(begin_step)
            self.simulation.context.setTime(self.inputs.dt * begin_step)

            end_step = self.inputs.nstep
            remaining_steps = end_step - begin_step

            print(
                f"\nMD run: begin {begin_step}, end {end_step}, total {remaining_steps}"
            )

            # Setup reporters (not appending for new simulation)
            self._setup_reporters(self.simulation, append=False)

            # Run simulation
            try:
                self.simulation.step(remaining_steps)
            except KeyboardInterrupt:
                print("\nSimulation interrupted!")
                report_time(start_time)
                cleanup(signal.SIGINT, self.simulation, self.inputs)

        # Write final output
        if self.inputs.output:
            from .base import load_structure

            write_output(self.inputs.output, self.simulation, self.inputs.input)
        if self.inputs.output_pdb:
            from .base import load_structure

            write_output(
                self.inputs.output_pdb, self.simulation, self.inputs.input
            )

        if self.inputs.ochk:
            write_checkpoint(self.simulation, self.inputs.ochk)

        print("\nSimulation finished!")
        report_time(start_time)

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
        if self.inputs.ochk and os.path.isfile(self.inputs.ochk):
            try:
                with open(self.inputs.ochk, "r") as f:
                    self.simulation.context.setState(
                        mm.XmlSerializer.deserialize(f.read())
                    )
                print(f"\nLoaded checkpoint for continuation: {self.inputs.ochk}")
            except Exception:
                with open(self.inputs.ochk, "rb") as f:
                    self.simulation.context.loadCheckpoint(f.read())
                print(f"\nLoaded binary checkpoint for continuation: {self.inputs.ochk}")

        # Determine steps to run
        current_step = self.simulation.context.getStepCount()
        target_step = self.inputs.nstep

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
            print("\nSimulation interrupted!")
            report_time(start_time)
            cleanup(signal.SIGINT, self.simulation, self.inputs)

        # Write final output
        if self.inputs.output:
            from .base import load_structure

            write_output(self.inputs.output, self.simulation, self.inputs.input)
        if self.inputs.output_pdb:
            from .base import load_structure

            write_output(
                self.inputs.output_pdb, self.simulation, self.inputs.input
            )

        if self.inputs.ochk:
            write_checkpoint(self.simulation, self.inputs.ochk)

        print("\nSimulation extension finished!")
        report_time(start_time)
