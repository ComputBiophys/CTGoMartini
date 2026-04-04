"""
XTC MultiState Reporter module for CTGoMartini.

Provides XTCMultiStateReporter class with two modes:
- Standard mode: fully compatible with MultiStateReporter
- XTC separate mode: stores trajectory to XTC files, checkpoint keeps only latest frame
"""

from __future__ import annotations

import os
from typing import Any

try:
    from openmm import unit
    from openmm.app import XTCFile
except ImportError:
    from simtk import unit
    from simtk.app import XTCFile

import numpy as np

# Import parent class from openmmtools
try:
    from openmmtools.multistate import MultiStateReporter
    from openmmtools import states
except ImportError:
    MultiStateReporter = None  # Will be handled on delayed import


class XTCMultiStateReporter(MultiStateReporter):
    """
    Dual-mode MultiStateReporter with XTC separate storage support.

    Mode 1 (remd_xtc_output = 'no'):
        Standard behavior, fully compatible with MultiStateReporter

    Mode 2 (remd_xtc_output = 'yes'):
        XTC separate storage + lightweight checkpoint
        - XTC files: store all historical trajectories (one file per replica)
        - NetCDF checkpoint: keeps only latest frame (for recovery)

    Attributes:
        _xtc_output: Whether to enable XTC separate storage
        _xtc_dir: Directory for XTC files
        _xtc_handles: Cache of XTCFile handles
        _is_lite_mode: Whether using lightweight checkpoint
    """

    def __init__(
        self,
        storage: str,
        open_mode: str | None = None,
        checkpoint_interval: int = 50,
        checkpoint_storage: str | None = None,
        analysis_particle_indices: tuple = (),
        position_interval: int = 1,
        velocity_interval: int = 1,
        xtc_output: str = 'no',
        xtc_dir: str = 'xtc_trajs',
        **kwargs: Any,
    ) -> None:
        """
        Initialize XTCMultiStateReporter.

        Args:
            storage: Path to main NetCDF file
            open_mode: File open mode ('r', 'w', 'a')
            checkpoint_interval: Checkpoint write frequency (iterations)
            checkpoint_storage: Checkpoint filename (default auto)
            analysis_particle_indices: Atom indices for additional analysis
            position_interval: Position write interval
            velocity_interval: Velocity write interval
            xtc_output: 'yes' to enable XTC separate storage, 'no' for standard mode
            xtc_dir: Directory for XTC files (relative to storage directory)
            **kwargs: Additional arguments passed to parent class
        """
        super().__init__(
            storage=storage,
            open_mode=open_mode,
            checkpoint_interval=checkpoint_interval,
            checkpoint_storage=checkpoint_storage,
            analysis_particle_indices=analysis_particle_indices,
            position_interval=position_interval,
            velocity_interval=velocity_interval,
        )

        self._xtc_output = (xtc_output.lower() == 'yes')
        self._xtc_dir = xtc_dir
        self._xtc_handles: dict[int, XTCFile] = {}
        self._xtc_timestep: float | None = None  # Inferred from mcmc_moves

        if self._xtc_output:
            # XTC mode: create directory
            storage_dir = os.path.dirname(storage) or '.'
            self._xtc_dir = os.path.join(storage_dir, xtc_dir)
            os.makedirs(self._xtc_dir, exist_ok=True)
            self._is_lite_mode = True
        else:
            self._is_lite_mode = False

    def write_sampler_states(
        self,
        sampler_states: list,
        iteration: int,
    ) -> None:
        """
        Write sampler states.

        Behavior depends on _xtc_output mode:
        - Standard mode: call parent class method
        - XTC mode: separate storage to XTC and lightweight checkpoint

        Args:
            sampler_states: List of sampler states for all replicas
            iteration: Current iteration number
        """
        if self._xtc_output:
            self._write_sampler_states_xtc_mode(sampler_states, iteration)
        else:
            super().write_sampler_states(sampler_states, iteration)

    def read_sampler_states(
        self,
        iteration: int | None = None,
        analysis_particles_only: bool = False,
    ) -> list | None:
        """
        Read sampler states.

        Args:
            iteration: Iteration number (ignored in XTC mode)
            analysis_particles_only: Whether to read only analysis particles

        Returns:
            List of sampler states, or None if not available
        """
        if self._xtc_output:
            return self._read_sampler_states_xtc_mode(iteration)
        else:
            return super().read_sampler_states(iteration, analysis_particles_only)

    # -------------------------------------------------------------------------
    # XTC mode private methods
    # -------------------------------------------------------------------------

    def _write_sampler_states_xtc_mode(
        self,
        sampler_states: list,
        iteration: int,
    ) -> None:
        """
        XTC mode write.

        1. Write to XTC files (all historical trajectories)
        2. Write to lightweight NetCDF checkpoint (only latest frame + velocities)
        """
        # 1. Write to XTC (only on checkpoint interval)
        if self._on_checkpoint_interval(iteration):
            for replica_idx, state in enumerate(sampler_states):
                xtc = self._get_xtc_handle(replica_idx, state)
                xtc.writeModel(
                    state.positions,
                    periodicBoxVectors=state.box_vectors,
                )

        # 2. Write to lightweight NetCDF checkpoint
        self._write_lite_checkpoint(sampler_states, iteration)

    def _write_lite_checkpoint(
        self,
        sampler_states: list,
        iteration: int,
    ) -> None:
        """
        Write lightweight checkpoint (only latest frame).

        Checkpoint file has fixed dimension of 1, overwritten each time.
        """
        storage = self._storage_checkpoint

        # Lazy initialization: create fixed-size variables on first write
        if 'positions' not in storage.variables:
            self._init_lite_checkpoint_storage(sampler_states)

        n_replicas = len(sampler_states)

        # Overwrite index 0 (latest frame)
        for i, state in enumerate(sampler_states):
            # Positions (nm)
            pos = state.positions.value_in_unit(unit.nanometers)
            storage.variables['positions'][0, i, :, :] = pos

            # Velocities (nm/ps) - for exact recovery
            if state.velocities is not None:
                vel = state.velocities.value_in_unit(
                    unit.nanometer / unit.picosecond
                )
                storage.variables['velocities'][0, i, :, :] = vel

            # Box vectors
            if state.box_vectors is not None:
                box = state.box_vectors.value_in_unit(unit.nanometers)
                storage.variables['box_vectors'][0, i, :, :] = box

        # Atomic marker: sync, write current_iteration, sync
        self.sync()
        storage.variables['current_iteration'][0] = iteration
        self.sync()

    def _init_lite_checkpoint_storage(
        self,
        sampler_states: list,
    ) -> None:
        """Initialize fixed-size checkpoint storage."""
        storage = self._storage_checkpoint

        n_replicas = len(sampler_states)
        n_atoms = sampler_states[0].n_particles

        # Key: fixed size of 1 (not unlimited)
        storage.createDimension('iteration_lite', 1)
        storage.createDimension('replica', n_replicas)
        storage.createDimension('atom', n_atoms)
        storage.createDimension('spatial', 3)

        # Positions
        var_pos = storage.createVariable(
            'positions', 'f4',
            ('iteration_lite', 'replica', 'atom', 'spatial'),
            zlib=True,
        )
        var_pos.units = 'nm'

        # Velocities
        var_vel = storage.createVariable(
            'velocities', 'f4',
            ('iteration_lite', 'replica', 'atom', 'spatial'),
            zlib=True,
        )
        var_vel.units = 'nm/ps'

        # Box vectors
        var_box = storage.createVariable(
            'box_vectors', 'f4',
            ('iteration_lite', 'replica', 'spatial', 'spatial'),
            zlib=False,
        )
        var_box.units = 'nm'

        # Current iteration marker
        storage.createVariable('current_iteration', 'i4', ('scalar',))

    def _read_sampler_states_xtc_mode(
        self,
        iteration: int | None = None,
    ) -> list | None:
        """
        Read sampler states from lightweight checkpoint.

        Always returns latest saved frame (index 0).

        Returns:
            List of sampler states
        """
        storage = self._storage_checkpoint

        # Check if initialized
        if 'positions' not in storage.variables:
            return None

        sampler_states = []
        n_replicas = storage.dimensions['replica'].size

        for i in range(n_replicas):
            # Positions
            pos_data = storage.variables['positions'][0, i, :, :]
            positions = unit.Quantity(pos_data.data, unit.nanometers)

            # Velocities
            vel_data = storage.variables['velocities'][0, i, :, :]
            velocities = unit.Quantity(
                vel_data.data,
                unit.nanometer / unit.picosecond,
            )

            # Box vectors
            box_data = storage.variables['box_vectors'][0, i, :, :]
            if box_data.any():
                box_vectors = unit.Quantity(box_data.data, unit.nanometers)
            else:
                box_vectors = None

            # Create SamplerState
            sampler_states.append(
                states.SamplerState(
                    positions=positions,
                    velocities=velocities,
                    box_vectors=box_vectors,
                )
            )

        return sampler_states

    def _get_xtc_handle(
        self,
        replica_idx: int,
        state: states.SamplerState,
    ) -> XTCFile:
        """
        Get or create XTCFile handle.

        Args:
            replica_idx: Replica index
            state: Sampler state (for topology)

        Returns:
            XTCFile handle
        """
        if replica_idx not in self._xtc_handles:
            path = os.path.join(self._xtc_dir, f"replica_{replica_idx}.xtc")

            # Check if append needed
            append = os.path.exists(path) and os.path.getsize(path) > 0

            # Infer timestep from mcmc_moves
            if self._xtc_timestep is None:
                self._xtc_timestep = self._infer_timestep()

            # Need topology - infer from state
            topology = self._create_topology_from_state(state)

            self._xtc_handles[replica_idx] = XTCFile(
                path,
                topology,
                dt=self._xtc_timestep * unit.picosecond,
                interval=self._checkpoint_interval,
                append=append,
            )

        return self._xtc_handles[replica_idx]

    def _on_checkpoint_interval(self, iteration: int) -> bool:
        """Check if current iteration is on checkpoint interval."""
        return iteration % self._checkpoint_interval == 0

    def _infer_timestep(self) -> float:
        """
        Infer timestep from mcmc_moves.

        Returns:
            Timestep in ps
        """
        try:
            # Try to read mcmc_moves from analysis file
            moves = self.read_mcmc_moves()
            if moves and len(moves) > 0:
                timestep = moves[0].timestep
                if unit.is_quantity(timestep):
                    return timestep.value_in_unit(unit.picosecond)
                return timestep
        except Exception:
            pass

        # Default value
        return 0.02  # 20 fs = 0.02 ps (Martini default)

    def _create_topology_from_state(
        self,
        state: states.SamplerState,
    ) -> Any:
        """
        Create OpenMM Topology from SamplerState.

        This is a simplified version, actual use may require more complete topology.
        """
        from openmm.app import Topology as OpenMMTopology, Element

        n_atoms = state.n_particles
        top = OpenMMTopology()
        chain = top.addChain()
        residue = top.addResidue('UNK', chain)

        for i in range(n_atoms):
            # Use carbon as placeholder
            element = Element.getBySymbol('C')
            top.addAtom(f'AT{i}', element, residue)

        # Set periodic box
        if state.box_vectors is not None:
            top.setPeriodicBoxVectors(state.box_vectors)

        return top

    def close(self) -> None:
        """Close reporter and clean up XTC handles."""
        self._xtc_handles.clear()
        super().close()

    def read_checkpoint_iterations(self) -> np.ndarray:
        """
        Read all checkpoint iterations.

        Lightweight mode has only one valid checkpoint.
        """
        if not self._xtc_output:
            return super().read_checkpoint_iterations()

        # XTC mode: return current_iteration
        try:
            storage = self._storage_checkpoint
            current = int(storage.variables['current_iteration'][0])
            return np.array([current])
        except (KeyError, AttributeError):
            return np.array([0])

    def read_last_iteration(self, last_checkpoint: bool = True) -> int:
        """
        Read last saved iteration.

        Args:
            last_checkpoint: Whether to return last checkpoint iteration

        Returns:
            Last iteration number
        """
        if not self._xtc_output:
            return super().read_last_iteration(last_checkpoint)

        # XTC mode: read from current_iteration
        try:
            storage = self._storage_checkpoint
            return int(storage.variables['current_iteration'][0])
        except (KeyError, AttributeError):
            return 0
