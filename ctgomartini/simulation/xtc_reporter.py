"""
XTC MultiState Reporter module for CTGoMartini.

Provides XTCMultiStateReporter class for XTC separate storage:
- XTC files: store all historical trajectories (one file per replica)
- NetCDF checkpoint: keeps only latest frame for recovery (fixed size)

This significantly reduces checkpoint file size while maintaining exact recovery.
"""

from __future__ import annotations

import os
import time
from typing import Any


from openmm import unit, Vec3
from openmm.app import XTCFile, Topology as OpenMMTopology, Element

import numpy as np

from openmmtools.multistate import MultiStateReporter
from openmmtools import states


class XTCMultiStateReporter(MultiStateReporter):
    """
    MultiStateReporter with XTC separate storage.

    Stores trajectory data to XTC files (one per replica) and keeps only
    the latest frame in NetCDF checkpoint for recovery.

    Attributes:
        _xtc_dir: Directory for XTC files
        _xtc_handles: Cache of XTCFile handles (replica_idx -> XTCFile)
        _xtc_timestep: Timestep in ps for XTC files
        _total_iterations: Total number of iterations for progress display
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
        xtc_dir: str = 'xtc_trajs',
        total_iterations: int | None = None,
        progress_interval: int | None = None,
        n_replicas: int = 1,
        exc_freq: int | None = None,
        dt: float | None = None,
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
            xtc_dir: Directory for XTC files (relative to storage directory)
            total_iterations: Total number of iterations for progress display
            progress_interval: Progress output interval (iterations, defaults to checkpoint_interval)
            n_replicas: Number of replicas for ns/day calculation
            exc_freq: Exchange frequency (MD steps per iteration) for ns/day calculation
            dt: Time step in ps for ns/day calculation
            **kwargs: Additional arguments passed to parent class
        """
        # Initialize XTC handles before super().__init__() 
        # because close() may be called during parent init
        self._xtc_handles: dict[int, XTCFile] = {}
        self._xtc_timestep: float | None = None
        self._xtc_interval: int | None = None  # MD steps between checkpoints
        
        # Progress tracking
        self._total_iterations: int | None = total_iterations
        self._progress_interval: int = progress_interval if progress_interval is not None else checkpoint_interval
        self._progress_header_printed: bool = False
        self._progress_start_time: float | None = None
        self._progress_start_iter: int | None = None
        
        # For ns/day calculation
        self._n_replicas: int = n_replicas
        self._exc_freq: int | None = exc_freq
        self._dt: float | None = dt
        
        # Setup XTC directory
        storage_dir = os.path.dirname(storage) or '.'
        self._xtc_dir = os.path.join(storage_dir, xtc_dir)
        os.makedirs(self._xtc_dir, exist_ok=True)

        super().__init__(
            storage=storage,
            open_mode=open_mode,
            checkpoint_interval=checkpoint_interval,
            checkpoint_storage=checkpoint_storage,
            analysis_particle_indices=analysis_particle_indices,
            position_interval=position_interval,
            velocity_interval=velocity_interval,
        )

    def write_sampler_states(
        self,
        sampler_states: list,
        iteration: int,
    ) -> None:
        """
        Write sampler states to XTC and lightweight checkpoint.

        Args:
            sampler_states: List of sampler states for all replicas
            iteration: Current iteration number
        """
        # 1. Write to XTC (only on checkpoint interval)
        if self._on_checkpoint_interval(iteration):
            # Try to infer timestep/interval again if not set (mcmc_moves may now be available)
            if self._xtc_timestep is None:
                self._xtc_timestep = self._infer_timestep()
            if self._xtc_interval is None:
                self._xtc_interval = self._infer_interval()
                # If we just got the correct interval, update existing XTCFile handles
                if self._xtc_interval is not None:
                    for xtc in self._xtc_handles.values():
                        xtc._interval = self._xtc_interval
            
            for replica_idx, state in enumerate(sampler_states):
                xtc = self._get_xtc_handle(replica_idx, state)
                # Convert box_vectors to tuple of Vec3 for XTCFile
                box_vectors = self._convert_box_vectors(state.box_vectors)
                xtc.writeModel(
                    state.positions,
                    periodicBoxVectors=box_vectors,
                )

        # 2. Write to lightweight NetCDF checkpoint (always overwrite index 0)
        self._write_lite_checkpoint(sampler_states, iteration)

    def read_sampler_states(
        self,
        iteration: int | None = None,
        analysis_particles_only: bool = False,
    ) -> list | None:
        """
        Read sampler states from lightweight checkpoint.

        Always returns the latest saved frame (index 0), iteration is ignored.

        Args:
            iteration: Ignored in XTC mode (for API compatibility)
            analysis_particles_only: Whether to read only analysis particles

        Returns:
            List of sampler states, or None if not available
        """
        if analysis_particles_only:
            # Fall back to parent for analysis particles
            return super().read_sampler_states(iteration, True)

        return self._read_lite_checkpoint()

    # -------------------------------------------------------------------------
    # Private methods
    # -------------------------------------------------------------------------

    def _write_lite_checkpoint(
        self,
        sampler_states: list,
        iteration: int,
    ) -> None:
        """
        Write lightweight checkpoint (only latest frame, fixed size).

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

        # Atomic write: sync, write current_iteration, sync
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

        # Fixed size of 1 (not unlimited) - key for lightweight checkpoint
        if 'iteration_lite' not in storage.dimensions:
            storage.createDimension('iteration_lite', 1)
        if 'replica' not in storage.dimensions:
            storage.createDimension('replica', n_replicas)
        if 'atom' not in storage.dimensions:
            storage.createDimension('atom', n_atoms)
        if 'spatial' not in storage.dimensions:
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

    def _read_lite_checkpoint(self) -> list | None:
        """
        Read sampler states from lightweight checkpoint.

        Always reads from index 0 (latest frame).

        Returns:
            List of sampler states
        """
        storage = self._storage_checkpoint

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
        Get or create XTCFile handle for a replica.

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

            # Infer interval (MD steps between checkpoints) from mcmc_moves
            if self._xtc_interval is None:
                self._xtc_interval = self._infer_interval()
            
            # If still None, use checkpoint_interval as fallback (will show iteration numbers)
            # This happens on first write before mcmc_moves is stored in NetCDF
            xtc_interval = self._xtc_interval if self._xtc_interval is not None else self._checkpoint_interval

            # Create topology from state
            topology = self._create_topology_from_state(state)

            self._xtc_handles[replica_idx] = XTCFile(
                path,
                topology,
                dt=self._xtc_timestep * unit.picosecond,
                interval=xtc_interval,
                append=append,
            )

        return self._xtc_handles[replica_idx]

    def _convert_box_vectors(self, box_vectors) -> tuple | None:
        """
        Convert box_vectors to tuple of Vec3 for XTCFile.
        
        Handles various input formats: Quantity with tuple, Quantity with ndarray,
        or raw tuple of Vec3.
        
        Args:
            box_vectors: Box vectors in various formats
            
        Returns:
            Tuple of 3 Vec3, or None
        """
        if box_vectors is None:
            return None
            
        # If it's a Quantity, extract the value
        if hasattr(box_vectors, 'value_in_unit'):
            box_vectors = box_vectors.value_in_unit(unit.nanometers)
        
        # If it's already a tuple of Vec3, return as-is
        if isinstance(box_vectors, tuple) and len(box_vectors) == 3:
            if hasattr(box_vectors[0], 'x'):  # Check if it's Vec3
                return box_vectors
        
        # If it's a numpy array or list, convert to Vec3 tuple
        if hasattr(box_vectors, '__len__') and len(box_vectors) == 3:
            vec_array = np.asarray(box_vectors)
            return (
                Vec3(vec_array[0][0], vec_array[0][1], vec_array[0][2]),
                Vec3(vec_array[1][0], vec_array[1][1], vec_array[1][2]),
                Vec3(vec_array[2][0], vec_array[2][1], vec_array[2][2]),
            )
        
        return None

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
            moves = self.read_mcmc_moves()
            if moves and len(moves) > 0:
                timestep = moves[0].timestep
                if unit.is_quantity(timestep):
                    return timestep.value_in_unit(unit.picosecond)
                return timestep
        except Exception:
            pass

        # Default for Martini
        return 0.02  # 20 fs = 0.02 ps

    def _infer_interval(self) -> int | None:
        """
        Infer MD step interval between checkpoints from mcmc_moves.

        In REMD, interval = n_steps_per_iteration * checkpoint_interval

        Returns:
            Number of MD steps between XTC writes, or None if cannot infer
        """
        try:
            moves = self.read_mcmc_moves()
            if moves and len(moves) > 0:
                n_steps = moves[0].n_steps
                # interval = steps per iteration * checkpoint interval (iterations)
                return n_steps * self._checkpoint_interval
        except Exception:
            pass

        # Return None to indicate failure - will retry later
        return None

    def _create_topology_from_state(
        self,
        state: states.SamplerState,
    ) -> OpenMMTopology:
        """
        Create OpenMM Topology from SamplerState.

        Creates a minimal topology with placeholder atoms.
        """
        n_atoms = state.n_particles
        top = OpenMMTopology()
        chain = top.addChain()
        residue = top.addResidue('UNK', chain)

        for i in range(n_atoms):
            element = Element.getBySymbol('C')
            top.addAtom(f'AT{i}', element, residue)

        if state.box_vectors is not None:
            top.setPeriodicBoxVectors(state.box_vectors)

        return top

    def close(self) -> None:
        """Close reporter and clean up XTC handles."""
        self._xtc_handles.clear()
        super().close()

    def read_checkpoint_iterations(self) -> np.ndarray:
        """Read checkpoint iterations (only one in lite mode)."""
        try:
            storage = self._storage_checkpoint
            current = int(storage.variables['current_iteration'][0])
            return np.array([current])
        except (KeyError, AttributeError):
            return np.array([0])

    def read_last_iteration(self, last_checkpoint: bool = True) -> int:
        """Read last saved iteration from current_iteration marker."""
        try:
            storage = self._storage_checkpoint
            return int(storage.variables['current_iteration'][0])
        except (KeyError, AttributeError):
            return 0

    def write_timestamp(self, iteration: int) -> None:
        """Write timestamp and print progress.
        
        Overrides parent method to add progress output.
        
        Args:
            iteration: Current iteration number
        """
        super().write_timestamp(iteration)
        
        # Initialize progress tracking on first call
        if not self._progress_initialized:
            self._progress_start_time = time.time()
            self._progress_start_iter = iteration
            self._progress_initialized = True
        
        # Print progress at specified intervals
        if iteration % self._progress_interval == 0:
            self._print_progress(iteration)
    
    def _print_progress(self, iteration: int) -> None:
        """Print progress information in StateDataReporter format.
        
        Args:
            iteration: Current iteration number
        """
        if self._progress_start_time is None or self._progress_start_iter is None:
            return
            
        elapsed = time.time() - self._progress_start_time
        it_done = iteration - self._progress_start_iter
        
        if it_done <= 0:
            return
        
        # Calculate speed (iterations per second)
        speed_it_s = it_done / elapsed
        
        # Calculate progress and ETA
        if self._total_iterations is not None and self._total_iterations > 0:
            progress = 100.0 * iteration / self._total_iterations
            remaining_iters = self._total_iterations - iteration
            eta_seconds = remaining_iters / speed_it_s if speed_it_s > 0 else float('inf')
            eta_str = self._format_time(eta_seconds)
        else:
            progress = 0.0
            eta_str = "--"
        
        # Calculate ns/day (total across all replicas)
        # ns/day = it/s × n_replicas × exc_freq × dt(ps) × 86400 / 1000
        if self._exc_freq is not None and self._dt is not None:
            speed_ns_day = speed_it_s * self._n_replicas * self._exc_freq * self._dt * 86400 / 1000
        else:
            speed_ns_day = 0.0
        
        # Print header on first call
        if not self._progress_header_printed:
            print('#"Progress (%)"\t"Iteration"\t"Speed (it/s)"\t"Speed (ns/day)"\t"Time Remaining"')
            self._progress_header_printed = True
        
        # Print data row (tab-separated)
        print(f"{progress:.1f}%\t{iteration}\t{speed_it_s:.2f}\t{speed_ns_day:.2f}\t{eta_str}")
    
    @staticmethod
    def _format_time(seconds: float) -> str:
        """Format time in human-readable format.
        
        Args:
            seconds: Time in seconds
            
        Returns:
            Formatted time string (e.g., "2:30:45" or "45:30")
        """
        if seconds == float('inf') or seconds < 0:
            return "--"
        
        seconds = int(seconds)
        days = seconds // 86400
        seconds -= days * 86400
        hours = seconds // 3600
        seconds -= hours * 3600
        minutes = seconds // 60
        seconds -= minutes * 60
        
        if days > 0:
            return f"{days}:{hours:02d}:{minutes:02d}:{seconds:02d}"
        elif hours > 0:
            return f"{hours}:{minutes:02d}:{seconds:02d}"
        else:
            return f"{minutes}:{seconds:02d}"
