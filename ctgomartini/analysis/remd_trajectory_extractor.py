#!/usr/bin/env python3
"""
REMD trajectory extraction tool.

Extracts replica and thermodynamic state trajectories from NetCDF files
with support for parallel processing and memory-efficient streaming.

Author: Song Yang
Date: 2025
"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import warnings
from pathlib import Path
from typing import Optional, Literal

import numpy as np

warnings.filterwarnings("ignore")

# Lazy imports
_MultiStateReporter = None

def _import_openmmtools():
    """Lazy import openmmtools."""
    global _MultiStateReporter
    if _MultiStateReporter is None:
        from openmmtools.multistate import MultiStateReporter
        _MultiStateReporter = MultiStateReporter
    return _MultiStateReporter


class REMDTrajectoryExtractor:
    """
    Extract trajectories from REMD NetCDF files.
    
    Supports:
    - Replica trajectory extraction (continuous)
    - State trajectory extraction (discontinuous, constant thermodynamic state)
    - Single frame extraction
    - Parallel processing
    - Memory-efficient streaming write
    
    Examples:
        >>> extractor = REMDTrajectoryExtractor("output.nc", "topology.pdb")
        >>> 
        >>> # Extract all replica trajectories
        >>> extractor.save_replica_trajectories("output_dir/", format="xtc")
        >>> 
        >>> # Extract specific state
        >>> extractor.save_state_trajectory(0, "state_0.xtc")
        >>> 
        >>> # Extract single frame
        >>> extractor.save_frame(frame_idx=100, replica_idx=0, output="frame_100.pdb")
    """
    
    def __init__(
        self,
        netcdf_file: str,
        topology_file: str,
        checkpoint_file: Optional[str] = None,
    ):
        """
        Initialize extractor.
        
        Args:
            netcdf_file: Path to NetCDF trajectory file
            topology_file: Path to topology file (pdb/gro)
            checkpoint_file: Path to checkpoint file (optional)
        """
        self.netcdf_file = netcdf_file
        self.checkpoint_file = checkpoint_file or "output_checkpoint.nc"
        self.topology_file = topology_file
        
        # Load topology using MDAnalysis (lightweight)
        self._load_topology()
        
        # Get metadata from NetCDF
        self._load_metadata()
    
    def _load_topology(self):
        """Load topology using MDAnalysis."""
        import MDAnalysis as mda
        self.u = mda.Universe(self.topology_file)
        self.n_atoms = len(self.u.atoms)
        
        # Store atom names and residues for writing
        self.atom_names = self.u.atoms.names
        self.resnames = self.u.atoms.resnames
        self.resids = self.u.atoms.resids
    
    def _load_metadata(self):
        """Load metadata from NetCDF without loading positions."""
        MultiStateReporter = _import_openmmtools()
        
        reporter = MultiStateReporter(
            self.netcdf_file,
            open_mode='r',
            checkpoint_storage=self.checkpoint_file
        )
        
        try:
            storage = reporter._storage_checkpoint
            positions_var = storage.variables['positions']
            
            self.n_frames = positions_var.shape[0]
            self.n_replicas = positions_var.shape[1]
            self.n_atoms_nc = positions_var.shape[2]
            
            if self.n_atoms_nc != self.n_atoms:
                raise ValueError(
                    f"Atom count mismatch: NetCDF has {self.n_atoms_nc}, "
                    f"topology has {self.n_atoms}"
                )
            
            # Get thermodynamic states information
            try:
                states_var = reporter._storage_analysis.variables['states']
                self.state_trajectories = states_var[:]
            except Exception:
                # If no state information, assume replica = state
                self.state_trajectories = np.arange(self.n_replicas).reshape(1, -1)
                self.state_trajectories = np.repeat(
                    self.state_trajectories, self.n_frames, axis=0
                )
            
            # Get time information
            try:
                from openmm import unit
                mcmove = reporter.read_mcmc_moves()[0]
                self.timestep = mcmove.timestep.value_in_unit(unit.picosecond)
                self.n_steps_per_frame = mcmove.n_steps
                self.dt = self.timestep * self.n_steps_per_frame
            except Exception:
                self.dt = 5.0  # Default 5 ps
            
        finally:
            reporter.close()
    
    def get_frame(
        self,
        frame_idx: int,
        replica_idx: Optional[int] = None,
        state_idx: Optional[int] = None,
    ) -> np.ndarray:
        """
        Extract single frame positions.
        
        Args:
            frame_idx: Frame index
            replica_idx: Replica index (mutually exclusive with state_idx)
            state_idx: Thermodynamic state index (mutually exclusive with replica_idx)
            
        Returns:
            Positions array (n_atoms x 3) in Angstroms
        """
        if replica_idx is not None and state_idx is not None:
            raise ValueError("Cannot specify both replica_idx and state_idx")
        if replica_idx is None and state_idx is None:
            raise ValueError("Must specify either replica_idx or state_idx")
        
        MultiStateReporter = _import_openmmtools()
        reporter = MultiStateReporter(
            self.netcdf_file,
            open_mode='r',
            checkpoint_storage=self.checkpoint_file
        )
        
        try:
            storage = reporter._storage_checkpoint
            
            if replica_idx is not None:
                # Direct replica extraction
                positions = storage.variables['positions'][
                    frame_idx, replica_idx, :, :
                ]
            else:
                # State extraction - find which replica is in this state
                replica_indices = self.state_trajectories[frame_idx]
                matching_replicas = np.where(replica_indices == state_idx)[0]
                
                if len(matching_replicas) == 0:
                    raise ValueError(
                        f"No replica found in state {state_idx} at frame {frame_idx}"
                    )
                
                # Take first matching replica
                replica_idx = matching_replicas[0]
                positions = storage.variables['positions'][
                    frame_idx, replica_idx, :, :
                ]
            
            return positions * 10.0  # nm to Angstrom
            
        finally:
            reporter.close()
    
    def save_frame(
        self,
        frame_idx: int,
        output: str,
        replica_idx: Optional[int] = None,
        state_idx: Optional[int] = None,
    ):
        """
        Save single frame to file.
        
        Args:
            frame_idx: Frame index
            output: Output file path
            replica_idx: Replica index
            state_idx: State index
        """
        positions = self.get_frame(frame_idx, replica_idx, state_idx)
        time_ps = frame_idx * self.dt
        
        # Write based on extension
        ext = Path(output).suffix.lower()
        
        if ext == '.pdb':
            self._write_pdb(positions, output, time_ps)
        elif ext == '.gro':
            self._write_gro(positions, output, time_ps)
        else:
            raise ValueError(f"Unsupported format: {ext}")
        
        print(f"Saved frame {frame_idx} to {output}")
    
    def _write_pdb(self, positions: np.ndarray, output: str, time_ps: float):
        """Write PDB format."""
        with open(output, 'w') as f:
            f.write(f"REMARK   Frame at {time_ps:.3f} ps\n")
            f.write(f"CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n")
            
            for i, (pos, name, resname, resid) in enumerate(
                zip(positions, self.atom_names, self.resnames, self.resids)
            ):
                f.write(
                    f"ATOM  {i+1:5d} {name:4s} {resname:3s} A{resid:4d}    "
                    f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}"
                    f"  1.00  0.00           {name[0]:1s}\n"
                )
            
            f.write("END\n")
    
    def _write_gro(self, positions: np.ndarray, output: str, time_ps: float):
        """Write GRO format."""
        # Convert to nm
        positions_nm = positions / 10.0
        
        with open(output, 'w') as f:
            f.write(f"Frame at t={time_ps:.3f} ps\n")
            f.write(f"{len(positions):5d}\n")
            
            for i, (pos, name, resname, resid) in enumerate(
                zip(positions_nm, self.atom_names, self.resnames, self.resids)
            ):
                f.write(
                    f"{resid:5d}{resname:5s}{name:5s}{i+1:5d}"
                    f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n"
                )
            
            f.write("  10.00000  10.00000  10.00000\n")
    
    def _extract_replica_worker(
        self,
        replica_idx: int,
        output: str,
        format: str,
        frame_begin: int,
        frame_end: int,
        frame_stride: int,
    ):
        """Worker function for extracting single replica trajectory."""
        MultiStateReporter = _import_openmmtools()
        reporter = MultiStateReporter(
            self.netcdf_file,
            open_mode='r',
            checkpoint_storage=self.checkpoint_file
        )
        
        try:
            storage = reporter._storage_checkpoint
            frame_indices = range(frame_begin, frame_end, frame_stride)
            n_frames = len(frame_indices)
            
            if format in ['dcd', 'xtc']:
                # Use MDTraj for binary formats
                from mdtraj import Trajectory, Topology as MDTrajTopology
                
                positions = np.zeros((n_frames, self.n_atoms, 3), dtype=np.float32)
                
                for i, frame_idx in enumerate(frame_indices):
                    pos = storage.variables['positions'][frame_idx, replica_idx, :, :]
                    positions[i] = pos.astype(np.float32)
                
                times = np.array([f * self.dt for f in frame_indices])
                
                traj = Trajectory(
                    positions,
                    MDTrajTopology.from_openmm(self._get_openmm_topology()),
                    time=times,
                )
                
                if format == 'dcd':
                    Trajectory.save_dcd(traj, output)
                else:
                    Trajectory.save_xtc(traj, output)
                    
            elif format == 'multi-pdb':
                # Stream write multi-model PDB
                with open(output, 'w') as f:
                    for i, frame_idx in enumerate(frame_indices):
                        pos = storage.variables['positions'][frame_idx, replica_idx, :, :]
                        pos_angstrom = pos * 10.0
                        time_ps = frame_idx * self.dt
                        
                        f.write(f"MODEL     {i+1:4d}\n")
                        f.write(f"REMARK   Frame {frame_idx} at {time_ps:.3f} ps\n")
                        
                        for j, (p, name, resname, resid) in enumerate(
                            zip(pos_angstrom, self.atom_names, self.resnames, self.resids)
                        ):
                            f.write(
                                f"ATOM  {j+1:5d} {name:4s} {resname:3s} A{resid:4d}    "
                                f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}"
                                f"  1.00  0.00           {name[0]:1s}\n"
                            )
                        
                        f.write("ENDMDL\n")
            
            print(f"  Replica {replica_idx}: {n_frames} frames -> {output}")
            
        finally:
            reporter.close()
    
    def _get_openmm_topology(self):
        """Convert MDAnalysis topology to OpenMM topology."""
        # Lazy import
        from openmm.app import Topology as OpenMMTopology, Element
        
        omm_top = OpenMMTopology()
        chain = omm_top.addChain()
        
        current_resid = None
        residue = None
        
        for atom in self.u.atoms:
            if atom.resid != current_resid:
                residue = omm_top.addResidue(atom.resname, chain)
                current_resid = atom.resid
            
            # Guess element from atom name
            element = Element.getBySymbol(atom.name[0]) if atom.name[0].isalpha() else Element.getBySymbol('C')
            omm_top.addAtom(atom.name, element, residue)
        
        return omm_top
    
    def save_replica_trajectories(
        self,
        output_dir: str,
        format: Literal['dcd', 'xtc', 'multi-pdb'] = 'xtc',
        frame_begin: int = 0,
        frame_end: Optional[int] = None,
        frame_stride: int = 1,
        num_workers: Optional[int] = None,
        replicas: Optional[list[int]] = None,
    ):
        """
        Extract and save replica trajectories.
        
        Args:
            output_dir: Output directory
            format: Output format ('dcd', 'xtc', 'multi-pdb')
            frame_begin: Start frame
            frame_end: End frame (None = all)
            frame_stride: Frame stride
            num_workers: Number of parallel workers
            replicas: List of replica indices to extract (None = all)
        """
        os.makedirs(output_dir, exist_ok=True)
        
        frame_end = frame_end or self.n_frames
        
        if replicas is None:
            replicas = list(range(self.n_replicas))
        
        print(f"Extracting {len(replicas)} replica trajectories...")
        print(f"  Frames: {frame_begin}:{frame_end}:{frame_stride}")
        print(f"  Format: {format}")
        
        # Prepare tasks
        tasks = []
        for replica_idx in replicas:
            output = os.path.join(output_dir, f"replica_{replica_idx}.{format}")
            tasks.append((
                replica_idx, output, format,
                frame_begin, frame_end, frame_stride
            ))
        
        # Parallel execution
        num_workers = num_workers or min(multiprocessing.cpu_count(), len(replicas))
        
        if num_workers == 1:
            for task in tasks:
                self._extract_replica_worker(*task)
        else:
            with multiprocessing.Pool(num_workers) as pool:
                pool.starmap(self._extract_replica_worker, tasks)
        
        print(f"Done. Saved to {output_dir}/")
    
    def save_state_trajectories(
        self,
        output_dir: str,
        format: Literal['dcd', 'xtc', 'multi-pdb'] = 'xtc',
        frame_begin: int = 0,
        frame_end: Optional[int] = None,
        frame_stride: int = 1,
        num_workers: Optional[int] = None,
        states: Optional[list[int]] = None,
    ):
        """
        Extract and save thermodynamic state trajectories.
        
        Note: State trajectories are discontinuous (replicas exchange states).
        
        Args:
            output_dir: Output directory
            format: Output format
            frame_begin: Start frame
            frame_end: End frame
            frame_stride: Frame stride
            num_workers: Number of parallel workers
            states: List of state indices to extract (None = all)
        """
        os.makedirs(output_dir, exist_ok=True)
        
        frame_end = frame_end or self.n_frames
        
        n_states = len(np.unique(self.state_trajectories[0]))
        if states is None:
            states = list(range(n_states))
        
        print(f"Extracting {len(states)} state trajectories...")
        print(f"  Frames: {frame_begin}:{frame_end}:{frame_stride}")
        print(f"  Format: {format}")
        
        # State extraction is more complex - each frame may have different replicas
        # For now, use sequential with streaming write
        for state_idx in states:
            output = os.path.join(output_dir, f"state_{state_idx}.{format}")
            self._extract_state_streaming(
                state_idx, output, format,
                frame_begin, frame_end, frame_stride
            )
        
        print(f"Done. Saved to {output_dir}/")
    
    def _extract_state_streaming(
        self,
        state_idx: int,
        output: str,
        format: str,
        frame_begin: int,
        frame_end: int,
        frame_stride: int,
    ):
        """Stream extract state trajectory (memory efficient)."""
        MultiStateReporter = _import_openmmtools()
        reporter = MultiStateReporter(
            self.netcdf_file,
            open_mode='r',
            checkpoint_storage=self.checkpoint_file
        )
        
        try:
            storage = reporter._storage_checkpoint
            frame_indices = list(range(frame_begin, frame_end, frame_stride))
            
            if format == 'multi-pdb':
                # Stream write PDB
                with open(output, 'w') as f:
                    model_idx = 1
                    for frame_idx in frame_indices:
                        # Find replica in this state
                        replica_indices = self.state_trajectories[frame_idx]
                        matching = np.where(replica_indices == state_idx)[0]
                        
                        if len(matching) == 0:
                            continue
                        
                        replica_idx = matching[0]
                        pos = storage.variables['positions'][frame_idx, replica_idx, :, :]
                        pos_angstrom = pos * 10.0
                        time_ps = frame_idx * self.dt
                        
                        f.write(f"MODEL     {model_idx:4d}\n")
                        f.write(f"REMARK   Frame {frame_idx} at {time_ps:.3f} ps\n")
                        
                        for j, (p, name, resname, resid) in enumerate(
                            zip(pos_angstrom, self.atom_names, self.resnames, self.resids)
                        ):
                            f.write(
                                f"ATOM  {j+1:5d} {name:4s} {resname:3s} A{resid:4d}    "
                                f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}"
                                f"  1.00  0.00           {name[0]:1s}\n"
                            )
                        
                        f.write("ENDMDL\n")
                        model_idx += 1
                
                print(f"  State {state_idx}: {model_idx-1} frames -> {output}")
                
            else:
                # For binary formats, collect positions then write
                positions_list = []
                times = []
                
                for frame_idx in frame_indices:
                    replica_indices = self.state_trajectories[frame_idx]
                    matching = np.where(replica_indices == state_idx)[0]
                    
                    if len(matching) == 0:
                        continue
                    
                    replica_idx = matching[0]
                    pos = storage.variables['positions'][frame_idx, replica_idx, :, :]
                    positions_list.append(pos.astype(np.float32))
                    times.append(frame_idx * self.dt)
                
                if len(positions_list) > 0:
                    from mdtraj import Trajectory, Topology as MDTrajTopology
                    
                    positions = np.array(positions_list)
                    traj = Trajectory(
                        positions,
                        MDTrajTopology.from_openmm(self._get_openmm_topology()),
                        time=np.array(times),
                    )
                    
                    if format == 'dcd':
                        Trajectory.save_dcd(traj, output)
                    else:
                        Trajectory.save_xtc(traj, output)
                
                print(f"  State {state_idx}: {len(positions_list)} frames -> {output}")
                
        finally:
            reporter.close()


def main():
    """Command-line interface for trajectory extraction."""
    parser = argparse.ArgumentParser(
        description="Extract replica and state trajectories from REMD simulations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    parser.add_argument("-nc", "--netcdf", required=True,
                        help="NetCDF trajectory file")
    parser.add_argument("-p", "--topology", required=True,
                        help="Topology file (pdb/gro)")
    parser.add_argument("-c", "--checkpoint", default="output_checkpoint.nc",
                        help="Checkpoint NetCDF file")
    
    # Extraction mode
    parser.add_argument("--mode", choices=["replica", "state", "frame"], required=True,
                        help="Extraction mode")
    
    # Output options
    parser.add_argument("-o", "--output", default="./extracted",
                        help="Output directory or file (for single frame)")
    parser.add_argument("--format", choices=["dcd", "xtc", "pdb"], default="xtc",
                        help="Output format")
    
    # Frame selection
    parser.add_argument("-b", "--begin", type=int, default=0,
                        help="Start frame")
    parser.add_argument("-e", "--end", type=int, default=None,
                        help="End frame")
    parser.add_argument("-s", "--stride", type=int, default=1,
                        help="Frame stride")
    
    # For single frame mode
    parser.add_argument("--frame", type=int, default=None,
                        help="Frame index (for --mode frame)")
    parser.add_argument("--replica", type=int, default=None,
                        help="Replica index")
    parser.add_argument("--state", type=int, default=None,
                        help="State index")
    
    # Parallel options
    parser.add_argument("-np", "--num-workers", type=int, default=None,
                        help="Number of parallel workers")
    parser.add_argument("--replicas", type=str, default=None,
                        help="Comma-separated replica indices (e.g., '0,1,2')")
    parser.add_argument("--states", type=str, default=None,
                        help="Comma-separated state indices")
    
    args = parser.parse_args()
    
    # Initialize extractor
    extractor = REMDTrajectoryExtractor(
        args.netcdf,
        args.topology,
        args.checkpoint
    )
    
    # Execute based on mode
    if args.mode == "frame":
        if args.frame is None:
            print("Error: --frame required for frame mode")
            return 1
        
        output = args.output
        if os.path.isdir(output):
            output = os.path.join(output, f"frame_{args.frame}.pdb")
        
        extractor.save_frame(
            frame_idx=args.frame,
            output=output,
            replica_idx=args.replica,
            state_idx=args.state,
        )
    
    elif args.mode == "replica":
        replicas = None
        if args.replicas:
            replicas = [int(x) for x in args.replicas.split(",")]
        
        extractor.save_replica_trajectories(
            output_dir=args.output,
            format=args.format,
            frame_begin=args.begin,
            frame_end=args.end,
            frame_stride=args.stride,
            num_workers=args.num_workers,
            replicas=replicas,
        )
    
    elif args.mode == "state":
        states = None
        if args.states:
            states = [int(x) for x in args.states.split(",")]
        
        extractor.save_state_trajectories(
            output_dir=args.output,
            format=args.format,
            frame_begin=args.begin,
            frame_end=args.end,
            frame_stride=args.stride,
            num_workers=args.num_workers,
            states=states,
        )
    
    return 0


if __name__ == "__main__":
    exit(main())
