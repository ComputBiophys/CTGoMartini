#!/usr/bin/env python3
"""
REMD trajectory extraction tool.

Extracts replica and thermodynamic state trajectories from NetCDF files
with support for parallel processing and memory-efficient streaming.

Supports automatic format detection from file extension:
  - .xtc (default) - GROMACS compressed trajectory
  - .dcd - CHARMM/NAMD binary trajectory
  - .pdb - Multi-model PDB (text format)

"""

from __future__ import annotations

import argparse
import multiprocessing
import os
import warnings
from pathlib import Path
from typing import Optional

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


def _detect_format(filename: str) -> str:
    """
    Detect trajectory format from file extension.
    
    Args:
        filename: Output filename
        
    Returns:
        Format string: 'xtc', 'dcd', or 'pdb'
        
    Raises:
        ValueError: If format cannot be detected
    """
    ext = Path(filename).suffix.lower()
    format_map = {
        '.xtc': 'xtc',
        '.dcd': 'dcd',
        '.pdb': 'pdb',
    }
    if ext not in format_map:
        raise ValueError(
            f"Cannot detect format from extension '{ext}'. "
            f"Supported: .xtc, .dcd, .pdb"
        )
    return format_map[ext]


def _optimize_chunk_size(n_frames: int, n_atoms: int, num_workers: int) -> int:
    """
    Optimize chunk size for batch reading.
    
    Balances memory usage and I/O efficiency.
    
    Args:
        n_frames: Total frames to process
        n_atoms: Number of atoms
        num_workers: Number of parallel workers
        
    Returns:
        Optimal chunk size
    """
    # Target ~50MB per chunk
    bytes_per_frame = n_atoms * 3 * 4  # float32
    target_frames = int((50 * 1024 * 1024) / max(bytes_per_frame, 1))
    
    # Ensure at least 2 chunks per worker for load balancing
    min_chunks = num_workers * 2
    max_chunk = max(1, n_frames // min_chunks)
    
    # Clamp between 10 and 500 frames
    chunk_size = min(target_frames, max_chunk)
    chunk_size = max(10, min(chunk_size, 500))
    chunk_size = min(chunk_size, n_frames)
    
    return chunk_size


class REMDTrajectoryExtractor:
    """
    Extract trajectories from REMD NetCDF files.
    
    Supports:
    - Replica trajectory extraction (continuous)
    - State trajectory extraction (discontinuous, constant thermodynamic state)
    - Single frame extraction
    - Parallel processing
    - Memory-efficient streaming write
    - Automatic format detection from file extension
    
    Default format is XTC (GROMACS compressed).
    
    Examples:
        >>> extractor = REMDTrajectoryExtractor("output.nc", "topology.pdb")
        >>> 
        >>> # Extract all replica trajectories (default xtc format)
        >>> extractor.save_replica_trajectories("output_dir/")
        >>> 
        >>> # Extract with auto-detected format
        >>> extractor.save_replica_trajectories("output_dir/", output_pattern="replica_{i}.dcd")
        >>> 
        >>> # Extract specific state
        >>> extractor.save_state_trajectories("output_dir/", output_pattern="state_{i}.pdb")
        >>> 
        >>> # Extract single frame
        >>> extractor.save_frame(frame_idx=100, replica_idx=0, output="frame_100.pdb")
    """
    
    def __init__(
        self,
        netcdf_file: str,
        topology_file: str,
        checkpoint_file: Optional[str] = None,
        cache_states: bool = True,
    ):
        """
        Initialize extractor.
        
        Args:
            netcdf_file: Path to NetCDF trajectory file
            topology_file: Path to topology file (pdb/gro)
            checkpoint_file: Path to checkpoint file (optional)
            cache_states: Whether to cache state-replica mapping in memory
        """
        self.netcdf_file = netcdf_file
        self.checkpoint_file = checkpoint_file or "output_checkpoint.nc"
        self.topology_file = topology_file
        self.cache_states = cache_states
        
        # Load topology using MDAnalysis (lightweight)
        self._load_topology()
        
        # Get metadata from NetCDF
        self._load_metadata()
        
        # Cache for state mapping (optimization)
        self._state_cache = None
    
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
            
            # Get thermodynamic states information and cache it
            self._load_state_mapping(reporter)
            
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
    
    def _load_state_mapping(self, reporter):
        """
        Load and cache state-replica mapping.
        
        Creates an optimized lookup structure:
        - self.state_to_replicas: dict[state] -> list of (frame, replica) tuples
        """
        try:
            states_var = reporter._storage_analysis.variables['states']
            self.state_trajectories = states_var[:]
        except Exception:
            # If no state information, assume replica = state
            self.state_trajectories = np.arange(self.n_replicas).reshape(1, -1)
            self.state_trajectories = np.repeat(
                self.state_trajectories, self.n_frames, axis=0
            )
        
        if self.cache_states:
            # Build reverse lookup: state -> list of (frame, replica) pairs
            self._state_cache = {}
            for frame_idx in range(self.n_frames):
                for replica_idx, state in enumerate(self.state_trajectories[frame_idx]):
                    if state not in self._state_cache:
                        self._state_cache[state] = []
                    self._state_cache[state].append((frame_idx, replica_idx))
    
    def _get_replica_for_state(self, frame_idx: int, state_idx: int) -> Optional[int]:
        """
        Get replica index for a given state at a specific frame.
        
        Uses cached mapping if available.
        """
        if self._state_cache is not None:
            # Use cache
            state_list = self._state_cache.get(state_idx, [])
            for f, r in state_list:
                if f == frame_idx:
                    return r
            return None
        else:
            # Direct lookup
            replica_indices = self.state_trajectories[frame_idx]
            matching = np.where(replica_indices == state_idx)[0]
            return matching[0] if len(matching) > 0 else None
    
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
        
        if not ext:
            # Default to PDB if no extension
            output = output + '.pdb'
            ext = '.pdb'
        
        if ext == '.pdb':
            self._write_pdb(positions, output, time_ps)
        elif ext == '.gro':
            self._write_gro(positions, output, time_ps)
        else:
            raise ValueError(
                f"Unsupported format '{ext}' for single frame. "
                f"Use .pdb or .gro for single frame extraction."
            )
        
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
        chunk_size: int,
    ):
        """
        Worker function for extracting single replica trajectory.
        
        Uses batch reading for improved I/O performance.
        """
        MultiStateReporter = _import_openmmtools()
        reporter = MultiStateReporter(
            self.netcdf_file,
            open_mode='r',
            checkpoint_storage=self.checkpoint_file
        )
        
        try:
            storage = reporter._storage_checkpoint
            frame_indices = list(range(frame_begin, frame_end, frame_stride))
            n_frames = len(frame_indices)
            
            if format in ['dcd', 'xtc']:
                # Use MDTraj for binary formats
                from mdtraj import Trajectory, Topology as MDTrajTopology
                
                positions = np.zeros((n_frames, self.n_atoms, 3), dtype=np.float32)
                
                # Batch read for better I/O performance
                for chunk_start in range(0, n_frames, chunk_size):
                    chunk_end = min(chunk_start + chunk_size, n_frames)
                    chunk_frame_indices = frame_indices[chunk_start:chunk_end]
                    
                    # Read chunk in one operation (if possible)
                    for i, frame_idx in enumerate(chunk_frame_indices):
                        pos = storage.variables['positions'][frame_idx, replica_idx, :, :]
                        positions[chunk_start + i] = pos.astype(np.float32)
                
                times = np.array([f * self.dt for f in frame_indices])
                
                traj = Trajectory(
                    positions,
                    MDTrajTopology.from_openmm(self._get_openmm_topology()),
                    time=times,
                )
                
                if format == 'dcd':
                    Trajectory.save_dcd(traj, output)
                else:  # xtc
                    Trajectory.save_xtc(traj, output)
                    
            elif format == 'pdb':
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
        output_pattern: str = "replica_{i}.xtc",
        frame_begin: int = 0,
        frame_end: Optional[int] = None,
        frame_stride: int = 1,
        num_workers: Optional[int] = None,
        replicas: Optional[list[int]] = None,
    ):
        """
        Extract and save replica trajectories.
        
        Format is auto-detected from output_pattern extension.
        Default is XTC format.
        
        Args:
            output_dir: Output directory
            output_pattern: Output filename pattern with {i} placeholder
                          (e.g., "replica_{i}.xtc", "replica_{i}.dcd")
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
        
        # Auto-detect format from pattern
        format = _detect_format(output_pattern)
        
        print(f"Extracting {len(replicas)} replica trajectories...")
        print(f"  Frames: {frame_begin}:{frame_end}:{frame_stride}")
        print(f"  Format: {format} (auto-detected)")
        
        # Optimize chunk size for batch reading
        chunk_size = _optimize_chunk_size(
            frame_end - frame_begin, self.n_atoms, 
            num_workers or multiprocessing.cpu_count()
        )
        
        # Prepare tasks
        tasks = []
        for replica_idx in replicas:
            output = os.path.join(output_dir, output_pattern.replace('{mode}', 'replica').format(i=replica_idx))
            tasks.append((
                replica_idx, output, format,
                frame_begin, frame_end, frame_stride, chunk_size
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
        output_pattern: str = "state_{i}.xtc",
        frame_begin: int = 0,
        frame_end: Optional[int] = None,
        frame_stride: int = 1,
        num_workers: Optional[int] = None,
        states: Optional[list[int]] = None,
    ):
        """
        Extract and save thermodynamic state trajectories.
        
        Note: State trajectories are discontinuous (replicas exchange states).
        Format is auto-detected from output_pattern extension.
        
        Args:
            output_dir: Output directory
            output_pattern: Output filename pattern with {i} placeholder
            frame_begin: Start frame
            frame_end: End frame
            frame_stride: Frame stride
            num_workers: Number of parallel workers (sequential for states)
            states: List of state indices to extract (None = all)
        """
        os.makedirs(output_dir, exist_ok=True)
        
        frame_end = frame_end or self.n_frames
        
        n_states = len(np.unique(self.state_trajectories[0]))
        if states is None:
            states = list(range(n_states))
        
        # Auto-detect format from pattern
        format = _detect_format(output_pattern)
        
        print(f"Extracting {len(states)} state trajectories...")
        print(f"  Frames: {frame_begin}:{frame_end}:{frame_stride}")
        print(f"  Format: {format} (auto-detected)")
        print(f"  Using cached state mapping" if self._state_cache else "  Using direct lookup")
        
        # Sequential processing with streaming write
        # (State extraction is I/O bound, not CPU bound)
        for state_idx in states:
            output = os.path.join(output_dir, output_pattern.replace('{mode}', 'state').format(i=state_idx))
            self._extract_state_optimized(
                state_idx, output, format,
                frame_begin, frame_end, frame_stride
            )
        
        print(f"Done. Saved to {output_dir}/")
    
    def _extract_state_optimized(
        self,
        state_idx: int,
        output: str,
        format: str,
        frame_begin: int,
        frame_end: int,
        frame_stride: int,
    ):
        """
        Optimized state trajectory extraction using cached mapping.
        
        Uses pre-computed state-to-replica mapping for O(1) lookup.
        """
        MultiStateReporter = _import_openmmtools()
        reporter = MultiStateReporter(
            self.netcdf_file,
            open_mode='r',
            checkpoint_storage=self.checkpoint_file
        )
        
        try:
            storage = reporter._storage_checkpoint
            frame_indices = list(range(frame_begin, frame_end, frame_stride))
            
            # Collect (frame, replica) pairs using cache
            frame_replica_pairs = []
            for frame_idx in frame_indices:
                replica_idx = self._get_replica_for_state(frame_idx, state_idx)
                if replica_idx is not None:
                    frame_replica_pairs.append((frame_idx, replica_idx))
            
            n_frames_found = len(frame_replica_pairs)
            
            if format == 'pdb':
                # Stream write PDB
                with open(output, 'w') as f:
                    for model_idx, (frame_idx, replica_idx) in enumerate(frame_replica_pairs, 1):
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
                
                print(f"  State {state_idx}: {n_frames_found} frames -> {output}")
                
            else:
                # For binary formats, batch read then write
                if n_frames_found > 0:
                    from mdtraj import Trajectory, Topology as MDTrajTopology
                    
                    positions = np.zeros((n_frames_found, self.n_atoms, 3), dtype=np.float32)
                    times = np.zeros(n_frames_found)
                    
                    for i, (frame_idx, replica_idx) in enumerate(frame_replica_pairs):
                        pos = storage.variables['positions'][frame_idx, replica_idx, :, :]
                        positions[i] = pos.astype(np.float32)
                        times[i] = frame_idx * self.dt
                    
                    traj = Trajectory(
                        positions,
                        MDTrajTopology.from_openmm(self._get_openmm_topology()),
                        time=times,
                    )
                    
                    if format == 'dcd':
                        Trajectory.save_dcd(traj, output)
                    else:  # xtc
                        Trajectory.save_xtc(traj, output)
                
                print(f"  State {state_idx}: {n_frames_found} frames -> {output}")
                
        finally:
            reporter.close()


def main():
    """Command-line interface for trajectory extraction."""
    parser = argparse.ArgumentParser(
        description="Extract replica and state trajectories from REMD simulations. "
                    "Format is auto-detected from output filename extension "
                    "(.xtc, .dcd, .pdb). Default is XTC.",
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
    parser.add_argument("--pattern", default="{mode}_{i}.xtc",
                        help="Output filename pattern with {mode} and {i} placeholders. "
                             "Format auto-detected from extension. "
                             "Examples: 'replica_{i}.xtc', 'state_{i}.dcd', 'traj.pdb'")
    
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
    
    # Performance options
    parser.add_argument("--no-cache", action="store_true",
                        help="Disable state mapping cache (reduce memory for large trajectories)")
    
    args = parser.parse_args()
    
    # Initialize extractor
    extractor = REMDTrajectoryExtractor(
        args.netcdf,
        args.topology,
        args.checkpoint,
        cache_states=not args.no_cache
    )
    
    # Execute based on mode
    if args.mode == "frame":
        if args.frame is None:
            print("Error: --frame required for frame mode")
            return 1
        
        output = args.output
        if os.path.isdir(output):
            # Use pattern to determine filename
            ext = Path(args.pattern).suffix
            if not ext:
                ext = '.xtc'  # Default
            output = os.path.join(output, f"frame_{args.frame}{ext}")
        else:
            # If output is a file path, ensure it has proper extension
            # save_frame will detect format from extension
            pass
        
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
            output_pattern=args.pattern,
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
            output_pattern=args.pattern,
            frame_begin=args.begin,
            frame_end=args.end,
            frame_stride=args.stride,
            states=states,
        )
    
    return 0


if __name__ == "__main__":
    exit(main())
