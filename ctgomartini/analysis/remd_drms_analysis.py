#!/usr/bin/env python3
"""
REMD dRMS trajectory analysis and plotting.

Analyzes distance root-mean-square (dRMS) trajectories from replica exchange
molecular dynamics simulations and generates publication-quality plots.
"""

from __future__ import annotations

import argparse
import multiprocessing
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Tuple

import h5py
import numpy as np

warnings.filterwarnings("ignore")


def _get_metadata(chk_file: str) -> Tuple[int, int, int]:
    """Read basic dimensions from checkpoint file via h5py.

    Args:
        chk_file: Path to the checkpoint NetCDF file.

    Returns:
        Tuple of (n_frames, n_replicas, n_atoms).
    """
    with h5py.File(chk_file, "r") as f:
        shape = f["positions"].shape
    return shape[0], shape[1], shape[2]


def _get_timestep(nc_file: str, chk_file: str) -> float:
    """Read the frame timestep from REMD output files.

    Uses openmmtools.MultiStateReporter to read MCMC move parameters,
    then computes the timestep per frame. Falls back to 5.0 ps on error.

    Args:
        nc_file: Path to the main output.nc file.
        chk_file: Path to the checkpoint file.

    Returns:
        Timestep in picoseconds per frame.
    """
    try:
        from openmmtools.multistate import MultiStateReporter
        from openmm import unit

        reporter = MultiStateReporter(
            nc_file, open_mode="r", checkpoint_storage=chk_file
        )
        try:
            mcmove = reporter.read_mcmc_moves()[0]
            dt = mcmove.timestep.value_in_unit(unit.picosecond) * mcmove.n_steps
        finally:
            reporter.close()
        return dt
    except Exception:
        return 5.0


def _format_bytes(size_bytes: int) -> str:
    """Format bytes to human readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def _estimate_chunk_size(
    n_process_frames: int,
    n_replicas: int,
    n_pairs: int,
    num_workers: int,
) -> int:
    """Automatically estimate optimal chunk_size based on worker count.

    Targets one task per worker since dRMS has uniform per-frame cost.
    Bounded by memory limit (500 MB) and absolute limits [10, 1000].

    Args:
        n_process_frames: Number of frames to process.
        n_replicas: Number of replicas.
        n_pairs: Number of atom pairs (for dRMS calculation).
        num_workers: Number of worker processes.

    Returns:
        Optimal chunk_size.
    """
    if n_process_frames <= 0:
        return 1

    chunk_size = max(10, n_process_frames // max(num_workers, 1))

    bytes_per_frame = n_replicas * n_pairs * 3 * 4
    if bytes_per_frame > 0:
        max_mem = 500 * 1024 * 1024
        max_cs_memory = max(10, int(max_mem / bytes_per_frame))
        chunk_size = min(chunk_size, max_cs_memory)

    chunk_size = min(chunk_size, 1000, n_process_frames)
    chunk_size = max(10, chunk_size)

    return chunk_size


def calculate_reference_distances(
    ref_files: list[str],
    selected_atom: str = "name BB",
    min_dist: float = 6.0,
    max_dist: float = 50.0,
    min_diff: float = 5.0,
    excl_res: int = 4
) -> list[np.ndarray]:
    """Calculate reference distances from multiple structure files.

    This function mirrors the behavior of drms_analysis.DRMSAnalyzer.
    It filters atom pairs based on distance differences across multiple reference states,
    only keeping pairs where the minimum difference between any two states is >= min_diff.
    Each state's distances are filtered independently for the distance range.

    Args:
        ref_files: List of reference structure files (pdb/gro).
        selected_atom: Atom selection string.
        min_dist: Minimum distance in Angstrom.
        max_dist: Maximum distance in Angstrom.
        min_diff: Minimum difference between states
            (pairs with smaller differences are excluded).
        excl_res: Exclude residues within this separation.

    Returns:
        List of reference distance arrays, one per state (n_pairs x 3):
        [atom_i, atom_j, distance].
    """
    import MDAnalysis as mda

    ref_universes = [mda.Universe(f) for f in ref_files]
    ref_atoms = [u.select_atoms(selected_atom) for u in ref_universes]

    base_count = len(ref_atoms[0])
    for i, atoms in enumerate(ref_atoms[1:], 1):
        if len(atoms) != base_count:
            raise ValueError(f"Reference structure {i+1} has different atom count!")
        for j in range(base_count):
            if (atoms[j].segid != ref_atoms[0][j].segid or
                atoms[j].resid != ref_atoms[0][j].resid):
                raise ValueError(f"Atom {j} segid/resid mismatch in references")

    n_atoms = base_count
    n_states = len(ref_atoms)
    results = [[] for _ in range(n_states)]

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            atom_i, atom_j = ref_atoms[0][i], ref_atoms[0][j]

            if atom_i.segid == atom_j.segid and atom_i.resid + excl_res > atom_j.resid:
                continue

            dists = []
            for atoms in ref_atoms:
                pos_i = atoms[i].position
                pos_j = atoms[j].position
                dist = np.linalg.norm(pos_i - pos_j)
                dists.append(dist)

            min_difference = min(abs(dists[k] - dists[l])
                                for k in range(n_states)
                                for l in range(k + 1, n_states))

            if min_difference >= min_diff:
                for state_idx in range(n_states):
                    dist = dists[state_idx]
                    if min_dist <= dist <= max_dist:
                        results[state_idx].append([atoms[i].ix, atoms[j].ix, dist])

    return [np.array(state_data) if state_data else np.array([]) for state_data in results]


# ---- Worker process globals (set once per process via initializer) ----

_CHK_HANDLE = None
_UNIQUE_ATOMS = None
_STATE_IDX_I = None
_STATE_IDX_J = None
_STATE_REF_DIST = None


def _worker_setup(chk_file, ref_distances_list):
    """Initialize worker process globals.

    Called once per worker process by ProcessPoolExecutor initializer.
    Opens an h5py handle (reused across all tasks in this worker) and
    pre-computes atom index mappings to avoid per-task pickle overhead.

    Args:
        chk_file: Path to the checkpoint HDF5 file.
        ref_distances_list: List of reference distance arrays (n_pairs x 3).
    """
    import atexit
    import numpy as np
    import h5py

    global _CHK_HANDLE, _UNIQUE_ATOMS
    global _STATE_IDX_I, _STATE_IDX_J, _STATE_REF_DIST

    _CHK_HANDLE = h5py.File(chk_file, "r")
    atexit.register(_CHK_HANDLE.close)

    atom_groups = []
    for ref_dist in ref_distances_list:
        atom_groups.append(ref_dist[:, 0].astype(int))
        atom_groups.append(ref_dist[:, 1].astype(int))

    _UNIQUE_ATOMS = np.unique(np.concatenate(atom_groups))
    orig_to_compact = {orig: compact for compact, orig in enumerate(_UNIQUE_ATOMS)}

    _STATE_IDX_I = []
    _STATE_IDX_J = []
    _STATE_REF_DIST = []
    for ref_dist in ref_distances_list:
        _STATE_IDX_I.append(
            np.array([orig_to_compact[int(idx)] for idx in ref_dist[:, 0]])
        )
        _STATE_IDX_J.append(
            np.array([orig_to_compact[int(idx)] for idx in ref_dist[:, 1]])
        )
        _STATE_REF_DIST.append(ref_dist[:, 2])


def _worker_chunk(args):
    """Process a chunk of frames using the process-global HDF5 handle.

    Reads positions for all frames in the chunk via a single HDF5 slice,
    then computes dRMS for each reference state. Reference data and
    atom index mappings are accessed from process globals (set by initializer).

    Args:
        args: Tuple of (frame_indices, replica_indices, dt).

    Returns:
        List of (times, drms) tuples, one per reference state.
    """
    import numpy as np

    frame_indices, replica_indices, dt = args

    pos_all = _CHK_HANDLE["positions"][
        list(frame_indices), :, :, :
    ][:, :, _UNIQUE_ATOMS, :]

    # Filter to requested replicas (read all from HDF5, slice in numpy)
    pos_all = pos_all[:, replica_indices, :, :]

    pos_all *= 10.0  # nm to Angstrom

    n_states = len(_STATE_IDX_I)
    results = []

    for state_idx in range(n_states):
        pos_i = pos_all[:, :, _STATE_IDX_I[state_idx], :]
        pos_j = pos_all[:, :, _STATE_IDX_J[state_idx], :]
        ref_dist = _STATE_REF_DIST[state_idx]

        distances = np.linalg.norm(pos_i - pos_j, axis=-1)
        drms = np.sqrt(np.mean((distances - ref_dist) ** 2, axis=-1))

        times = np.array(frame_indices, dtype=np.float32) * dt
        results.append((times, drms))

    return results


def calculate_trajectory_drms(
    netcdf_file: str,
    checkpoint_file: str,
    ref_distances: np.ndarray,
    replica_indices: Optional[np.ndarray] = None,
    skip: int = 1,
    num_workers: Optional[int] = None,
    chunk_size: Optional[int] = None,
    auto_chunk: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate dRMS for entire trajectory using parallel processing.

    Uses direct h5py reads and ProcessPoolExecutor with per-worker
    HDF5 handle reuse. Reference data is loaded once per worker via
    initializer, eliminating per-task pickle overhead.

    Args:
        netcdf_file: Path to NetCDF file (for timestep metadata).
        checkpoint_file: Path to checkpoint file (for positions).
        ref_distances: Reference distances array (n_pairs x 3).
        replica_indices: Replica indices to analyze (None for all).
        skip: Process every N-th frame.
        num_workers: Number of parallel workers (default: CPU count).
        chunk_size: Frames per worker task (auto-optimized if None).
        auto_chunk: Whether to automatically optimize chunk_size.

    Returns:
        Tuple of (times, drms_array) where drms_array is (n_frames, n_replicas).
    """
    n_frames, n_replicas, _n_atoms = _get_metadata(checkpoint_file)
    dt = _get_timestep(netcdf_file, checkpoint_file)

    if replica_indices is None:
        replica_indices = np.arange(n_replicas)

    frame_indices = np.arange(0, n_frames, skip)
    n_process_frames = len(frame_indices)

    if n_process_frames == 0:
        return np.array([]), np.array([])

    print(f"Processing {n_process_frames} frames (every {skip} of {n_frames})...")
    print(f"Analyzing {len(replica_indices)} replicas")

    if num_workers is None:
        num_workers = multiprocessing.cpu_count()
    num_workers = min(num_workers, n_process_frames)

    n_pairs = len(ref_distances)
    if chunk_size is None and auto_chunk:
        chunk_size = _estimate_chunk_size(
            n_process_frames=n_process_frames,
            n_replicas=len(replica_indices),
            n_pairs=n_pairs,
            num_workers=num_workers,
        )
        mem_per_chunk = len(replica_indices) * n_pairs * 3 * 4 * chunk_size
        print(f"Auto-optimized chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    elif chunk_size is None:
        chunk_size = 100
        print(f"Using default chunk_size={chunk_size}")
    else:
        mem_per_chunk = len(replica_indices) * n_pairs * 3 * 4 * chunk_size
        print(f"Using specified chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")

    tasks = []
    for i in range(0, n_process_frames, chunk_size):
        chunk_end = min(i + chunk_size, n_process_frames)
        chunk_frames = frame_indices[i:chunk_end]
        tasks.append((chunk_frames, replica_indices, dt))

    print(f"Using {num_workers} workers for {len(tasks)} chunks")

    start_time = time.time()
    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_worker_setup,
        initargs=(checkpoint_file, [ref_distances]),
    ) as executor:
        chunk_results = list(executor.map(_worker_chunk, tasks))
    elapsed = time.time() - start_time

    all_times = np.concatenate([r[0][0] for r in chunk_results])
    all_drms = np.vstack([r[0][1] for r in chunk_results])

    print(f"Total time: {elapsed:.1f}s ({n_process_frames/elapsed:.1f} fps)")

    return all_times, all_drms


def save_drms_results(
    times: np.ndarray,
    drms_data: np.ndarray,
    output_file: str,
    netcdf_file: str = "",
    checkpoint_file: str = "",
    replica_indices: Optional[np.ndarray] = None
):
    """Save dRMS results to file.

    Args:
        times: Time array (n_frames,).
        drms_data: dRMS array (n_frames, n_replicas).
        output_file: Output file path.
        netcdf_file: NetCDF path for header comment.
        checkpoint_file: Checkpoint path for header comment.
        replica_indices: Replica indices for header comment.
    """
    header_lines = [
        f"# time(ps), dRMS values for each replica",
        f"# NetCDF: {netcdf_file}",
        f"# Checkpoint: {checkpoint_file}",
    ]
    if replica_indices is not None:
        header_lines.append(f"# Replicas: {','.join(map(str, replica_indices))}")
    header = "\n".join(header_lines)

    data = np.column_stack([times, drms_data])
    fmt = ['%.1f'] + ['%.4f'] * drms_data.shape[1]
    np.savetxt(output_file, data, header=header, fmt=fmt)
    print(f"Saved dRMS results to: {output_file}")


def main():
    """Command-line interface for REMD dRMS analysis."""
    parser = argparse.ArgumentParser(
        description="Calculate REMD dRMS trajectories",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-nc", "--netcdf", required=True,
                        help="NetCDF trajectory file")
    parser.add_argument("-c", "--checkpoint", default="output_checkpoint.nc",
                        help="Checkpoint NetCDF file")
    parser.add_argument("-ref", "--ref-files", nargs="+", required=True,
                        help="Reference structure files (pdb/gro)")

    parser.add_argument("-o", "--output", default="dRMSTraj_nc",
                        help="Output prefix")

    parser.add_argument("--sel", "--selected-atom", default="name BB",
                        help="Atom selection string")
    parser.add_argument("--min-dist", type=float, default=6.0,
                        help="Minimum distance in Angstrom")
    parser.add_argument("--max-dist", type=float, default=50.0,
                        help="Maximum distance in Angstrom")
    parser.add_argument("--min-diff", type=float, default=5.0,
                        help="Minimum difference between references")
    parser.add_argument("--excl-res", type=int, default=4,
                        help="Residue exclusion distance")
    parser.add_argument("--skip", type=int, default=1,
                        help="Process every N-th frame")
    parser.add_argument("--num-workers", type=int, default=None,
                        help="Number of parallel workers (default: CPU count)")
    parser.add_argument("--chunk-size", type=int, default=None,
                        help="Frames per chunk (auto-optimized if not specified)")
    parser.add_argument("--auto-chunk", action="store_true", default=True,
                        help="Automatically optimize chunk size (enabled by default)")
    parser.add_argument("--no-auto-chunk", action="store_false", dest="auto_chunk",
                        help="Disable automatic chunk size optimization")
    parser.add_argument("--replicas", type=str, default="all",
                        help="Comma-separated replica indices or 'all'")

    args = parser.parse_args()

    print("\nCalculating reference distances...")
    states_str = 'ABCDEFGHIJKLMN'

    all_ref_distances = calculate_reference_distances(
        args.ref_files,
        selected_atom=args.sel,
        min_dist=args.min_dist,
        max_dist=args.max_dist,
        min_diff=args.min_diff,
        excl_res=args.excl_res
    )

    n_states = len(all_ref_distances)
    for i, ref_dist in enumerate(all_ref_distances):
        print(f"  State {states_str[i]}: {len(ref_dist)} distance pairs")

    if args.replicas == "all":
        replica_indices = None
    else:
        replica_indices = np.array([int(x) for x in args.replicas.split(",")])

    n_frames, n_replicas, _n_atoms = _get_metadata(args.checkpoint)
    dt = _get_timestep(args.netcdf, args.checkpoint)

    if replica_indices is None:
        replica_indices = np.arange(n_replicas)

    frame_indices = np.arange(0, n_frames, args.skip)
    n_process_frames = len(frame_indices)

    num_workers = args.num_workers if args.num_workers else multiprocessing.cpu_count()

    print(f"\nProcessing {n_process_frames} frames (every {args.skip} of {n_frames})...")

    max_pairs = max(len(ref_dist) for ref_dist in all_ref_distances)
    if args.chunk_size is None and args.auto_chunk:
        chunk_size = _estimate_chunk_size(
            n_process_frames=n_process_frames,
            n_replicas=len(replica_indices),
            n_pairs=max_pairs,
            num_workers=num_workers,
        )
        mem_per_chunk = len(replica_indices) * max_pairs * 3 * 4 * chunk_size
        print(f"Auto-optimized chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    elif args.chunk_size is not None:
        chunk_size = args.chunk_size
        mem_per_chunk = len(replica_indices) * max_pairs * 3 * 4 * chunk_size
        print(f"Using specified chunk_size={chunk_size} "
              f"(~{_format_bytes(mem_per_chunk)} per chunk)")
    else:
        chunk_size = 100
        print(f"Using default chunk_size={chunk_size}")

    print(f"Using {num_workers} workers")

    tasks = []
    for i in range(0, n_process_frames, chunk_size):
        chunk_end = min(i + chunk_size, n_process_frames)
        chunk_frames = frame_indices[i:chunk_end]
        tasks.append((chunk_frames, replica_indices, dt))

    start_time = time.time()
    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_worker_setup,
        initargs=(args.checkpoint, all_ref_distances),
    ) as executor:
        all_chunk_results = list(executor.map(_worker_chunk, tasks))

    elapsed = time.time() - start_time
    print(f"\nComputation time: {elapsed:.1f}s ({n_process_frames/elapsed:.1f} fps)")

    for state_idx in range(n_states):
        print(f"\nSaving State {states_str[state_idx]}...")
        state_times = np.concatenate(
            [chunk_res[state_idx][0] for chunk_res in all_chunk_results]
        )
        state_drms = np.vstack(
            [chunk_res[state_idx][1] for chunk_res in all_chunk_results]
        )

        output_file = f"{args.output}_State{states_str[state_idx]}.dat"
        save_drms_results(
            state_times, state_drms, output_file,
            args.netcdf, args.checkpoint, replica_indices
        )


if __name__ == "__main__":
    main()
