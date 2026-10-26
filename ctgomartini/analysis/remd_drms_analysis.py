#!/usr/bin/env python3
"""REMD dRMS trajectory analysis.

Analyzes distance root-mean-square (dRMS) trajectories from the per-replica
XTC trajectories (``replica_<i>.xtc``) of an OpenMM >= 8.6 built-in REMD
output directory. The XTC trajectories carry no topology, so a topology file
(pdb/gro, ``--top``) is required to interpret the atom indices.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np


def _get_timestep(output_dir: str | Path) -> float:
    """Read the frame timestep from a built-in REMD output directory.

    Computes ``dt * stepsPerIteration * reportInterval`` from metadata.json
    (one log/trajectory row every reportInterval iterations of
    stepsPerIteration MD steps each). Falls back to 5.0 ps on error.

    Args:
        output_dir: Path to the built-in REMD output directory.

    Returns:
        Timestep in picoseconds per frame.
    """
    try:
        metadata = json.loads((Path(output_dir) / "metadata.json").read_text())
        dt = float(metadata["dt"])
        steps_per_iteration = float(
            metadata.get("stepsPerIteration", metadata["exc_freq"])
        )
        report_interval = float(metadata["reportInterval"])
        return dt * steps_per_iteration * report_interval
    except Exception:
        return 5.0


def _get_n_replicas(output_dir: str | Path) -> int:
    """Read the replica count from metadata.json."""
    metadata = json.loads((Path(output_dir) / "metadata.json").read_text())
    return int(metadata["replica_count"])


def _count_xtc_frames(xtc_path: str | Path, topology_file: str | Path) -> int:
    """Return the number of frames of an XTC trajectory.

    Uses the XTC frame offsets (no trajectory decoding).

    Raises:
        FileNotFoundError: If the trajectory does not exist.
    """
    import MDAnalysis as mda  # noqa: N813

    if not Path(xtc_path).is_file():
        raise FileNotFoundError(
            f"Replica trajectory {xtc_path} does not exist. Run REMD with "
            f"'remd_xtc_output = replica' (or 'both') to write per-replica "
            f"XTC trajectories."
        )
    return len(mda.Universe(str(topology_file), str(xtc_path)).trajectory)


def _format_bytes(size_bytes: int) -> str:
    """Format bytes to human readable string."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} TB"


def _run_with_progress(executor, tasks):
    """Execute tasks via ProcessPoolExecutor with progress and ETA reporting.

    Submits all tasks, then collects results as they complete,
    printing progress at most once every 2 seconds.

    Args:
        executor: ProcessPoolExecutor instance.
        tasks: List of task tuples.

    Returns:
        List of results in submission order.
    """
    n_total = len(tasks)
    futures = {executor.submit(_worker_replica, t): i for i, t in enumerate(tasks)}
    results = [None] * n_total
    t_start = time.perf_counter()
    last_print = 0.0
    print_interval = 2.0  # seconds

    for n_done, future in enumerate(as_completed(futures), start=1):
        idx = futures[future]
        results[idx] = future.result()

        now = time.perf_counter()
        if now - last_print >= print_interval or n_done == n_total:
            elapsed = now - t_start
            fps = n_done / elapsed if elapsed > 0 else 0
            eta = (n_total - n_done) / fps if fps > 0 else 0
            pct = n_done * 100 / n_total
            print(
                f"  [{n_done}/{n_total} replicas, {pct:.0f}%] "
                f"elapsed: {elapsed:.0f}s, ETA: {eta:.0f}s, "
                f"speed: {fps:.1f} replicas/s"
            )
            last_print = now

    return results


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
    import MDAnalysis as mda  # noqa: N813

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

            min_difference = min(abs(dists[k] - dists[m])
                                for k in range(n_states)
                                for m in range(k + 1, n_states))

            if min_difference >= min_diff:
                for state_idx in range(n_states):
                    dist = dists[state_idx]
                    if min_dist <= dist <= max_dist:
                        results[state_idx].append([atoms[i].ix, atoms[j].ix, dist])

    return [np.array(state_data) if state_data else np.array([]) for state_data in results]


# ---- Worker process globals (set once per process via initializer) ----

_TOP_FILE = None
_UNIQUE_ATOMS = None
_STATE_IDX_I = None
_STATE_IDX_J = None
_STATE_REF_DIST = None


def _worker_setup(topology_file, ref_distances_list):
    """Initialize worker process globals.

    Called once per worker process by ProcessPoolExecutor initializer.
    Pre-computes atom index mappings to avoid per-task pickle overhead.

    Args:
        topology_file: Topology file (pdb/gro) for reading the XTC
            trajectories with MDAnalysis.
        ref_distances_list: List of reference distance arrays (n_pairs x 3).
    """
    global _TOP_FILE, _UNIQUE_ATOMS
    global _STATE_IDX_I, _STATE_IDX_J, _STATE_REF_DIST

    _TOP_FILE = str(topology_file)

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


def _worker_replica(args):
    """Process all frames of one per-replica XTC trajectory.

    Loads the trajectory with MDAnalysis (subsetting to the atoms used by the
    reference pairs), then computes dRMS for each reference state. Reference
    data and atom index mappings are accessed from process globals (set by
    initializer).

    Args:
        args: Tuple of (replica_index, xtc_path, frame_indices, dt).

    Returns:
        List of (times, drms) tuples, one per reference state; drms has
        shape (n_frames,) for this replica.
    """
    import MDAnalysis as mda  # noqa: N813

    _replica_index, xtc_path, frame_indices, dt = args

    universe = mda.Universe(_TOP_FILE, xtc_path)
    atoms = universe.atoms[_UNIQUE_ATOMS]
    pos = np.empty((len(frame_indices), len(_UNIQUE_ATOMS), 3))
    for k, frame_index in enumerate(frame_indices):
        universe.trajectory[int(frame_index)]
        pos[k] = atoms.positions  # already in Angstrom

    n_states = len(_STATE_IDX_I)
    results = []

    for state_idx in range(n_states):
        pos_i = pos[:, _STATE_IDX_I[state_idx], :]
        pos_j = pos[:, _STATE_IDX_J[state_idx], :]
        ref_dist = _STATE_REF_DIST[state_idx]

        distances = np.linalg.norm(pos_i - pos_j, axis=-1)
        drms = np.sqrt(np.mean((distances - ref_dist) ** 2, axis=-1))

        times = np.array(frame_indices, dtype=np.float64) * dt
        results.append((times, drms))

    return results


def _calculate_all_states(
    output_dir: str | Path,
    topology_file: str,
    ref_distances_list: list[np.ndarray],
    replica_indices: np.ndarray,
    frame_indices: np.ndarray,
    dt: float,
    num_workers: int,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Compute dRMS against every reference state for the given replicas.

    One worker task per replica (each loads its own ``replica_<i>.xtc``).

    Args:
        output_dir: Built-in REMD output directory with replica_<i>.xtc.
        topology_file: Topology file (pdb/gro) for MDAnalysis.
        ref_distances_list: Reference distance arrays, one per state.
        replica_indices: Replica indices to analyze.
        frame_indices: Frame indices to process (already begin/end/skip).
        dt: Timestep per frame in ps.
        num_workers: Number of parallel workers.

    Returns:
        List of (times, drms) per state, drms with shape
        (n_frames, n_replicas).
    """
    output_dir = Path(output_dir)
    tasks = []
    for replica_index in replica_indices:
        xtc_path = output_dir / f"replica_{replica_index}.xtc"
        if not xtc_path.is_file():
            raise FileNotFoundError(
                f"Replica trajectory {xtc_path} does not exist. Run REMD with "
                f"'remd_xtc_output = replica' (or 'both') to write per-replica "
                f"XTC trajectories."
            )
        tasks.append((int(replica_index), str(xtc_path), frame_indices, dt))

    num_workers = min(num_workers, len(tasks))
    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=_worker_setup,
        initargs=(str(topology_file), ref_distances_list),
    ) as executor:
        per_replica = _run_with_progress(executor, tasks)

    n_states = len(ref_distances_list)
    results = []
    for state_idx in range(n_states):
        times = per_replica[0][state_idx][0]
        drms = np.column_stack(
            [per_replica[k][state_idx][1] for k in range(len(replica_indices))]
        )
        results.append((times, drms))
    return results


def calculate_trajectory_drms(
    output_dir: str | Path,
    topology_file: str,
    ref_distances: np.ndarray,
    replica_indices: np.ndarray | None = None,
    skip: int = 1,
    num_workers: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Calculate dRMS for entire trajectory using parallel processing.

    Reads the per-replica XTC trajectories of a built-in REMD output
    directory with MDAnalysis, one worker task per replica. Reference data
    is loaded once per worker via initializer, eliminating per-task pickle
    overhead.

    Args:
        output_dir: Path to the built-in REMD output directory.
        topology_file: Topology file (pdb/gro) for reading the XTC
            trajectories with MDAnalysis.
        ref_distances: Reference distances array (n_pairs x 3).
        replica_indices: Replica indices to analyze (None for all).
        skip: Process every N-th frame.
        num_workers: Number of parallel workers (default: CPU count).

    Returns:
        Tuple of (times, drms_array) where drms_array is (n_frames, n_replicas).
    """
    output_dir = Path(output_dir)
    if not output_dir.is_dir():
        raise ValueError(
            f"calculate_trajectory_drms requires a built-in REMD output "
            f"directory (containing replica_<i>.xtc); got {output_dir}. "
            f"openmmtools NetCDF (.nc) output is no longer supported."
        )
    if replica_indices is None:
        replica_indices = np.arange(_get_n_replicas(output_dir))
    n_frames = _count_xtc_frames(
        output_dir / f"replica_{replica_indices[0]}.xtc", topology_file
    )
    dt = _get_timestep(output_dir)

    frame_indices = np.arange(0, n_frames, skip)
    n_process_frames = len(frame_indices)

    if n_process_frames == 0:
        return np.array([]), np.array([])

    print(f"Processing {n_process_frames} frames (every {skip} of {n_frames})...")
    print(f"Analyzing {len(replica_indices)} replicas")

    if num_workers is None:
        num_workers = multiprocessing.cpu_count()

    start_time = time.time()
    all_times, all_drms = _calculate_all_states(
        output_dir,
        topology_file,
        [ref_distances],
        replica_indices,
        frame_indices,
        dt,
        num_workers,
    )[0]
    elapsed = time.time() - start_time

    print(f"Total time: {elapsed:.1f}s ({n_process_frames/elapsed:.1f} fps)")

    return all_times, all_drms


def save_drms_results(
    times: np.ndarray,
    drms_data: np.ndarray,
    output_file: str,
    output_dir: str | Path = "",
    replica_indices: np.ndarray | None = None
):
    """Save dRMS results to file.

    Args:
        times: Time array (n_frames,).
        drms_data: dRMS array (n_frames, n_replicas).
        output_file: Output file path.
        output_dir: REMD output directory path for header comment.
        replica_indices: Replica indices for header comment.
    """
    header_lines = [
        "# time(ps), dRMS values for each replica",
        f"# Output directory: {output_dir}",
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

    parser.add_argument("-d", "--dir", default="output_remd",
                        help="Built-in REMD output directory")
    parser.add_argument("--top", "--topology", required=True,
                        help="Topology file (pdb/gro) for reading the XTC trajectories")
    parser.add_argument("-ref", "--ref-files", nargs="+", required=True,
                        help="Reference structure files (pdb/gro)")

    parser.add_argument("-o", "--output", default="dRMSTraj",
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
    parser.add_argument(
        "-b", "--begin", type=int, default=0,
        help="Start frame (inclusive)",
    )
    parser.add_argument(
        "-e", "--end", type=int, default=None,
        help="End frame (exclusive, default: all)",
    )
    parser.add_argument("--skip", type=int, default=1,
                        help="Process every N-th frame")
    parser.add_argument("-n", "--num-workers", type=int, default=10,
                        help="Number of parallel workers")
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

    output_dir = Path(args.dir)
    if not output_dir.is_dir():
        print(f"Error: {output_dir} is not a built-in REMD output directory")
        return

    if args.replicas == "all":
        replica_indices = np.arange(_get_n_replicas(output_dir))
    else:
        replica_indices = np.array([int(x) for x in args.replicas.split(",")])

    n_frames = _count_xtc_frames(
        output_dir / f"replica_{replica_indices[0]}.xtc", args.top
    )
    dt = _get_timestep(output_dir)

    frame_end = args.end if args.end is not None else n_frames
    frame_indices = np.arange(args.begin, frame_end, args.skip)
    n_process_frames = len(frame_indices)

    if n_process_frames == 0:
        print("Error: zero frames to process (check --begin, --end, --skip)")
        return

    print(f"\nProcessing {n_process_frames} frames (every {args.skip} of {n_frames})...")
    print(f"Using {args.num_workers} workers")

    start_time = time.time()
    per_state_results = _calculate_all_states(
        output_dir,
        args.top,
        all_ref_distances,
        replica_indices,
        frame_indices,
        dt,
        args.num_workers,
    )

    elapsed = time.time() - start_time
    print(f"\nComputation time: {elapsed:.1f}s ({n_process_frames/elapsed:.1f} fps)")

    for state_idx in range(n_states):
        print(f"\nSaving State {states_str[state_idx]}...")
        state_times, state_drms = per_state_results[state_idx]

        output_file = f"{args.output}_State{states_str[state_idx]}.dat"
        save_drms_results(
            state_times, state_drms, output_file,
            output_dir=output_dir, replica_indices=replica_indices
        )


if __name__ == "__main__":
    main()
