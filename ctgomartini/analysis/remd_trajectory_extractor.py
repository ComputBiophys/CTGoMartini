#!/usr/bin/env python3
"""
REMD trajectory extraction tool.

Extracts replica trajectories from NetCDF checkpoint files using block-level
parallel reads (multi-process HDF5) and parallel writes (ctypes xdrfile append).

Architecture:
  1. Split total frames into blocks that each fit in RAM
  2. Each block: PP×N h5py readers → shared memory → PP×M xdrfile writers
  3. First block: 'w' mode (create XTC). Subsequent blocks: 'a' mode (append).

For datasets that fit entirely in RAM, set block_frames >= total frames for
maximum throughput (single-block path, zero append overhead).
"""

from __future__ import annotations

import argparse
import ctypes
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import shared_memory
from typing import Optional

import h5py
import numpy as np

warnings.filterwarnings("ignore")

# ---- xdrfile C library setup ----
# The xdrfile C library (bundled with MDAnalysis) supports 'a' (append) mode,
# which neither MDTraj nor MDAnalysis expose in their Python wrappers.
# We call it directly via ctypes.

_xdr = None

def _xdr_setup():
    """Load the xdrfile shared library and configure ctypes signatures."""
    global _xdr
    if _xdr is not None:
        return

    import MDAnalysis.lib.formats.libmdaxdr as _mod

    so_dir = os.path.dirname(_mod.__file__)
    for f in os.listdir(so_dir):
        if f.startswith("libmdaxdr") and f.endswith(".so"):
            _xdr = ctypes.CDLL(os.path.join(so_dir, f))
            break
    else:
        raise RuntimeError("Cannot find libmdaxdr shared library")

    _xdr.xdrfile_open.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    _xdr.xdrfile_open.restype = ctypes.c_void_p
    _xdr.write_xtc.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
        ctypes.c_float, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_float,
    ]
    _xdr.write_xtc.restype = ctypes.c_int
    _xdr.xdrfile_close.argtypes = [ctypes.c_void_p]
    _xdr.xdrfile_close.restype = ctypes.c_int


# Shared memory globals for writer worker processes
_SHM_DTYPE = np.float32
_SHM_NAME = ""
_SHM_SHAPE = ()


def _get_timestep(nc_file: str, chk_file: str) -> float:
    """Read the frame timestep from an openmmtools NetCDF file.

    Args:
        nc_file: Path to the main output.nc file.
        chk_file: Path to the checkpoint file.

    Returns:
        Timestep in picoseconds per frame. Defaults to 5.0 ps on failure.
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


def _reader_worker(args):
    """Read a frame range from HDF5 directly into a shared memory slice.

    Args:
        args: Tuple of (chk_file, shm_name, shm_shape, r_start_local,
              r_end_global, r_start_global, _unused, stride).
    """
    chk_file, shm_name, shm_shape, r_start_local, r_end_global, r_start_global, _unused, stride = args

    shm = shared_memory.SharedMemory(name=shm_name)
    try:
        buf = np.ndarray(shm_shape, dtype=_SHM_DTYPE, buffer=shm.buf)
        with h5py.File(chk_file, "r") as f:
            n_local = (r_end_global - r_start_global) // stride
            buf[r_start_local : r_start_local + n_local, :, :, :] = (
                f["positions"][r_start_global:r_end_global:stride, :, :, :]
            )
    finally:
        shm.close()


def _writer_setup(shm_name, shm_shape):
    """Initialize writer worker globals with shared memory info.

    Called once per worker process by ProcessPoolExecutor initializer.

    Args:
        shm_name: Name of the shared memory block.
        shm_shape: Shape of the numpy array in shared memory.
    """
    global _SHM_NAME, _SHM_SHAPE
    _SHM_NAME = shm_name
    _SHM_SHAPE = shm_shape


def _writer_worker(args):
    """Write one replica from shared memory to XTC via ctypes xdrfile.

    Uses 'w' mode for the first block (creates file) and 'a' mode for
    subsequent blocks (appends frames).

    Args:
        args: Tuple of (rep_idx, _unused, block_frames, times, steps,
              path, is_first_block).

    Returns:
        Tuple of (path, file_size_bytes).
    """
    rep_idx, _unused, block_frames, times, steps, path, is_first_block = args

    _xdr_setup()

    shm = shared_memory.SharedMemory(name=_SHM_NAME)
    try:
        data = np.ndarray(_SHM_SHAPE, dtype=_SHM_DTYPE, buffer=shm.buf)
        rep_pos = np.ascontiguousarray(data[:, rep_idx, :, :], dtype=np.float32)

        mode = b"w" if is_first_block else b"a"
        xfp = _xdr.xdrfile_open(path.encode(), mode)
        if not xfp:
            raise RuntimeError(f"xdrfile_open({path}) failed")

        box = np.zeros((3, 3), dtype=np.float32)
        try:
            for i in range(block_frames):
                ret = _xdr.write_xtc(
                    xfp,
                    rep_pos.shape[1],
                    int(steps[i]),
                    float(times[i]),
                    box.ctypes.data,
                    rep_pos[i].ctypes.data,
                    1000.0,
                )
                if ret != 0:
                    raise RuntimeError(
                        f"write_xtc failed at frame {i}: ret={ret}"
                    )
        finally:
            _xdr.xdrfile_close(xfp)
    finally:
        shm.close()
    return path, os.path.getsize(path)


def extract_replicas(
    nc_file: str,
    chk_file: str,
    pdb_file: str,
    output_dir: str,
    output_pattern: str = "replica_{i}.xtc",
    frame_begin: int = 0,
    frame_end: Optional[int] = None,
    frame_stride: int = 1,
    replicas: Optional[list[int]] = None,
    block_frames: int = 5000,
    num_readers: int = 8,
    num_writers: int = 21,
):
    """Extract replica trajectories from a REMD checkpoint.

    Splits the frame range into blocks that fit in RAM. Each block is
    processed with parallel HDF5 reads and parallel xdrfile writes.
    Blocks are appended to the same XTC files via ctypes xdrfile 'a' mode.

    For datasets that fit entirely in RAM, set block_frames >= total frames
    to use a single block (maximum throughput).

    Args:
        nc_file: Path to the main output.nc file (for timestep metadata).
        chk_file: Path to the checkpoint NetCDF file.
        pdb_file: Path to the topology file (PDB or GRO).
        output_dir: Directory for output XTC files.
        output_pattern: Output filename pattern with {i} placeholder.
        frame_begin: First frame index (inclusive).
        frame_end: Last frame index (exclusive, None = all).
        frame_stride: Frame stride.
        replicas: List of replica indices to extract (None = all).
        block_frames: Frames per block (tune for available RAM).
        num_readers: Number of parallel HDF5 reader processes.
        num_writers: Number of parallel XTC writer processes.
    """
    import MDAnalysis as mda

    os.makedirs(output_dir, exist_ok=True)

    # Validate atom count
    u = mda.Universe(pdb_file)
    n_atoms_topo = len(u.atoms)

    with h5py.File(chk_file, "r") as f:
        pos_var = f["positions"]
        n_frames_total, n_replicas_total, n_atoms_nc, _ = pos_var.shape

    if n_atoms_nc != n_atoms_topo:
        raise ValueError(
            f"Atom count mismatch: topology has {n_atoms_topo}, "
            f"checkpoint has {n_atoms_nc}"
        )

    frame_end = frame_end or n_frames_total
    if replicas is None:
        replicas = list(range(n_replicas_total))

    frame_indices = list(range(frame_begin, frame_end, frame_stride))
    n_frames_out = len(frame_indices)
    frame_dt = _get_timestep(nc_file, chk_file)
    block_frames = min(block_frames, n_frames_out)

    n_replicas_sel = len(replicas)
    block_gb = block_frames * n_replicas_total * n_atoms_nc * 3 * 4 / 1e9
    n_blocks = (n_frames_out + block_frames - 1) // block_frames

    print(
        f"Extracting {n_frames_out} frames × {n_replicas_sel} replicas "
        f"× {n_atoms_nc} atoms"
    )
    print(
        f"  Blocks: {n_blocks} × {block_frames} frames "
        f"({block_gb:.1f} GB/block)"
    )
    print(f"  Readers: {num_readers}, Writers: {num_writers}")

    t_start = time.perf_counter()
    t_read_total = 0.0
    t_write_total = 0.0

    for block_idx in range(n_blocks):
        b_start = block_idx * block_frames
        b_end = min(b_start + block_frames, n_frames_out)
        block_indices = frame_indices[b_start:b_end]
        n_block_frames = b_end - b_start
        is_first_block = block_idx == 0

        # Allocate shared memory for all replicas.
        # HDF5 chunks contain all replicas per frame, so we always read
        # the full replica dimension even when extracting a subset.
        shm_shape = (n_block_frames, n_replicas_total, n_atoms_nc, 3)
        shm_size = int(np.prod(shm_shape) * 4)
        shm = shared_memory.SharedMemory(create=True, size=shm_size)

        # Phase 1: Read block from HDF5 into shared memory.
        t0 = time.perf_counter()

        if n_block_frames <= 100:
            # Small block: read in the main process.
            r_start_global = frame_indices[b_start]
            r_end_global = frame_indices[b_end - 1] + frame_stride
            with h5py.File(chk_file, "r") as f:
                block_data = f["positions"][
                    r_start_global:r_end_global:frame_stride, :, :, :
                ]
            shm_buf = np.ndarray(shm_shape, dtype=_SHM_DTYPE, buffer=shm.buf)
            np.copyto(shm_buf, block_data)
            del block_data
        else:
            # Large block: parallel reads via ProcessPool.
            reader_args = []
            chunk_f = n_block_frames // num_readers
            for r_idx in range(num_readers):
                r_start_local = r_idx * chunk_f
                if r_idx == num_readers - 1:
                    r_end_local = n_block_frames
                else:
                    r_end_local = r_start_local + chunk_f
                if r_start_local >= n_block_frames:
                    break
                r_start_global = frame_indices[b_start + r_start_local]
                r_end_global = (
                    frame_indices[b_start + r_end_local - 1] + frame_stride
                )
                reader_args.append((
                    chk_file, shm.name, shm_shape,
                    r_start_local, r_end_local,
                    r_start_global, r_end_global, frame_stride,
                ))

            with ProcessPoolExecutor(max_workers=num_readers) as rpool:
                list(rpool.map(_reader_worker, reader_args))
        t_read_total += time.perf_counter() - t0

        # Phase 2: Build write tasks (metadata only, data in shared memory).
        times = np.array(
            [
                frame_indices[b_start + i] * frame_dt
                for i in range(n_block_frames)
            ],
            dtype=np.float32,
        )
        steps = np.array(block_indices, dtype=np.int32)

        write_tasks = []
        for local_idx, rep_idx in enumerate(replicas):
            path = os.path.join(output_dir, output_pattern.format(i=rep_idx))
            write_tasks.append((
                rep_idx, local_idx, n_block_frames,
                times, steps, path, is_first_block,
            ))

        # Phase 3: Parallel XTC writes via ctypes xdrfile.
        t0 = time.perf_counter()
        with ProcessPoolExecutor(
            max_workers=num_writers,
            initializer=_writer_setup,
            initargs=(shm.name, shm_shape),
        ) as wpool:
            list(wpool.map(_writer_worker, write_tasks))
        t_write_total += time.perf_counter() - t0

        shm.close()
        shm.unlink()

    elapsed = time.perf_counter() - t_start
    total_gb = sum(
        os.path.getsize(
            os.path.join(output_dir, output_pattern.format(i=rep_idx))
        )
        for rep_idx in replicas
    ) / 1e9
    print(
        f"Done: {n_replicas_sel} replicas in {elapsed:.1f}s "
        f"(read {t_read_total:.1f}s, write {t_write_total:.1f}s, "
        f"{total_gb:.2f} GB, {total_gb / elapsed * 1000:.0f} MB/s)"
    )


def extract_states(
    nc_file: str,
    chk_file: str,
    pdb_file: str,
    output_dir: str,
    output_pattern: str = "state_{i}.xtc",
    frame_begin: int = 0,
    frame_end: Optional[int] = None,
    frame_stride: int = 1,
    states: Optional[list[int]] = None,
    block_frames: int = 5000,
    num_readers: int = 8,
    num_writers: int = 8,
):
    """Extract thermodynamic state trajectories from a REMD simulation.

    Unlike replica trajectories, state trajectories are discontinuous:
    replicas exchange states during REMD, so a given state may be visited
    by different replicas at different frames. This function collects all
    frames belonging to each thermodynamic state.

    Uses block-level parallelism: each block of frames is read once from
    HDF5 (all replicas), then frames are distributed to per-state writers.

    Args:
        nc_file: Path to the main output.nc file (contains states mapping).
        chk_file: Path to the checkpoint NetCDF file (contains positions).
        pdb_file: Path to the topology file (PDB or GRO).
        output_dir: Directory for output XTC files.
        output_pattern: Output filename pattern with {i} placeholder.
        frame_begin: First frame index (inclusive).
        frame_end: Last frame index (exclusive, None = all).
        frame_stride: Frame stride.
        states: List of state indices to extract (None = all).
        block_frames: Frames per block (tune for available RAM).
        num_readers: Number of parallel HDF5 reader processes.
        num_writers: Number of parallel XTC writer processes.
    """
    import MDAnalysis as mda

    os.makedirs(output_dir, exist_ok=True)

    u = mda.Universe(pdb_file)
    n_atoms_topo = len(u.atoms)

    # Read state mapping from output.nc
    with h5py.File(nc_file, "r") as f:
        if "states" not in f:
            raise ValueError(
                "No 'states' variable found in output.nc. "
                "State extraction requires REMD with state information."
            )
        states_all = f["states"][:]
        n_frames_total, n_replicas_total = states_all.shape

    # Read positions metadata from checkpoint
    with h5py.File(chk_file, "r") as f:
        pos_var = f["positions"]
        n_frames_chk, n_replicas_chk, n_atoms_nc, _ = pos_var.shape

    if n_frames_chk != n_frames_total:
        print(
            f"Warning: output.nc has {n_frames_total} state records, "
            f"checkpoint has {n_frames_chk} frames. "
            f"Using min({n_frames_total}, {n_frames_chk})."
        )
        n_frames_total = min(n_frames_total, n_frames_chk)
    if n_replicas_chk != n_replicas_total:
        raise ValueError(
            f"Replica count mismatch: output.nc has {n_replicas_total}, "
            f"checkpoint has {n_replicas_chk}"
        )
    if n_atoms_nc != n_atoms_topo:
        raise ValueError(
            f"Atom count mismatch: topology has {n_atoms_topo}, "
            f"checkpoint has {n_atoms_nc}"
        )

    frame_end = frame_end or n_frames_total
    frame_indices = list(range(frame_begin, frame_end, frame_stride))
    n_frames_out = len(frame_indices)
    frame_dt = _get_timestep(nc_file, chk_file)
    block_frames = min(block_frames, n_frames_out)

    # Determine which states to extract
    unique_states = np.unique(states_all[frame_begin:frame_end:frame_stride])
    if states is not None:
        unique_states = np.array([s for s in states if s in unique_states])
    n_states_sel = len(unique_states)

    block_gb = block_frames * n_replicas_total * n_atoms_nc * 3 * 4 / 1e9
    n_blocks = (n_frames_out + block_frames - 1) // block_frames

    print(
        f"Extracting {n_states_sel} states from {n_frames_out} frames "
        f"× {n_replicas_total} replicas × {n_atoms_nc} atoms"
    )
    print(
        f"  Blocks: {n_blocks} × {block_frames} frames "
        f"({block_gb:.1f} GB/block)"
    )

    t_start = time.perf_counter()
    t_read_total = 0.0
    t_write_total = 0.0

    for block_idx in range(n_blocks):
        b_start = block_idx * block_frames
        b_end = min(b_start + block_frames, n_frames_out)
        block_indices = frame_indices[b_start:b_end]
        n_block_frames = b_end - b_start
        is_first_block = block_idx == 0

        # Read state mapping for this block (small, in memory)
        block_states = states_all[
            block_indices[0] : block_indices[-1] + frame_stride : frame_stride
        ]

        # Allocate shared memory and read ALL replicas for this block
        shm_shape = (n_block_frames, n_replicas_total, n_atoms_nc, 3)
        shm_size = int(np.prod(shm_shape) * 4)
        shm = shared_memory.SharedMemory(create=True, size=shm_size)

        # Phase 1: Read block from HDF5 (same as replica extraction).
        t0 = time.perf_counter()
        if n_block_frames <= 100:
            r_start_global = block_indices[0]
            r_end_global = block_indices[-1] + frame_stride
            with h5py.File(chk_file, "r") as f:
                block_data = f["positions"][
                    r_start_global:r_end_global:frame_stride, :, :, :
                ]
            shm_buf = np.ndarray(shm_shape, dtype=_SHM_DTYPE, buffer=shm.buf)
            np.copyto(shm_buf, block_data)
            del block_data
        else:
            reader_args = []
            chunk_f = n_block_frames // num_readers
            for r_idx in range(num_readers):
                r_start_local = r_idx * chunk_f
                if r_idx == num_readers - 1:
                    r_end_local = n_block_frames
                else:
                    r_end_local = r_start_local + chunk_f
                if r_start_local >= n_block_frames:
                    break
                r_start_global = block_indices[r_start_local]
                r_end_global = (
                    block_indices[r_end_local - 1] + frame_stride
                )
                reader_args.append((
                    chk_file, shm.name, shm_shape,
                    r_start_local, r_end_local,
                    r_start_global, r_end_global, frame_stride,
                ))
            with ProcessPoolExecutor(max_workers=num_readers) as rpool:
                list(rpool.map(_reader_worker, reader_args))
        t_read_total += time.perf_counter() - t0

        # Phase 2: Write each state's frames from shared memory.
        t0 = time.perf_counter()
        write_tasks = []
        for state_idx in unique_states:
            # Find (local_frame, replica) pairs for this state in this block
            state_pairs = []
            for local_f in range(n_block_frames):
                global_f = block_indices[local_f]
                matching = np.where(block_states[local_f] == state_idx)[0]
                for rep in matching:
                    state_pairs.append((local_f, rep, global_f))

            if state_pairs:
                path = os.path.join(output_dir, output_pattern.format(i=state_idx))
                write_tasks.append((state_idx, state_pairs, path,
                                   is_first_block, frame_dt))

        with ProcessPoolExecutor(
            max_workers=num_writers,
            initializer=_writer_setup,
            initargs=(shm.name, shm_shape),
        ) as wpool:
            list(wpool.map(_state_writer_worker, write_tasks))
        t_write_total += time.perf_counter() - t0

        shm.close()
        shm.unlink()

    elapsed = time.perf_counter() - t_start
    print(
        f"Done: {n_states_sel} states in {elapsed:.1f}s "
        f"(read {t_read_total:.1f}s, write {t_write_total:.1f}s)"
    )


def _state_writer_worker(args):
    """Write one state's block frames from shared memory to XTC.

    Args:
        args: Tuple of (state_idx, state_pairs, path, is_first_block).
            state_pairs is a list of (local_frame, replica, global_frame) tuples.
    """
    state_idx, state_pairs, path, is_first_block, frame_dt = args

    _xdr_setup()

    shm = shared_memory.SharedMemory(name=_SHM_NAME)
    try:
        data = np.ndarray(_SHM_SHAPE, dtype=_SHM_DTYPE, buffer=shm.buf)

        mode = b"w" if is_first_block else b"a"
        xfp = _xdr.xdrfile_open(path.encode(), mode)
        if not xfp:
            raise RuntimeError(f"xdrfile_open({path}) failed")

        box = np.zeros((3, 3), dtype=np.float32)
        try:
            for local_f, rep, global_f in state_pairs:
                coords = np.ascontiguousarray(
                    data[local_f, rep, :, :], dtype=np.float32
                )
                ret = _xdr.write_xtc(
                    xfp, coords.shape[0], global_f,
                    float(global_f * frame_dt),
                    box.ctypes.data, coords.ctypes.data, 1000.0,
                )
                if ret != 0:
                    raise RuntimeError(
                        f"write_xtc failed for state {state_idx}, "
                        f"frame {global_f}: ret={ret}"
                    )
        finally:
            _xdr.xdrfile_close(xfp)
    finally:
        shm.close()
    return path, os.path.getsize(path)


def extract_frame(
    chk_file: str,
    output: str,
    frame_idx: int,
    replica_idx: Optional[int] = None,
    state_idx: Optional[int] = None,
    nc_file: Optional[str] = None,
):
    """Extract a single frame from a REMD checkpoint.

    Args:
        chk_file: Path to the checkpoint NetCDF file.
        output: Output file path (extension determines format: .gro, .pdb, .xtc).
        frame_idx: Global frame index.
        replica_idx: Replica index (mutually exclusive with state_idx).
        state_idx: Thermodynamic state index (mutually exclusive with replica_idx).
        nc_file: Path to output.nc, required when using state_idx.
    """
    import mdtraj as md

    if replica_idx is not None and state_idx is not None:
        raise ValueError("Cannot specify both replica_idx and state_idx")
    if replica_idx is None and state_idx is None:
        raise ValueError("Must specify either replica_idx or state_idx")

    # Read positions
    with h5py.File(chk_file, "r") as f:
        pos_var = f["positions"]
        n_atoms = pos_var.shape[2]

        if replica_idx is not None:
            xyz = pos_var[frame_idx, replica_idx, :, :]
        else:
            # Find which replica is in the requested state
            if nc_file is None:
                raise ValueError("nc_file required for state_idx lookup")
            with h5py.File(nc_file, "r") as nf:
                states = nf["states"][:]
            row = states[frame_idx]
            matching = np.where(row == state_idx)[0]
            if len(matching) == 0:
                raise ValueError(
                    f"No replica found in state {state_idx} at frame {frame_idx}"
                )
            xyz = pos_var[frame_idx, matching[0], :, :]

    xyz = np.ascontiguousarray(xyz, dtype=np.float32).reshape(1, n_atoms, 3)

    # Write output
    ext = os.path.splitext(output)[1].lower()
    if ext == ".gro":
        _write_gro(output, xyz[0])
    elif ext == ".pdb":
        _write_pdb(output, xyz[0])
    elif ext in (".xtc", ".dcd"):
        traj = md.Trajectory(xyz, None, time=np.array([frame_idx], dtype=np.float32))
        if ext == ".xtc":
            traj.save_xtc(output)
        else:
            traj.save_dcd(output)
    else:
        raise ValueError(f"Unsupported format: {ext}")

    print(f"Frame {frame_idx} → {output}")


def _write_pdb(path: str, xyz: np.ndarray):
    """Write a single frame as PDB."""
    with open(path, "w") as f:
        f.write("CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n")
        for i, pos in enumerate(xyz):
            f.write(
                f"ATOM  {i+1:5d}  CA  ALA A   1    "
                f"{pos[0]*10:8.3f}{pos[1]*10:8.3f}{pos[2]*10:8.3f}"
                f"  1.00  0.00           C  \n"
            )
        f.write("END\n")


def _write_gro(path: str, xyz: np.ndarray):
    """Write a single frame as GRO."""
    n = len(xyz)
    with open(path, "w") as f:
        f.write("Single frame\n")
        f.write(f"{n:5d}\n")
        for i, pos in enumerate(xyz):
            f.write(
                f"{1:5d}ALA  CA  {i+1:5d}"
                f"{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n"
            )
        f.write("  10.00000  10.00000  10.00000\n")


def main():
    """Command-line interface for REMD trajectory extraction."""
    parser = argparse.ArgumentParser(
        description="Extract replica trajectories from REMD checkpoint files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-nc", "--netcdf", required=True,
        help="NetCDF trajectory file (output.nc)",
    )
    parser.add_argument(
        "-p", "--topology", required=True,
        help="Topology file (PDB or GRO)",
    )
    parser.add_argument(
        "-c", "--checkpoint", default="output_checkpoint.nc",
        help="Checkpoint NetCDF file",
    )
    parser.add_argument(
        "-o", "--output", default="./extracted",
        help="Output directory or file (for --mode frame)",
    )
    parser.add_argument(
        "--pattern", default="{mode}_{i}.xtc",
        help="Output filename pattern with {mode} and {i} placeholders",
    )
    parser.add_argument(
        "-b", "--begin", type=int, default=0,
        help="Start frame (inclusive)",
    )
    parser.add_argument(
        "-e", "--end", type=int, default=None,
        help="End frame (exclusive, default: all)",
    )
    parser.add_argument(
        "-s", "--stride", type=int, default=1,
        help="Frame stride",
    )
    parser.add_argument(
        "--replicas", type=str, default=None,
        help="Comma-separated replica indices, e.g. '0,1,2' (default: all)",
    )
    parser.add_argument(
        "--states", type=str, default=None,
        help="Comma-separated state indices, e.g. '0,1,2' (default: all)",
    )
    parser.add_argument(
        "--frame", type=int, default=None,
        help="Frame index (for --mode frame)",
    )
    parser.add_argument(
        "--replica", type=int, default=None,
        help="Replica index (for --mode frame)",
    )
    parser.add_argument(
        "--state", type=int, default=None,
        help="State index (for --mode frame)",
    )
    parser.add_argument(
        "--block-frames", type=int, default=5000,
        help="Frames per block; tune for available RAM",
    )
    parser.add_argument(
        "--num-readers", type=int, default=8,
        help="Number of parallel HDF5 reader processes",
    )
    parser.add_argument(
        "--num-writers", type=int, default=11,
        help="Number of parallel XTC writer processes",
    )
    parser.add_argument(
        "-np", "--num-workers", type=int, default=None,
        help="Set both --num-readers and --num-writers to this value",
    )
    parser.add_argument(
        "--mode", choices=["replica", "state", "frame"], required=True,
        help="Extraction mode: replica, state, or frame",
    )
    args = parser.parse_args()

    if args.end is not None and args.end == 0:
        args.end = None
    if args.num_workers is not None:
        args.num_readers = args.num_workers
        args.num_writers = args.num_workers

    pattern = args.pattern.replace("{mode}", args.mode)

    if args.mode == "replica":
        replicas = None
        if args.replicas:
            replicas = [int(x) for x in args.replicas.split(",")]
        extract_replicas(
            nc_file=args.netcdf,
            chk_file=args.checkpoint,
            pdb_file=args.topology,
            output_dir=args.output,
            output_pattern=pattern,
            frame_begin=args.begin,
            frame_end=args.end,
            frame_stride=args.stride,
            replicas=replicas,
            block_frames=args.block_frames,
            num_readers=args.num_readers,
            num_writers=args.num_writers,
        )
    elif args.mode == "state":
        states = None
        if args.states:
            states = [int(x) for x in args.states.split(",")]
        extract_states(
            nc_file=args.netcdf,
            chk_file=args.checkpoint,
            pdb_file=args.topology,
            output_dir=args.output,
            output_pattern=pattern,
            frame_begin=args.begin,
            frame_end=args.end,
            frame_stride=args.stride,
            states=states,
            block_frames=args.block_frames,
            num_readers=args.num_readers,
            num_writers=args.num_writers,
        )
    elif args.mode == "frame":
        if args.frame is None:
            print("Error: --frame required for --mode frame")
            return 1
        output = args.output
        if os.path.isdir(output):
            ext = os.path.splitext(pattern)[1] or ".pdb"
            output = os.path.join(output, f"frame_{args.frame}{ext}")
        extract_frame(
            chk_file=args.checkpoint,
            output=output,
            frame_idx=args.frame,
            replica_idx=args.replica,
            state_idx=args.state,
            nc_file=args.netcdf,
        )


if __name__ == "__main__":
    main()
