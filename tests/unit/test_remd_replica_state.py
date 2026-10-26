"""
Unit tests for REMD replica state analysis module.

Tests cover state loading, dt auto-detection, state saving, and occupancy
calculation using synthetic built-in REMD output directories (log.csv +
metadata.json), so no OpenMM 8.6 or openmmtools is required.
"""

import json

import numpy as np
import pytest

from ctgomartini.analysis.remd_replica_state import (
    _get_dt_from_metadata,
    compute_state_occupancies,
    load_replica_states,
    save_replica_states,
)


def _make_remd_dir(
    tmp_path,
    n_replicas: int = 3,
    n_frames: int = 4,
    dt: float = 0.02,
    steps_per_iteration: int = 250,
    report_interval: int = 5,
    seed: int = 0,
):
    """Create a synthetic built-in REMD output directory (log.csv + metadata.json)."""
    directory = tmp_path / "output_remd"
    directory.mkdir()

    rng = np.random.default_rng(seed)
    lines = ["Iteration,Step," + ",".join(f"Replica_{i}_State" for i in range(n_replicas))]
    states_written = []
    for frame in range(n_frames):
        iteration = frame * report_interval
        step = iteration * steps_per_iteration
        states = rng.integers(0, n_replicas, size=n_replicas)
        states_written.append(states)
        lines.append(f"{iteration},{step}," + ",".join(str(s) for s in states))
    (directory / "log.csv").write_text("\n".join(lines) + "\n")

    metadata = {
        "dt": dt,
        "exc_freq": steps_per_iteration,
        "stepsPerIteration": steps_per_iteration,
        "reportInterval": report_interval,
        "replica_count": n_replicas,
    }
    (directory / "metadata.json").write_text(json.dumps(metadata))
    return directory, np.array(states_written)


class TestLoadReplicaStates:
    """Tests for load_replica_states function."""

    def test_load_replica_states(self, tmp_path):
        """Test loading replica states from a synthetic REMD directory."""
        directory, states_written = _make_remd_dir(tmp_path)

        states = load_replica_states(directory)

        # Should return 2D array (n_replicas, n_timesteps)
        assert states.ndim == 2
        assert states.shape == (3, 4)

        # States should be integers (state indices)
        assert np.issubdtype(states.dtype, np.integer)

        # All state indices should be non-negative
        assert np.all(states >= 0)

        # Rows are replicas, columns are frames (transpose of log.csv columns)
        np.testing.assert_array_equal(states, states_written.T)

    def test_load_replica_states_rejects_netcdf(self, tmp_path):
        """Test that a non-directory input raises ValueError."""
        nc_file = tmp_path / "output.nc"
        nc_file.write_text("not a directory")

        with pytest.raises(ValueError, match="no longer supported"):
            load_replica_states(nc_file)

    def test_load_replica_states_missing_log(self, tmp_path):
        """Test that a directory without log.csv raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError, match="log.csv"):
            load_replica_states(tmp_path)


class TestGetDtFromMetadata:
    """Tests for _get_dt_from_metadata function."""

    def test_auto_detect_dt(self, tmp_path):
        """Test dt auto-detection: dt * stepsPerIteration * reportInterval."""
        directory, _ = _make_remd_dir(tmp_path)

        dt = _get_dt_from_metadata(directory)

        # 0.02 ps * 250 steps * 5 report interval = 25.0 ps
        assert dt is not None
        assert isinstance(dt, float)
        assert dt == pytest.approx(25.0)

    def test_auto_detect_dt_failure(self, tmp_path):
        """Test that a directory without metadata.json returns None."""
        (tmp_path / "log.csv").write_text("Iteration,Step,Replica_0_State\n0,0,0\n")

        assert _get_dt_from_metadata(tmp_path) is None


class TestSaveReplicaStates:
    """Tests for save_replica_states function."""

    def test_save_replica_states_auto_dt(self, tmp_path):
        """Test saving replica states with auto-detected dt."""
        directory, _ = _make_remd_dir(tmp_path)
        output_file = tmp_path / "replica_states.dat"

        save_replica_states(directory, output_data=output_file)

        # Check file was created
        assert output_file.exists()

        # Check file content
        data = np.loadtxt(output_file, skiprows=1)
        assert data.ndim == 2
        assert data.shape == (4, 4)  # 4 frames, time + 3 replicas

        # First column should be time (starting from 0)
        assert data[0, 0] == 0.0

        # Time should be increasing with the auto-detected step
        assert np.all(np.diff(data[:, 0]) > 0)
        assert data[1, 0] - data[0, 0] == pytest.approx(25.0)

    def test_save_replica_states_custom_dt(self, tmp_path):
        """Test saving replica states with custom dt."""
        directory, _ = _make_remd_dir(tmp_path)
        output_file = tmp_path / "replica_states.dat"
        custom_dt = 10.0

        save_replica_states(directory, output_data=output_file, dt=custom_dt)

        assert output_file.exists()

        # Load and check time step
        data = np.loadtxt(output_file, skiprows=1)
        time_diff = data[1, 0] - data[0, 0]
        assert time_diff == pytest.approx(custom_dt, rel=0.01)

    def test_save_replica_states_header(self, tmp_path):
        """Test that output file has correct header."""
        directory, _ = _make_remd_dir(tmp_path)
        output_file = tmp_path / "replica_states.dat"

        save_replica_states(directory, output_data=output_file)

        with open(output_file) as f:
            header = f.readline()

        # Header should contain Time(ps) and Replica columns
        assert "Time(ps)" in header
        assert "Replica_0" in header


class TestComputeStateOccupancies:
    """Tests for compute_state_occupancies function."""

    def test_compute_state_occupancies(self, tmp_path):
        """Test computing state occupancies."""
        directory, _ = _make_remd_dir(tmp_path)

        occupancies = compute_state_occupancies(directory)

        # Should return dictionary with replica keys
        assert isinstance(occupancies, dict)
        assert len(occupancies) == 3

        # Check each replica has required keys
        for replica_key, stats in occupancies.items():
            assert replica_key.startswith("replica_")
            assert "states" in stats
            assert "counts" in stats
            assert "fractions" in stats

            # Fractions should sum to 1
            assert pytest.approx(sum(stats["fractions"]), abs=0.01) == 1.0

    def test_compute_state_occupancies_with_output(self, tmp_path):
        """Test computing and saving state occupancies."""
        directory, _ = _make_remd_dir(tmp_path)
        output_file = tmp_path / "occupancies.dat"

        compute_state_occupancies(directory, output_file_path=output_file)

        # Check file was created
        assert output_file.exists()

        # Check file content
        with open(output_file) as f:
            lines = f.readlines()

        # Should have header comments and data lines
        header_lines = [ln for ln in lines if ln.startswith("#")]
        data_lines = [ln for ln in lines if not ln.startswith("#")]

        assert len(header_lines) > 0
        assert len(data_lines) > 0
