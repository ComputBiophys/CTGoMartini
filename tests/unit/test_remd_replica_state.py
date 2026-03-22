"""
Unit tests for REMD replica state analysis module.

Tests cover state loading, dt auto-detection, state saving, and occupancy calculation.
"""

import numpy as np
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch

from ctgomartini.analysis.remd_replica_state import (
    load_replica_states,
    save_replica_states,
    compute_state_occupancies,
    _get_dt_from_netcdf,
)


class TestLoadReplicaStates:
    """Tests for load_replica_states function."""

    def test_load_replica_states_from_fixture(self, remd_output_nc):
        """Test loading replica states from actual REMD fixture."""
        states = load_replica_states(remd_output_nc)
        
        # Should return 2D array (n_replicas, n_timesteps)
        assert states.ndim == 2
        assert states.shape[0] > 0  # At least one replica
        assert states.shape[1] > 0  # At least one timestep
        
        # States should be integers (state indices)
        assert np.issubdtype(states.dtype, np.number)
        
        # All state indices should be non-negative
        assert np.all(states >= 0)


class TestGetDtFromNetcdf:
    """Tests for _get_dt_from_netcdf function."""

    def test_auto_detect_dt_from_fixture(self, remd_output_nc):
        """Test auto-detection of dt from actual NetCDF file."""
        dt = _get_dt_from_netcdf(remd_output_nc)
        
        # Should return a positive float
        assert dt is not None
        assert isinstance(dt, float)
        assert dt > 0
        # Typical Martini timestep with 250 steps * 0.02 ps = 5 ps
        assert 0.1 < dt < 100

    def test_auto_detect_dt_failure(self, tmp_path):
        """Test that invalid file returns None."""
        invalid_file = tmp_path / "invalid.nc"
        invalid_file.write_text("not a netcdf file")
        
        dt = _get_dt_from_netcdf(invalid_file)
        assert dt is None


class TestSaveReplicaStates:
    """Tests for save_replica_states function."""

    def test_save_replica_states_auto_dt(self, remd_output_nc, tmp_path):
        """Test saving replica states with auto-detected dt."""
        output_file = tmp_path / "replica_states.dat"
        
        save_replica_states(remd_output_nc, output_data=output_file)
        
        # Check file was created
        assert output_file.exists()
        
        # Check file content
        data = np.loadtxt(output_file, skiprows=1)
        assert data.ndim == 2
        assert data.shape[1] >= 2  # Time column + at least one replica
        
        # First column should be time (starting from 0)
        assert data[0, 0] == 0.0
        
        # Time should be increasing
        assert np.all(np.diff(data[:, 0]) > 0)

    def test_save_replica_states_custom_dt(self, remd_output_nc, tmp_path):
        """Test saving replica states with custom dt."""
        output_file = tmp_path / "replica_states.dat"
        custom_dt = 10.0
        
        save_replica_states(remd_output_nc, output_data=output_file, dt=custom_dt)
        
        assert output_file.exists()
        
        # Load and check time step
        data = np.loadtxt(output_file, skiprows=1)
        time_diff = data[1, 0] - data[0, 0]
        assert time_diff == pytest.approx(custom_dt, rel=0.01)

    def test_save_replica_states_header(self, remd_output_nc, tmp_path):
        """Test that output file has correct header."""
        output_file = tmp_path / "replica_states.dat"
        
        save_replica_states(remd_output_nc, output_data=output_file)
        
        with open(output_file) as f:
            header = f.readline()
        
        # Header should contain Time(ps) and Replica columns
        assert "Time(ps)" in header
        assert "Replica_0" in header


class TestComputeStateOccupancies:
    """Tests for compute_state_occupancies function."""

    def test_compute_state_occupancies(self, remd_output_nc):
        """Test computing state occupancies."""
        occupancies = compute_state_occupancies(remd_output_nc)
        
        # Should return dictionary with replica keys
        assert isinstance(occupancies, dict)
        assert len(occupancies) > 0
        
        # Check each replica has required keys
        for replica_key, stats in occupancies.items():
            assert replica_key.startswith("replica_")
            assert "states" in stats
            assert "counts" in stats
            assert "fractions" in stats
            
            # Fractions should sum to 1
            assert pytest.approx(sum(stats["fractions"]), abs=0.01) == 1.0

    def test_compute_state_occupancies_with_output(self, remd_output_nc, tmp_path):
        """Test computing and saving state occupancies."""
        output_file = tmp_path / "occupancies.dat"
        
        occupancies = compute_state_occupancies(
            remd_output_nc, 
            output_file_path=output_file
        )
        
        # Check file was created
        assert output_file.exists()
        
        # Check file content
        with open(output_file) as f:
            lines = f.readlines()
        
        # Should have header comments and data lines
        header_lines = [l for l in lines if l.startswith("#")]
        data_lines = [l for l in lines if not l.startswith("#")]
        
        assert len(header_lines) > 0
        assert len(data_lines) > 0
