"""
Unit tests for REMD MBAR FES analysis module.

Tests cover FES initialization, single-state analysis, mixing parameter analysis,
and parameter sweeps.
"""

import numpy as np
import pytest
from pathlib import Path

from ctgomartini.analysis.remd_mbar import FESAnalyzer


class TestFESAnalyzerInitialization:
    """Tests for FESAnalyzer initialization."""

    def test_init_with_fixture(self, remd_output_nc, tmp_path):
        """Test initialization with actual REMD fixture."""
        # Create a mock CV file
        cv_file = tmp_path / "cv.dat"
        # Create CV data: time + 11 replicas x 4 frames
        cv_data = np.random.rand(4, 12)
        cv_data[:, 0] = np.arange(4) * 5.0  # time column
        np.savetxt(cv_file, cv_data)
        
        analyzer = FESAnalyzer(remd_output_nc, cv_file, interval=5)
        
        assert analyzer.n_states == 11
        assert len(analyzer.temperatures_k) == 11
        assert analyzer.cv_values_replica.shape[0] == 11

    def test_init_file_not_found(self, tmp_path):
        """Test that FileNotFoundError is raised for missing files."""
        with pytest.raises(FileNotFoundError):
            FESAnalyzer(tmp_path / "nonexistent.nc", tmp_path / "cv.dat")


class TestFESAnalyzerFunctionality:
    """Tests for FESAnalyzer main functionality with actual data."""

    def test_initialize_fes(self, remd_output_nc, tmp_path):
        """Test FES initialization."""
        # Create mock CV file with sufficient frames
        cv_file = tmp_path / "cv.dat"
        n_frames = 100  # Need more frames for MBAR
        cv_data = np.random.rand(n_frames, 12)
        cv_data[:, 0] = np.arange(n_frames) * 5.0
        np.savetxt(cv_file, cv_data)
        
        analyzer = FESAnalyzer(remd_output_nc, cv_file, interval=5)
        analyzer.initialize_fes(g=1.0, start_ratio=0.0)
        
        assert analyzer.fes is not None
        assert analyzer.subsampled_indices is not None
        assert analyzer.u_kn is not None
        assert analyzer.concatenated_cv is not None

    def test_analyze_invalid_state(self, remd_output_nc, tmp_path):
        """Test that invalid state raises ValueError."""
        cv_file = tmp_path / "cv.dat"
        cv_data = np.random.rand(10, 12)
        cv_data[:, 0] = np.arange(10) * 5.0
        np.savetxt(cv_file, cv_data)
        
        analyzer = FESAnalyzer(remd_output_nc, cv_file, interval=5)
        analyzer.initialize_fes(g=1.0)
        
        with pytest.raises(ValueError):
            analyzer.analyze_onestate(selected_state=20)  # Invalid state


class TestMixingMethods:
    """Tests for EXP and HAM mixing methods."""

    def test_invalid_mixing_method(self, remd_output_nc, tmp_path):
        """Test that invalid mixing method raises ValueError."""
        cv_file = tmp_path / "cv.dat"
        cv_data = np.random.rand(10, 12)
        cv_data[:, 0] = np.arange(10) * 5.0
        np.savetxt(cv_file, cv_data)
        
        analyzer = FESAnalyzer(remd_output_nc, cv_file, interval=5)
        analyzer.initialize_fes(g=1.0)
        
        with pytest.raises(ValueError):
            analyzer.analyze_one_mixing_parameters(
                {}, temperature=310, method="INVALID"
            )


class TestParameterSweep:
    """Tests for parameter sweep functionality."""
    pass  # Parameter sweep tests require large datasets, tested manually


class TestSaveLoadResults:
    """Tests for save/load functionality."""

    def test_save_and_load_results(self, tmp_path):
        """Test saving and loading analysis results."""
        # Create mock results
        results = {
            "cv_values": np.array([1.0, 2.0, 3.0]),
            "pmf": np.array([0.0, 1.0, 2.0]),
            "pmf_uncertainty": np.array([0.1, 0.1, 0.1]),
            "metrics": {"barrier": 5.0, "keq": 1.0}
        }
        
        # Save results
        save_file = tmp_path / "results.pkl"
        FESAnalyzer.save_results(results, save_file)
        assert save_file.exists()
        
        # Load results
        loaded = FESAnalyzer.load_results(save_file)
        assert "cv_values" in loaded
        assert "pmf" in loaded
        assert "metrics" in loaded
        
        # Check data consistency
        np.testing.assert_array_almost_equal(loaded["cv_values"], results["cv_values"])


class TestCalculateMetrics:
    """Tests for metric calculation methods."""

    def test_calculate_barrier(self):
        """Test barrier calculation."""
        cv_values = np.linspace(0, 10, 100)
        pmf = -(cv_values - 5) ** 2  # Parabolic with max at 5
        pmf = pmf - pmf.min()  # Normalize
        
        barrier_pos, barrier_height = FESAnalyzer.calculate_barrier(
            cv_values, pmf, left_bound=3, right_bound=7
        )
        
        assert 4.5 < barrier_pos < 5.5
        assert barrier_height > 0

    def test_calculate_basins(self):
        """Test basin calculation."""
        cv_values = np.linspace(0, 10, 100)
        # Double well potential
        pmf = np.where(cv_values < 5, (cv_values - 2) ** 2, (cv_values - 8) ** 2)
        pmf = pmf - pmf.min()
        
        b1_pos, b2_pos, b1_pmf, b2_pmf = FESAnalyzer.calculate_basins(
            cv_values, pmf, barrier_position=5
        )
        
        assert b1_pos < 5
        assert b2_pos > 5
        assert b1_pmf >= 0
        assert b2_pmf >= 0

    def test_calculate_equilibrium_constant(self):
        """Test equilibrium constant calculation."""
        cv_values = np.linspace(0, 10, 100)
        pmf = np.where(cv_values < 5, (cv_values - 2) ** 2, (cv_values - 8) ** 2)
        pmf = pmf - pmf.min()
        
        keq = FESAnalyzer.calculate_equilibrium_constant(
            cv_values, pmf, barrier_position=5, temperature=310
        )
        
        assert keq > 0
        assert np.isfinite(keq)
