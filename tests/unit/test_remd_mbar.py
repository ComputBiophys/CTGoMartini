"""
Unit tests for REMD MBAR FES analysis module.

Tests cover FESAnalyzer initialization on synthetic built-in REMD output
directories (metadata.json + energy.csv + log.csv + cv_energies.csv), mixing
method validation, and the pure save/load/metrics helpers. Tests that need
pymbar use importorskip so they run only where pymbar is installed.
"""

import json

import numpy as np
import pytest

from ctgomartini.analysis.remd_mbar import FESAnalyzer

TEMPERATURES = [300.0, 310.0]


def _make_openmm_dir(tmp_path, n_frames: int = 100, seed: int = 0):
    """Create a synthetic OpenMM >= 8.6 built-in REMD output directory.

    Reduced energies are built from a shared base plus small state offsets so
    the two states overlap strongly (MBAR converges on the synthetic data).
    Returns (directory, cv_file).
    """
    directory = tmp_path / "output_remd"
    directory.mkdir()
    n_states = len(TEMPERATURES)

    metadata = {
        "backend": "openmm",
        "method": "exp",
        "molname": "gbp",
        "replica_count": n_states,
        "n_basins": 2,
        "dt": 0.02,
        "exc_freq": 250,
        "stepsPerIteration": 250,
        "reportInterval": 1,
        "states": [
            {"temperature": t, "ctg_gbp_beta": 1.0, "ctg_gbp_C1": 0.0, "ctg_gbp_C2": 0.0}
            for t in TEMPERATURES
        ],
    }
    (directory / "metadata.json").write_text(json.dumps(metadata))

    rng = np.random.default_rng(seed)
    base = rng.normal(0.0, 1.0, size=(n_frames, n_states))
    # energy.csv rows: [replica][state] flattened per frame
    energy = np.empty((n_frames, n_states * n_states))
    for state in range(n_states):
        energy[:, state::n_states] = base + 0.05 * state + rng.normal(
            0.0, 0.01, size=(n_frames, n_states)
        )
    np.savetxt(directory / "energy.csv", energy, delimiter=",")

    # log.csv: replicas stay in their own states
    log_lines = ["Iteration,Step,Replica_0_State,Replica_1_State"]
    for frame in range(n_frames):
        log_lines.append(f"{frame},{frame * 250},0,1")
    (directory / "log.csv").write_text("\n".join(log_lines) + "\n")

    # cv_energies.csv: one row per replica per frame
    cv_lines = ["Iteration,Replica,Molecule,E1,E2"]
    for frame in range(n_frames):
        for replica in range(n_states):
            e1 = 1000.0 + replica + 0.01 * frame
            e2 = 2000.0 + replica + 0.01 * frame
            cv_lines.append(f"{frame},{replica},gbp,{e1},{e2}")
    (directory / "cv_energies.csv").write_text("\n".join(cv_lines) + "\n")

    # External CV file: time + one column per state
    cv_data = np.empty((n_frames, n_states + 1))
    cv_data[:, 0] = np.arange(n_frames) * 5.0
    cv_data[:, 1:] = rng.normal(5.0, 1.0, size=(n_frames, n_states))
    cv_file = tmp_path / "cv.dat"
    np.savetxt(cv_file, cv_data)

    return directory, cv_file


class TestFESAnalyzerInitialization:
    """Tests for FESAnalyzer initialization."""

    def test_init_with_openmm_directory(self, tmp_path):
        """Test initialization with a synthetic built-in REMD directory."""
        directory, cv_file = _make_openmm_dir(tmp_path, n_frames=4)

        analyzer = FESAnalyzer(directory, cv_file)

        assert analyzer.n_states == 2
        assert len(analyzer.temperatures_k) == 2
        assert analyzer.temperatures_k.tolist() == TEMPERATURES
        assert analyzer.cv_values_replica.shape[0] == 2
        assert analyzer.replica_energies.shape == (2, 2, 4)
        assert analyzer.unsampled_energies.shape == (2, 2, 4)
        assert analyzer._basin_energies_kj.shape == (2, 2, 4)

    def test_init_file_not_found(self, tmp_path):
        """Test that FileNotFoundError is raised for missing files."""
        with pytest.raises(FileNotFoundError):
            FESAnalyzer(tmp_path / "nonexistent_remd", tmp_path / "cv.dat")

    def test_init_rejects_netcdf(self, tmp_path):
        """Test that a non-directory (.nc) input raises a clear error."""
        nc_file = tmp_path / "output.nc"
        nc_file.write_text("not a directory")
        cv_file = tmp_path / "cv.dat"
        np.savetxt(cv_file, np.zeros((4, 3)))

        with pytest.raises(ValueError, match="no longer supported"):
            FESAnalyzer(nc_file, cv_file)


class TestFESAnalyzerFunctionality:
    """Tests for FESAnalyzer main functionality with synthetic data."""

    def test_initialize_fes(self, tmp_path):
        """Test FES initialization (requires pymbar)."""
        pytest.importorskip("pymbar")
        directory, cv_file = _make_openmm_dir(tmp_path, n_frames=100)

        analyzer = FESAnalyzer(directory, cv_file)
        analyzer.initialize_fes(g=1.0, start_ratio=0.0)

        assert analyzer.fes is not None
        assert analyzer.subsampled_indices is not None
        assert analyzer.u_kn is not None
        assert analyzer.concatenated_cv is not None

    def test_analyze_invalid_state(self, tmp_path):
        """Test that invalid state raises ValueError."""
        directory, cv_file = _make_openmm_dir(tmp_path, n_frames=10)

        analyzer = FESAnalyzer(directory, cv_file)

        with pytest.raises(ValueError):
            analyzer.analyze_onestate(selected_state=20)  # Invalid state


class TestMixingMethods:
    """Tests for EXP and HAM mixing methods."""

    def test_invalid_mixing_method(self, tmp_path):
        """Test that invalid mixing method raises ValueError."""
        directory, cv_file = _make_openmm_dir(tmp_path, n_frames=10)

        analyzer = FESAnalyzer(directory, cv_file)

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
