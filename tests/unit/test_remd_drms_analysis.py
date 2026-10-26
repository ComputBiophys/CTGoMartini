"""
Unit tests for REMD dRMS analysis module.

Tests cover reference distance calculation (MDAnalysis fixtures) and dRMS
computation over per-replica XTC trajectories of synthetic built-in REMD
output directories (replica_<i>.xtc + metadata.json), so no OpenMM 8.6 or
openmmtools is required.
"""

import json
from pathlib import Path

import numpy as np
import pytest

# Fixture paths
MD_RUN_DIR = Path(__file__).parent.parent / "fixtures" / "md_run" / "EXP" / "template"
TOP_FILE = MD_RUN_DIR / "GBP_open_cg.pdb"
# insane-generated full system GRO: 4435 atoms (protein = first 732, matching
# GBP_*_cg.pdb index order), 9-value triclinic box line, trailing blank line.
IONS_GRO = MD_RUN_DIR / "ions.gro"


def _make_remd_dir(
    tmp_path,
    n_replicas: int = 2,
    n_frames: int = 6,
    dt: float = 0.02,
    steps_per_iteration: int = 250,
    report_interval: int = 5,
    seed: int = 1,
):
    """Create a synthetic built-in REMD output directory with per-replica XTCs.

    Coordinates are random-walk-like noise around the open-state structure so
    dRMS values are finite and non-negative.
    """
    import MDAnalysis as mda  # noqa: N813

    directory = tmp_path / "output_remd"
    directory.mkdir()

    universe = mda.Universe(str(TOP_FILE))
    rng = np.random.default_rng(seed)
    for replica in range(n_replicas):
        with mda.Writer(
            str(directory / f"replica_{replica}.xtc"), universe.atoms.n_atoms
        ) as writer:
            for _ in range(n_frames):
                universe.atoms.positions = rng.normal(
                    loc=20.0, scale=2.0, size=(universe.atoms.n_atoms, 3)
                )
                writer.write(universe.atoms)

    metadata = {
        "dt": dt,
        "exc_freq": steps_per_iteration,
        "stepsPerIteration": steps_per_iteration,
        "reportInterval": report_interval,
        "replica_count": n_replicas,
    }
    (directory / "metadata.json").write_text(json.dumps(metadata))
    return directory


@pytest.fixture
def ref_files():
    """Provide reference file paths."""
    return [
        str(MD_RUN_DIR / "GBP_open_cg.pdb"),
        str(MD_RUN_DIR / "GBP_closed_cg.pdb"),
    ]


class TestCalculateReferenceDistances:
    """Tests for calculate_reference_distances function."""

    def test_calculate_reference_distances_basic(self, ref_files):
        """Test basic reference distance calculation."""
        from ctgomartini.analysis.remd_drms_analysis import calculate_reference_distances

        ref_distances = calculate_reference_distances(ref_files)

        # Should return list of 2 arrays (one per state)
        assert len(ref_distances) == 2

        # Each should be a numpy array with shape (n_pairs, 3)
        for dist_array in ref_distances:
            assert isinstance(dist_array, np.ndarray)
            assert dist_array.ndim == 2
            assert dist_array.shape[1] == 3  # [atom_i, atom_j, distance]

    def test_calculate_reference_distances_with_custom_params(self, ref_files):
        """Test reference distance calculation with custom parameters."""
        from ctgomartini.analysis.remd_drms_analysis import calculate_reference_distances

        ref_distances = calculate_reference_distances(
            ref_files,
            selected_atom="name BB",
            min_dist=5.0,
            max_dist=40.0,
            min_diff=3.0,
            excl_res=2
        )

        assert len(ref_distances) == 2

        # With less strict parameters, should get more pairs
        ref_distances_strict = calculate_reference_distances(
            ref_files,
            min_dist=6.0,
            max_dist=50.0,
            min_diff=5.0,
            excl_res=4
        )

        # Strict params should give fewer or equal pairs
        for i in range(2):
            assert len(ref_distances[i]) >= len(ref_distances_strict[i])

    def test_calculate_reference_distances_different_atoms(self, ref_files):
        """Test with different atom selections."""
        from ctgomartini.analysis.remd_drms_analysis import calculate_reference_distances

        # BB only
        ref_bb = calculate_reference_distances(ref_files, selected_atom="name BB")

        # BB and SC1
        ref_bb_sc = calculate_reference_distances(
            ref_files,
            selected_atom="name BB or name SC1"
        )

        # More atoms should give more pairs
        for i in range(2):
            assert len(ref_bb_sc[i]) > len(ref_bb[i])

    def test_invalid_reference_count(self):
        """Test that single reference raises error."""
        from ctgomartini.analysis.remd_drms_analysis import calculate_reference_distances

        single_ref = [str(MD_RUN_DIR / "GBP_open_cg.pdb")]

        # Should raise ValueError when trying to calculate min() of empty iterable
        with pytest.raises((ValueError, IndexError)):
            calculate_reference_distances(single_ref)


class TestCalculateTrajectoryDrms:
    """Tests for calculate_trajectory_drms function on synthetic XTC output."""

    def test_calculate_trajectory_drms_basic(self, tmp_path, ref_files):
        """Test basic dRMS calculation."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory = _make_remd_dir(tmp_path)
        ref_distances = calculate_reference_distances(ref_files)[0]  # State A

        times, drms = calculate_trajectory_drms(
            directory,
            str(TOP_FILE),
            ref_distances,
            replica_indices=np.array([0]),
            skip=1,
            num_workers=1,
        )

        assert isinstance(times, np.ndarray)
        assert isinstance(drms, np.ndarray)
        assert times.shape[0] == drms.shape[0] == 6  # 6 frames
        assert drms.shape[1] == 1  # One replica
        assert np.all(drms >= 0)  # dRMS should be non-negative
        # dt frame = 0.02 * 250 * 5 = 25 ps
        assert times[1] - times[0] == pytest.approx(25.0)

    def test_calculate_trajectory_drms_multiple_replicas(self, tmp_path, ref_files):
        """Test dRMS calculation for multiple replicas."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory = _make_remd_dir(tmp_path, n_replicas=3)
        ref_distances = calculate_reference_distances(ref_files)[0]

        times, drms = calculate_trajectory_drms(
            directory,
            str(TOP_FILE),
            ref_distances,
            replica_indices=np.array([0, 1, 2]),
            skip=1,
            num_workers=2,
        )

        assert drms.shape == (6, 3)  # 6 frames, three replicas

    def test_calculate_trajectory_drms_with_skip(self, tmp_path, ref_files):
        """Test dRMS calculation with frame skipping."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory = _make_remd_dir(tmp_path)
        ref_distances = calculate_reference_distances(ref_files)[0]

        times1, drms1 = calculate_trajectory_drms(
            directory,
            str(TOP_FILE),
            ref_distances,
            replica_indices=np.array([0]),
            skip=1,
            num_workers=1,
        )

        times2, drms2 = calculate_trajectory_drms(
            directory,
            str(TOP_FILE),
            ref_distances,
            replica_indices=np.array([0]),
            skip=2,
            num_workers=1,
        )

        # Skip=2 should give about half the frames
        assert len(times2) <= len(times1)
        assert len(times2) == 3

    def test_rejects_netcdf_input(self, tmp_path, ref_files):
        """Test that a non-directory input raises ValueError."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        nc_file = tmp_path / "output.nc"
        nc_file.write_text("not a directory")
        ref_distances = calculate_reference_distances(ref_files)[0]

        with pytest.raises(ValueError, match="no longer supported"):
            calculate_trajectory_drms(
                nc_file, str(TOP_FILE), ref_distances, num_workers=1
            )

    def test_missing_xtc_raises(self, tmp_path, ref_files):
        """Test that a directory without replica XTCs raises FileNotFoundError."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory = tmp_path / "output_remd"
        directory.mkdir()
        (directory / "metadata.json").write_text(json.dumps({"replica_count": 2}))
        ref_distances = calculate_reference_distances(ref_files)[0]

        with pytest.raises(FileNotFoundError, match="remd_xtc_output"):
            calculate_trajectory_drms(
                directory, str(TOP_FILE), ref_distances, num_workers=1
            )


class TestIntegrationWithSyntheticOutput:
    """Integration tests on synthetic built-in REMD output directories."""

    def test_full_analysis_pipeline(self, tmp_path, ref_files):
        """Test full analysis pipeline: reference -> dRMS -> save."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
            save_drms_results,
        )

        directory = _make_remd_dir(tmp_path)

        # Calculate reference distances for both states
        all_ref_distances = calculate_reference_distances(ref_files)

        # Process each state
        for state_idx, ref_dist in enumerate(all_ref_distances):
            times, drms = calculate_trajectory_drms(
                directory,
                str(TOP_FILE),
                ref_dist,
                replica_indices=np.array([0, 1]),
                skip=1,
                num_workers=1,
            )

            # Save results
            output_file = str(tmp_path / f"drms_state{state_idx}.dat")
            save_drms_results(
                times,
                drms,
                output_file,
                output_dir=directory,
                replica_indices=np.array([0, 1]),
            )

            # Verify file was created
            assert Path(output_file).exists()

            # Read and verify content
            data = np.loadtxt(output_file, comments='#')
            assert data.shape[0] == len(times)
            assert data.shape[1] == 3  # time + 2 replicas

    def test_different_replica_selections(self, tmp_path, ref_files):
        """Test analysis with different replica selections."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory = _make_remd_dir(tmp_path)
        ref_distances = calculate_reference_distances(ref_files)[0]

        # Test with single replica
        times1, drms1 = calculate_trajectory_drms(
            directory,
            str(TOP_FILE),
            ref_distances,
            replica_indices=np.array([0]),
            skip=1,
            num_workers=1,
        )
        assert drms1.shape[1] == 1

        # Test with all replicas (None reads replica_count from metadata.json)
        times_all, drms_all = calculate_trajectory_drms(
            directory,
            str(TOP_FILE),
            ref_distances,
            replica_indices=None,  # All replicas
            skip=1,
            num_workers=1,
        )
        assert drms_all.shape[1] == 2


class TestFullSystemTopologyFormats:
    """dRMS on a realistic full-system XTC with GRO vs PDB topology.

    Uses the insane-generated ions.gro fixture (4435 atoms, 9-value triclinic
    box line, trailing blank line) whose format mdtraj's GRO parser cannot
    read. MDAnalysis must read it, and a PDB topology with the same atom
    order must give identical dRMS values.
    """

    N_FRAMES = 4

    def _make_full_system_remd_dir(self, tmp_path):
        """Write replica_0.xtc for the full ions.gro system + metadata.json.

        Returns the output directory and a PDB topology with the same atom
        order as ions.gro.
        """
        import MDAnalysis as mda  # noqa: N813

        directory = tmp_path / "output_remd"
        directory.mkdir()

        universe = mda.Universe(str(IONS_GRO))
        rng = np.random.default_rng(7)
        with mda.Writer(
            str(directory / "replica_0.xtc"), universe.atoms.n_atoms
        ) as writer:
            for _ in range(self.N_FRAMES):
                universe.atoms.positions += rng.normal(
                    0.0, 1.0, size=(universe.atoms.n_atoms, 3)
                )
                writer.write(universe.atoms)

        metadata = {
            "dt": 0.02,
            "exc_freq": 250,
            "stepsPerIteration": 250,
            "reportInterval": 5,
            "replica_count": 1,
        }
        (directory / "metadata.json").write_text(json.dumps(metadata))

        pdb_top = tmp_path / "system_full.pdb"
        universe.atoms.write(str(pdb_top))
        return directory, str(pdb_top)

    def test_gro_topology_accepted(self, tmp_path, ref_files):
        """The insane-format GRO works as --top (regression for mdtraj failure)."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory, _ = self._make_full_system_remd_dir(tmp_path)
        ref_distances = calculate_reference_distances(ref_files)[0]

        times, drms = calculate_trajectory_drms(
            directory,
            str(IONS_GRO),
            ref_distances,
            replica_indices=np.array([0]),
            num_workers=1,
        )

        assert times.shape[0] == self.N_FRAMES
        assert drms.shape == (self.N_FRAMES, 1)
        assert np.all(drms >= 0)

    def test_gro_and_pdb_topologies_give_identical_drms(self, tmp_path, ref_files):
        """GRO and PDB topologies of the same system give identical dRMS."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_reference_distances,
            calculate_trajectory_drms,
        )

        directory, pdb_top = self._make_full_system_remd_dir(tmp_path)
        ref_distances = calculate_reference_distances(ref_files)[0]

        times_gro, drms_gro = calculate_trajectory_drms(
            directory,
            str(IONS_GRO),
            ref_distances,
            replica_indices=np.array([0]),
            num_workers=1,
        )
        times_pdb, drms_pdb = calculate_trajectory_drms(
            directory,
            pdb_top,
            ref_distances,
            replica_indices=np.array([0]),
            num_workers=1,
        )

        np.testing.assert_allclose(times_gro, times_pdb)
        np.testing.assert_allclose(drms_gro, drms_pdb, rtol=1e-6, atol=1e-8)
