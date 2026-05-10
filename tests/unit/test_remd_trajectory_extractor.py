"""
Unit tests for REMD trajectory extractor.

Tests cover the functional API: extract_replicas, extract_states, extract_frame,
and CLI argument parsing.
"""

import os
from pathlib import Path
from unittest.mock import patch

from ctgomartini.analysis.remd_trajectory_extractor import (
    extract_replicas,
    extract_states,
    extract_frame,
    main,
)


class TestExtractReplicas:
    """Tests for extract_replicas function using fixture data."""

    def test_extract_replicas_xtc(
        self,
        remd_netcdf_file: Path,
        remd_topology_file: Path,
        remd_checkpoint_file: Path,
        tmp_path: Path,
    ):
        """Test extracting replica trajectories to XTC format."""
        output_dir = tmp_path / "replicas"
        extract_replicas(
            nc_file=str(remd_netcdf_file),
            chk_file=str(remd_checkpoint_file),
            pdb_file=str(remd_topology_file),
            output_dir=str(output_dir),
            output_pattern="replica_{i}.xtc",
            replicas=[0, 1],
            num_readers=1,
            num_writers=1,
        )
        assert (output_dir / "replica_0.xtc").exists()
        assert (output_dir / "replica_1.xtc").exists()

    def test_extract_replicas_with_stride(
        self,
        remd_netcdf_file: Path,
        remd_topology_file: Path,
        remd_checkpoint_file: Path,
        tmp_path: Path,
    ):
        """Test extracting replica trajectories with frame stride."""
        output_dir = tmp_path / "replicas_stride"
        extract_replicas(
            nc_file=str(remd_netcdf_file),
            chk_file=str(remd_checkpoint_file),
            pdb_file=str(remd_topology_file),
            output_dir=str(output_dir),
            output_pattern="replica_{i}.xtc",
            replicas=[0],
            frame_stride=2,
            num_readers=1,
            num_writers=1,
        )
        assert (output_dir / "replica_0.xtc").exists()


class TestExtractStates:
    """Tests for extract_states function using fixture data."""

    def test_extract_states(
        self,
        remd_netcdf_file: Path,
        remd_topology_file: Path,
        remd_checkpoint_file: Path,
        tmp_path: Path,
    ):
        """Test extracting state trajectories to XTC format."""
        output_dir = tmp_path / "states"
        extract_states(
            nc_file=str(remd_netcdf_file),
            chk_file=str(remd_checkpoint_file),
            pdb_file=str(remd_topology_file),
            output_dir=str(output_dir),
            output_pattern="state_{i}.xtc",
            states=[0, 1],
            num_readers=1,
            num_writers=1,
        )
        assert (output_dir / "state_0.xtc").exists()
        assert (output_dir / "state_1.xtc").exists()


class TestExtractFrame:
    """Tests for extract_frame function."""

    def test_extract_frame_pdb(
        self,
        remd_checkpoint_file: Path,
        tmp_path: Path,
    ):
        """Test extracting a single frame to PDB format."""
        output = str(tmp_path / "frame.pdb")
        extract_frame(
            chk_file=str(remd_checkpoint_file),
            output=output,
            frame_idx=0,
            replica_idx=0,
        )
        assert os.path.exists(output)

    def test_extract_frame_gro(
        self,
        remd_checkpoint_file: Path,
        tmp_path: Path,
    ):
        """Test extracting a single frame to GRO format."""
        output = str(tmp_path / "frame.gro")
        extract_frame(
            chk_file=str(remd_checkpoint_file),
            output=output,
            frame_idx=0,
            replica_idx=0,
        )
        assert os.path.exists(output)

    def test_extract_frame_xtc(
        self,
        remd_checkpoint_file: Path,
        tmp_path: Path,
    ):
        """Test extracting a single frame to XTC format."""
        output = str(tmp_path / "frame.xtc")
        extract_frame(
            chk_file=str(remd_checkpoint_file),
            output=output,
            frame_idx=0,
            replica_idx=0,
        )
        assert os.path.exists(output)


class TestCLIArguments:
    """Tests for CLI argument parsing."""

    def test_cli_replica_mode(self):
        """Test CLI arg parsing for replica mode dispatches correctly."""
        test_args = [
            "-nc", "output.nc",
            "-p", "topology.pdb",
            "--mode", "replica",
            "-o", "./output",
            "--replicas", "0,1,2",
        ]
        with patch("sys.argv", ["remd_trajectory_extractor"] + test_args):
            with patch("ctgomartini.analysis.remd_trajectory_extractor.extract_replicas") as mock_extract:
                main()
                mock_extract.assert_called_once()

    def test_cli_state_mode(self):
        """Test CLI arg parsing for state mode dispatches correctly."""
        test_args = [
            "-nc", "output.nc",
            "-p", "topology.pdb",
            "--mode", "state",
            "-o", "./output",
            "--states", "0,1",
        ]
        with patch("sys.argv", ["remd_trajectory_extractor"] + test_args):
            with patch("ctgomartini.analysis.remd_trajectory_extractor.extract_states") as mock_extract:
                main()
                mock_extract.assert_called_once()

    def test_cli_frame_mode_requires_frame(self):
        """Test that frame mode requires --frame argument."""
        test_args = [
            "-nc", "output.nc",
            "-p", "topology.pdb",
            "--mode", "frame",
            "-o", "frame.pdb",
        ]
        with patch("sys.argv", ["remd_trajectory_extractor"] + test_args):
            result = main()
            assert result == 1

    def test_cli_num_workers_sets_both(self):
        """Test that -np sets both num_readers and num_writers."""
        test_args = [
            "-nc", "output.nc",
            "-p", "topology.pdb",
            "--mode", "replica",
            "-o", "./output",
            "-np", "4",
        ]
        with patch("sys.argv", ["remd_trajectory_extractor"] + test_args):
            with patch("ctgomartini.analysis.remd_trajectory_extractor.extract_replicas") as mock_extract:
                main()
                call_kwargs = mock_extract.call_args.kwargs
                assert call_kwargs["num_readers"] == 4
                assert call_kwargs["num_writers"] == 4
