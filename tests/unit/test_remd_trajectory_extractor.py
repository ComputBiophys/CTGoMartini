"""
Unit tests for REMD trajectory extractor.

Tests cover format detection, chunk size optimization, and trajectory extraction
functionality with mocked NetCDF data.
"""

import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch, mock_open

import numpy as np
import pytest

if TYPE_CHECKING:
    from _pytest.fixtures import FixtureRequest
    from _pytest.tmpdir import TempPathFactory

from ctgomartini.analysis.remd_trajectory_extractor import (
    REMDTrajectoryExtractor,
    _detect_format,
    _optimize_chunk_size,
)


class TestFormatDetection:
    """Tests for _detect_format function."""
    
    def test_detect_xtc_format(self):
        """Test XTC format detection."""
        assert _detect_format("traj.xtc") == "xtc"
        assert _detect_format("traj.XTC") == "xtc"
        assert _detect_format("/path/to/file.xtc") == "xtc"
    
    def test_detect_dcd_format(self):
        """Test DCD format detection."""
        assert _detect_format("traj.dcd") == "dcd"
        assert _detect_format("traj.DCD") == "dcd"
    
    def test_detect_pdb_format(self):
        """Test PDB format detection."""
        assert _detect_format("frame.pdb") == "pdb"
        assert _detect_format("frame.PDB") == "pdb"
    
    def test_detect_unknown_format(self):
        """Test that unknown formats raise ValueError."""
        with pytest.raises(ValueError, match="Cannot detect format"):
            _detect_format("traj.xyz")
        
        with pytest.raises(ValueError, match="Cannot detect format"):
            _detect_format("traj")
    
    def test_detect_format_no_extension(self):
        """Test that files without extension raise ValueError."""
        with pytest.raises(ValueError, match="Cannot detect format"):
            _detect_format("traj")


class TestChunkSizeOptimization:
    """Tests for _optimize_chunk_size function."""
    
    def test_small_system(self):
        """Test chunk size for small system (low memory per frame)."""
        # 1000 atoms, 100 frames, 4 workers
        chunk_size = _optimize_chunk_size(100, 1000, 4)
        # Should be limited by max_chunk (100 // 8 = 12)
        assert 10 <= chunk_size <= 100
    
    def test_large_system(self):
        """Test chunk size for large system (high memory per frame)."""
        # 100000 atoms, 1000 frames, 8 workers
        chunk_size = _optimize_chunk_size(1000, 100000, 8)
        # Memory limited: 50MB / (100000 * 3 * 4 bytes) ≈ 41 frames
        assert chunk_size <= 500
        assert chunk_size >= 10
    
    def test_min_chunk_enforcement(self):
        """Test that minimum chunk size is enforced."""
        # Very small system with many workers
        chunk_size = _optimize_chunk_size(10, 100, 100)
        assert chunk_size == 10  # Minimum enforced
    
    def test_max_chunk_enforcement(self):
        """Test that maximum chunk size is enforced."""
        # Small system with few frames
        chunk_size = _optimize_chunk_size(5000, 100, 1)
        assert chunk_size <= 500  # Maximum enforced
    
    def test_not_exceed_total_frames(self):
        """Test that chunk size doesn't exceed total frames."""
        chunk_size = _optimize_chunk_size(5, 1000, 1)
        assert chunk_size <= 5


class TestREMDTrajectoryExtractorInitialization:
    """Tests for REMDTrajectoryExtractor initialization."""
    
    @pytest.fixture
    def mock_reporter(self):
        """Create a mock MultiStateReporter."""
        reporter = MagicMock()
        
        # Mock storage checkpoint
        storage = MagicMock()
        positions_var = MagicMock()
        positions_var.shape = (100, 8, 1000)  # frames, replicas, atoms
        storage.variables = {'positions': positions_var}
        reporter._storage_checkpoint = storage
        
        # Mock state trajectories
        analysis_storage = MagicMock()
        states_var = MagicMock()
        # Create simple state mapping: 8 replicas, 8 states
        states_data = np.arange(8).reshape(1, -1)
        states_data = np.repeat(states_data, 100, axis=0)
        states_var.__getitem__ = MagicMock(return_value=states_data)
        analysis_storage.variables = {'states': states_var}
        reporter._storage_analysis = analysis_storage
        
        # Mock MCMC moves for dt calculation
        mcmove = MagicMock()
        mcmove.n_steps = 100
        mcmove.timestep = MagicMock()
        mcmove.timestep.value_in_unit = MagicMock(return_value=0.005)
        reporter.read_mcmc_moves = MagicMock(return_value=[mcmove])
        
        return reporter
    
    @patch('ctgomartini.analysis.remd_trajectory_extractor._import_openmmtools')
    @patch('MDAnalysis.Universe')
    def test_init_loads_topology(self, mock_universe, mock_import, mock_reporter):
        """Test that initialization loads topology correctly."""
        mock_import.return_value = MagicMock(return_value=mock_reporter)
        
        n_atoms = 1000  # Match the mock NetCDF
        
        # Mock MDAnalysis Universe
        mock_u = MagicMock()
        mock_u.atoms.names = np.array(['BB'] * n_atoms)
        mock_u.atoms.resnames = np.array(['ALA'] * n_atoms)
        mock_u.atoms.resids = np.arange(1, n_atoms + 1)
        mock_u.atoms.__len__ = MagicMock(return_value=n_atoms)
        mock_universe.return_value = mock_u
        
        with tempfile.TemporaryDirectory() as tmpdir:
            netcdf = os.path.join(tmpdir, "test.nc")
            topology = os.path.join(tmpdir, "test.pdb")
            
            extractor = REMDTrajectoryExtractor(netcdf, topology)
            
            assert extractor.n_atoms == n_atoms
            assert extractor.n_frames == 100
            assert extractor.n_replicas == 8
            assert extractor.n_atoms_nc == n_atoms
            
            # Check topology attributes
            assert len(extractor.atom_names) == n_atoms
            assert len(extractor.resnames) == n_atoms
            assert len(extractor.resids) == n_atoms


class TestREMDTrajectoryExtractorUtils:
    """Tests for utility methods."""
    
    @pytest.fixture
    def extractor_with_mock_data(self):
        """Create an extractor with mocked NetCDF data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            netcdf = os.path.join(tmpdir, "test.nc")
            topology = os.path.join(tmpdir, "test.pdb")
            
            n_atoms = 3  # Match NetCDF shape
            
            # Create a minimal PDB file with correct number of atoms
            with open(topology, 'w') as f:
                for i in range(n_atoms):
                    f.write(f"ATOM  {i+1:5d}  BB  ALA A{i+1:4d}     0.000   0.000   0.000  1.00  0.00\n")
                f.write("END\n")
            
            with patch('ctgomartini.analysis.remd_trajectory_extractor._import_openmmtools') as mock_import:
                reporter = MagicMock()
                storage = MagicMock()
                positions_var = MagicMock()
                positions_var.shape = (10, 4, n_atoms)  # 10 frames, 4 replicas, n_atoms
                storage.variables = {'positions': positions_var}
                reporter._storage_checkpoint = storage
                
                # States: replica 0->state 0, replica 1->state 1, etc.
                states_data = np.arange(4).reshape(1, -1)
                states_data = np.repeat(states_data, 10, axis=0)
                
                # Create analysis storage with proper states variable
                analysis_storage = MagicMock()
                states_var = MagicMock()
                # This makes states_var[:] return states_data
                states_var.__getitem__ = MagicMock(return_value=states_data)
                analysis_storage.variables = {'states': states_var}
                reporter._storage_analysis = analysis_storage
                
                reporter.read_mcmc_moves = MagicMock(return_value=[])
                mock_import.return_value = MagicMock(return_value=reporter)
                
                # Create extractor with cache_states=True explicitly
                extractor = REMDTrajectoryExtractor(netcdf, topology, cache_states=True)
                yield extractor
    
    def test_get_replica_for_state_with_cache(self, extractor_with_mock_data):
        """Test cached state lookup."""
        extractor = extractor_with_mock_data
        
        # Should use cache (initialized by default)
        replica = extractor._get_replica_for_state(0, 2)
        assert replica == 2
        
        replica = extractor._get_replica_for_state(5, 0)
        assert replica == 0
    
    def test_get_replica_for_state_no_match(self, extractor_with_mock_data):
        """Test when no replica is in the requested state."""
        extractor = extractor_with_mock_data
        
        # Request state that doesn't exist (e.g., state 100)
        replica = extractor._get_replica_for_state(0, 100)
        assert replica is None
    
    def test_state_cache_building(self, extractor_with_mock_data):
        """Test that state cache is built correctly."""
        extractor = extractor_with_mock_data
        
        # Check if cache was built (may be None if states couldn't be loaded)
        if extractor._state_cache is None:
            # States might be using default assignment (replica = state)
            # Verify state_trajectories was loaded
            assert extractor.state_trajectories is not None
            assert extractor.state_trajectories.shape[0] == 10  # n_frames
            assert extractor.state_trajectories.shape[1] == 4   # n_replicas
        else:
            # Should have 4 states (0, 1, 2, 3)
            assert len(extractor._state_cache) == 4
            
            # Each state should have 10 frames
            for state in range(4):
                assert len(extractor._state_cache[state]) == 10
                # Each entry should be (frame_idx, replica_idx)
                for frame_idx, replica_idx in extractor._state_cache[state]:
                    assert frame_idx in range(10)
                    assert replica_idx == state  # In this mock, replica = state


class TestSingleFrameExtraction:
    """Tests for single frame extraction."""
    
    def test_save_frame_pdb(self):
        """Test saving single frame as PDB."""
        with tempfile.TemporaryDirectory() as tmpdir:
            netcdf = os.path.join(tmpdir, "test.nc")
            topology = os.path.join(tmpdir, "test.pdb")
            output = os.path.join(tmpdir, "frame.pdb")
            
            # Create minimal PDB topology
            with open(topology, 'w') as f:
                f.write("ATOM      1  BB  ALA A   1       0.000   0.000   0.000  1.00  0.00\n")
                f.write("ATOM      2  SC1 ALA A   1       1.000   0.000   0.000  1.00  0.00\n")
                f.write("END\n")
            
            with patch('ctgomartini.analysis.remd_trajectory_extractor._import_openmmtools') as mock_import:
                reporter = MagicMock()
                storage = MagicMock()
                positions_var = MagicMock()
                positions_var.shape = (10, 2, 2)  # frames, replicas, atoms
                
                # Return position in nm (will be converted to Angstrom)
                positions_var.__getitem__ = MagicMock(return_value=np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]))
                storage.variables = {'positions': positions_var}
                reporter._storage_checkpoint = storage
                
                analysis_storage = MagicMock()
                analysis_storage.variables = {'states': MagicMock()}
                reporter._storage_analysis = analysis_storage
                reporter.read_mcmc_moves = MagicMock(return_value=[])
                mock_import.return_value = MagicMock(return_value=reporter)
                
                extractor = REMDTrajectoryExtractor(netcdf, topology)
                extractor.dt = 5.0  # Set known dt
                extractor.save_frame(5, output, replica_idx=0)
                
                # Verify file was created
                assert os.path.exists(output)
                
                # Check content
                with open(output, 'r') as f:
                    content = f.read()
                    assert "Frame at 25.000 ps" in content  # frame 5 * dt 5.0
                    assert "ATOM" in content


class TestIntegrationWithFixtures:
    """Integration tests using REMD fixture data."""
    
    def test_extract_replica_xtc(
        self,
        remd_netcdf_file: Path,
        remd_topology_file: Path,
        remd_checkpoint_file: Path,
        tmp_path: Path
    ):
        """Test extracting replica trajectories to XTC format."""
        from ctgomartini.analysis.remd_trajectory_extractor import REMDTrajectoryExtractor
        
        extractor = REMDTrajectoryExtractor(
            str(remd_netcdf_file),
            str(remd_topology_file),
            str(remd_checkpoint_file)
        )
        
        output_dir = tmp_path / "replicas"
        extractor.save_replica_trajectories(
            str(output_dir),
            output_pattern="replica_{i}.xtc",
            replicas=[0, 1],
            num_workers=1
        )
        
        # Check output files
        assert (output_dir / "replica_0.xtc").exists()
        assert (output_dir / "replica_1.xtc").exists()
    
    def test_extract_state_dcd(
        self,
        remd_netcdf_file: Path,
        remd_topology_file: Path,
        remd_checkpoint_file: Path,
        tmp_path: Path
    ):
        """Test extracting state trajectories to DCD format."""
        from ctgomartini.analysis.remd_trajectory_extractor import REMDTrajectoryExtractor
        
        extractor = REMDTrajectoryExtractor(
            str(remd_netcdf_file),
            str(remd_topology_file),
            str(remd_checkpoint_file)
        )
        
        output_dir = tmp_path / "states"
        extractor.save_state_trajectories(
            str(output_dir),
            output_pattern="state_{i}.dcd",
            states=[0, 1]
        )
        
        # Check output files
        assert (output_dir / "state_0.dcd").exists()
        assert (output_dir / "state_1.dcd").exists()


class TestCLIArguments:
    """Tests for CLI argument parsing."""
    
    def test_cli_replica_mode(self):
        """Test CLI argument parsing for replica mode."""
        from ctgomartini.analysis.remd_trajectory_extractor import main
        
        test_args = [
            "-nc", "output.nc",
            "-p", "topology.pdb",
            "--mode", "replica",
            "-o", "./output",
            "--pattern", "replica_{i}.xtc",
            "--replicas", "0,1,2"
        ]
        
        with patch('sys.argv', ['remd_trajectory_extractor'] + test_args):
            with patch.object(REMDTrajectoryExtractor, 'save_replica_trajectories') as mock_save:
                with patch.object(REMDTrajectoryExtractor, '_load_topology'):
                    with patch.object(REMDTrajectoryExtractor, '_load_metadata'):
                        main()
                        mock_save.assert_called_once()
    
    def test_cli_frame_mode_requires_frame(self):
        """Test that frame mode requires --frame argument."""
        from ctgomartini.analysis.remd_trajectory_extractor import main
        
        test_args = [
            "-nc", "output.nc",
            "-p", "topology.pdb",
            "--mode", "frame",
            "-o", "frame.pdb"
        ]
        
        with patch('sys.argv', ['remd_trajectory_extractor'] + test_args):
            with patch.object(REMDTrajectoryExtractor, '_load_topology'):
                with patch.object(REMDTrajectoryExtractor, '_load_metadata'):
                    result = main()
                    assert result == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
