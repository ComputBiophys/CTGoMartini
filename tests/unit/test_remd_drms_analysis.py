"""
Unit tests for REMD dRMS analysis module.

Tests cover reference distance calculation, position extraction, and dRMS computation
using fixtures from tests/fixtures/remd/.
"""

import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from _pytest.fixtures import FixtureRequest

# Fixture paths
FIXTURES_DIR = Path(__file__).parent.parent / "fixtures" / "remd"
MD_RUN_DIR = Path(__file__).parent.parent / "fixtures" / "md_run" / "EXP" / "test"


class TestCalculateReferenceDistances:
    """Tests for calculate_reference_distances function."""
    
    @pytest.fixture
    def ref_files(self):
        """Provide reference file paths."""
        return [
            str(MD_RUN_DIR / "GBP_open_cg.pdb"),
            str(MD_RUN_DIR / "GBP_closed_cg.pdb"),
        ]
    
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
    """Tests for calculate_trajectory_drms function."""
    
    @pytest.fixture
    def test_files(self):
        """Provide test file paths."""
        return {
            'nc': str(FIXTURES_DIR / "output.nc"),
            'checkpoint': str(FIXTURES_DIR / "output_checkpoint.nc"),
            'ref_open': str(MD_RUN_DIR / "GBP_open_cg.pdb"),
            'ref_closed': str(MD_RUN_DIR / "GBP_closed_cg.pdb"),
        }
    
    def test_calculate_trajectory_drms_basic(self, test_files):
        """Test basic dRMS calculation."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_trajectory_drms,
            calculate_reference_distances,
        )
        
        # Calculate reference distances
        ref_distances = calculate_reference_distances([
            test_files['ref_open'],
            test_files['ref_closed'],
        ])[0]  # Use state A
        
        # Calculate dRMS for replica 0 only
        times, drms = calculate_trajectory_drms(
            test_files['nc'],
            test_files['checkpoint'],
            ref_distances,
            replica_indices=np.array([0]),
            skip=1,
            num_workers=1,
            chunk_size=2,
        )
        
        assert isinstance(times, np.ndarray)
        assert isinstance(drms, np.ndarray)
        assert times.shape[0] == drms.shape[0]  # Same number of frames
        assert drms.shape[1] == 1  # One replica
        assert np.all(drms >= 0)  # dRMS should be non-negative
    
    def test_calculate_trajectory_drms_multiple_replicas(self, test_files):
        """Test dRMS calculation for multiple replicas."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_trajectory_drms,
            calculate_reference_distances,
        )
        
        ref_distances = calculate_reference_distances([
            test_files['ref_open'],
            test_files['ref_closed'],
        ])[0]
        
        times, drms = calculate_trajectory_drms(
            test_files['nc'],
            test_files['checkpoint'],
            ref_distances,
            replica_indices=np.array([0, 1, 2]),
            skip=1,
            num_workers=2,
            chunk_size=2,
        )
        
        assert drms.shape[1] == 3  # Three replicas
    
    def test_calculate_trajectory_drms_with_skip(self, test_files):
        """Test dRMS calculation with frame skipping."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_trajectory_drms,
            calculate_reference_distances,
        )
        
        ref_distances = calculate_reference_distances([
            test_files['ref_open'],
            test_files['ref_closed'],
        ])[0]
        
        times1, drms1 = calculate_trajectory_drms(
            test_files['nc'],
            test_files['checkpoint'],
            ref_distances,
            replica_indices=np.array([0]),
            skip=1,
            num_workers=1,
            chunk_size=4,
        )
        
        times2, drms2 = calculate_trajectory_drms(
            test_files['nc'],
            test_files['checkpoint'],
            ref_distances,
            replica_indices=np.array([0]),
            skip=2,
            num_workers=1,
            chunk_size=2,
        )
        
        # Skip=2 should give about half the frames
        assert len(times2) <= len(times1)


class TestIntegrationWithFixtures:
    """Integration tests using actual fixture data."""
    
    @pytest.fixture
    def test_files(self):
        """Provide test file paths."""
        return {
            'nc': str(FIXTURES_DIR / "output.nc"),
            'checkpoint': str(FIXTURES_DIR / "output_checkpoint.nc"),
            'ref_open': str(MD_RUN_DIR / "GBP_open_cg.pdb"),
            'ref_closed': str(MD_RUN_DIR / "GBP_closed_cg.pdb"),
        }
    
    def test_full_analysis_pipeline(self, test_files, tmp_path):
        """Test full analysis pipeline."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_trajectory_drms,
            calculate_reference_distances,
            save_drms_results,
        )
        
        # Calculate reference distances for both states
        all_ref_distances = calculate_reference_distances([
            test_files['ref_open'],
            test_files['ref_closed'],
        ])
        
        # Process each state
        for state_idx, ref_dist in enumerate(all_ref_distances):
            times, drms = calculate_trajectory_drms(
                test_files['nc'],
                test_files['checkpoint'],
                ref_dist,
                replica_indices=np.array([0, 1]),
                skip=1,
                num_workers=1,
                chunk_size=4,
            )
            
            # Save results
            output_file = str(tmp_path / f"drms_state{state_idx}.dat")
            save_drms_results(
                times,
                drms,
                output_file,
                netcdf_file=test_files['nc'],
                checkpoint_file=test_files['checkpoint'],
                replica_indices=np.array([0, 1]),
            )
            
            # Verify file was created
            assert Path(output_file).exists()
            
            # Read and verify content
            data = np.loadtxt(output_file, comments='#')
            assert data.shape[0] == len(times)
            assert data.shape[1] == 3  # time + 2 replicas
    
    def test_different_replica_selections(self, test_files):
        """Test analysis with different replica selections."""
        from ctgomartini.analysis.remd_drms_analysis import (
            calculate_trajectory_drms,
            calculate_reference_distances,
        )
        
        ref_distances = calculate_reference_distances([
            test_files['ref_open'],
            test_files['ref_closed'],
        ])[0]
        
        # Test with single replica
        times1, drms1 = calculate_trajectory_drms(
            test_files['nc'],
            test_files['checkpoint'],
            ref_distances,
            replica_indices=np.array([0]),
            skip=1,
            num_workers=1,
            chunk_size=4,
        )
        assert drms1.shape[1] == 1
        
        # Test with all replicas
        times_all, drms_all = calculate_trajectory_drms(
            test_files['nc'],
            test_files['checkpoint'],
            ref_distances,
            replica_indices=None,  # All replicas
            skip=1,
            num_workers=1,
            chunk_size=4,
        )
        assert drms_all.shape[1] >= 2  # At least 2 replicas
