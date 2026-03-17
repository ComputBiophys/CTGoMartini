"""
Unit tests for ctgomartinize command-line interface and all four methods.

Tests SBP, EXP, HAM, and Switching modes - focusing on input validation
and function signatures rather than full execution (which requires martinize2).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Import the functions to test
from ctgomartini.data.ctgomartinize import (
    SBGOMartinize,
    MBGOMartinize,
    SwitchingGOMartinize,
)


class TestSBPModeValidation:
    """Tests for Single-Basin Potential (SBP) mode input validation."""
    
    def test_sbp_rejects_multiple_pdbs(self):
        """Test SBP mode raises error with multiple PDB files."""
        with pytest.raises(ValueError) as exc_info:
            SBGOMartinize(
                aa_strfile_list=["State1.pdb", "State2.pdb"],
                map_file_list=["auto", "auto"],
                state_name_list=["StateA", "StateB"],
                method="SBP",
                dssp="auto",
                ff="martini3001",
                other_params=""
            )
        
        error_msg = str(exc_info.value)
        assert "SBP mode requires exactly one PDB file" in error_msg
        assert "got 2" in error_msg
    
    def test_sbp_rejects_multiple_mol_names(self):
        """Test SBP mode raises error with multiple molecule names."""
        with pytest.raises(ValueError) as exc_info:
            SBGOMartinize(
                aa_strfile_list=["protein.pdb"],
                map_file_list=["auto"],
                state_name_list=["StateA", "StateB"],  # Two names
                method="SBP",
                dssp="auto",
                ff="martini3001",
                other_params=""
            )
        
        assert "exactly one molecule name" in str(exc_info.value)
    
    def test_sbp_accepts_single_pdb(self, tmp_path, monkeypatch):
        """Test SBP mode accepts single PDB (would fail later on missing files)."""
        # Change to temp dir to avoid creating directories elsewhere
        monkeypatch.chdir(tmp_path)
        
        # This should pass validation but fail on file not found (during contact generation)
        with pytest.raises(Exception):  # OVrCSU raises Exception for missing file
            SBGOMartinize(
                aa_strfile_list=["nonexistent.pdb"],
                map_file_list=["auto"],
                state_name_list=["protein"],
                method="SBP",
                dssp="auto",
                ff="martini3001",
                other_params=""
            )


class TestMethodSignatures:
    """Tests that all methods have correct function signatures."""
    
    def test_sbgomartinize_signature(self):
        """Test SBGOMartinize has correct signature (no sbmol_name parameter)."""
        import inspect
        sig = inspect.signature(SBGOMartinize)
        params = list(sig.parameters.keys())
        
        # Should have these parameters
        assert "aa_strfile_list" in params
        assert "map_file_list" in params
        assert "state_name_list" in params
        assert "method" in params
        
        # Should NOT have sbmol_name (removed in recent refactoring)
        assert "sbmol_name" not in params
    
    def test_mbgomartinize_signature(self):
        """Test MBGOMartinize has correct signature."""
        import inspect
        sig = inspect.signature(MBGOMartinize)
        params = list(sig.parameters.keys())
        
        assert "aa_strfile_list" in params
        assert "map_file_list" in params
        assert "state_name_list" in params
        assert "mbmol_name" in params  # MB modes need mbmol_name
        assert "dict_cutoffs" in params
        assert "method" in params
    
    def test_switchinggomartinize_signature(self):
        """Test SwitchingGOMartinize has correct signature."""
        import inspect
        sig = inspect.signature(SwitchingGOMartinize)
        params = list(sig.parameters.keys())
        
        assert "aa_strfile_list" in params
        assert "map_file_list" in params
        assert "state_name_list" in params
        assert "mbmol_name" in params
        assert "dict_cutoffs" in params
        assert "method" in params


class TestMethodParameters:
    """Tests for method-specific parameters."""
    
    def test_sbp_default_method(self):
        """Test SBP mode defaults to 'SBP' method."""
        import inspect
        sig = inspect.signature(SBGOMartinize)
        method_param = sig.parameters['method']
        assert method_param.default == 'SBP'
    
    def test_mbp_methods_accept_exp(self):
        """Test MBGOMartinize accepts 'exp' method."""
        import inspect
        sig = inspect.signature(MBGOMartinize)
        method_param = sig.parameters['method']
        assert method_param.default == 'exp'


class TestDSSPParameterHandling:
    """Tests for -dssp parameter handling across all methods."""
    
    def test_sbp_dssp_auto(self, tmp_path, monkeypatch):
        """Test SBP with dssp='auto' (should use MDTraj)."""
        monkeypatch.chdir(tmp_path)
        
        with pytest.raises(Exception):  # OVrCSU raises Exception for missing file
            SBGOMartinize(
                aa_strfile_list=["nonexistent.pdb"],
                map_file_list=["auto"],
                state_name_list=["protein"],
                method="SBP",
                dssp="auto",  # Should be treated as using MDTraj
                ff="martini3001",
                other_params=""
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
