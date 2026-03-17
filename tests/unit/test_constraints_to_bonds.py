"""
Tests for constraints_to_bonds utility module.

This module tests the conversion of GROMACS constraints to harmonic bonds.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from ctgomartini.utils.constraints_to_bonds import (
    convert_constraints_to_bonds_in_molecule,
    convert_constraints_to_bonds,
)


class MockMolecule:
    """Mock Molecule class for testing."""

    def __init__(self, topology: dict | None = None):
        self._topology = topology or {}


class TestConvertConstraintsToBondsInMolecule:
    """Tests for convert_constraints_to_bonds_in_molecule function."""

    def test_no_constraints_section(self):
        """Test molecule without constraints section returns 0."""
        mol = MockMolecule({"bonds": [["1", "2", "1", "0.5"]]})
        result = convert_constraints_to_bonds_in_molecule(mol, fc=50000.0)
        assert result == 0
        assert "constraints" not in mol._topology

    def test_empty_constraints(self):
        """Test molecule with empty constraints returns 0."""
        mol = MockMolecule({"constraints": [], "bonds": []})
        result = convert_constraints_to_bonds_in_molecule(mol, fc=50000.0)
        assert result == 0
        assert "constraints" not in mol._topology

    def test_single_constraint(self):
        """Test converting a single constraint."""
        mol = MockMolecule({
            "constraints": [["1", "2", "1", "0.147"]],
            "bonds": []
        })
        result = convert_constraints_to_bonds_in_molecule(mol, fc=50000.0)
        assert result == 1
        assert "constraints" not in mol._topology
        assert len(mol._topology["bonds"]) == 1
        assert mol._topology["bonds"][0] == ["1", "2", "1", "0.147", "50000.0"]

    def test_multiple_constraints(self):
        """Test converting multiple constraints."""
        mol = MockMolecule({
            "constraints": [
                ["1", "2", "1", "0.147"],
                ["3", "4", "1", "0.153"],
                ["5", "6", "2", "0.134"],
            ],
            "bonds": [["7", "8", "1", "0.1", "1000.0"]]
        })
        result = convert_constraints_to_bonds_in_molecule(mol, fc=2000.0)
        assert result == 3
        assert "constraints" not in mol._topology
        assert len(mol._topology["bonds"]) == 4  # 1 existing + 3 converted
        # Check converted bonds have correct FC
        assert mol._topology["bonds"][1] == ["1", "2", "1", "0.147", "2000.0"]
        assert mol._topology["bonds"][2] == ["3", "4", "1", "0.153", "2000.0"]
        assert mol._topology["bonds"][3] == ["5", "6", "2", "0.134", "2000.0"]

    def test_custom_force_constant(self):
        """Test using custom force constant."""
        mol = MockMolecule({
            "constraints": [["1", "2", "1", "0.147"]],
            "bonds": []
        })
        result = convert_constraints_to_bonds_in_molecule(mol, fc=12345.6)
        assert result == 1
        assert mol._topology["bonds"][0][-1] == "12345.6"

    def test_original_bonds_preserved(self):
        """Test that original bonds are preserved."""
        original_bonds = [
            ["10", "11", "1", "0.2", "5000.0"],
            ["12", "13", "1", "0.3", "6000.0"],
        ]
        mol = MockMolecule({
            "constraints": [["1", "2", "1", "0.147"]],
            "bonds": original_bonds.copy()
        })
        convert_constraints_to_bonds_in_molecule(mol, fc=50000.0)
        assert mol._topology["bonds"][0:2] == original_bonds


class TestConvertConstraintsToBondsIntegration:
    """Integration tests for convert_constraints_to_bonds with real files."""

    @pytest.fixture
    def sample_topology_with_constraints(self):
        """Sample topology content with constraints."""
        return """;
; Sample topology for testing
;

[ moleculetype ]
; name       nrexcl
Protein      3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
   1   P1     1   ALA     BB   1   0.000  72.000
   2   P1     2   GLY     BB   2   0.000  72.000
   3   P1     3   LEU     BB   3   0.000  72.000
   4   P1     4   VAL     BB   4   0.000  72.000

[ bonds ]
; ai  aj  funct   r0  fc
   1   2      1  0.35  5000
   2   3      1  0.35  5000

[ constraints ]
; ai  aj  funct   r0
   3   4      1  0.40

[ angles ]
; ai  aj  ak  funct   theta0  fc
   1   2   3      2     127.0  50.0
"""

    @pytest.fixture
    def sample_topology_no_constraints(self):
        """Sample topology without constraints."""
        return """[ moleculetype ]
Protein      3

[ atoms ]
   1   P1     1   ALA     BB   1   0.000  72.000
   2   P1     2   GLY     BB   2   0.000  72.000

[ bonds ]
   1   2      1  0.35  5000
"""

    def test_convert_constraints_in_file(self, sample_topology_with_constraints, tmp_path):
        """Test converting constraints in a real file."""
        # Create a temporary topology file
        top_file = tmp_path / "test.itp"
        top_file.write_text(sample_topology_with_constraints)

        # Convert constraints
        result = convert_constraints_to_bonds(top_file, "Protein", fc=50000.0)
        assert result == 1

        # Read and verify the modified file
        modified_content = top_file.read_text()
        
        # Should not have constraints section anymore
        assert "[ constraints ]" not in modified_content
        
        # Should have the constraint converted to bond
        assert "50000" in modified_content  # Force constant added

    def test_no_constraints_in_file(self, sample_topology_no_constraints, tmp_path):
        """Test file without constraints returns 0."""
        top_file = tmp_path / "test.itp"
        top_file.write_text(sample_topology_no_constraints)

        result = convert_constraints_to_bonds(top_file, "Protein", fc=50000.0)
        assert result == 0

    def test_file_not_found(self, tmp_path):
        """Test FileNotFoundError for non-existent file."""
        non_existent = tmp_path / "does_not_exist.itp"
        
        with pytest.raises(FileNotFoundError, match="does_not_exist.itp"):
            convert_constraints_to_bonds(non_existent, "Protein")

    def test_molecule_not_found(self, sample_topology_with_constraints, tmp_path):
        """Test ValueError for non-existent molecule."""
        top_file = tmp_path / "test.itp"
        top_file.write_text(sample_topology_with_constraints)

        with pytest.raises(ValueError, match="'NonExistent' not found"):
            convert_constraints_to_bonds(top_file, "NonExistent")

    def test_custom_fc_in_file(self, sample_topology_with_constraints, tmp_path):
        """Test custom force constant in file conversion."""
        top_file = tmp_path / "test.itp"
        top_file.write_text(sample_topology_with_constraints)

        convert_constraints_to_bonds(top_file, "Protein", fc=12345.6)
        
        modified_content = top_file.read_text()
        assert "12345.6" in modified_content


class TestArgumentParsing:
    """Tests for command-line argument parsing of -constraints2bonds."""

    def test_argument_default_none(self):
        """Test that default is None (flag not provided)."""
        import argparse
        
        parser = argparse.ArgumentParser()
        parser.add_argument('-constraints2bonds', dest='constraints2bonds',
                           nargs='?', const=50000.0, default=None, type=float)
        
        # No flag
        args = parser.parse_args([])
        assert args.constraints2bonds is None

    def test_argument_flag_only(self):
        """Test -constraints2bonds without value uses const."""
        import argparse
        
        parser = argparse.ArgumentParser()
        parser.add_argument('-constraints2bonds', dest='constraints2bonds',
                           nargs='?', const=50000.0, default=None, type=float)
        
        # Flag only
        args = parser.parse_args(['-constraints2bonds'])
        assert args.constraints2bonds == 50000.0

    def test_argument_with_value(self):
        """Test -constraints2bonds with custom value."""
        import argparse
        
        parser = argparse.ArgumentParser()
        parser.add_argument('-constraints2bonds', dest='constraints2bonds',
                           nargs='?', const=50000.0, default=None, type=float)
        
        # Flag with value
        args = parser.parse_args(['-constraints2bonds', '2000'])
        assert args.constraints2bonds == 2000.0

    def test_argument_with_float_value(self):
        """Test -constraints2bonds with float value."""
        import argparse
        
        parser = argparse.ArgumentParser()
        parser.add_argument('-constraints2bonds', dest='constraints2bonds',
                           nargs='?', const=50000.0, default=None, type=float)
        
        # Flag with float value
        args = parser.parse_args(['-constraints2bonds', '12345.67'])
        assert args.constraints2bonds == 12345.67
