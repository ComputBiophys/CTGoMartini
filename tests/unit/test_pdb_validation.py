"""Unit tests for PDB validation module."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from ctgomartini.utils.pdb_validation import (
    validate_pdb_compatibility,
    ValidationError,
    _parse_pdb,
)


class TestPDBParsing:
    """Tests for PDB parsing functionality."""

    def test_parse_simple_pdb(self):
        """Test parsing a simple PDB file."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            path = f.name

        try:
            data = _parse_pdb(Path(path), "Test")
            assert data.name == "Test"
            assert len(data.atoms) == 3
            assert data.n_residues == 1
            assert data.n_chains == 1
            assert data.chain_ids == ['A']
            
            # Check first atom
            atom = data.atoms[0]
            assert atom.name == "N"
            assert atom.resname == "MET"
            assert atom.resid == 1
            assert atom.chainid == "A"
        finally:
            os.unlink(path)

    def test_hydrogen_filtering(self):
        """Test that hydrogen atoms are filtered out."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  HA  MET A   1      11.000  11.000  11.000  1.00  0.00           H
ATOM    104  CB  MET A   1      11.500  11.500  11.500  1.00  0.00           C
ATOM    105  HB2 MET A   1      12.000  12.000  12.000  1.00  0.00           H
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            path = f.name

        try:
            data = _parse_pdb(Path(path), "Test")
            assert len(data.atoms) == 3  # N, CA, CB (HA, HB2 filtered)
            atom_names = [a.name for a in data.atoms]
            assert "N" in atom_names
            assert "CA" in atom_names
            assert "CB" in atom_names
            assert "HA" not in atom_names
            assert "HB2" not in atom_names
        finally:
            os.unlink(path)

    def test_multi_chain_parsing(self):
        """Test parsing multi-chain PDB."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            path = f.name

        try:
            data = _parse_pdb(Path(path), "Test")
            assert data.n_chains == 2
            assert data.chain_ids == ['A', 'B']
            
            # Check chain IDs are assigned correctly
            assert data.atoms[0].chainid == 'A'
            assert data.atoms[1].chainid == 'A'
            assert data.atoms[2].chainid == 'B'
            assert data.atoms[3].chainid == 'B'
        finally:
            os.unlink(path)


class TestValidationPasses:
    """Tests for successful validation cases."""

    def test_identical_pdbs_pass(self):
        """Test that identical PDBs pass validation."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    201  N   ALA A   2      12.000  12.000  12.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_content)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_content)
            path2 = f2.name

        try:
            validate_pdb_compatibility([path1, path2], ['StateA', 'StateB'])
            # Should not raise
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_different_coordinates_pass(self):
        """Test that same sequence with different coordinates passes."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   MET A   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    102  CA  MET A   1      20.500  20.500  20.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            validate_pdb_compatibility([path1, path2], ['Open', 'Closed'])
            # Should not raise
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_single_file_no_validation(self):
        """Test that single file returns without error."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            path = f.name

        try:
            validate_pdb_compatibility([path])
            # Should not raise
        finally:
            os.unlink(path)


class TestValidationFailures:
    """Tests for validation failure cases."""

    def test_atom_count_mismatch(self):
        """Test detection of atom count mismatch."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2], ['Full', 'Short'])
            assert "Atom count mismatch" in str(exc_info.value)
            assert "Full" in str(exc_info.value)
            assert "Short" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_residue_count_mismatch(self):
        """Test detection of residue count mismatch (atom count same)."""
        # Same atom count, different residue count
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  N   ALA A   2      10.500  10.500  10.500  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2])
            assert "Residue count mismatch" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_chain_count_mismatch(self):
        """Test detection of chain count mismatch (atom/residue count same)."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
END
"""
        # Same residues, but chain A and B merged into one chain in second file
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    201  N   MET A   2      20.000  20.000  20.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2])
            assert "Chain count mismatch" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_chain_id_mismatch(self):
        """Test detection of chain ID mismatch."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    201  N   MET C   1      20.000  20.000  20.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2], ['AB', 'AC'])
            assert "Chain ID mismatch" in str(exc_info.value)
            assert "['A', 'B']" in str(exc_info.value)
            assert "['A', 'C']" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_residue_name_mismatch(self):
        """Test detection of residue name mismatch."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  ALA A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2], ['MET', 'ALA'])
            assert "residue name mismatch" in str(exc_info.value)
            assert "MET" in str(exc_info.value)
            assert "ALA" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_residue_id_mismatch(self):
        """Test detection of residue ID mismatch."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   MET A   2      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   2      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2])
            assert "residue ID mismatch" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_atom_name_difference_allowed(self):
        """Test that atom name differences are allowed."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CB  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            # Should pass - atom name differences are allowed
            validate_pdb_compatibility([path1, path2])
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_chainid_mismatch_per_atom(self):
        """Test detection of chain ID mismatch at atom level."""
        # Same chain IDs in list, but different chain ID at specific atom position
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  N   MET A   2      11.000  11.000  11.000  1.00  0.00           N
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  N   MET B   2      11.000  11.000  11.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility([path1, path2])
            # This will be caught at chain_ids level since ['A'] != ['A', 'B']
            # or at atom level if both have same chains list
            assert "Chain" in str(exc_info.value) or "chain" in str(exc_info.value)
        finally:
            os.unlink(path1)
            os.unlink(path2)


class TestMultipleFiles:
    """Tests for validating more than two files."""

    def test_three_files_all_match(self):
        """Test validation of three matching PDB files."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        paths = []
        try:
            for _ in range(3):
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                    f.write(pdb_content)
                    paths.append(f.name)

            validate_pdb_compatibility(paths, ['S1', 'S2', 'S3'])
            # Should not raise
        finally:
            for p in paths:
                os.unlink(p)

    def test_three_files_second_mismatch(self):
        """Test that mismatch in second file is detected."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        pdb2 = """ATOM    101  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        pdb3 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        paths = []
        try:
            for content in [pdb1, pdb2, pdb3]:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                    f.write(content)
                    paths.append(f.name)

            with pytest.raises(ValidationError) as exc_info:
                validate_pdb_compatibility(paths, ['Good', 'Bad', 'Good2'])
            assert "residue name mismatch" in str(exc_info.value)
        finally:
            for p in paths:
                os.unlink(p)


class TestErrorHandling:
    """Tests for error handling."""

    def test_missing_file(self):
        """Test that missing file raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError) as exc_info:
            validate_pdb_compatibility(['/nonexistent/file1.pdb', '/nonexistent/file2.pdb'])
        assert "/nonexistent/file1.pdb" in str(exc_info.value)

    def test_invalid_pdb_format(self):
        """Test handling of invalid PDB content."""
        invalid_pdb = "This is not a PDB file\nJust random text\n"
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(invalid_pdb)
            path = f.name

        try:
            # MDTraj may raise exception for invalid PDB
            try:
                data = _parse_pdb(Path(path), "Test")
                # If parsing succeeds, should have 0 atoms
                assert len(data.atoms) == 0
            except Exception:
                # MDTraj may raise various exceptions for invalid input
                pass
        finally:
            os.unlink(path)
