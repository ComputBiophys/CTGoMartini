"""
Unit tests for PDB validation module.

Tests the PDBCompatibilityValidator and related functionality
for checking PDB file compatibility in multiple-basin topology generation.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from ctgomartini.utils.pdb_validation import (
    PDBCompatibilityValidator,
    validate_pdb_compatibility,
    ResidueInfo,
    StructureSummary,
    ValidationReport,
    MismatchDetail,
)


class TestResidueInfo:
    """Tests for ResidueInfo dataclass."""

    def test_residue_info_creation(self):
        """Test basic ResidueInfo creation."""
        residue = ResidueInfo(
            index=0,
            res_seq=1,
            res_name="MET",
            chain_id="A",
            atoms=["N", "CA", "C", "O", "CB", "CG", "SD", "CE"]
        )
        assert residue.atom_count == 8
        assert str(residue) == "MET1(A):8atoms"

    def test_residue_info_empty_atoms(self):
        """Test ResidueInfo with empty atom list."""
        residue = ResidueInfo(
            index=1,
            res_seq=2,
            res_name="ALA",
            chain_id="A",
            atoms=[]
        )
        assert residue.atom_count == 0


class TestPDBCompatibilityValidatorBasic:
    """Basic tests for PDBCompatibilityValidator."""

    def test_validator_creation(self):
        """Test validator initialization."""
        validator = PDBCompatibilityValidator(verbose=True)
        assert validator.verbose is True
        assert validator.PROTONATION_VARIANTS is not None

    def test_single_file_no_validation_needed(self):
        """Test that single file returns empty list."""
        validator = PDBCompatibilityValidator()

        # Create a simple PDB
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            pdb_path = f.name

        try:
            reports = validator.validate([pdb_path])
            assert reports == []
        finally:
            os.unlink(pdb_path)


class TestPDBCompatibilityMatching:
    """Tests for matching/compatible PDB files."""

    def test_identical_pdbs_pass(self):
        """Test that identical PDBs pass validation."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   ALA A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  ALA A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   ALA A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   ALA A   2      13.500  13.500  13.500  1.00  0.00           O
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_content)
            pdb1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_content)
            pdb2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate([pdb1, pdb2], ['StateA', 'StateB'])

            assert len(reports) == 1
            assert reports[0].is_valid is True
            assert reports[0].ref_residue_count == 2
            assert reports[0].other_residue_count == 2
        finally:
            os.unlink(pdb1)
            os.unlink(pdb2)

    def test_pdbs_with_different_coordinates_pass(self):
        """Test that PDBs with same sequence but different coordinates pass."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
END
"""
        # Same sequence, different coordinates (different conformation)
        pdb2 = """ATOM    101  N   MET A   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    102  CA  MET A   1      20.500  20.500  20.500  1.00  0.00           C
ATOM    103  C   MET A   1      21.000  21.000  21.000  1.00  0.00           C
ATOM    104  O   MET A   1      21.500  21.500  21.500  1.00  0.00           O
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            validate_pdb_compatibility([path1, path2], ['Open', 'Closed'], verbose=False)
            # Should not raise
        finally:
            os.unlink(path1)
            os.unlink(path2)


class TestPDBCompatibilityResidueCountMismatch:
    """Tests for residue count mismatch detection."""

    def test_missing_terminal_residue(self):
        """Test detection of missing terminal residue."""
        pdb_full = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   ALA A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  ALA A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   ALA A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   ALA A   2      13.500  13.500  13.500  1.00  0.00           O
END
"""
        pdb_truncated = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_full)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_truncated)
            path2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)

            with pytest.raises(ValueError) as exc_info:
                validator.validate([path1, path2], ['Full', 'Truncated'])

            error_msg = str(exc_info.value)
            assert "RESIDUE COUNT MISMATCH" in error_msg
            assert "difference: 1" in error_msg
            assert "ALA2" in error_msg  # The missing residue
            assert "missing in 'Truncated'" in error_msg
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_extra_residue_in_second_structure(self):
        """Test detection of extra residue in second structure."""
        pdb_short = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
END
"""
        pdb_long = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   LYS A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  LYS A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   LYS A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   LYS A   2      13.500  13.500  13.500  1.00  0.00           O
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_short)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_long)
            path2 = f2.name

        try:
            with pytest.raises(ValueError) as exc_info:
                validate_pdb_compatibility([path1, path2], ['Short', 'Long'], verbose=False)

            error_msg = str(exc_info.value)
            assert "RESIDUE COUNT MISMATCH" in error_msg
            assert "LYS2" in error_msg
            assert "present in 'Long' but not in 'Short'" in error_msg
        finally:
            os.unlink(path1)
            os.unlink(path2)


class TestPDBCompatibilitySequenceMismatch:
    """Tests for residue sequence mismatch detection."""

    def test_point_mutation_detection(self):
        """Test detection of point mutation (ALA vs GLY)."""
        pdb_wt = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   ALA A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  ALA A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   ALA A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   ALA A   2      13.500  13.500  13.500  1.00  0.00           O
END
"""
        pdb_mutant = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   GLY A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  GLY A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   GLY A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   GLY A   2      13.500  13.500  13.500  1.00  0.00           O
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_wt)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_mutant)
            path2 = f2.name

        try:
            with pytest.raises(ValueError) as exc_info:
                validate_pdb_compatibility([path1, path2], ['WT', 'Mutant'], verbose=False)

            error_msg = str(exc_info.value)
            assert "RESIDUE SEQUENCE MISMATCH" in error_msg
            assert "Position 2" in error_msg
            assert "ALA" in error_msg
            assert "GLY" in error_msg
            assert "Different residue types" in error_msg
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_multiple_sequence_differences(self):
        """Test detection of multiple sequence differences."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   ALA A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  ALA A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   ALA A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   ALA A   2      13.500  13.500  13.500  1.00  0.00           O
ATOM    301  N   VAL A   3      14.000  14.000  14.000  1.00  0.00           N
ATOM    302  CA  VAL A   3      14.500  14.500  14.500  1.00  0.00           C
ATOM    303  C   VAL A   3      15.000  15.000  15.000  1.00  0.00           C
ATOM    304  O   VAL A   3      15.500  15.500  15.500  1.00  0.00           O
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
ATOM    201  N   GLY A   2      12.000  12.000  12.000  1.00  0.00           N
ATOM    202  CA  GLY A   2      12.500  12.500  12.500  1.00  0.00           C
ATOM    203  C   GLY A   2      13.000  13.000  13.000  1.00  0.00           C
ATOM    204  O   GLY A   2      13.500  13.500  13.500  1.00  0.00           O
ATOM    301  N   LEU A   3      14.000  14.000  14.000  1.00  0.00           N
ATOM    302  CA  LEU A   3      14.500  14.500  14.500  1.00  0.00           C
ATOM    303  C   LEU A   3      15.000  15.000  15.000  1.00  0.00           C
ATOM    304  O   LEU A   3      15.500  15.500  15.500  1.00  0.00           O
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb1)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb2)
            path2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)

            with pytest.raises(ValueError) as exc_info:
                validator.validate([path1, path2], ['Structure1', 'Structure2'])

            error_msg = str(exc_info.value)
            # Should detect both mismatches
            assert "ALA" in error_msg and "GLY" in error_msg
            assert "VAL" in error_msg and "LEU" in error_msg
        finally:
            os.unlink(path1)
            os.unlink(path2)


class TestPDBCompatibilityMultiFile:
    """Tests for validating multiple files at once."""

    def test_three_file_validation(self):
        """Test validation of three PDB files."""
        pdb_base = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    103  C   MET A   1      11.000  11.000  11.000  1.00  0.00           C
ATOM    104  O   MET A   1      11.500  11.500  11.500  1.00  0.00           O
END
"""
        paths = []
        try:
            for i in range(3):
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                    f.write(pdb_base)
                    paths.append(f.name)

            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate(paths, ['State1', 'State2', 'State3'])

            # Should have 2 reports (State2 vs State1, State3 vs State1)
            assert len(reports) == 2
            for report in reports:
                assert report.is_valid is True
        finally:
            for p in paths:
                os.unlink(p)

    def test_three_files_with_mismatch(self):
        """Test that mismatch is caught in multi-file validation."""
        pdb1 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        pdb3 = """ATOM    101  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  ALA A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""  # Different residue name

        paths = []
        try:
            for content in [pdb1, pdb2, pdb3]:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                    f.write(content)
                    paths.append(f.name)

            with pytest.raises(ValueError) as exc_info:
                validate_pdb_compatibility(paths, ['Good1', 'Good2', 'Bad'], verbose=False)

            error_msg = str(exc_info.value)
            assert "MET" in error_msg and "ALA" in error_msg
        finally:
            for p in paths:
                os.unlink(p)


class TestPDBCompatibilityErrors:
    """Tests for error conditions."""

    def test_missing_file_raises_file_not_found(self):
        """Test that missing file raises FileNotFoundError."""
        validator = PDBCompatibilityValidator()

        with pytest.raises(FileNotFoundError) as exc_info:
            validator.validate(['/nonexistent/path.pdb', '/another/missing.pdb'])

        assert "/nonexistent/path.pdb" in str(exc_info.value)

    def test_invalid_pdb_format(self):
        """Test handling of invalid PDB content."""
        invalid_pdb = """This is not a valid PDB file
Just some random text
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(invalid_pdb)
            path = f.name

        try:
            validator = PDBCompatibilityValidator()
            # Should not crash, just find no atoms
            with pytest.raises(ValueError) as exc_info:
                validator.validate([path, path])
            # Error message will indicate 0 residues found
        finally:
            os.unlink(path)


class TestValidationReport:
    """Tests for ValidationReport formatting."""

    def test_valid_report_format(self):
        """Test formatting of valid report."""
        report = ValidationReport(
            is_valid=True,
            ref_name="Ref",
            other_name="Other",
            ref_file=Path("ref.pdb"),
            other_file=Path("other.pdb"),
            ref_residue_count=100,
            other_residue_count=100
        )

        error_msg = report.format_error()
        assert error_msg == ""  # Valid report has no error

    def test_report_with_warnings(self):
        """Test report with atom count warnings."""
        report = ValidationReport(
            is_valid=True,
            ref_name="Ref",
            other_name="Other",
            ref_file=Path("ref.pdb"),
            other_file=Path("other.pdb"),
            warnings=["Residue X has significant atom count difference"]
        )

        assert len(report.warnings) == 1


class TestConvenienceFunction:
    """Tests for the validate_pdb_compatibility convenience function."""

    def test_convenience_function_passes(self):
        """Test that convenience function works for valid case."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_content)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_content)
            path2 = f2.name

        try:
            # Should not raise
            validate_pdb_compatibility([path1, path2], ['A', 'B'], verbose=False)
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_convenience_function_single_file(self):
        """Test that single file returns immediately."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            path = f.name

        try:
            # Single file should not raise
            validate_pdb_compatibility([path], verbose=False)
        finally:
            os.unlink(path)


class TestChainOrderValidation:
    """Tests for chain order validation."""

    def test_correct_chain_order_no_warning(self):
        """Test that correct alphabetical chain order produces no warnings."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_content)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_content)
            path2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate([path1, path2], ['StateA', 'StateB'])
            
            # Check no chain order warnings
            for report in reports:
                chain_warnings = [w for w in report.warnings if 'Chain order' in w]
                assert len(chain_warnings) == 0, f"Unexpected chain order warning: {chain_warnings}"
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_non_alphabetical_chain_order_warning(self):
        """Test that non-alphabetical chain order produces warning."""
        pdb_correct = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
END
"""
        # B chain before A chain - incorrect order
        pdb_incorrect = """ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_correct)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_incorrect)
            path2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate([path1, path2], ['StateA', 'StateB'])
            
            # The second structure should have chain order warning
            state_b_warnings = [w for w in reports[0].warnings if 'not alphabetical' in w]
            assert len(state_b_warnings) > 0, "Expected chain order warning for StateB"
            assert 'B -> A' in state_b_warnings[0], f"Warning should mention B -> A order: {state_b_warnings[0]}"
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_chain_order_mismatch_between_structures(self):
        """Test that different chain orders between structures produce warning."""
        pdb_ab = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
END
"""
        # Same chains but different order
        pdb_ba = """ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_ab)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_ba)
            path2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate([path1, path2], ['StateAB', 'StateBA'])
            
            # Should have chain order mismatch warning
            mismatch_warnings = [w for w in reports[0].warnings if 'Chain order mismatch' in w]
            assert len(mismatch_warnings) > 0, "Expected chain order mismatch warning"
            assert 'A -> B' in mismatch_warnings[0], f"Warning should mention chain order: {mismatch_warnings[0]}"
            assert 'B -> A' in mismatch_warnings[0], f"Warning should mention both orders: {mismatch_warnings[0]}"
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_three_chains_correct_order(self):
        """Test correct order with three chains A -> B -> C."""
        pdb_content = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    102  CA  MET A   1      10.500  10.500  10.500  1.00  0.00           C
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    202  CA  MET B   1      20.500  20.500  20.500  1.00  0.00           C
ATOM    301  N   MET C   1      30.000  30.000  30.000  1.00  0.00           N
ATOM    302  CA  MET C   1      30.500  30.500  30.500  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f1:
            f1.write(pdb_content)
            path1 = f1.name

        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f2:
            f2.write(pdb_content)
            path2 = f2.name

        try:
            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate([path1, path2], ['State1', 'State2'])
            
            for report in reports:
                chain_warnings = [w for w in report.warnings if 'Chain order' in w]
                assert len(chain_warnings) == 0, f"Unexpected chain order warning: {chain_warnings}"
        finally:
            os.unlink(path1)
            os.unlink(path2)

    def test_multiple_files_chain_order_consistency(self):
        """Test chain order consistency across three files."""
        pdb_abc = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    301  N   MET C   1      30.000  30.000  30.000  1.00  0.00           N
END
"""
        pdb_acb = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    301  N   MET C   1      30.000  30.000  30.000  1.00  0.00           N
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
END
"""
        pdb_abc2 = """ATOM    101  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM    201  N   MET B   1      20.000  20.000  20.000  1.00  0.00           N
ATOM    301  N   MET C   1      30.000  30.000  30.000  1.00  0.00           N
END
"""
        paths = []
        try:
            for content in [pdb_abc, pdb_acb, pdb_abc2]:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                    f.write(content)
                    paths.append(f.name)

            validator = PDBCompatibilityValidator(verbose=False)
            reports = validator.validate(paths, ['Ref', 'ACB', 'ABC2'])
            
            # Second file (ACB) should have:
            # 1. Internal chain order warning (not alphabetical)
            # 2. Chain order mismatch with reference
            acb_warnings = reports[0].warnings
            internal_warnings = [w for w in acb_warnings if 'not alphabetical' in w]
            mismatch_warnings = [w for w in acb_warnings if 'Chain order mismatch' in w]
            
            assert len(internal_warnings) > 0, "Expected internal chain order warning for ACB"
            assert len(mismatch_warnings) > 0, "Expected chain order mismatch warning"
            
            # Third file (ABC2) should have no warnings (same order as reference)
            abc2_warnings = reports[1].warnings
            abc2_chain_warnings = [w for w in abc2_warnings if 'Chain order' in w]
            assert len(abc2_chain_warnings) == 0, f"ABC2 should have no chain warnings: {abc2_chain_warnings}"
        finally:
            for p in paths:
                os.unlink(p)
