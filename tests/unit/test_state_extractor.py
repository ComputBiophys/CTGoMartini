"""
Unit tests for state extraction utilities.

Tests the extraction of single-state topologies from multiple-basin
topologies for use in REMD unsampled states.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from ctgomartini.utils import extract_state_topology, extract_all_states, write_state_itp
from ctgomartini.topology import MartiniTopFile


class TestStateExtractor:
    """Test cases for state extraction functionality."""

    @pytest.fixture
    def sample_mb_topology(self):
        """Create a sample multiple-basin topology file for testing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a minimal system.top with multi-basin molecule
            top_content = """#include "martini_v3.0.0.itp"

[ moleculetype ]
TestMB     1

[ atoms ]
   1   Qd     1   ASP   BB     1   1.000
   2   Qd     1   ASP   SC1    2   0.000
   3   P4     2   LEU   BB     3   0.000
   4   AC2    2   LEU   SC1    4   0.000

[ bonds ]
    1     2     1   0.27000  20000.0
    3     4     1   0.33000  20000.0

[ multiple_basin ]
True exp 2 beta C1 C2

[ multi_contacts ]
    1     2     2     1   6   0.3000  12.00000
    3     4     2     2   6   0.4000  12.00000

[ contacts ]

[ angles ]
    1     2     3     2   100.0   25.0

[ multi_angles ]
    1     2     3     2     1   2   100.0   25.0
    1     2     3     2     2   2   110.0   25.0

[ dihedrals ]

[ multi_dihedrals ]
    1     2     3     4     1   2     1   180.0   25.0   1
    1     2     3     4     2   2     1   190.0   25.0   1

[ system ]
Test System

[ molecules ]
TestMB     1
"""
            top_file = Path(tmpdir) / "test_system.top"
            ff_file = Path(tmpdir) / "martini_v3.0.0.itp"
            
            # Create minimal force field file
            ff_content = """[ defaults ]
1   1   no   1.0   1.0

[ atomtypes ]
Qd   72.0   0.0000   A   0.0   0.0
P4   72.0   0.0000   A   0.0   0.0
AC2  72.0   0.0000   A   0.0   0.0
"""
            
            with open(top_file, "w") as f:
                f.write(top_content)
            with open(ff_file, "w") as f:
                f.write(ff_content)
            
            yield str(top_file), str(ff_file)

    def test_extract_state_topology_state1(self, sample_mb_topology):
        """Test extracting state 1 from multi-basin topology."""
        top_file, _ = sample_mb_topology
        
        mol = extract_state_topology(top_file, "TestMB", state_id=1)
        
        # Verify molecule name is preserved (keep_original_name=True)
        assert mol.name == "TestMB"
        
        # Verify multi-basin sections are removed
        assert "multiple_basin" not in mol._topology
        assert "multi_contacts" not in mol._topology
        assert "multi_angles" not in mol._topology
        assert "multi_dihedrals" not in mol._topology
        
        # Verify state 1 contacts are extracted
        contacts = mol._topology.get("contacts", [])
        assert len(contacts) > 0
        # First contact should be state 1 (atom 1-2)
        assert contacts[0][0] == "1"
        assert contacts[0][1] == "2"

    def test_extract_state_topology_state2(self, sample_mb_topology):
        """Test extracting state 2 from multi-basin topology."""
        top_file, _ = sample_mb_topology
        
        mol = extract_state_topology(top_file, "TestMB", state_id=2)
        
        # Verify molecule name is preserved
        assert mol.name == "TestMB"
        
        # Verify state 2 contacts are extracted
        contacts = mol._topology.get("contacts", [])
        assert len(contacts) > 0
        # Contact should be state 2 (atom 3-4)
        assert contacts[0][0] == "3"
        assert contacts[0][1] == "4"

    def test_extract_state_topology_invalid_state(self, sample_mb_topology):
        """Test that invalid state_id raises ValueError."""
        top_file, _ = sample_mb_topology
        
        with pytest.raises(ValueError, match="Invalid state_id"):
            extract_state_topology(top_file, "TestMB", state_id=3)

    def test_extract_state_topology_invalid_molecule(self, sample_mb_topology):
        """Test that invalid molecule name raises ValueError."""
        top_file, _ = sample_mb_topology
        
        with pytest.raises(ValueError, match="not found"):
            extract_state_topology(top_file, "NonExistent", state_id=1)

    def test_extract_state_topology_not_multibasin(self):
        """Test that non-multi-basin molecule raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a single-basin topology
            top_content = """#include "martini_v3.0.0.itp"

[ moleculetype ]
SingleMol     1

[ atoms ]
   1   Qd     1   ASP   BB     1   1.000

[ system ]
Test System

[ molecules ]
SingleMol     1
"""
            ff_content = """[ defaults ]
1   1   no   1.0   1.0

[ atomtypes ]
Qd   72.0   0.0000   A   0.0   0.0
"""
            top_file = Path(tmpdir) / "test.top"
            ff_file = Path(tmpdir) / "martini_v3.0.0.itp"
            
            with open(top_file, "w") as f:
                f.write(top_content)
            with open(ff_file, "w") as f:
                f.write(ff_content)
            
            with pytest.raises(ValueError, match="not a multi-basin"):
                extract_state_topology(str(top_file), "SingleMol", state_id=1)

    def test_extract_all_states(self, sample_mb_topology):
        """Test extracting all states at once."""
        top_file, _ = sample_mb_topology
        
        with tempfile.TemporaryDirectory() as tmpdir:
            output_prefix = str(Path(tmpdir) / "TestMB")
            output_files = extract_all_states(
                top_file, 
                "TestMB", 
                output_prefix=output_prefix
            )
            
            # Should generate 2 files for 2-state system
            assert len(output_files) == 2
            assert output_files[0] == f"{output_prefix}_stateA.itp"
            assert output_files[1] == f"{output_prefix}_stateB.itp"
            
            # Verify files exist
            assert os.path.exists(output_files[0])
            assert os.path.exists(output_files[1])

    def test_write_state_itp(self, sample_mb_topology):
        """Test writing extracted state to ITP file."""
        top_file, _ = sample_mb_topology
        
        with tempfile.TemporaryDirectory() as tmpdir:
            mol = extract_state_topology(top_file, "TestMB", state_id=1)
            output_file = Path(tmpdir) / "TestMB_state.itp"
            
            write_state_itp(mol, output_file)
            
            assert output_file.exists()
            
            # Verify content by reading back
            content = output_file.read_text()
            assert "[ moleculetype ]" in content
            assert "TestMB" in content

    def test_extract_all_states_default_prefix(self, sample_mb_topology):
        """Test extract_all_states with default prefix."""
        top_file, _ = sample_mb_topology
        
        with tempfile.TemporaryDirectory() as tmpdir:
            original_dir = os.getcwd()
            os.chdir(tmpdir)
            try:
                output_files = extract_all_states(top_file, "TestMB")
                
                # Should use molecule name as prefix
                assert len(output_files) == 2
                assert "TestMB_stateA.itp" in output_files[0]
                assert "TestMB_stateB.itp" in output_files[1]
            finally:
                os.chdir(original_dir)

    def test_extract_state_keep_original_name_false(self, sample_mb_topology):
        """Test extracting state with keep_original_name=False."""
        top_file, _ = sample_mb_topology
        
        mol = extract_state_topology(
            top_file, 
            "TestMB", 
            state_id=1, 
            keep_original_name=False
        )
        
        # Molecule name should have _stateA suffix
        assert mol.name == "TestMB_stateA"
