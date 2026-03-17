"""
Unit tests for MD simulation runner.

Tests the MDRunner class for standard molecular dynamics simulations.
"""

from __future__ import annotations

import os
import tempfile
from unittest.mock import Mock, patch, MagicMock

import pytest

from ctgomartini.core.simulation import MDRunner


class TestMDRunner:
    """Test cases for MDRunner class."""

    @pytest.fixture
    def temp_dir(self):
        """Provide a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    @pytest.fixture
    def mock_inp_file(self, temp_dir):
        """Create a mock input file for testing."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
platform = CPU
input = npt.gro
topol = system.top
output = md.gro
output_pdb = md.pdb
odcd = md.dcd
ochk = md.chk
nstout = 100
nstdcd = 100
append = no
mini_nstep = 0
gen_vel = no
rest = no
pcouple = no
nonbonded_cutoff = 1.1
epsilon_r = 15.0
fric_coeff = 1.0
"""
        inp_file = os.path.join(temp_dir, "md.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)
        return inp_file

    def test_md_runner_init(self, mock_inp_file):
        """Test MDRunner initialization."""
        runner = MDRunner(mock_inp_file)

        assert runner.inputs.nstep == 1000
        assert runner.inputs.dt == 0.020
        assert runner.inputs.temp == 310.0
        assert runner.inputs.platform == "CPU"
        assert runner.simulation is None

    def test_detect_existing_output_no_files(self, mock_inp_file, temp_dir):
        """Test detect_existing_output when no files exist."""
        runner = MDRunner(mock_inp_file)

        # No files should exist yet
        result = runner.detect_existing_output()
        assert result is False

    def test_detect_existing_output_with_dcd(self, mock_inp_file, temp_dir):
        """Test detect_existing_output with existing DCD file."""
        runner = MDRunner(mock_inp_file)
        
        # Create the DCD file in the correct location (relative to working dir)
        dcd_file = runner.inputs.odcd
        with open(dcd_file, "w") as f:
            f.write("")

        try:
            result = runner.detect_existing_output()
            assert result is True
        finally:
            # Cleanup
            if os.path.exists(dcd_file):
                os.remove(dcd_file)

    def test_detect_existing_output_with_chk(self, mock_inp_file, temp_dir):
        """Test detect_existing_output with existing checkpoint file."""
        runner = MDRunner(mock_inp_file)
        
        # Create the checkpoint file in the correct location
        chk_file = runner.inputs.ochk
        with open(chk_file, "w") as f:
            f.write("")

        try:
            result = runner.detect_existing_output()
            assert result is True
        finally:
            # Cleanup
            if os.path.exists(chk_file):
                os.remove(chk_file)


class TestMDRunnerRemdDetection:
    """Test that MDRunner correctly handles REMD inputs."""

    @pytest.fixture
    def temp_dir(self):
        """Provide a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_md_runner_with_remd_no_exc_freq(self, temp_dir):
        """Test that MDRunner raises error when REMD=yes but exc_freq missing."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
platform = CPU
REMD = yes
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        with pytest.raises(ValueError, match="exc_freq"):
            MDRunner(inp_file)
