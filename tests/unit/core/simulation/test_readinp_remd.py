"""
Unit tests for REMD parameter parsing in ReadInp module.

Tests the REMD and exc_freq parameter parsing functionality
added to the _OpenMMReadInputs class.
"""

from __future__ import annotations

import os
import tempfile

import pytest

from ctgomartini.simulation import load_config as read_inputs


class TestRemdParameters:
    """Test cases for REMD parameter parsing."""

    @pytest.fixture
    def temp_dir(self):
        """Provide a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_remd_default_disabled(self, temp_dir):
        """Test that REMD is disabled by default."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
"""
        inp_file = os.path.join(temp_dir, "md.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)

        assert inputs.remd == "no"
        assert inputs.exc_freq is None

    def test_remd_enabled_with_exc_freq(self, temp_dir):
        """Test REMD enabled with exc_freq specified."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
REMD = yes
exc_freq = 1000
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)

        assert inputs.remd == "yes"
        assert inputs.exc_freq == 1000

    def test_remd_explicitly_disabled(self, temp_dir):
        """Test explicit REMD disabling."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
REMD = no
exc_freq = 1000
"""
        inp_file = os.path.join(temp_dir, "md.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)

        # When REMD = no, exc_freq should still be parsed
        assert inputs.remd == "no"
        assert inputs.exc_freq == 1000

    def test_remd_enabled_missing_exc_freq(self, temp_dir):
        """Test that REMD enabled without exc_freq raises ValueError."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
REMD = yes
"""
        inp_file = os.path.join(temp_dir, "invalid.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        with pytest.raises(ValueError, match="exc_freq"):
            read_inputs(inp_file)

    def test_remd_case_insensitive(self, temp_dir):
        """Test that REMD parameter is case insensitive."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
REMD = YES
exc_freq = 500
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)

        assert inputs.remd == "yes"
        assert inputs.exc_freq == 500

    def test_all_standard_parameters_plus_remd(self, temp_dir):
        """Test that all standard parameters work with REMD."""
        inp_content = """
mini_nstep = 0
mini_Tol = 1000.0
gen_vel = no
gen_temp = 310
nstep = 50000000
dt = 0.020
b_step = 0
append = no
input = npt.gro
topol = system.top
ichk = npt.chk
nstout = 5000
nstdcd = 5000
output = md.gro
output_pdb = md.pdb
odcd = 
oxtc = md.xtc
ochk = md.chk
defines = 
rest = no
rest_ref = npt.gro
rest_file = restraints.txt
gen_rest = no
atomname = BB
fc = 1000.0
gen_rest_file = restraints.txt
plumed = no
plumed_file = 
platform = CUDA
precision = single
GPU_id = 
temp = 310
fric_coeff = 0.1
nonbonded_cutoff = 1.1
epsilon_r = 15.0
const_tol = 
pcouple = yes
p_ref = 1.0
p_type = isotropic
p_freq = 100
REMD = yes
exc_freq = 1000
"""
        inp_file = os.path.join(temp_dir, "full.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)

        # Verify standard parameters
        assert inputs.nstep == 50000000
        assert inputs.dt == 0.020
        assert inputs.temp == 310.0
        assert inputs.platform == "CUDA"

        # Verify REMD parameters
        assert inputs.remd == "yes"
        assert inputs.exc_freq == 1000
