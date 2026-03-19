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


class TestRemdExtendedParameters:
    """Test cases for extended REMD parameters (replica_*)."""

    @pytest.fixture
    def temp_dir(self):
        """Provide a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_replica_count_default(self, temp_dir):
        """Test default replica_count."""
        inp_content = """
REMD = yes
exc_freq = 1000
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.replica_count == 1  # Default value

    def test_replica_count_custom(self, temp_dir):
        """Test custom replica_count."""
        inp_content = """
REMD = yes
exc_freq = 1000
replica_count = 11
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.replica_count == 11

    def test_replica_c1_single_value_expansion(self, temp_dir):
        """Test that single replica_c1 value is expanded to all replicas."""
        inp_content = """
REMD = yes
exc_freq = 1000
replica_count = 5
replica_c1 = 0
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert len(inputs.replica_c1) == 5
        assert all(c == 0.0 for c in inputs.replica_c1)

    def test_replica_c1_multiple_values(self, temp_dir):
        """Test explicit replica_c1 values."""
        inp_content = """
REMD = yes
exc_freq = 1000
replica_count = 3
replica_c1 = -100 0 100
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.replica_c1 == [-100.0, 0.0, 100.0]

    def test_replica_c1_wrong_count_raises_error(self, temp_dir):
        """Test that wrong number of replica_c1 values raises error."""
        inp_content = """
REMD = yes
exc_freq = 1000
replica_count = 5
replica_c1 = -100 0 100
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        with pytest.raises(ValueError, match="replica_c1"):
            read_inputs(inp_file)

    def test_replica_coupling_math_expression(self, temp_dir):
        """Test replica_coupling with math expression (1/300)."""
        inp_content = """
REMD = yes
exc_freq = 1000
replica_count = 3
replica_coupling = 1/300
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        expected = 1.0 / 300.0
        assert len(inputs.replica_coupling) == 3
        assert abs(inputs.replica_coupling[0] - expected) < 1e-6

    def test_remd_output_default(self, temp_dir):
        """Test default remd_output."""
        inp_content = """
REMD = yes
exc_freq = 1000
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.remd_output == "output.nc"

    def test_remd_output_custom(self, temp_dir):
        """Test custom remd_output."""
        inp_content = """
REMD = yes
exc_freq = 1000
remd_output = my_remd.nc
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.remd_output == "my_remd.nc"

    def test_remd_advanced_parameters(self, temp_dir):
        """Test advanced REMD parameters."""
        inp_content = """
REMD = yes
exc_freq = 1000
replica_molname = MOL
replica_method = exp
remd_checkpoint_interval = 10
remd_online_analysis_interval = 5000
remd_mixing_scheme = swap-all
remd_reassign_velocities = yes
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.replica_molname == "MOL"
        assert inputs.replica_method == "exp"
        assert inputs.remd_checkpoint_interval == 10
        assert inputs.remd_online_analysis_interval == 5000
        assert inputs.remd_mixing_scheme == "swap-all"
        assert inputs.remd_reassign_velocities == True

    def test_remd_unsampled_topfiles(self, temp_dir):
        """Test remd_unsampled_topfiles parsing."""
        inp_content = """
REMD = yes
exc_freq = 1000
remd_unsampled_topfiles = stateA.top stateB.top stateC.top
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        assert inputs.remd_unsampled_topfiles == ["stateA.top", "stateB.top", "stateC.top"]

    def test_full_remd_configuration(self, temp_dir):
        """Test a complete REMD configuration like GlnBP tutorial."""
        inp_content = """
REMD = yes
exc_freq = 250
replica_count = 11
replica_molname = GlnBP
replica_method = exp
replica_temp = 310.0
replica_coupling = 1/300
replica_c1 = -480 -440 -400 -360 -320 -280 -240 -200 -160 -120 -80
replica_c2 = 0
remd_output = output.nc
remd_checkpoint_interval = 5
remd_online_analysis_interval = 20000
remd_mixing_scheme = swap-neighbors
remd_reassign_velocities = no
remd_unsampled_topfiles = system_stateA.top system_stateB.top
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)

        inputs = read_inputs(inp_file)
        
        # Verify all parameters
        assert inputs.remd == "yes"
        assert inputs.exc_freq == 250
        assert inputs.replica_count == 11
        assert inputs.replica_molname == "GlnBP"
        assert inputs.replica_method == "exp"
        assert inputs.replica_temp == [310.0] * 11
        assert len(inputs.replica_coupling) == 11
        assert inputs.replica_c1 == [-480.0, -440.0, -400.0, -360.0, -320.0, 
                                      -280.0, -240.0, -200.0, -160.0, -120.0, -80.0]
        assert inputs.replica_c2 == [0.0] * 11
        assert inputs.remd_output == "output.nc"
        assert inputs.remd_checkpoint_interval == 5
        assert inputs.remd_online_analysis_interval == 20000
        assert inputs.remd_mixing_scheme == "swap-neighbors"
        assert inputs.remd_reassign_velocities == False
        assert inputs.remd_unsampled_topfiles == ["system_stateA.top", "system_stateB.top"]
