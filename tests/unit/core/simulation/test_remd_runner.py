"""
Unit tests for REMD simulation runner.

Tests the REMDRunner class, focusing on extend/resume behavior
with and without XTC separate storage.
"""

from __future__ import annotations

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from ctgomartini.simulation.remd import REMDRunner
from ctgomartini.simulation.xtc_reporter import XTCMultiStateReporter


class TestREMDRunnerExtend:
    """Test cases for REMDRunner.extend() method."""

    @pytest.fixture
    def temp_dir(self):
        """Provide a temporary directory for test files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    @pytest.fixture
    def mock_inp_file_no_xtc(self, temp_dir):
        """Create a mock input file without XTC output."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
platform = CPU
input = npt.gro
topol = system.top
REMD = yes
exc_freq = 250
replica_count = 2
replica_molname = MOL
replica_method = exp
replica_temp = 310.0
replica_coupling = 1/300
replica_c1 = -100 0
replica_c2 = 0
remd_output = output.nc
remd_checkpoint_interval = 5
remd_xtc_output = no
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)
        return inp_file

    @pytest.fixture
    def mock_inp_file_xtc(self, temp_dir):
        """Create a mock input file with XTC output enabled."""
        inp_content = """
nstep = 1000
dt = 0.020
temp = 310
platform = CPU
input = npt.gro
topol = system.top
REMD = yes
exc_freq = 250
replica_count = 2
replica_molname = MOL
replica_method = exp
replica_temp = 310.0
replica_coupling = 1/300
replica_c1 = -100 0
replica_c2 = 0
remd_output = output.nc
remd_checkpoint_interval = 5
remd_xtc_output = yes
remd_xtc_dir = xtc_trajs
"""
        inp_file = os.path.join(temp_dir, "remd.inp")
        with open(inp_file, "w") as f:
            f.write(inp_content)
        return inp_file

    @pytest.fixture
    def mock_openmmtools(self):
        """Return mock objects for openmmtools imports."""
        mock_cache = MagicMock()
        mock_reporter = MagicMock()
        mock_sampler_cls = MagicMock()
        mock_sampler = MagicMock()
        mock_sampler.iteration = 0
        mock_sampler.number_of_iterations = 100
        mock_sampler_cls.from_storage.return_value = mock_sampler
        return mock_cache, mock_reporter, mock_sampler_cls

    def test_extend_without_xtc_passes_string_to_from_storage(
        self, mock_inp_file_no_xtc, mock_openmmtools
    ):
        """When XTC is disabled, from_storage should receive the path string."""
        mock_cache, mock_reporter, mock_sampler_cls = mock_openmmtools

        with patch(
            "ctgomartini.simulation.remd._import_openmmtools",
            return_value=(mock_cache, mock_reporter, mock_sampler_cls),
        ), patch("ctgomartini.simulation.remd.os.path.exists", return_value=True):
            runner = REMDRunner(mock_inp_file_no_xtc)
            runner.extend()

        mock_sampler_cls.from_storage.assert_called_once()
        call_args = mock_sampler_cls.from_storage.call_args[0]
        assert isinstance(call_args[0], str)
        assert os.path.basename(call_args[0]) == "output.nc"

    def test_extend_with_xtc_passes_reporter_instance_to_from_storage(
        self, mock_inp_file_xtc, mock_openmmtools
    ):
        """When XTC is enabled, from_storage should receive an XTCMultiStateReporter instance."""
        mock_cache, mock_reporter, mock_sampler_cls = mock_openmmtools

        with patch(
            "ctgomartini.simulation.remd._import_openmmtools",
            return_value=(mock_cache, mock_reporter, mock_sampler_cls),
        ), patch("ctgomartini.simulation.remd.os.path.exists", return_value=True):
            runner = REMDRunner(mock_inp_file_xtc)
            runner.extend()

        mock_sampler_cls.from_storage.assert_called_once()
        call_args = mock_sampler_cls.from_storage.call_args[0]
        assert isinstance(call_args[0], XTCMultiStateReporter)

    def test_extend_raises_when_output_missing(self, mock_inp_file_no_xtc):
        """extend() should raise FileNotFoundError when output file is missing."""
        mock_cache = MagicMock()
        mock_reporter = MagicMock()
        mock_sampler_cls = MagicMock()

        with patch(
            "ctgomartini.simulation.remd._import_openmmtools",
            return_value=(mock_cache, mock_reporter, mock_sampler_cls),
        ):
            runner = REMDRunner(mock_inp_file_no_xtc)
            with pytest.raises(FileNotFoundError, match="output file"):
                runner.extend()
