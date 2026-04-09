"""
Unit tests for XTCMultiStateReporter.

Tests the XTC separate storage reporter for REMD simulations,
including lite checkpoint read/write and iteration tracking.
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest
from openmm import unit, Vec3
from openmmtools import states

from ctgomartini.simulation.xtc_reporter import XTCMultiStateReporter


class TestXTCMultiStateReporter:
    """Test cases for XTCMultiStateReporter class."""

    @staticmethod
    def _create_test_sampler_states(
        n_replicas: int = 2,
        n_atoms: int = 10,
        box_size: float = 5.0,
    ) -> list:
        """Create test sampler states with random positions/velocities."""
        sampler_states = []
        for i in range(n_replicas):
            positions = np.random.rand(n_atoms, 3).astype(np.float32) * box_size + i
            velocities = np.random.randn(n_atoms, 3).astype(np.float32) * 0.1
            box_vectors = (
                unit.Quantity(Vec3(box_size, 0, 0), unit.nanometers),
                unit.Quantity(Vec3(0, box_size, 0), unit.nanometers),
                unit.Quantity(Vec3(0, 0, box_size), unit.nanometers),
            )
            sampler_states.append(
                states.SamplerState(
                    positions=unit.Quantity(positions, unit.nanometers),
                    velocities=unit.Quantity(
                        velocities, unit.nanometer / unit.picosecond
                    ),
                    box_vectors=box_vectors,
                )
            )
        return sampler_states

    def test_write_and_read_sampler_states_roundtrip(self):
        """Test that written sampler states can be read back correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            storage = os.path.join(tmpdir, "output.nc")
            xtc_dir = os.path.join(tmpdir, "xtc_trajs")

            reporter = XTCMultiStateReporter(
                storage=storage,
                open_mode="w",
                checkpoint_interval=5,
                xtc_dir=xtc_dir,
            )
            sampler_states = self._create_test_sampler_states(
                n_replicas=2, n_atoms=10
            )

            for iteration in [0, 5, 10]:
                reporter.write_sampler_states(sampler_states, iteration)
            reporter.close()

            # Read back via a new reporter instance
            reporter_read = XTCMultiStateReporter(
                storage=storage,
                open_mode="r",
                checkpoint_interval=5,
                xtc_dir=xtc_dir,
            )
            read_states = reporter_read.read_sampler_states()

            assert len(read_states) == 2
            assert read_states[0].positions.shape == (10, 3)
            assert read_states[0].velocities is not None
            assert read_states[0].box_vectors is not None

            # Verify positions match the last write
            expected_pos = sampler_states[0].positions.value_in_unit(
                unit.nanometers
            )
            actual_pos = read_states[0].positions.value_in_unit(unit.nanometers)
            np.testing.assert_allclose(actual_pos, expected_pos, rtol=1e-4)

            # Verify velocities match
            expected_vel = sampler_states[0].velocities.value_in_unit(
                unit.nanometer / unit.picosecond
            )
            actual_vel = read_states[0].velocities.value_in_unit(
                unit.nanometer / unit.picosecond
            )
            np.testing.assert_allclose(actual_vel, expected_vel, rtol=1e-4)

            reporter_read.close()

            # Verify XTC files were created
            for i in range(2):
                assert os.path.exists(os.path.join(xtc_dir, f"replica_{i}.xtc"))

    def test_read_last_iteration(self):
        """Test read_last_iteration returns the current_iteration marker."""
        with tempfile.TemporaryDirectory() as tmpdir:
            storage = os.path.join(tmpdir, "output.nc")
            reporter = XTCMultiStateReporter(
                storage=storage,
                open_mode="w",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
            )
            sampler_states = self._create_test_sampler_states(
                n_replicas=1, n_atoms=5
            )
            reporter.write_sampler_states(sampler_states, iteration=15)
            reporter.close()

            reporter_read = XTCMultiStateReporter(
                storage=storage,
                open_mode="r",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
            )
            assert reporter_read.read_last_iteration() == 15
            reporter_read.close()

    def test_read_checkpoint_iterations(self):
        """Test read_checkpoint_iterations in lite mode."""
        with tempfile.TemporaryDirectory() as tmpdir:
            storage = os.path.join(tmpdir, "output.nc")
            reporter = XTCMultiStateReporter(
                storage=storage,
                open_mode="w",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
            )
            sampler_states = self._create_test_sampler_states(
                n_replicas=1, n_atoms=5
            )
            reporter.write_sampler_states(sampler_states, iteration=20)
            reporter.close()

            reporter_read = XTCMultiStateReporter(
                storage=storage,
                open_mode="r",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
            )
            checkpoints = reporter_read.read_checkpoint_iterations()
            assert checkpoints.shape == (1,)
            assert checkpoints[0] == 20
            reporter_read.close()

    def test_lite_checkpoint_overwrites_index_zero(self):
        """Test that multiple writes only keep the latest frame in checkpoint."""
        with tempfile.TemporaryDirectory() as tmpdir:
            storage = os.path.join(tmpdir, "output.nc")
            reporter = XTCMultiStateReporter(
                storage=storage,
                open_mode="w",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
            )
            states1 = self._create_test_sampler_states(n_replicas=1, n_atoms=5)
            states2 = self._create_test_sampler_states(n_replicas=1, n_atoms=5)

            reporter.write_sampler_states(states1, iteration=5)
            reporter.write_sampler_states(states2, iteration=10)
            reporter.close()

            reporter_read = XTCMultiStateReporter(
                storage=storage,
                open_mode="r",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
            )
            read_states = reporter_read.read_sampler_states()

            # Should return states from iteration 10, not 5
            expected_pos = states2[0].positions.value_in_unit(unit.nanometers)
            actual_pos = read_states[0].positions.value_in_unit(unit.nanometers)
            np.testing.assert_allclose(actual_pos, expected_pos, rtol=1e-4)

            assert reporter_read.read_last_iteration() == 10
            reporter_read.close()

    def test_progress_initialization(self):
        """Test progress tracking attributes are initialized correctly."""
        with tempfile.TemporaryDirectory() as tmpdir:
            storage = os.path.join(tmpdir, "output.nc")
            reporter = XTCMultiStateReporter(
                storage=storage,
                open_mode="w",
                checkpoint_interval=5,
                xtc_dir=os.path.join(tmpdir, "xtc"),
                total_iterations=100,
                n_replicas=4,
                exc_freq=100,
                dt=0.02,
            )
            
            # Check progress attributes are initialized
            assert reporter._total_iterations == 100
            assert reporter._n_replicas == 4
            assert reporter._exc_freq == 100
            assert reporter._dt == 0.02
            assert reporter._progress_interval == 5
            assert reporter._progress_initialized is False
            assert reporter._progress_header_printed is False
            assert reporter._progress_start_time is None
            assert reporter._progress_start_iter is None
            reporter.close()

    def test_format_time(self):
        """Test time formatting utility."""
        # Test various durations
        assert XTCMultiStateReporter._format_time(45) == "0:45"
        assert XTCMultiStateReporter._format_time(125) == "2:05"
        assert XTCMultiStateReporter._format_time(3665) == "1:01:05"
        assert XTCMultiStateReporter._format_time(90061) == "1:01:01:01"
        
        # Test edge cases
        assert XTCMultiStateReporter._format_time(0) == "0:00"
        assert XTCMultiStateReporter._format_time(float('inf')) == "--"
        assert XTCMultiStateReporter._format_time(-1) == "--"
