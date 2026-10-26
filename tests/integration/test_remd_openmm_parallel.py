"""Integration tests for the parallel multi-Context REMD driver (OpenMM >= 8.6).

Runs real (tiny) replica exchange simulations with remd_devices set, driving
one worker Context per replica in threads (CPU platform, so no GPU is needed).
Checks energy self-consistency against the mixing formulas, output structure
(identical to the serial driver), worker seed assignment, and extend/resume.
The whole module is skipped when OpenMM < 8.6.
"""

import json
import os
import shutil

import numpy as np
import pytest

pytest.importorskip("openmm.app.replicaexchangesampler")

from ctgomartini.simulation.remd_openmm import (
    BuiltinREMDRunner,
    read_exchange_statistics,
)
from tests.conftest import WorkingDirectoryContext
from tests.integration.test_remd_openmm_builtin import (
    FIXTURE_BASE,
    _check_energy_consistency,
    _minimize_fixture,
    _write_remd_inp,
)


@pytest.fixture(scope="module")
def exp_dir(tmp_path_factory):
    test_dir = tmp_path_factory.mktemp("remd_par_exp")
    shutil.copytree(
        os.path.join(FIXTURE_BASE, "EXP", "template"), test_dir, dirs_exist_ok=True
    )
    _minimize_fixture(test_dir)
    return str(test_dir)


@pytest.fixture(scope="module")
def ham_dir(tmp_path_factory):
    test_dir = tmp_path_factory.mktemp("remd_par_ham")
    shutil.copytree(
        os.path.join(FIXTURE_BASE, "HAM", "template"), test_dir, dirs_exist_ok=True
    )
    _minimize_fixture(test_dir)
    return str(test_dir)


class TestParallelREMDEXP:
    def test_parallel_run_energy_consistency(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_par")
        with WorkingDirectoryContext(exp_dir):
            # gen_vel + gen_seed make the run fully deterministic: seeded
            # initial velocities and worker integrator seeds.
            inp = _write_remd_inp(
                exp_dir, nstep=15, method="exp", devices="0 1",
                gen_seed=1234, gen_vel="yes",
            )
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()

        assert runner.sampler.currentIteration == 3
        assert runner.sampler.n_proposed[0][1] == 12  # 3 iterations x 2^2 attempts

        # Same output structure as the serial driver.
        for name in ("log.csv", "energy.csv", "cv_energies.csv", "exchange.csv",
                     "metadata.json", "checkpoint_0.xml", "checkpoint_1.xml",
                     "replica_0.xtc", "replica_1.xtc"):
            assert os.path.isfile(os.path.join(outdir, name)), name

        # Round-robin device assignment and per-worker seeds (gen_seed + i).
        assert runner.sampler.worker_devices == [0, 1]
        assert runner.sampler.worker_seeds == [1234, 1235]
        for worker, seed in zip(
            runner.sampler._workers, runner.sampler.worker_seeds, strict=True
        ):
            assert worker.integrator.getRandomNumberSeed() == seed

        with open(os.path.join(outdir, "metadata.json")) as f:
            metadata = json.load(f)
        assert metadata["parallel"] is True
        assert metadata["devices"] == [0, 1]
        assert metadata["worker_seeds"] == [1234, 1235]

        # log.csv Step follows the serial semantics: the last replica's
        # conformation step count (stepsPerIteration x iteration), even
        # though the workers stepped their own Contexts.
        log = np.loadtxt(os.path.join(outdir, "log.csv"), delimiter=",", skiprows=1)
        assert log[:, 1].tolist() == [5, 10, 15]

        max_dev = _check_energy_consistency(outdir, "exp", n_iterations=3)
        print(f"\nParallel EXP max |DU - D(E_mixing)| = {max_dev:.3e} kJ/mol")
        assert max_dev < 1e-3

    def test_parallel_extend_continues_output(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_par_extend")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp", devices="0 1",
                                  gen_seed=1234, gen_vel="yes")
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()
            assert runner.sampler.currentIteration == 2

            resumed_inp = _write_remd_inp(
                exp_dir, nstep=10, method="exp", devices="1 0",
                gen_seed=1234, gen_vel="yes",
            )
            # A different device list is allowed on resume (e.g. new machine).
            resumed = BuiltinREMDRunner(resumed_inp, output_data=outdir)
            resumed.extend(n_iterations=2)
            assert resumed.sampler.currentIteration == 4
            # 2 more iterations on top of the restored 4 proposals
            assert resumed.sampler.n_proposed[0][1] == 16

        log = np.loadtxt(os.path.join(outdir, "log.csv"), delimiter=",", skiprows=1)
        assert log[:, 0].tolist() == [1, 2, 3, 4]
        # Step keeps counting per-replica steps across the resume.
        assert log[:, 1].tolist() == [5, 10, 15, 20]
        energy = np.loadtxt(os.path.join(outdir, "energy.csv"), delimiter=",", ndmin=2)
        assert energy.shape == (4, 4)
        stats = read_exchange_statistics(os.path.join(outdir, "exchange.csv"))
        assert stats[(0, 1)][1] == 16
        # Energy consistency also holds for the resumed frames.
        max_dev = _check_energy_consistency(outdir, "exp", n_iterations=4)
        assert max_dev < 1e-3


class TestParallelREMDHAM:
    def test_parallel_ham_energy_consistency(self, ham_dir):
        outdir = os.path.join(ham_dir, "output_par")
        with WorkingDirectoryContext(ham_dir):
            inp = _write_remd_inp(ham_dir, nstep=15, method="ham", devices="0 1",
                                  gen_seed=5678, gen_vel="yes")
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()

        assert runner.sampler.currentIteration == 3
        assert runner.sampler.worker_seeds == [5678, 5679]

        max_dev = _check_energy_consistency(outdir, "ham", n_iterations=3)
        print(f"\nParallel HAM max |DU - D(E_mixing)| = {max_dev:.3e} kJ/mol")
        assert max_dev < 1e-3


class TestParallelREMDNPT:
    def test_parallel_npt_smoke(self, exp_dir):
        """The parallel driver inherits the NPT exchange correction unchanged."""
        outdir = os.path.join(exp_dir, "output_par_npt")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(
                exp_dir, nstep=15, method="exp", devices="0 1", pcouple="yes",
                gen_seed=1234, gen_vel="yes",
            )
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()

        assert runner.sampler.currentIteration == 3
        # NPT kwargs pass through the parallel sampler's __init__.
        assert runner.sampler.pressure_bar == 1.0
        assert runner.sampler.surface_tension_bar_nm == 0.0

        volume = np.loadtxt(os.path.join(outdir, "volume.csv"), delimiter=",")
        assert volume.shape == (3, 2)
        with open(os.path.join(outdir, "metadata.json")) as f:
            metadata = json.load(f)
        assert metadata["npt"]["enabled"] is True
        assert metadata["parallel"] is True

        max_dev = _check_energy_consistency(outdir, "exp", n_iterations=3)
        print(f"\nParallel NPT max |DU - D(E_mixing)| = {max_dev:.3e} kJ/mol")
        assert max_dev < 1e-3
