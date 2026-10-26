"""Integration tests for the OpenMM >= 8.6 built-in REMD backend.

Runs real (tiny) replica exchange simulations on the GlnBP MBP fixtures and
checks that the reduced energy matrix written by the built-in reporter is
consistent with the basin energies recorded by BasinEnergyReporter through the
EXP/HAM mixing formulas. The whole module is skipped when OpenMM < 8.6.
"""

import csv
import json
import math
import os
import shutil

import numpy as np
import pytest

pytest.importorskip("openmm.app.replicaexchangesampler")

import openmm.unit as unit

from ctgomartini.analysis.remd_mbar import FESAnalyzer, kB
from ctgomartini.simulation.remd_openmm import (
    BuiltinREMDRunner,
    read_exchange_statistics,
)
from tests.conftest import WorkingDirectoryContext

FIXTURE_BASE = os.path.join(os.path.dirname(__file__), "../fixtures/md_run")

EXP_COUPLINGS = [1 / 1000, 1 / 500]
HAM_COUPLINGS = [1000.0, 800.0]
C1_LIST = [-50.0, 0.0]
C2_LIST = [0.0, 0.0]
TEMPERATURES = [300.0, 310.0]


def exp_mixing(e1, e2, beta, c1, c2):
    return -math.log(math.exp(-beta * (e1 + c1)) + math.exp(-beta * (e2 + c2))) / beta


def ham_mixing(e1, e2, delta, c1, c2):
    dv = c2 - c1
    return (e1 + e2 + dv) / 2 - math.sqrt(((e1 - e2 - dv) / 2) ** 2 + delta**2)


def _write_remd_inp(test_dir, nstep, method, pcouple="no", c1="-50 0"):
    inp_path = os.path.join(test_dir, "remd.inp")
    coupling = "1/1000 1/500" if method == "exp" else "1000 800"
    with open(inp_path, "w") as f:
        f.write(
            f"""
nstep       = {nstep}
dt          = 0.002
temp        = 310
fric_coeff  = 1
platform    = CPU
nonbonded_cutoff = 1.1
epsilon_r   = 15.0
input       = ions.gro
topol       = system.top
ichk        = mini.chk
mini_nstep  = 0
gen_vel     = no
rest        = no
pcouple     = {pcouple}

REMD = yes
remd_backend = openmm
exc_freq = 5
replica_count = 2
replica_molname = gbp
replica_method = {method}
replica_temp = 300 310
replica_coupling = {coupling}
replica_c1 = {c1}
replica_c2 = 0 0
remd_checkpoint_interval = 1
remd_xtc_output = yes
"""
        )
    return inp_path


def _minimize_fixture(test_dir):
    """Minimize the insane-solvated starting structure once per fixture copy.

    The REMD runners (both backends) expect a clash-free input structure; the
    raw ions.gro from insane needs a short minimization first. The minimized
    conformation is written as an XML checkpoint (mini.chk), which also
    exercises the runner's ichk path.
    """
    import openmm as mm
    from openmm.app import GromacsGroFile, Simulation

    from ctgomartini.topology import MartiniTopFile

    with WorkingDirectoryContext(test_dir):
        conf = GromacsGroFile("ions.gro")
        top = MartiniTopFile(
            "system.top", periodicBoxVectors=conf.getPeriodicBoxVectors()
        )
        system = top.create_system()
        integrator = mm.LangevinIntegrator(
            310 * unit.kelvin, 1.0 / unit.picosecond, 0.001 * unit.picosecond
        )
        simulation = Simulation(
            top.topology, system, integrator, mm.Platform.getPlatformByName("CPU")
        )
        simulation.context.setPositions(conf.getPositions())
        simulation.minimizeEnergy(maxIterations=200)
        state = simulation.context.getState(getPositions=True, getParameters=True)
        with open("mini.chk", "w") as f:
            f.write(mm.XmlSerializer.serialize(state))


@pytest.fixture(scope="module")
def exp_dir(tmp_path_factory):
    test_dir = tmp_path_factory.mktemp("remd_exp")
    shutil.copytree(
        os.path.join(FIXTURE_BASE, "EXP", "template"), test_dir, dirs_exist_ok=True
    )
    _minimize_fixture(test_dir)
    return str(test_dir)


@pytest.fixture(scope="module")
def ham_dir(tmp_path_factory):
    test_dir = tmp_path_factory.mktemp("remd_ham")
    shutil.copytree(
        os.path.join(FIXTURE_BASE, "HAM", "template"), test_dir, dirs_exist_ok=True
    )
    _minimize_fixture(test_dir)
    return str(test_dir)


def _read_cv_energies(outdir):
    """Return {(iteration, replica): (E1, E2)} in kJ/mol."""
    table = {}
    with open(os.path.join(outdir, "cv_energies.csv")) as f:
        for row in csv.DictReader(f):
            table[(int(row["Iteration"]), int(row["Replica"]))] = (
                float(row["E1"]),
                float(row["E2"]),
            )
    return table


def _check_energy_consistency(outdir, method, n_iterations):
    """energy.csv (reduced) must match basin energies via the mixing formula.

    The total energy in state s is E_other(replica) + E_mixing,s(E1, E2), so
    the state-to-state difference eliminates all state-independent terms:

        U(r, s1) - U(r, s0) == E_mixing,s1 - E_mixing,s0

    where U(r, s) is the unreduced energy of replica r's conformation in
    state s. This simultaneously validates the global-parameter switching,
    the sampler's energy matrix, and the CV reporter.
    """
    couplings = EXP_COUPLINGS if method == "exp" else HAM_COUPLINGS
    mixing = exp_mixing if method == "exp" else ham_mixing

    energy = np.loadtxt(os.path.join(outdir, "energy.csv"), delimiter=",", ndmin=2)
    assert energy.shape == (n_iterations, 4)  # [frame][replica x state]
    cv_table = _read_cv_energies(outdir)

    max_dev = 0.0
    for frame in range(n_iterations):
        iteration = frame + 1
        for replica in range(2):
            e1, e2 = cv_table[(iteration, replica)]
            unreduced = [
                energy[frame][replica * 2 + s] * kB * TEMPERATURES[s] for s in range(2)
            ]
            modeled = [
                mixing(e1, e2, couplings[s], C1_LIST[s], C2_LIST[s]) for s in range(2)
            ]
            deviation = abs(
                (unreduced[1] - unreduced[0]) - (modeled[1] - modeled[0])
            )
            max_dev = max(max_dev, deviation)
    return max_dev


class TestBuiltinREMDEXP:
    def test_run_energy_consistency(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_run")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=15, method="exp")
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()

        # 3 iterations, 2 replicas
        assert runner.sampler.currentIteration == 3
        assert runner.sampler.n_proposed[0][1] == 12  # 3 iterations x 2^2 attempts
        assert runner.sampler.n_proposed[1][0] == 0  # only unordered i<j counted
        assert runner.sampler.n_accepted[0][1] <= 12

        # Output files
        for name in ("log.csv", "energy.csv", "cv_energies.csv", "exchange.csv",
                     "metadata.json", "checkpoint_0.xml", "checkpoint_1.xml",
                     "replica_0.xtc", "replica_1.xtc"):
            assert os.path.isfile(os.path.join(outdir, name)), name

        # log.csv: 3 rows, valid state indices
        log = np.loadtxt(os.path.join(outdir, "log.csv"), delimiter=",", skiprows=1)
        assert log.shape == (3, 4)
        assert set(log[:, 2:].flatten().astype(int)) <= {0, 1}
        assert log[:, 0].tolist() == [1, 2, 3]

        # exchange.csv rows read "counters after completing iteration X":
        # 2^2 random pair attempts per iteration over the single pair 0-1.
        exchange = np.loadtxt(
            os.path.join(outdir, "exchange.csv"),
            delimiter=",", skiprows=1, ndmin=2, usecols=(0, 2, 3),
        )
        assert exchange[:, 0].tolist() == [0, 1, 2, 3]
        assert exchange[:, 2].tolist() == [0, 4, 8, 12]

        # metadata.json
        with open(os.path.join(outdir, "metadata.json")) as f:
            metadata = json.load(f)
        assert metadata["method"] == "exp"
        assert metadata["molname"] == "gbp"
        assert metadata["parameter_names"]["coupling"] == "ctg_gbp_beta"
        assert [s["temperature"] for s in metadata["states"]] == TEMPERATURES
        assert [s["ctg_gbp_beta"] for s in metadata["states"]] == pytest.approx(
            EXP_COUPLINGS
        )

        max_dev = _check_energy_consistency(outdir, "exp", n_iterations=3)
        print(f"\nEXP max |DU - D(E_mixing)| = {max_dev:.3e} kJ/mol")
        assert max_dev < 1e-3

        assert runner.detect_existing_output()

    def test_extend_continues_output(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_extend")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp")
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()
            assert runner.sampler.currentIteration == 2

            # nstep may change between run and extend; the state definition
            # must not.
            inp = _write_remd_inp(exp_dir, nstep=5000, method="exp")
            resumed = BuiltinREMDRunner(inp, output_data=outdir)
            resumed.extend(n_iterations=2)
            assert resumed.sampler.currentIteration == 4
            # 2 more iterations were counted on top of the restored 4 proposals
            assert resumed.sampler.n_proposed[0][1] == 16

        log = np.loadtxt(os.path.join(outdir, "log.csv"), delimiter=",", skiprows=1)
        assert log[:, 0].tolist() == [1, 2, 3, 4]
        energy = np.loadtxt(os.path.join(outdir, "energy.csv"), delimiter=",", ndmin=2)
        assert energy.shape == (4, 4)
        # The totals row restores the counters exactly: 4 iterations x 2^2
        # random pair attempts over the single pair 0-1.
        stats = read_exchange_statistics(os.path.join(outdir, "exchange.csv"))
        assert stats[(0, 1)][1] == 16
        # Energy consistency also holds for the resumed frames
        max_dev = _check_energy_consistency(outdir, "exp", n_iterations=4)
        assert max_dev < 1e-3

    def test_run_refuses_nonempty_directory(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_refuse")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp")
            BuiltinREMDRunner(inp, output_data=outdir).run()
            with pytest.raises(ValueError, match="not empty"):
                BuiltinREMDRunner(inp, output_data=outdir).run()
            with pytest.raises(FileNotFoundError):
                BuiltinREMDRunner(inp, output_data=os.path.join(exp_dir, "missing")).extend()

    def test_pcouple_rejected(self, exp_dir):
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp", pcouple="yes")
            with pytest.raises(ValueError, match="pressure coupling"):
                BuiltinREMDRunner(inp, output_data=os.path.join(exp_dir, "output_npt"))

    def test_run_zero_iterations_rejected(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_zero")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=3, method="exp")  # < exc_freq = 5
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            with pytest.raises(ValueError, match="nstep"):
                runner.run()
            assert not os.path.exists(outdir)

    def test_extend_mismatched_c1_rejected(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_c1mismatch")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp")
            BuiltinREMDRunner(inp, output_data=outdir).run()

            inp = _write_remd_inp(exp_dir, nstep=10, method="exp", c1="-30 0")
            resumed = BuiltinREMDRunner(inp, output_data=outdir)
            with pytest.raises(ValueError, match="ctg_gbp_C1"):
                resumed.extend()

    def test_extend_missing_metadata_rejected(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_nometa")
        os.makedirs(outdir)
        with open(os.path.join(outdir, "log.csv"), "w") as f:
            f.write("Iteration,Step,Replica_0_State,Replica_1_State\n")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp")
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            with pytest.raises(FileNotFoundError, match="metadata.json"):
                runner.extend()

    def test_extend_already_complete(self, exp_dir):
        outdir = os.path.join(exp_dir, "output_complete")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=10, method="exp")
            BuiltinREMDRunner(inp, output_data=outdir).run()
            with open(os.path.join(outdir, "exchange.csv"), "rb") as f:
                exchange_before = f.read()

            resumed = BuiltinREMDRunner(inp, output_data=outdir)
            resumed.extend()  # already at the nstep-implied iteration count
            assert resumed.sampler.currentIteration == 2
            # The completed path must not append another totals row.
            with open(os.path.join(outdir, "exchange.csv"), "rb") as f:
                assert f.read() == exchange_before

    def test_fes_analyzer_reads_run_directory(self, exp_dir, tmp_path):
        outdir = os.path.join(exp_dir, "output_fes")
        with WorkingDirectoryContext(exp_dir):
            inp = _write_remd_inp(exp_dir, nstep=15, method="exp")
            BuiltinREMDRunner(inp, output_data=outdir).run()

        cv_file = tmp_path / "cv.dat"
        np.savetxt(cv_file, np.arange(9, dtype=float).reshape(3, 3))
        analyzer = FESAnalyzer(outdir, cv_file)
        assert analyzer.n_states == 2
        assert analyzer.temperatures_k.tolist() == TEMPERATURES
        assert analyzer.replica_energies.shape == (2, 2, 3)

        # unsampled energies are the basin energies reduced by the temperature
        # of the state the replica occupies at that frame; the layout is
        # [replica][basin][frame], same as the openmmtools path. The raw
        # kJ/mol channel mirrors the cv_energies.csv values exactly.
        cv_table = _read_cv_energies(outdir)
        log = np.loadtxt(os.path.join(outdir, "log.csv"), delimiter=",", skiprows=1)
        for frame in range(3):
            for replica in range(2):
                state = int(log[frame][2 + replica])
                e1, e2 = cv_table[(frame + 1, replica)]
                kt = kB * TEMPERATURES[state]
                assert analyzer.unsampled_energies[replica, 0, frame] == pytest.approx(
                    e1 / kt
                )
                assert analyzer.unsampled_energies[replica, 1, frame] == pytest.approx(
                    e2 / kt
                )
                assert analyzer._basin_energies_kj[replica, 0, frame] == pytest.approx(
                    e1
                )
                assert analyzer._basin_energies_kj[replica, 1, frame] == pytest.approx(
                    e2
                )


class TestBuiltinREMDHAM:
    def test_run_energy_consistency(self, ham_dir):
        outdir = os.path.join(ham_dir, "output_run")
        with WorkingDirectoryContext(ham_dir):
            inp = _write_remd_inp(ham_dir, nstep=15, method="ham")
            runner = BuiltinREMDRunner(inp, output_data=outdir)
            runner.run()

        assert runner.sampler.currentIteration == 3
        with open(os.path.join(outdir, "metadata.json")) as f:
            metadata = json.load(f)
        assert metadata["method"] == "ham"
        assert metadata["parameter_names"]["coupling"] == "ctg_gbp_delta"

        max_dev = _check_energy_consistency(outdir, "ham", n_iterations=3)
        print(f"\nHAM max |DU - D(E_mixing)| = {max_dev:.3e} kJ/mol")
        assert max_dev < 1e-3
