"""Unit tests for the OpenMM >= 8.6 built-in REMD backend helpers.

These tests cover the parts of ctgomartini.simulation.remd_openmm that do not
require OpenMM 8.6: state construction, MBP force discovery, parameter
validation, the custom reporters (with a mock sampler), exchange statistics
parsing, config parsing, and CLI backend dispatch.
"""

import argparse
import csv
import types

import openmm as mm
import openmm.unit as unit
import pytest

import ctgomartini.cli.run_ctgomartini as cli
from ctgomartini.simulation.config import load_config
from ctgomartini.simulation.remd_openmm import (
    BasinEnergyReporter,
    ExchangeStatisticsReporter,
    basin_cv_groups,
    build_replica_states,
    find_mbp_forces,
    read_exchange_statistics,
    validate_state_parameters,
)
from ctgomartini.topology.interactions import EXPInteraction, HAMInteraction


def _bond_force(k: float) -> mm.HarmonicBondForce:
    """Harmonic bond between particles 0 and 1, r0 = 1.0 nm."""
    force = mm.HarmonicBondForce()
    force.addBond(0, 1, 1.0, k)
    return force


def _exp_force_dict() -> dict[str, list[mm.Force]]:
    return {'1': [_bond_force(100.0)], '2': [_bond_force(200.0)]}


def _make_system_with_mbp(molname: str = 'gbp', method: str = 'exp'):
    """Two-particle System with one MBP CustomCVForce; returns (system, force)."""
    if method == 'exp':
        mbp_force = EXPInteraction().addForce(
            _exp_force_dict(), 5.0, [10.0, 20.0], molname=molname
        )
    else:
        mbp_force = HAMInteraction().addForce(
            _exp_force_dict(), 5.0, [10.0, 20.0], molname=molname
        )
    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)
    system.addForce(mbp_force)
    return system, mbp_force


def _make_context(system: mm.System) -> mm.Context:
    context = mm.Context(
        system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference")
    )
    context.setPositions([mm.Vec3(0, 0, 0), mm.Vec3(0.3, 0, 0)] * unit.nanometer)
    return context


# Bond energies for k=100/200 at 0.3 nm: 0.5 * k * (0.3 - 1.0)^2
STATE1_ENERGY = 24.5
STATE2_ENERGY = 49.0


class TestBuildReplicaStates:
    def test_exp_states(self):
        states = build_replica_states(
            'gbp', 'exp', [300.0, 310.0], [0.001, 0.002], [-50.0, 0.0], [0.0, 5.0]
        )
        assert states == [
            {
                'temperature': 300.0,
                'ctg_gbp_beta': 0.001,
                'ctg_gbp_C1': -50.0,
                'ctg_gbp_C2': 0.0,
            },
            {
                'temperature': 310.0,
                'ctg_gbp_beta': 0.002,
                'ctg_gbp_C1': 0.0,
                'ctg_gbp_C2': 5.0,
            },
        ]

    def test_ham_states(self):
        states = build_replica_states(
            'gbp', 'ham', [310.0, 310.0], [1000.0, 800.0], [-50.0, -50.0], [0.0, 0.0]
        )
        assert states[0]['ctg_gbp_delta'] == 1000.0
        assert states[1]['ctg_gbp_delta'] == 800.0
        assert 'ctg_gbp_beta' not in states[0]

    def test_molname_sanitized(self):
        states = build_replica_states('Mol-A', 'exp', [300.0], [1.0], [0.0], [0.0])
        assert 'ctg_Mol_A_beta' in states[0]

    def test_invalid_method(self):
        with pytest.raises(ValueError, match="replica_method"):
            build_replica_states('gbp', 'lambda', [300.0], [1.0], [0.0], [0.0])

    def test_length_mismatch(self):
        with pytest.raises(ValueError, match="coupling"):
            build_replica_states('gbp', 'exp', [300.0, 310.0], [1.0], [0.0, 0.0], [0.0, 0.0])


class TestConfigRemdBackend:
    def _write_inp(self, tmp_path, extra: str = "") -> str:
        inp = tmp_path / "remd.inp"
        inp.write_text(
            "REMD = yes\n"
            "exc_freq = 100\n"
            "replica_count = 2\n"
            "replica_temp = 300 310\n"
            "platform = CPU\n"
            f"{extra}"
        )
        return str(inp)

    def test_default_backend(self, tmp_path):
        config = load_config(self._write_inp(tmp_path))
        assert config.remd_backend == 'openmmtools'

    def test_openmm_backend(self, tmp_path):
        config = load_config(self._write_inp(tmp_path, "remd_backend = openmm\n"))
        assert config.remd_backend == 'openmm'

    def test_invalid_backend(self, tmp_path):
        with pytest.raises(ValueError, match="remd_backend"):
            load_config(self._write_inp(tmp_path, "remd_backend = lammps\n"))


class TestFindAndValidateMBP:
    def test_find_mbp_forces(self):
        system, mbp_force = _make_system_with_mbp()
        decoy = mm.CustomCVForce("0.0")
        system.addForce(decoy)

        found = find_mbp_forces(system, 'gbp')
        assert len(found) == 1
        assert found[0].getGlobalParameterName(0) == 'ctg_gbp_beta'
        assert find_mbp_forces(system, 'other') == []

    def test_find_mbp_forces_ham(self):
        system, _ = _make_system_with_mbp(method='ham')
        found = find_mbp_forces(system, 'gbp')
        assert len(found) == 1
        assert found[0].getGlobalParameterName(0) == 'ctg_gbp_delta'

    def test_basin_cv_groups(self):
        _, mbp_force = _make_system_with_mbp()
        groups = basin_cv_groups(mbp_force)
        names = [
            mbp_force.getCollectiveVariableName(k)
            for k in range(mbp_force.getNumCollectiveVariables())
        ]
        assert [names[k] for k in groups[0]] == ['state1_force1']
        assert [names[k] for k in groups[1]] == ['state2_force1']

    def test_basin_cv_groups_requires_basin_cvs(self):
        with pytest.raises(ValueError, match="collective variables"):
            basin_cv_groups(mm.CustomCVForce("0.0"))

    def test_validate_state_parameters_ok(self):
        system, _ = _make_system_with_mbp()
        context = _make_context(system)
        states = build_replica_states(
            'gbp', 'exp', [300.0, 310.0], [0.001, 0.002], [-50.0, 0.0], [0.0, 0.0]
        )
        validate_state_parameters(context, states, 'gbp', 'exp')

    def test_validate_state_parameters_missing(self):
        system, _ = _make_system_with_mbp()
        context = _make_context(system)
        states = build_replica_states(
            'typo', 'exp', [300.0], [0.001], [0.0], [0.0]
        )
        with pytest.raises(ValueError, match="ctg_typo_beta"):
            validate_state_parameters(context, states, 'typo', 'exp')


def _mock_sampler(context, conformations, iteration=1):
    return types.SimpleNamespace(
        currentIteration=iteration,
        replicaConformation=conformations,
        simulation=types.SimpleNamespace(context=context),
        states=[{}],
    )


class TestBasinEnergyReporter:
    def test_writes_basin_energies(self, tmp_path):
        system, mbp_force = _make_system_with_mbp()
        context = _make_context(system)
        conformation = context.getState(
            getPositions=True, getVelocities=True, getParameters=True
        )
        sampler = _mock_sampler(context, [conformation])

        out = tmp_path / "cv_energies.csv"
        reporter = BasinEnergyReporter(out, 1, [('gbp', mbp_force)])
        reporter(sampler)
        reporter.close()

        with open(out) as f:
            rows = list(csv.reader(f))
        assert rows[0] == ['Iteration', 'Replica', 'Molecule', 'E1', 'E2']
        assert rows[1][0:3] == ['1', '0', 'gbp']
        assert float(rows[1][3]) == pytest.approx(STATE1_ENERGY)
        assert float(rows[1][4]) == pytest.approx(STATE2_ENERGY)

    def test_respects_report_interval_and_resume(self, tmp_path):
        system, mbp_force = _make_system_with_mbp()
        context = _make_context(system)
        conformation = context.getState(
            getPositions=True, getVelocities=True, getParameters=True
        )

        out = tmp_path / "cv_energies.csv"
        reporter = BasinEnergyReporter(out, 2, [('gbp', mbp_force)])
        reporter(_mock_sampler(context, [conformation], iteration=1))
        reporter(_mock_sampler(context, [conformation], iteration=2))
        reporter.close()
        with open(out) as f:
            rows = list(csv.reader(f))
        assert len(rows) == 2  # header + iteration 2 only

        reporter = BasinEnergyReporter(out, 2, [('gbp', mbp_force)], resume=True)
        reporter(_mock_sampler(context, [conformation], iteration=4))
        reporter.close()
        with open(out) as f:
            rows = list(csv.reader(f))
        assert len(rows) == 3  # no duplicated header on resume
        assert rows[2][0] == '4'


class TestExchangeStatistics:
    def _sampler(self, iteration, accepted, proposed):
        return types.SimpleNamespace(
            currentIteration=iteration,
            states=[{}, {}],
            n_accepted=[[0, accepted], [accepted, 0]],
            n_proposed=[[0, proposed], [proposed, 0]],
        )

    def test_reporter_and_readback(self, tmp_path):
        out = tmp_path / "exchange.csv"
        reporter = ExchangeStatisticsReporter(out, 1)
        # __call__ during iteration 2 reports the counters as of iteration 1
        reporter(self._sampler(2, 3, 8))
        # finalize writes the totals after the last completed iteration
        reporter.finalize(self._sampler(2, 3, 8))
        reporter.close()
        assert out.read_text().splitlines() == [
            "Iteration,Pair,Accepted,Proposed",
            "1,0-1,3,8",
            "2,0-1,3,8",
        ]

        # Resume appends without a new header; readback returns the last totals.
        reporter = ExchangeStatisticsReporter(out, 1, resume=True)
        reporter(self._sampler(4, 15, 24))
        reporter.finalize(self._sampler(4, 15, 24))
        reporter.close()
        assert out.read_text().splitlines() == [
            "Iteration,Pair,Accepted,Proposed",
            "1,0-1,3,8",
            "2,0-1,3,8",
            "3,0-1,15,24",
            "4,0-1,15,24",
        ]
        assert read_exchange_statistics(out) == {(0, 1): [15, 24]}

    def test_read_exchange_statistics_multi_pair(self, tmp_path):
        out = tmp_path / "exchange.csv"
        out.write_text(
            "Iteration,Pair,Accepted,Proposed\n"
            "1,0-1,1,2\n"
            "1,0-2,0,2\n"
            "1,1-2,2,2\n"
            "2,0-1,3,4\n"
            "2,0-2,1,4\n"
            "2,1-2,4,4\n"
        )
        assert read_exchange_statistics(out) == {
            (0, 1): [3, 4],
            (0, 2): [1, 4],
            (1, 2): [4, 4],
        }


class TestCLIDispatch:
    def _remd_inp(self, tmp_path, extra: str = "") -> str:
        inp = tmp_path / "remd.inp"
        inp.write_text(
            "REMD = yes\n"
            "exc_freq = 100\n"
            "replica_count = 2\n"
            "replica_molname = gbp\n"
            "replica_method = exp\n"
            "replica_temp = 300 310\n"
            "replica_coupling = 0.001 0.002\n"
            "replica_c1 = -50 0\n"
            "replica_c2 = 0 0\n"
            "platform = CPU\n"
            f"{extra}"
        )
        return str(inp)

    def test_openmm_backend_selected(self, tmp_path):
        from unittest.mock import MagicMock, patch

        inp = self._remd_inp(tmp_path, "remd_backend = openmm\n")
        config = load_config(inp)
        args = argparse.Namespace(inpfile=inp, append=False)
        mock_runner = MagicMock()
        with patch(
            "ctgomartini.simulation.BuiltinREMDRunner", return_value=mock_runner
        ) as runner_class:
            assert cli._run_remd(args, config) == 0
        runner_class.assert_called_once()
        # The openmmtools default file name becomes a directory name.
        assert runner_class.call_args.kwargs["output_data"] == "output_remd"
        mock_runner.run.assert_called_once()

    def test_openmm_backend_keeps_custom_output(self, tmp_path):
        from unittest.mock import MagicMock, patch

        inp = self._remd_inp(
            tmp_path, "remd_backend = openmm\nremd_output = my_dir\n"
        )
        config = load_config(inp)
        args = argparse.Namespace(inpfile=inp, append=True)
        mock_runner = MagicMock()
        with patch(
            "ctgomartini.simulation.BuiltinREMDRunner", return_value=mock_runner
        ) as runner_class:
            assert cli._run_remd(args, config) == 0
        assert runner_class.call_args.kwargs["output_data"] == "my_dir"
        mock_runner.extend.assert_called_once()

    def test_openmmtools_backend_default(self, tmp_path):
        from unittest.mock import MagicMock, patch

        inp = self._remd_inp(tmp_path, "remd_output = my.nc\n")
        config = load_config(inp)
        args = argparse.Namespace(inpfile=inp, append=False)
        mock_runner = MagicMock()
        with patch.object(cli, "REMDRunner", return_value=mock_runner) as runner_class:
            assert cli._run_remd(args, config) == 0
        assert runner_class.call_args.kwargs["output_data"] == "my.nc"
        mock_runner.run.assert_called_once()
