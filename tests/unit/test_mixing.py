import math

import numpy as np
import openmm as mm
import openmm.unit as unit
import pytest

from ctgomartini.topology.builder import _register_mbp_global_params
from ctgomartini.topology.interactions import EXPInteraction, HAMInteraction, MixingError


def _bond_force(k: float) -> mm.HarmonicBondForce:
    """Harmonic bond between particles 0 and 1, r0 = 1.0 nm."""
    force = mm.HarmonicBondForce()
    force.addBond(0, 1, 1.0, k)
    return force


def _make_context(*forces: mm.Force, distance_nm: float = 0.3) -> mm.Context:
    """Two-particle System on the Reference platform with particles distance_nm apart."""
    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)
    for force in forces:
        system.addForce(force)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions([mm.Vec3(0, 0, 0), mm.Vec3(distance_nm, 0, 0)] * unit.nanometer)
    return context


def _energy(context: mm.Context) -> float:
    """Potential energy in kJ/mol."""
    state = context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def _forces(context: mm.Context) -> np.ndarray:
    """Forces in kJ/mol/nm as an (N, 3) array."""
    state = context.getState(getForces=True)
    return state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)


def _exp_energy(e1: float, e2: float, beta: float, c1: float, c2: float) -> float:
    return -math.log(math.exp(-beta * (e1 + c1)) + math.exp(-beta * (e2 + c2))) / beta


def _ham_energy(e1: float, e2: float, delta: float, c1: float, c2: float) -> float:
    dv = c2 - c1
    return (e1 + e2 + dv) / 2 - math.sqrt(((e1 - e2 - dv) / 2) ** 2 + delta**2)


def _exp_weights(e1: float, e2: float, beta: float, c1: float, c2: float) -> tuple[float, float]:
    """Basin weights of the EXP mixing energy: dE/dE_i = exp(-beta*(E_i+C_i)) / sum_j."""
    a1 = -beta * (e1 + c1)
    a2 = -beta * (e2 + c2)
    m = max(a1, a2)
    p1, p2 = math.exp(a1 - m), math.exp(a2 - m)
    return p1 / (p1 + p2), p2 / (p1 + p2)


def _ham_weights(e1: float, e2: float, delta: float, c1: float, c2: float) -> tuple[float, float]:
    """Basin weights of the HAM mixing energy: dE/dE1 = (1 - alpha/sqrt(alpha^2+delta^2))/2."""
    dv = c2 - c1
    alpha = (e1 - e2 - dv) / 2
    w1 = 0.5 * (1 - alpha / math.sqrt(alpha**2 + delta**2))
    return w1, 1.0 - w1


# With r0=1.0 nm and particles 0.3 nm apart, a bond with force constant k has
# energy 0.5 * k * (0.3 - 1.0)^2 = 0.245 * k kJ/mol. The two states use
# different force constants so E1 != E2, which catches state/Ci misalignment
# in the mixing expressions.
STATE1_ENERGY = 0.245 * 100.0
STATE2_ENERGY = 0.245 * 200.0


def _exp_force_dict() -> dict[str, list[mm.Force]]:
    return {'1': [_bond_force(100.0)], '2': [_bond_force(200.0)]}


def _ham_force_dict() -> dict[str, list[mm.Force]]:
    return {'1': [_bond_force(100.0)], '2': [_bond_force(200.0)]}


def test_exp_registers_namespaced_global_parameters() -> None:
    interaction = EXPInteraction()
    force = interaction.addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    names = [force.getGlobalParameterName(i) for i in range(force.getNumGlobalParameters())]
    defaults = [
        force.getGlobalParameterDefaultValue(i) for i in range(force.getNumGlobalParameters())
    ]

    assert names == ['ctg_GBP_beta', 'ctg_GBP_C1', 'ctg_GBP_C2']
    assert defaults == pytest.approx([5.0, 10.0, 20.0])


def test_exp_energy_matches_analytic_and_follows_setparameter() -> None:
    interaction = EXPInteraction()
    mbp_force = interaction.addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    context = _make_context(mbp_force)
    expected = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 5.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(expected)

    context.setParameter('ctg_GBP_beta', 2.0)
    context.setParameter('ctg_GBP_C1', 30.0)
    expected = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 2.0, 30.0, 20.0)
    assert _energy(context) == pytest.approx(expected)


def test_ham_registers_namespaced_global_parameters() -> None:
    interaction = HAMInteraction()
    force = interaction.addForce(_ham_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    names = [force.getGlobalParameterName(i) for i in range(force.getNumGlobalParameters())]
    defaults = [
        force.getGlobalParameterDefaultValue(i) for i in range(force.getNumGlobalParameters())
    ]

    assert names == ['ctg_GBP_delta', 'ctg_GBP_C1', 'ctg_GBP_C2']
    assert defaults == pytest.approx([5.0, 10.0, 20.0])


def test_ham_energy_matches_analytic_and_follows_setparameter() -> None:
    interaction = HAMInteraction()
    mbp_force = interaction.addForce(_ham_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    context = _make_context(mbp_force)
    expected = _ham_energy(STATE1_ENERGY, STATE2_ENERGY, 5.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(expected)

    context.setParameter('ctg_GBP_delta', 2.0)
    context.setParameter('ctg_GBP_C2', 30.0)
    expected = _ham_energy(STATE1_ENERGY, STATE2_ENERGY, 2.0, 10.0, 30.0)
    assert _energy(context) == pytest.approx(expected)


def test_distinct_molnames_are_independent() -> None:
    mbp_a = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBPA')
    mbp_b = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBPB')

    context = _make_context(mbp_a, mbp_b)
    expected_a = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 5.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(2 * expected_a)

    context.setParameter('ctg_GBPA_beta', 2.0)
    shifted_a = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 2.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(shifted_a + expected_a)


def test_copies_of_same_molecule_share_parameters() -> None:
    mbp_1 = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBP')
    mbp_2 = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    # Same name + same default across forces is allowed by OpenMM.
    context = _make_context(mbp_1, mbp_2)
    expected = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 5.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(2 * expected)

    # One setParameter call drives both copies.
    context.setParameter('ctg_GBP_beta', 2.0)
    shifted = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 2.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(2 * shifted)


def test_same_name_different_default_rejected_by_registry() -> None:
    mbp_1 = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBP')
    mbp_2 = EXPInteraction().addForce(_exp_force_dict(), 5.0, [30.0, 20.0], molname='GBP')

    registry: dict[str, tuple[float, str]] = {}
    _register_mbp_global_params(registry, mbp_1, 'GBP', 'EXP', 2)
    with pytest.raises(MixingError, match='ctg_GBP_C1'):
        _register_mbp_global_params(registry, mbp_2, 'GBP', 'EXP', 2)

    # Without the registry check, OpenMM itself refuses to create a Context for
    # same-name parameters with different defaults.
    with pytest.raises(mm.OpenMMException):
        _make_context(mbp_1, mbp_2)


def test_sanitized_namespace_collision_rejected() -> None:
    # 'Mol-A' and 'Mol_A' sanitize to the same namespace; even with identical
    # defaults they must not silently share parameters.
    mbp_1 = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='Mol-A')
    mbp_2 = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='Mol_A')

    registry: dict[str, tuple[float, str]] = {}
    _register_mbp_global_params(registry, mbp_1, 'Mol-A', 'EXP', 2)
    with pytest.raises(MixingError) as excinfo:
        _register_mbp_global_params(registry, mbp_2, 'Mol_A', 'EXP', 2)

    assert 'Mol-A' in str(excinfo.value)
    assert 'Mol_A' in str(excinfo.value)


def test_ham_rejects_non_two_states() -> None:
    force_dict = {
        '1': [_bond_force(100.0)],
        '2': [_bond_force(200.0)],
        '3': [_bond_force(300.0)],
    }
    with pytest.raises(MixingError, match='exactly 2 states'):
        HAMInteraction().addForce(force_dict, 5.0, [10.0, 20.0, 30.0], molname='GBP')

    with pytest.raises(MixingError, match='exactly 2 basin energies'):
        HAMInteraction().addForce(_ham_force_dict(), 5.0, [10.0], molname='GBP')


def test_exp_rejects_invalid_input() -> None:
    with pytest.raises(MixingError, match='Empty force dictionary'):
        EXPInteraction().addForce({}, 5.0, [], molname='GBP')

    with pytest.raises(MixingError, match="doesn't match energy count"):
        EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0], molname='GBP')


def test_molname_sanitized_in_parameter_names() -> None:
    interaction = EXPInteraction()
    mbp_force = interaction.addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='Mol-A')

    names = [
        mbp_force.getGlobalParameterName(i) for i in range(mbp_force.getNumGlobalParameters())
    ]
    assert names == ['ctg_Mol_A_beta', 'ctg_Mol_A_C1', 'ctg_Mol_A_C2']

    context = _make_context(mbp_force)
    context.setParameter('ctg_Mol_A_beta', 2.0)
    expected = _exp_energy(STATE1_ENERGY, STATE2_ENERGY, 2.0, 10.0, 20.0)
    assert _energy(context) == pytest.approx(expected)


def test_parameterized_exp_matches_legacy_hardcoded_expression() -> None:
    # The legacy implementation hardcoded the parameters into the expression
    # string; the parameterized force must reproduce its default energy.
    legacy = mm.CustomCVForce(
        "-1/beta * log(exp(-beta * (energy1 + C1)) + exp(-beta * (energy2 + C2)));\n"
        "energy1 = state1_force1;\n"
        "energy2 = state2_force1;\n"
        "beta=5.0;\n"
        "C1=10.0;\n"
        "C2=20.0;"
    )
    legacy.addCollectiveVariable("state1_force1", _bond_force(100.0))
    legacy.addCollectiveVariable("state2_force1", _bond_force(200.0))

    mbp_new = EXPInteraction().addForce(_exp_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    assert _energy(_make_context(mbp_new)) == pytest.approx(_energy(_make_context(legacy)))


def test_parameterized_ham_matches_legacy_hardcoded_expression() -> None:
    legacy = mm.CustomCVForce(
        "(energy1+energy2+deltaV)/2 - sqrt(((energy1-energy2-deltaV)/2)^2+delta^2);\n"
        "energy1 = state1_force1;\n"
        "energy2 = state2_force1;\n"
        "delta=5.0;\n"
        "deltaV=C2-C1;\n"
        "C1=10.0;\n"
        "C2=20.0;"
    )
    legacy.addCollectiveVariable("state1_force1", _bond_force(100.0))
    legacy.addCollectiveVariable("state2_force1", _bond_force(200.0))

    mbp_new = HAMInteraction().addForce(_ham_force_dict(), 5.0, [10.0, 20.0], molname='GBP')

    assert _energy(_make_context(mbp_new)) == pytest.approx(_energy(_make_context(legacy)))


def test_exp_forces_are_weighted_combination_of_basin_forces() -> None:
    # The energy tests above use beta=5, where w2 underflows to zero for these
    # basin energies; a force test would then be insensitive to state 2.
    # beta=0.05 keeps both weights well away from 0/1.
    mbp_force = EXPInteraction().addForce(_exp_force_dict(), 0.05, [10.0, 20.0], molname='GBP')
    context = _make_context(mbp_force)
    forces1 = _forces(_make_context(_bond_force(100.0)))
    forces2 = _forces(_make_context(_bond_force(200.0)))

    def check(beta: float, c1: float, c2: float) -> np.ndarray:
        w1, w2 = _exp_weights(STATE1_ENERGY, STATE2_ENERGY, beta, c1, c2)
        assert 0.01 < w2 < 0.99  # both basins must contribute
        mixed = _forces(context)
        np.testing.assert_allclose(mixed, w1 * forces1 + w2 * forces2, rtol=1e-9, atol=1e-9)
        return mixed

    default_forces = check(0.05, 10.0, 20.0)

    context.setParameter('ctg_GBP_beta', 0.06)
    context.setParameter('ctg_GBP_C1', 8.0)
    context.setParameter('ctg_GBP_C2', 22.0)
    shifted_forces = check(0.06, 8.0, 22.0)

    assert np.abs(shifted_forces - default_forces).max() > 1.0


def test_ham_forces_are_weighted_combination_of_basin_forces() -> None:
    mbp_force = HAMInteraction().addForce(_ham_force_dict(), 50.0, [10.0, 20.0], molname='GBP')
    context = _make_context(mbp_force)
    forces1 = _forces(_make_context(_bond_force(100.0)))
    forces2 = _forces(_make_context(_bond_force(200.0)))

    def check(delta: float, c1: float, c2: float) -> np.ndarray:
        w1, w2 = _ham_weights(STATE1_ENERGY, STATE2_ENERGY, delta, c1, c2)
        assert 0.01 < w2 < 0.99  # both basins must contribute
        mixed = _forces(context)
        np.testing.assert_allclose(mixed, w1 * forces1 + w2 * forces2, rtol=1e-9, atol=1e-9)
        return mixed

    default_forces = check(50.0, 10.0, 20.0)

    context.setParameter('ctg_GBP_delta', 10.0)
    context.setParameter('ctg_GBP_C1', 30.0)
    context.setParameter('ctg_GBP_C2', -40.0)
    shifted_forces = check(10.0, 30.0, -40.0)

    assert np.abs(shifted_forces - default_forces).max() > 1.0
