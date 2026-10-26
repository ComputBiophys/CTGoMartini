"""Force-level consistency tests for the Multiple Basin Potential on a real system.

Complements test_mb_go_martini.py (OpenMM vs precomputed GROMACS references)
with checks that need OpenMM only, on the GlnBP EXP/HAM fixtures
(4435 particles, Reference platform, single frame Strfiles/GlnBP_No0.gro):

1. Finite-difference check: for representative non-virtual-site atoms the
   analytic forces match the central-difference gradient of the total
   potential energy, at both the fixture-default and shifted mixing
   parameters.
2. Analytic mixing check: the MBP CustomCVForce contribution (isolated via a
   force group) equals the EXP/HAM mixing formula applied to the collective
   variable values, before and after Context.setParameter.
3. Legacy equivalence: a System rebuilt with the pre-parameterization
   hardcoded expression (same inner forces, XML-copied) reproduces the
   parameterized system's energy and forces to (near) machine precision.

Only single-point energy/force evaluations are performed (no MD). The
fixtures are used read-only; mixing parameters are changed at runtime via
Context.setParameter instead of rewriting the itp file.
"""

import math
import os
import re

import numpy as np
import openmm as mm
import openmm.unit as unit
import pytest
from openmm import app

from ctgomartini.topology import MartiniTopFile
from tests.conftest import WorkingDirectoryContext

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "../fixtures/mb_go_martini/GlnBP")
STRFILE = os.path.join(FIXTURE_DIR, "Strfiles/GlnBP_No0.gro")

MBP_FORCE_GROUP = 31
FD_STEP_NM = 1e-5
# Protein BB, protein SC1, mid-protein, first W, first NA, last CL.
FD_ATOM_CANDIDATES = (0, 1, 365, 732, 4358, 4434)

# Measured worst central-difference deviation on the Reference platform is
# 8.7e-6 kJ/mol/nm (both methods); tolerances carry a >10x margin.
FD_RTOL = 1e-4
FD_ATOL = 1e-4

# Shifted mixing parameters that keep both basin weights in (0.01, 0.99) for
# the GlnBP_No0 frame, so the checks stay sensitive to both basins.
SHIFTED_PARAMETERS = {
    "EXP": {"beta": 1 / 1000, "C1": -550.0, "C2": -50.0},
    "HAM": {"delta": 400.0, "C1": -800.0, "C2": -50.0},
}


def _build_mbp_case(method: str) -> dict:
    """Build the GlnBP MBP system from the fixtures (read-only)."""
    conf = app.GromacsGroFile(STRFILE)
    box_vectors = conf.getPeriodicBoxVectors()
    with WorkingDirectoryContext(os.path.join(FIXTURE_DIR, method)):
        top = MartiniTopFile("system.top", periodicBoxVectors=box_vectors)
        system = top.create_system(nonbonded_cutoff=1.1 * unit.nanometer, epsilon_r=15.0)

    mbp_force = None
    mbp_index = -1
    for index, force in enumerate(system.getForces()):
        if isinstance(force, mm.CustomCVForce) and "ctg_" in force.getEnergyFunction():
            assert mbp_force is None, "expected exactly one MBP CustomCVForce"
            mbp_force = force
            mbp_index = index
    assert mbp_force is not None, "no MBP CustomCVForce found"
    # Isolate the mixing contribution so it can be compared against the
    # analytic mixing formula; must happen before Context creation.
    mbp_force.setForceGroup(MBP_FORCE_GROUP)

    return {
        "method": method,
        "system": system,
        "mbp_force": mbp_force,
        "mbp_index": mbp_index,
        "positions": conf.getPositions(),
        "box_vectors": box_vectors,
    }


@pytest.fixture(scope="module", params=["EXP", "HAM"])
def mbp_case(request: pytest.FixtureRequest) -> dict:
    return _build_mbp_case(request.param)


def _make_reference_context(case: dict) -> mm.Context:
    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(
        case["system"], integrator, mm.Platform.getPlatformByName("Reference")
    )
    context.setPositions(case["positions"])
    context.setPeriodicBoxVectors(*case["box_vectors"])
    context.computeVirtualSites()
    return context


def _energy(context: mm.Context, groups: set[int] | None = None) -> float:
    kwargs = {"getEnergy": True}
    if groups is not None:
        kwargs["groups"] = groups
    return context.getState(**kwargs).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )


def _forces(context: mm.Context) -> np.ndarray:
    return context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )


def _global_parameters(mbp_force: mm.CustomCVForce) -> dict[str, float]:
    """Global parameter defaults keyed by suffix (beta/delta/C1/C2/...)."""
    return {
        mbp_force.getGlobalParameterName(i).rsplit("_", 1)[-1]:
            mbp_force.getGlobalParameterDefaultValue(i)
        for i in range(mbp_force.getNumGlobalParameters())
    }


def _set_mixing_parameters(
    context: mm.Context, mbp_force: mm.CustomCVForce, overrides: dict[str, float]
) -> None:
    for i in range(mbp_force.getNumGlobalParameters()):
        name = mbp_force.getGlobalParameterName(i)
        suffix = name.rsplit("_", 1)[-1]
        if suffix in overrides:
            context.setParameter(name, overrides[suffix])


def _basin_energies(mbp_force: mm.CustomCVForce, context: mm.Context) -> dict[int, float]:
    """Per-basin energy sums E_i = sum_j CV(state{i}_force{j})."""
    names = [
        mbp_force.getCollectiveVariableName(i)
        for i in range(mbp_force.getNumCollectiveVariables())
    ]
    values = mbp_force.getCollectiveVariableValues(context)
    energies: dict[int, float] = {}
    for name, value in zip(names, values, strict=True):
        match = re.fullmatch(r"state(\d+)_force\d+", name)
        assert match is not None, f"unexpected CV name: {name}"
        energies[int(match.group(1))] = energies.get(int(match.group(1)), 0.0) + value
    return energies


def _analytic_mixing(method: str, e1: float, e2: float, params: dict[str, float]) -> float:
    if method == "EXP":
        beta = params["beta"]
        a1 = -beta * (e1 + params["C1"])
        a2 = -beta * (e2 + params["C2"])
        m = max(a1, a2)
        return -(math.log(math.exp(a1 - m) + math.exp(a2 - m)) + m) / beta
    dv = params["C2"] - params["C1"]
    alpha = (e1 - e2 - dv) / 2
    return (e1 + e2 + dv) / 2 - math.sqrt(alpha**2 + params["delta"] ** 2)


def _basin_weights(
    method: str, e1: float, e2: float, params: dict[str, float]
) -> tuple[float, float]:
    """dE/dE1 and dE/dE2 of the mixing formula."""
    if method == "EXP":
        beta = params["beta"]
        a1 = -beta * (e1 + params["C1"])
        a2 = -beta * (e2 + params["C2"])
        m = max(a1, a2)
        p1, p2 = math.exp(a1 - m), math.exp(a2 - m)
        return p1 / (p1 + p2), p2 / (p1 + p2)
    dv = params["C2"] - params["C1"]
    alpha = (e1 - e2 - dv) / 2
    w1 = 0.5 * (1 - alpha / math.sqrt(alpha**2 + params["delta"] ** 2))
    return w1, 1.0 - w1


def _pick_fd_atoms(system: mm.System) -> list[int]:
    """Representative atoms for the finite-difference check (virtual sites skipped:
    their positions are slaved to the constructing atoms, so their gradient is
    redistributed rather than local)."""
    picked = []
    for candidate in FD_ATOM_CANDIDATES:
        atom = candidate
        while atom < system.getNumParticles() and system.isVirtualSite(atom):
            atom += 1
        assert atom < system.getNumParticles()
        picked.append(atom)
    return picked


def _fd_gradient(context: mm.Context, positions_nm: np.ndarray, atom: int, dim: int) -> float:
    """Central-difference dU/dx for one coordinate of one atom.

    Virtual site positions must be recomputed after every setPositions:
    OpenMM does not update them implicitly, and the energies of interactions
    involving virtual sites would otherwise be evaluated at stale positions.
    """
    plus = positions_nm.copy()
    plus[atom, dim] += FD_STEP_NM
    context.setPositions(plus * unit.nanometer)
    context.computeVirtualSites()
    e_plus = _energy(context)

    minus = positions_nm.copy()
    minus[atom, dim] -= FD_STEP_NM
    context.setPositions(minus * unit.nanometer)
    context.computeVirtualSites()
    e_minus = _energy(context)

    context.setPositions(positions_nm * unit.nanometer)
    context.computeVirtualSites()
    return (e_plus - e_minus) / (2 * FD_STEP_NM)


def _xml_force_copy(force: mm.Force) -> mm.Force:
    return mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(force))


def _copy_virtual_site(vsite: mm.VirtualSite) -> mm.VirtualSite:
    """Deep-copy a virtual site definition (getVirtualSite returns a borrowed
    reference, which System.setVirtualSite refuses to take ownership of)."""
    particles = [vsite.getParticle(i) for i in range(vsite.getNumParticles())]
    if isinstance(vsite, mm.LocalCoordinatesSite):
        return mm.LocalCoordinatesSite(
            particles,
            vsite.getOriginWeights(),
            vsite.getXWeights(),
            vsite.getYWeights(),
            vsite.getLocalPosition(),
        )
    if isinstance(vsite, mm.TwoParticleAverageSite):
        return mm.TwoParticleAverageSite(
            particles[0], particles[1], vsite.getWeight(0), vsite.getWeight(1)
        )
    if isinstance(vsite, mm.ThreeParticleAverageSite):
        return mm.ThreeParticleAverageSite(
            particles[0],
            particles[1],
            particles[2],
            vsite.getWeight(0),
            vsite.getWeight(1),
            vsite.getWeight(2),
        )
    if isinstance(vsite, mm.OutOfPlaneSite):
        return mm.OutOfPlaneSite(
            particles[0],
            particles[1],
            particles[2],
            vsite.getWeight12(),
            vsite.getWeight13(),
            vsite.getWeightCross(),
        )
    raise TypeError(f"Unsupported virtual site type: {type(vsite).__name__}")


def _legacy_expression(mbp_force: mm.CustomCVForce) -> str:
    """Reconstruct the pre-parameterization hardcoded energy expression.

    Mirrors the legacy implementation (parameters baked into the expression
    string, deltaV and delta^2 precomputed for HAM) around the same collective
    variable definitions as the parameterized force.
    """
    params = _global_parameters(mbp_force)
    groups: dict[int, list[str]] = {}
    for i in range(mbp_force.getNumCollectiveVariables()):
        name = mbp_force.getCollectiveVariableName(i)
        match = re.fullmatch(r"state(\d+)_force\d+", name)
        assert match is not None, f"unexpected CV name: {name}"
        groups.setdefault(int(match.group(1)), []).append(name)
    energy_defs = "\n".join(
        f"energy{i} = {' + '.join(groups[i])};" for i in sorted(groups)
    )

    if "beta" in params:
        beta = params["beta"]
        terms = " + ".join(
            f"exp(-{beta!r} * (energy{i} + {params[f'C{i}']!r}))" for i in sorted(groups)
        )
        return f"-1/{beta!r} * log({terms});\n{energy_defs}"

    delta_v = params["C2"] - params["C1"]
    delta_sq = params["delta"] ** 2
    return (
        f"(energy1+energy2+{delta_v!r})/2 - "
        f"sqrt(((energy1-energy2-{delta_v!r})/2)^2+{delta_sq!r});\n{energy_defs}"
    )


def _build_legacy_system(system: mm.System, mbp_index: int) -> mm.System:
    """Rebuild the system with the MBP force replaced by the legacy hardcoded
    expression. All other forces (and the MBP inner forces) are XML deep
    copies: Force objects cannot be shared between two Systems."""
    mbp_force = system.getForce(mbp_index)
    legacy_mbp = mm.CustomCVForce(_legacy_expression(mbp_force))
    for i in range(mbp_force.getNumCollectiveVariables()):
        legacy_mbp.addCollectiveVariable(
            mbp_force.getCollectiveVariableName(i),
            _xml_force_copy(mbp_force.getCollectiveVariable(i)),
        )

    legacy_system = mm.System()
    for i in range(system.getNumParticles()):
        legacy_system.addParticle(system.getParticleMass(i))
    legacy_system.setDefaultPeriodicBoxVectors(*system.getDefaultPeriodicBoxVectors())
    for i in range(system.getNumConstraints()):
        legacy_system.addConstraint(*system.getConstraintParameters(i))
    for i in range(system.getNumParticles()):
        if system.isVirtualSite(i):
            legacy_system.setVirtualSite(i, _copy_virtual_site(system.getVirtualSite(i)))
    for i in range(system.getNumForces()):
        legacy_system.addForce(
            legacy_mbp if i == mbp_index else _xml_force_copy(system.getForce(i))
        )
    return legacy_system


def test_mbp_forces_match_finite_difference(mbp_case: dict) -> None:
    system = mbp_case["system"]
    context = _make_reference_context(mbp_case)
    positions_nm = np.array(mbp_case["positions"].value_in_unit(unit.nanometer), dtype=float)
    atoms = _pick_fd_atoms(system)

    parameter_sets = [("default", None), ("shifted", SHIFTED_PARAMETERS[mbp_case["method"]])]
    previous_forces = None
    for label, overrides in parameter_sets:
        if overrides is not None:
            _set_mixing_parameters(context, mbp_case["mbp_force"], overrides)

        basins = _basin_energies(mbp_case["mbp_force"], context)
        params = _global_parameters(mbp_case["mbp_force"])
        if overrides is not None:
            params.update(overrides)
        w1, w2 = _basin_weights(mbp_case["method"], basins[1], basins[2], params)
        assert 0.01 < w2 < 0.99, f"{label} weights saturated: w2={w2}"

        forces = _forces(context)
        worst_abs = 0.0
        worst_rel = 0.0
        for atom in atoms:
            for dim in range(3):
                fd = _fd_gradient(context, positions_nm, atom, dim)
                abs_err = abs(fd - (-forces[atom, dim]))
                rel_err = abs_err / max(abs(forces[atom, dim]), 1e-30)
                worst_abs = max(worst_abs, abs_err)
                worst_rel = max(worst_rel, rel_err)
                assert np.isclose(fd, -forces[atom, dim], rtol=FD_RTOL, atol=FD_ATOL), (
                    f"{mbp_case['method']} {label}: atom {atom} dim {dim}: "
                    f"F={-forces[atom, dim]:.6f}, FD={fd:.6f}"
                )
        print(
            f"[{mbp_case['method']}:{label}] w1={w1:.4f} w2={w2:.4f} "
            f"FD worst abs={worst_abs:.3e} rel={worst_rel:.3e} kJ/mol/nm "
            f"(h={FD_STEP_NM} nm, atoms={atoms})"
        )
        if previous_forces is not None:
            force_shift = np.abs(forces - previous_forces).max()
            assert force_shift > 1e-2, "forces did not respond to setParameter"
            print(f"[{mbp_case['method']}:{label}] max force shift vs default: {force_shift:.3e}")
        previous_forces = forces


def test_mbp_energy_matches_analytic_mixing(mbp_case: dict) -> None:
    context = _make_reference_context(mbp_case)

    for label, overrides in [
        ("default", None),
        ("shifted", SHIFTED_PARAMETERS[mbp_case["method"]]),
    ]:
        if overrides is not None:
            _set_mixing_parameters(context, mbp_case["mbp_force"], overrides)

        basins = _basin_energies(mbp_case["mbp_force"], context)
        params = _global_parameters(mbp_case["mbp_force"])
        if overrides is not None:
            params.update(overrides)
        expected = _analytic_mixing(mbp_case["method"], basins[1], basins[2], params)
        actual = _energy(context, groups={MBP_FORCE_GROUP})
        print(
            f"[{mbp_case['method']}:{label}] E1={basins[1]:.6f} E2={basins[2]:.6f} "
            f"analytic={expected:.6f} group={actual:.6f} diff={abs(actual - expected):.3e}"
        )
        assert actual == pytest.approx(expected, rel=1e-12, abs=1e-9)


def test_parameterized_matches_legacy_hardcoded_expression(mbp_case: dict) -> None:
    system = mbp_case["system"]
    legacy_system = _build_legacy_system(system, mbp_case["mbp_index"])
    assert legacy_system.getNumParticles() == system.getNumParticles()
    assert legacy_system.getNumForces() == system.getNumForces()

    parameterized_context = _make_reference_context(mbp_case)
    legacy_context = _make_reference_context(
        {**mbp_case, "system": legacy_system}
    )

    e_parameterized = _energy(parameterized_context)
    e_legacy = _energy(legacy_context)
    f_parameterized = _forces(parameterized_context)
    f_legacy = _forces(legacy_context)

    energy_diff = abs(e_parameterized - e_legacy)
    force_diff = np.abs(f_parameterized - f_legacy).max()
    max_force = np.abs(f_parameterized).max()
    print(
        f"[{mbp_case['method']}] legacy vs parameterized: "
        f"dE={energy_diff:.3e} kJ/mol (E={e_parameterized:.6f}), "
        f"max dF={force_diff:.3e} kJ/mol/nm (max|F|={max_force:.6f})"
    )
    # Identical arithmetic on the Reference platform: differences are zero in
    # practice; the bounds still allow for last-ulp reordering.
    assert energy_diff <= 1e-9 * max(1.0, abs(e_parameterized))
    assert force_diff <= 1e-9 * max(1.0, max_force)
