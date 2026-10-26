"""Analytic validation of adjustable-cutoff Go contact potentials.

Covers the runtime contact cutoff (``-go-up`` / ``contact_rcut``) at the
interaction level for rcut in {0.9, 1.1, 1.3, 1.5} nm:

- ``Contacts6_12`` (functype 2) and ``Contacts10_12`` (functype 3):
  V(r) = step(rcut - r) * (LJ(r) - LJ(rcut)) evaluated on 2-particle systems
  against the independent analytic formula, including the long tail
  r in (1.1, 1.5) that only exists for rcut = 1.5, exact zeros beyond rcut,
  and V(rcut) = 0 (the force is discontinuous at rcut by construction).
- A 3-particle system with per-bond different (sigma, epsilon, rcut) to
  validate multi-bond accumulation and per-bond cutoff independence.
- ``MultiAllBonds`` multi_contacts functype 2 and 3 (per-state 6-12 / 10-12
  contacts), including a regression test that collected functype-3 bonds
  actually reach the built force.
- ``MultiAllBonds`` mixing several contact functypes (1+2, 2+3, 1+2+3) in one
  force: per-bond parameter names are disambiguated per handler on collision
  (single-type builds keep the canonical names), and the mixed force matches
  both the independent analytic sum and the sum of single-functype builds.
- ``ContactsLJ`` (functype 1) fixed 1.1 nm baseline.

All evaluations use the Reference platform (double precision).
"""

import re

import numpy as np
import openmm as mm
import openmm.unit as unit
import pytest

from ctgomartini.topology.interactions.bonds import Contacts6_12, Contacts10_12, ContactsLJ
from ctgomartini.topology.interactions.multi_state import MultiAllBonds

RCUTS = [0.9, 1.1, 1.3, 1.5]
# Fixed distance grid: for rcut=1.5 the points 1.2/1.4 probe the long tail
# beyond the legacy 1.1 cutoff; the same points must give exactly 0 for
# rcut <= 1.1, which is what makes the cutoff variable.
DISTANCES = [0.45, 0.55, 0.7, 0.9, 1.0, 1.2, 1.4, 1.6]
SIGMA = 0.5
EPSILON = 12.0

ENERGY_RTOL = 1e-12
ENERGY_ATOL = 1e-9
FORCE_RTOL = 1e-9
FORCE_ATOL = 1e-9


def lj_6_12(r: float, c6: float, c12: float, rcut: float) -> tuple[float, float]:
    """Shifted 6-12 contact potential and dU/dr; both exactly 0 beyond rcut."""
    if r > rcut:
        return 0.0, 0.0
    u = (c12 / r**12 - c6 / r**6) - (c12 / rcut**12 - c6 / rcut**6)
    dudr = -12 * c12 / r**13 + 6 * c6 / r**7
    return u, dudr


def lj_10_12(r: float, c10: float, c12: float, rcut: float) -> tuple[float, float]:
    """Shifted 10-12 contact potential and dU/dr; both exactly 0 beyond rcut."""
    if r > rcut:
        return 0.0, 0.0
    u = (c12 / r**12 - c10 / r**10) - (c12 / rcut**12 - c10 / rcut**10)
    dudr = -12 * c12 / r**13 + 10 * c10 / r**11
    return u, dudr


def make_system(force: mm.Force, positions: list[tuple[float, float, float]]) -> mm.Context:
    system = mm.System()
    for _ in positions:
        system.addParticle(1.0)
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions([mm.Vec3(*p) for p in positions] * unit.nanometer)
    return context


def evaluate(context: mm.Context) -> tuple[float, np.ndarray]:
    state = context.getState(getEnergy=True, getForces=True)
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    )
    return energy, np.array(forces, dtype=float)


def c6_c12_from_sigma_eps(sigma: float, eps: float) -> tuple[float, float]:
    return 4 * eps * sigma**6, 4 * eps * sigma**12


def assert_pair_matches_analytic(
    template: mm.CustomBondForce,
    potential,
    rcut: float,
) -> None:
    """2-particle system along x: energy/forces vs analytic over the full grid.

    The bond acts between particles 0 (at origin) and 1 (at r, 0, 0), so
    F0_x = +dU/dr and F1_x = -dU/dr, with transverse components exactly 0.
    """
    for r in DISTANCES + [rcut - 1e-9, rcut, rcut + 1e-9]:
        context = make_system(clone_bond_force(template), [(0.0, 0.0, 0.0), (r, 0.0, 0.0)])
        energy, forces = evaluate(context)
        u_ref, dudr_ref = potential(r, rcut)

        assert energy == pytest.approx(u_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL), (
            f"energy mismatch at r={r}, rcut={rcut}: {energy} vs {u_ref}"
        )
        np.testing.assert_allclose(
            forces[0], [dudr_ref, 0.0, 0.0], rtol=FORCE_RTOL, atol=FORCE_ATOL,
            err_msg=f"force on particle 0 at r={r}, rcut={rcut}",
        )
        np.testing.assert_allclose(
            forces[1], [-dudr_ref, 0.0, 0.0], rtol=FORCE_RTOL, atol=FORCE_ATOL,
            err_msg=f"force on particle 1 at r={r}, rcut={rcut}",
        )

        if r > rcut:
            # step(rcut-r) == 0: exactly zero energy and force outside the cutoff
            assert energy == 0.0
            assert np.abs(forces).max() == 0.0
        elif r == rcut:
            # shifted potential is continuous: zero energy at the boundary
            assert energy == pytest.approx(0.0, abs=1e-12)


def clone_bond_force(template: mm.CustomBondForce) -> mm.CustomBondForce:
    """Fresh CustomBondForce with identical bonds: a Force object cannot be
    owned by two Systems, so each evaluation gets its own copy."""
    fresh = mm.CustomBondForce(template.getEnergyFunction())
    for i in range(template.getNumPerBondParameters()):
        fresh.addPerBondParameter(template.getPerBondParameterName(i))
    for i in range(template.getNumBonds()):
        a1, a2, params = template.getBondParameters(i)
        fresh.addBond(a1, a2, params)
    return fresh


@pytest.mark.parametrize("rcut", RCUTS)
def test_contacts_6_12_matches_analytic(rcut: float) -> None:
    interaction = Contacts6_12()
    interaction.add_interaction(["1", "2", "2", str(SIGMA), str(EPSILON), str(rcut)])
    c6, c12 = c6_c12_from_sigma_eps(SIGMA, EPSILON)
    assert_pair_matches_analytic(
        interaction.mm_force,
        lambda r, rc: lj_6_12(r, c6, c12, rc),
        rcut,
    )


@pytest.mark.parametrize("rcut", RCUTS)
def test_contacts_10_12_matches_analytic(rcut: float) -> None:
    interaction = Contacts10_12()
    interaction.add_interaction(["1", "2", "3", str(SIGMA), str(EPSILON), str(rcut)])
    c10 = 4 * EPSILON * SIGMA**10
    c12 = 4 * EPSILON * SIGMA**12
    assert_pair_matches_analytic(
        interaction.mm_force,
        lambda r, rc: lj_10_12(r, c10, c12, rc),
        rcut,
    )


def test_contacts_6_12_rcut_1_5_extends_beyond_legacy_1_1() -> None:
    """The same distance r=1.2 nm is inside the cutoff for rcut=1.5 (nonzero)
    and outside for rcut=1.1 (exactly zero)."""
    c6, c12 = c6_c12_from_sigma_eps(SIGMA, EPSILON)
    energies = {}
    for rcut in (1.1, 1.5):
        interaction = Contacts6_12()
        interaction.add_interaction(["1", "2", "2", str(SIGMA), str(EPSILON), str(rcut)])
        context = make_system(interaction.mm_force, [(0.0, 0.0, 0.0), (1.2, 0.0, 0.0)])
        energies[rcut], _ = evaluate(context)
    assert energies[1.1] == 0.0
    assert energies[1.5] == pytest.approx(
        lj_6_12(1.2, c6, c12, 1.5)[0], rel=ENERGY_RTOL, abs=ENERGY_ATOL
    )
    assert energies[1.5] < 0.0  # attractive tail


def test_contacts_6_12_multi_bond_accumulation_with_mixed_rcuts() -> None:
    """3 particles, two bonds sharing particle 1 with different (sigma, eps, rcut):
    validates force accumulation and per-bond cutoff independence."""
    positions = [(0.0, 0.0, 0.0), (0.6, 0.0, 0.0), (1.1, 0.4, 0.0)]
    bond_specs = [
        # (atom i, atom j, sigma, epsilon, rcut)
        (0, 1, 0.50, 12.0, 0.9),
        (1, 2, 0.60, 8.0, 1.5),
    ]
    interaction = Contacts6_12()
    for i, j, sigma, eps, rcut in bond_specs:
        interaction.add_interaction([str(i + 1), str(j + 1), "2", str(sigma), str(eps), str(rcut)])

    context = make_system(interaction.mm_force, positions)
    energy, forces = evaluate(context)

    pos = np.array(positions, dtype=float)
    e_ref = 0.0
    f_ref = np.zeros_like(pos)
    for i, j, sigma, eps, rcut in bond_specs:
        c6, c12 = c6_c12_from_sigma_eps(sigma, eps)
        d = pos[j] - pos[i]
        r = float(np.linalg.norm(d))
        u, dudr = lj_6_12(r, c6, c12, rcut)
        e_ref += u
        f_ref[i] += dudr * d / r
        f_ref[j] -= dudr * d / r

    assert energy == pytest.approx(e_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL)
    np.testing.assert_allclose(forces, f_ref, rtol=FORCE_RTOL, atol=FORCE_ATOL)
    # Bond (1,2) at r=1.166 is outside rcut=0.9's reach only because it is the
    # second bond's rcut=1.5 that applies; both bonds contribute here.
    assert e_ref != 0.0


def test_contacts_lj_fixed_cutoff_baseline() -> None:
    """Legacy functype 1 contacts keep the fixed 1.1 nm cutoff."""
    interaction = ContactsLJ()
    interaction.add_interaction(["1", "2", "1", str(SIGMA), str(EPSILON)])
    c6, c12 = c6_c12_from_sigma_eps(SIGMA, EPSILON)
    assert_pair_matches_analytic(
        interaction.mm_force,
        lambda r, rc: lj_6_12(r, c6, c12, 1.1),
        1.1,
    )


def _multi_contacts_force(state: str, entries: list[list[str]]) -> mm.CustomCompoundBondForce:
    interaction = MultiAllBonds()
    for fields in entries:
        assert interaction.add_interaction(state, "multi_contacts", fields)
    return interaction.mm_force


def _eval_compound_force(
    force: mm.CustomCompoundBondForce, positions: list[tuple[float, float, float]]
) -> tuple[float, np.ndarray]:
    return evaluate(make_system(force, positions))


@pytest.mark.parametrize("rcut", RCUTS)
def test_multi_allbonds_contacts_6_12_matches_analytic(rcut: float) -> None:
    """MultiAllBonds multi_contacts functype 2: same shifted 6-12 form."""
    fields = ["1", "2", "2", "1", "2", str(SIGMA), str(EPSILON), str(rcut)]
    c6, c12 = c6_c12_from_sigma_eps(SIGMA, EPSILON)

    for r in DISTANCES + [rcut - 1e-9, rcut, rcut + 1e-9]:
        # Rebuilt per evaluation: a Force object cannot be owned by two Systems.
        force = _multi_contacts_force("1", [fields])
        energy, forces = _eval_compound_force(force, [(0.0, 0.0, 0.0), (r, 0.0, 0.0)])
        u_ref, dudr_ref = lj_6_12(r, c6, c12, rcut)
        assert energy == pytest.approx(u_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL)
        np.testing.assert_allclose(
            forces[0], [dudr_ref, 0.0, 0.0], rtol=FORCE_RTOL, atol=FORCE_ATOL
        )
        np.testing.assert_allclose(
            forces[1], [-dudr_ref, 0.0, 0.0], rtol=FORCE_RTOL, atol=FORCE_ATOL
        )
        if r > rcut:
            assert energy == 0.0
            assert np.abs(forces).max() == 0.0


def test_multi_allbonds_state_filtering() -> None:
    """Fields for state 2 are not collected by a state-1 MultiAllBonds."""
    interaction = MultiAllBonds()
    fields = ["1", "2", "2", "2", "2", str(SIGMA), str(EPSILON), "1.5"]
    assert not interaction.add_interaction("1", "multi_contacts", fields)


@pytest.mark.parametrize("rcut", RCUTS)
def test_multi_allbonds_contacts_10_12_matches_analytic(rcut: float) -> None:
    """MultiAllBonds multi_contacts functype 3: same shifted 10-12 form as
    Contacts10_12, evaluated through the per-state CustomCompoundBondForce."""
    fields = ["1", "2", "2", "1", "3", str(SIGMA), str(EPSILON), str(rcut)]
    c10 = 4 * EPSILON * SIGMA**10
    c12 = 4 * EPSILON * SIGMA**12

    for r in DISTANCES + [rcut - 1e-9, rcut, rcut + 1e-9]:
        # Rebuilt per evaluation: a Force object cannot be owned by two Systems.
        force = _multi_contacts_force("1", [fields])
        energy, forces = _eval_compound_force(force, [(0.0, 0.0, 0.0), (r, 0.0, 0.0)])
        u_ref, dudr_ref = lj_10_12(r, c10, c12, rcut)
        assert energy == pytest.approx(u_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL)
        np.testing.assert_allclose(
            forces[0], [dudr_ref, 0.0, 0.0], rtol=FORCE_RTOL, atol=FORCE_ATOL
        )
        np.testing.assert_allclose(
            forces[1], [-dudr_ref, 0.0, 0.0], rtol=FORCE_RTOL, atol=FORCE_ATOL
        )
        if r > rcut:
            assert energy == 0.0
            assert np.abs(forces).max() == 0.0
        elif r == rcut:
            assert energy == pytest.approx(0.0, abs=1e-12)


def test_multi_allbonds_10_12_bond_reaches_built_force() -> None:
    """Regression: a collected functype-3 contact must appear in the built force
    (previously 'contact_1012' was missing from MultiAllBonds._handler_order, so
    the bond was collected but silently dropped, yielding an empty '0;' force)."""
    interaction = MultiAllBonds()
    fields = ["1", "2", "2", "1", "3", str(SIGMA), str(EPSILON), "1.5"]
    assert interaction.add_interaction("1", "multi_contacts", fields)

    force = interaction.mm_force
    assert force.getEnergyFunction() != "0;"
    assert "delta_contact_1012" in force.getEnergyFunction()
    assert force.getNumBonds() == 1
    param_names = [
        force.getPerBondParameterName(i) for i in range(force.getNumPerBondParameters())
    ]
    assert "C10" in param_names


# ---------------------------------------------------------------------------
# Mixed contact functypes in one MultiAllBonds (parameter name disambiguation)
# ---------------------------------------------------------------------------

# Per-bond parameter names are shared between the contact handlers
# (C12/C6/rcut). When several functypes mix in one MultiAllBonds the first
# handler in fixed order keeps the canonical names and later colliding
# handlers are suffixed with their handler name; single-type builds must keep
# the canonical names exactly (runtime-cutoff consumers rely on "rcut").
EXPECTED_SINGLE_PARAMS = {
    "1": ["delta_contact", "C12", "C6"],
    "2": ["delta_contact_adj", "C12", "C6", "rcut"],
    "3": ["delta_contact_1012", "C12", "C10", "rcut"],
}
EXPECTED_MIXED_PARAMS = {
    "1+2": [
        "delta_contact", "C12", "C6",
        "delta_contact_adj", "C12_contact_lj_adj", "C6_contact_lj_adj", "rcut",
    ],
    "2+3": [
        "delta_contact_adj", "C12", "C6", "rcut",
        "delta_contact_1012", "C12_contact_1012", "C10", "rcut_contact_1012",
    ],
    "1+2+3": [
        "delta_contact", "C12", "C6",
        "delta_contact_adj", "C12_contact_lj_adj", "C6_contact_lj_adj", "rcut",
        "delta_contact_1012", "C12_contact_1012", "C10", "rcut_contact_1012",
    ],
}
# bond_specs: (0-based i, 0-based j, functype, sigma, epsilon, rcut or None).
# The (0, 2) pair sits beyond every cutoff in both geometries (exact zero);
# the other pairs probe inside-cutoff contributions for every functype.
MIXED_COMBOS = {
    "1+2": [(0, 1, "1", 0.60, 10.0, None), (1, 2, "2", 0.55, 15.0, 1.5)],
    "2+3": [(0, 1, "2", 0.55, 15.0, 1.5), (1, 2, "3", 0.50, 12.0, 1.3)],
    "1+2+3": [
        (0, 1, "1", 0.60, 10.0, None),
        (1, 2, "2", 0.55, 15.0, 1.5),
        (0, 2, "3", 0.50, 12.0, 1.3),
    ],
}
MIXED_GEOMETRIES = [
    [(0.0, 0.0, 0.0), (0.7, 0.2, 0.0), (1.3, 0.9, 0.4)],
    [(0.0, 0.0, 0.0), (0.9, 0.3, -0.2), (1.6, 0.5, 0.1)],
]
DELTA_OF_FUNCTYPE = {"1": "delta_contact", "2": "delta_contact_adj", "3": "delta_contact_1012"}
HANDLER_OF_FUNCTYPE = {"1": "contact_lj", "2": "contact_lj_adj", "3": "contact_1012"}


def _contact_fields(
    i: int, j: int, functype: str, sigma: float, eps: float, rcut: float | None
) -> list[str]:
    """Topology fields for a multi_contacts entry (i, j are 0-based)."""
    fields = [str(i + 1), str(j + 1), "2", "1", functype, str(sigma), str(eps)]
    if rcut is not None:
        fields.append(str(rcut))
    return fields


def _mixed_contacts_force(bond_specs: list[tuple]) -> mm.CustomCompoundBondForce:
    interaction = MultiAllBonds()
    for i, j, functype, sigma, eps, rcut in bond_specs:
        assert interaction.add_interaction(
            "1", "multi_contacts", _contact_fields(i, j, functype, sigma, eps, rcut)
        )
    return interaction.mm_force


def _analytic_mixed_contacts(
    bond_specs: list[tuple], positions: list[tuple[float, float, float]]
) -> tuple[float, np.ndarray]:
    """Independent analytic sum over the mixed contact bonds (energy, forces)."""
    pos = np.array(positions, dtype=float)
    e_ref = 0.0
    f_ref = np.zeros_like(pos)
    for i, j, functype, sigma, eps, rcut in bond_specs:
        c6, c12 = c6_c12_from_sigma_eps(sigma, eps)
        d = pos[j] - pos[i]
        r = float(np.linalg.norm(d))
        if functype == "3":
            u, dudr = lj_10_12(r, 4 * eps * sigma**10, c12, rcut)
        elif functype == "2":
            u, dudr = lj_6_12(r, c6, c12, rcut)
        else:
            u, dudr = lj_6_12(r, c6, c12, 1.1)
        e_ref += u
        f_ref[i] += dudr * d / r
        f_ref[j] -= dudr * d / r
    return e_ref, f_ref


def _param_segment(names: list[str], delta_name: str) -> list[str]:
    """Per-bond parameter names of one handler: its delta flag up to the next."""
    start = names.index(delta_name)
    stop = len(names)
    for k in range(start + 1, len(names)):
        if names[k].startswith("delta_"):
            stop = k
            break
    return names[start:stop]


def _registered_name(segment: list[str], handler: str, canonical: str) -> str:
    """Registered name of a handler's canonical parameter within its segment."""
    for candidate in (canonical, f"{canonical}_{handler}"):
        if candidate in segment:
            return candidate
    raise AssertionError(f"{canonical} not found in segment {segment}")


@pytest.mark.parametrize("functype", ["1", "2", "3"])
def test_multi_allbonds_single_contact_type_param_names_unchanged(functype: str) -> None:
    """Single-contact-type builds keep the canonical parameter names (consumers
    such as the runtime cutoff look up 'rcut' by name)."""
    rcut = None if functype == "1" else 1.3
    force = _mixed_contacts_force([(0, 1, functype, SIGMA, EPSILON, rcut)])
    names = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]
    assert names == EXPECTED_SINGLE_PARAMS[functype]


@pytest.mark.parametrize("combo", ["1+2", "2+3", "1+2+3"])
def test_multi_allbonds_mixed_contacts_param_names(combo: str) -> None:
    """Mixed builds: parameter names are unique, the first contact handler in
    fixed order keeps canonical names, colliding later handlers are suffixed,
    and every bond carries its own handler's values in the right slots."""
    bond_specs = MIXED_COMBOS[combo]
    force = _mixed_contacts_force(bond_specs)
    names = [force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())]

    assert len(set(names)) == len(names)
    assert names == EXPECTED_MIXED_PARAMS[combo]

    expr = force.getEnergyFunction()
    for name in names:
        # every registered parameter is actually referenced by the expression
        assert re.search(rf"\b{re.escape(name)}\b", expr), f"{name} unused in {expr}"

    for k, (_i, _j, functype, sigma, eps, rcut) in enumerate(bond_specs):
        _, params = force.getBondParameters(k)
        values = dict(zip(names, params, strict=True))
        handler = HANDLER_OF_FUNCTYPE[functype]
        segment = _param_segment(names, DELTA_OF_FUNCTYPE[functype])

        assert values[DELTA_OF_FUNCTYPE[functype]] == 1.0
        for other in set(DELTA_OF_FUNCTYPE.values()) - {DELTA_OF_FUNCTYPE[functype]}:
            if other in values:
                assert values[other] == 0.0
        c6, c12 = c6_c12_from_sigma_eps(sigma, eps)
        assert values[_registered_name(segment, handler, "C12")] == pytest.approx(c12)
        if functype == "3":
            assert values[_registered_name(segment, handler, "C10")] == pytest.approx(
                4 * eps * sigma**10
            )
            assert values[_registered_name(segment, handler, "rcut")] == pytest.approx(rcut)
        else:
            assert values[_registered_name(segment, handler, "C6")] == pytest.approx(c6)
        if functype == "2":
            assert values[_registered_name(segment, handler, "rcut")] == pytest.approx(rcut)


@pytest.mark.parametrize("combo", ["1+2", "2+3", "1+2+3"])
def test_multi_allbonds_mixed_contacts_matches_analytic(combo: str) -> None:
    """Mixed force vs independent analytic sum over both geometries: finite
    (no NaN from 0/0 on inactive zero-filled slots) and pointwise correct."""
    bond_specs = MIXED_COMBOS[combo]
    for positions in MIXED_GEOMETRIES:
        force = _mixed_contacts_force(bond_specs)
        energy, forces = _eval_compound_force(force, positions)
        e_ref, f_ref = _analytic_mixed_contacts(bond_specs, positions)

        assert np.isfinite(energy)
        assert np.isfinite(forces).all()
        assert energy == pytest.approx(e_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL)
        np.testing.assert_allclose(forces, f_ref, rtol=FORCE_RTOL, atol=FORCE_ATOL)
        # guards against a vacuous all-zero comparison
        assert e_ref != 0.0
        assert np.abs(f_ref).max() > 0.0


@pytest.mark.parametrize("combo", ["1+2", "2+3", "1+2+3"])
def test_multi_allbonds_mixed_contacts_equals_sum_of_single_builds(combo: str) -> None:
    """The mixed force equals the sum of separate single-functype builds."""
    bond_specs = MIXED_COMBOS[combo]
    for positions in MIXED_GEOMETRIES:
        energy, forces = _eval_compound_force(_mixed_contacts_force(bond_specs), positions)

        e_sum = 0.0
        f_sum = np.zeros((len(positions), 3))
        for spec in bond_specs:
            e_k, f_k = _eval_compound_force(_mixed_contacts_force([spec]), positions)
            e_sum += e_k
            f_sum += f_k

        assert energy == pytest.approx(e_sum, rel=ENERGY_RTOL, abs=ENERGY_ATOL)
        np.testing.assert_allclose(forces, f_sum, rtol=FORCE_RTOL, atol=FORCE_ATOL)
