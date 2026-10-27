"""Full-chain validation of the adjustable Go contact cutoff (``-go-up``).

Starting from the real, contact-map-derived GlnBP topologies (the .map files
are converted to Go contacts by martinize2 upstream; the fixtures keep those
topology products), this module runs the generator -> itp -> create_system ->
OpenMM chain at contact_rcut in {0.9, 1.1, 1.3, 1.5} nm and compares the
resulting contact energies and per-atom forces against an independent analytic
summation that re-parses the written itp file:

- SBP: ``create_sb_topology`` on the GlnBP_Open contacts fixture, both for the
  full solvated 4435-particle system (energy; virtual sites present) and for a
  protein-only no-vsite variant (energy + forces, exact per-atom comparison).
  The legacy fixed-cutoff format (functype 1, no contact_rcut) is checked to
  reproduce the adjustable format at the default 1.1 nm.
- MB: ``create_mb_topology`` on the genuine matched state pair
  (mb_go_martini StateA/Open + StateB/Closed), contacts-only molecule
  (multi_angles/multi_dihedrals/vsites dropped after generation so the MBP
  force is purely the mixed contact potential). Both the per-state
  multi_contacts (functype 2) and the shared contacts paths are checked for
  EXP and HAM mixing, with the analytic mixed forces w1*F1 + w2*F2.

All evaluations are single-point on the Reference platform (no MD).
"""

import math
import os
import re
import shutil

import numpy as np
import openmm as mm
import openmm.unit as unit
import pytest
from openmm import app

from ctgomartini.topology import MartiniTopFile
from ctgomartini.topology.generator.multi_basin import create_mb_topology
from ctgomartini.topology.generator.single_basin import create_sb_topology
from ctgomartini.utils import write_itp
from tests.conftest import WorkingDirectoryContext

FIXTURES_CONTACTS = os.path.join(os.path.dirname(__file__), "../fixtures/contacts")
FIXTURES_MB = os.path.join(os.path.dirname(__file__), "../fixtures/mb_go_martini/GlnBP")

RCUTS = [0.9, 1.1, 1.3, 1.5]
N_PROTEIN_ATOMS = 732

GROUP_CONTACTS = 30
GROUP_MBP = 31

# Measured deviations on the Reference platform: dE <= 9.1e-13 kJ/mol,
# max|dF| <= 5.3e-13 kJ/mol/nm; tolerances carry a >100x margin.
ENERGY_RTOL = 1e-12
ENERGY_ATOL = 1e-9
FORCE_ATOL = 1e-9
FORCE_RTOL = 1e-12

# Mixing parameters written into the MB molecule; chosen so both basin weights
# stay in (0.01, 0.99) for the evaluated conformation at every rcut.
MB_PARAMETERS = {
    "EXP": ["True", "exp", "2", "1/300", "-600", "0"],
    "HAM": ["True", "ham", "2", "400", "-400", "0"],
}

# Contacts forces are identified by this fragment of their energy expression;
# the ES self/exclusion forces also contain "rcut=..." but not this fragment.
CONTACT_EXPR_FRAGMENT = "corr = (C12/rcut^12 - C6/rcut^6)"


# --------------------------------------------------------------------------- helpers
def parse_itp_sections(path: str) -> dict[str, list[list[str]]]:
    """Minimal itp parser: section name -> list of whitespace-split rows."""
    sections: dict[str, list[list[str]]] = {}
    current = None
    with open(path) as handle:
        for line in handle:
            line = line.split(";")[0].strip()
            if not line:
                continue
            match = re.match(r"\[\s*(\S+)\s*\]", line)
            if match:
                current = match.group(1)
                sections.setdefault(current, [])
                continue
            if current is not None:
                sections[current].append(line.split())
    return sections


def contact_energy_force_6_12(r: float, c6: float, c12: float, rcut: float) -> tuple[float, float]:
    """Shifted 6-12 contact potential U and dU/dr; both exactly 0 beyond rcut."""
    if r > rcut:
        return 0.0, 0.0
    u = (c12 / r**12 - c6 / r**6) - (c12 / rcut**12 - c6 / rcut**6)
    dudr = -12 * c12 / r**13 + 6 * c6 / r**7
    return u, dudr


def analytic_contacts(
    entries: list[tuple[int, int, float, float, float]], positions: np.ndarray
) -> tuple[float, np.ndarray]:
    """Sum shifted 6-12 contacts: entries are (i, j, sigma, epsilon, rcut),
    0-based atom indices. Returns (energy, per-atom forces)."""
    energy = 0.0
    forces = np.zeros_like(positions)
    for i, j, sigma, eps, rcut in entries:
        c6 = 4 * eps * sigma**6
        c12 = 4 * eps * sigma**12
        delta = positions[j] - positions[i]
        r = float(np.linalg.norm(delta))
        u, dudr = contact_energy_force_6_12(r, c6, c12, rcut)
        energy += u
        fvec = dudr * delta / r  # F_i = dU/dr * (x_j - x_i) / r
        forces[i] += fvec
        forces[j] -= fvec
    return energy, forces


def make_reference_context(
    system: mm.System, positions, box_vectors
) -> mm.Context:
    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    context.setPeriodicBoxVectors(*box_vectors)
    context.computeVirtualSites()
    return context


def group_energy_forces(context: mm.Context, group: int) -> tuple[float, np.ndarray]:
    state = context.getState(getEnergy=True, getForces=True, groups={group})
    energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = np.array(
        state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer),
        dtype=float,
    )
    return energy, forces


def contact_forces_in_system(system: mm.System, adjustable: bool) -> list[mm.CustomBondForce]:
    """Contact CustomBondForces: adjustable ones carry a per-bond 'rcut'
    parameter, the legacy functype-1 force has 'rcut=1.1;' baked in."""
    found = []
    for force in system.getForces():
        if not isinstance(force, mm.CustomBondForce):
            continue
        if CONTACT_EXPR_FRAGMENT not in force.getEnergyFunction():
            continue
        names = [
            force.getPerBondParameterName(k) for k in range(force.getNumPerBondParameters())
        ]
        if ("rcut" in names) == adjustable:
            found.append(force)
    return found


def entries_from_contacts_section(
    rows: list[list[str]],
) -> list[tuple[int, int, float, float, float]]:
    entries = []
    for a1, a2, functype, sigma, eps, rcut in rows:
        assert functype == "2", f"expected adjustable functype 2, got {rows}"
        entries.append((int(a1) - 1, int(a2) - 1, float(sigma), float(eps), float(rcut)))
    return entries


def mixing_from_fields(method: str, fields: list[str]):
    """Return (mix_energy, mix_weights) closures for the parsed multiple_basin row."""
    coupling = eval(fields[3])  # same semantics as the builder ("1/300" etc.)
    basin_offsets = [float(v) for v in fields[4:6]]
    c1, c2 = basin_offsets

    if method == "EXP":
        beta = coupling

        def mix(e1, e2):
            a1 = -beta * (e1 + c1)
            a2 = -beta * (e2 + c2)
            m = max(a1, a2)
            p1, p2 = math.exp(a1 - m), math.exp(a2 - m)
            w1, w2 = p1 / (p1 + p2), p2 / (p1 + p2)
            return -(math.log(p1 + p2) + m) / beta, w1, w2

    else:
        delta = coupling
        dv = c2 - c1

        def mix(e1, e2):
            alpha = (e1 - e2 - dv) / 2
            root = math.sqrt(alpha**2 + delta**2)
            w1 = 0.5 * (1 - alpha / root)
            return (e1 + e2 + dv) / 2 - root, w1, 1.0 - w1

    return mix


def assert_matches(tag, e_omm, e_ref, f_omm, f_ref):
    """Report and assert energy/force agreement."""
    de = abs(e_omm - e_ref)
    df = float(np.abs(f_omm - f_ref).max())
    fmax = float(np.abs(f_ref).max())
    print(f"[{tag}] E_omm={e_omm:.9f} E_ref={e_ref:.9f} dE={de:.2e} "
          f"maxdF={df:.2e} (maxF={fmax:.3e})")
    assert e_omm == pytest.approx(e_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL), tag
    assert df <= FORCE_ATOL + FORCE_RTOL * fmax, tag


# --------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def sbp_work(tmp_path_factory) -> str:
    """Scratch copy of the GlnBP_Open contacts fixture."""
    work = tmp_path_factory.mktemp("sbp") / "Contacts"
    shutil.copytree(os.path.join(FIXTURES_CONTACTS, "GlnBP_Open/Contacts"), work)
    os.rename(os.path.join(work, "system.top"), os.path.join(work, "open.top"))
    return str(work)


@pytest.fixture(scope="module")
def mb_work(tmp_path_factory) -> str:
    """Flat scratch dir with the genuine MB state pair (Open + Closed)."""
    work = str(tmp_path_factory.mktemp("mb"))
    shutil.copytree(os.path.join(FIXTURES_MB, "StateA/Open"), os.path.join(work, "Open"))
    shutil.copytree(os.path.join(FIXTURES_MB, "StateB/Closed"), os.path.join(work, "Closed"))
    for name in ["open", "closed"]:
        mol = name.capitalize()
        with open(os.path.join(work, f"{name}.top"), "w") as handle:
            handle.write(
                "#define GO_VIRT\n"
                f'#include "{mol}/martini_v3.0.0.itp"\n#include "{mol}/{mol}.itp"\n'
                f"[ system ]\nTitle\n[ molecules ]\n{mol} 1\n"
            )
    return work


def _full_gro_data(sbp_work: str):
    gro = app.GromacsGroFile(os.path.join(sbp_work, "minimized.gro"))
    positions = np.array(gro.getPositions().value_in_unit(unit.nanometer), dtype=float)
    return gro, positions


def build_sb_system(
    work: str, rcut: float | None, protein_only: bool, contact_type: int = 2
) -> tuple[mm.System, dict, np.ndarray, object]:
    """create_sb_topology -> write_itp -> MartiniTopFile -> create_system.

    Returns (system, parsed_itp_sections, positions_nm, box_vectors).
    Positions are taken from the fixture minimized.gro (protein block for the
    protein-only variant, which drops virtual_sitesn so forces can be compared
    per atom without virtual-site redistribution).
    """
    gro, positions = _full_gro_data(work)
    box = gro.getPeriodicBoxVectors()
    with WorkingDirectoryContext(work):
        sbmol = create_sb_topology(
            [["open.top", "open"]], "sbrcut", contact_rcut=rcut, contact_type=contact_type
        )
        if protein_only:
            sbmol._topology.pop("virtual_sitesn", None)
        write_itp(sbmol)
        solvent_lines = "" if protein_only else "W 3626\nNA 40\nCL 37\n"
        with open("run.top", "w") as handle:
            handle.write(
                '#include "FF/martini_v3.0.0.itp"\n'
                '#include "sbrcut_params.itp"\n#include "sbrcut.itp"\n'
            )
            if not protein_only:
                handle.write(
                    '#include "FF/martini_v3.0.0_phospholipids_v1.itp"\n'
                    '#include "FF/martini_v3.0.0_solvents_v1.itp"\n'
                    '#include "FF/martini_v3.0.0_ions_v1.itp"\n'
                )
            handle.write(f"[ system ]\nTitle\n[ molecules ]\nsbrcut 1\n{solvent_lines}")
        top = MartiniTopFile("run.top", periodicBoxVectors=box)
        system = top.create_system(nonbonded_cutoff=1.1 * unit.nanometer, epsilon_r=15.0)
        sections = parse_itp_sections("sbrcut.itp")

    positions = positions[:N_PROTEIN_ATOMS] if protein_only else positions
    return system, sections, positions, box


def build_mb_system(
    work: str, method: str, rcut: float, cutoff_contacts: float
) -> tuple[mm.System, dict, np.ndarray, object]:
    """create_mb_topology -> contacts-only surgery -> write_itp -> create_system.

    The surgery keeps only the sections needed for the contact chain
    (atoms/multiple_basin/contacts/multi_contacts), so the MBP CustomCVForce
    mixes exactly the per-state contact energies and an independent analytic
    summation over the written itp is tractable.
    """
    pdb = app.PDBFile(os.path.join(work, "Open/Open_cg.pdb"))
    positions = np.array(
        pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float
    )
    assert len(positions) == N_PROTEIN_ATOMS
    # Any box comfortably larger than the protein works; reuse the SBP one.
    box = app.GromacsGroFile(
        os.path.join(FIXTURES_CONTACTS, "GlnBP_Open/Contacts/minimized.gro")
    ).getPeriodicBoxVectors()

    keep = {"moleculetype", "atoms", "multiple_basin", "contacts", "multi_contacts"}
    with WorkingDirectoryContext(work):
        mbmol = create_mb_topology(
            [["open.top", "Open"], ["closed.top", "Closed"]],
            "mbmol",
            dict_cutoffs={
                "cutoff_BBB_angles": 0.0,
                "cutoff_BBBB_dihedrals": 0.0,
                "cutoff_SBBS_dihedrals": 0.0,
                "cutoff_contacts": cutoff_contacts,
            },
            contact_rcut=rcut,
        )
        for key in list(mbmol._topology.keys()):
            if key not in keep:
                del mbmol._topology[key]
        mbmol._topology["multiple_basin"] = [MB_PARAMETERS[method]]
        write_itp(mbmol)
        with open("run.top", "w") as handle:
            handle.write(
                '#include "Open/martini_v3.0.0.itp"\n'
                '#include "mbmol_params.itp"\n#include "mbmol.itp"\n'
                "[ system ]\nTitle\n[ molecules ]\nmbmol 1\n"
            )
        top = MartiniTopFile("run.top", periodicBoxVectors=box)
        system = top.create_system(nonbonded_cutoff=1.1 * unit.nanometer, epsilon_r=15.0)
        sections = parse_itp_sections("mbmol.itp")

    assert system.getNumParticles() == N_PROTEIN_ATOMS
    return system, sections, positions, box


# --------------------------------------------------------------------------- SBP tests
@pytest.mark.parametrize("rcut", RCUTS)
def test_sb_chain_full_system_energy(sbp_work: str, rcut: float) -> None:
    """Full solvated system (4435 particles, virtual sites present): contact
    group energy vs analytic summation over the written itp."""
    system, sections, positions, box = build_sb_system(sbp_work, rcut, protein_only=False)
    contacts = contact_forces_in_system(system, adjustable=True)
    assert len(contacts) == 1
    contacts[0].setForceGroup(GROUP_CONTACTS)

    context = make_reference_context(system, positions * unit.nanometer, box)
    e_omm, _ = group_energy_forces(context, GROUP_CONTACTS)

    entries = entries_from_contacts_section(sections["contacts"])
    assert entries, "contacts section unexpectedly empty"
    assert len(entries) == contacts[0].getNumBonds()
    e_ref, _ = analytic_contacts(entries, positions)

    de = abs(e_omm - e_ref)
    print(f"[SBP-full rcut={rcut}] E_omm={e_omm:.9f} E_ref={e_ref:.9f} dE={de:.2e}")
    assert e_omm == pytest.approx(e_ref, rel=ENERGY_RTOL, abs=ENERGY_ATOL)


@pytest.mark.parametrize("rcut", RCUTS)
def test_sb_chain_protein_energy_forces(sbp_work: str, rcut: float) -> None:
    """Protein-only system without virtual sites: exact per-atom force check."""
    system, sections, positions, box = build_sb_system(sbp_work, rcut, protein_only=True)
    assert system.getNumParticles() == N_PROTEIN_ATOMS
    contacts = contact_forces_in_system(system, adjustable=True)
    assert len(contacts) == 1
    contacts[0].setForceGroup(GROUP_CONTACTS)

    context = make_reference_context(system, positions * unit.nanometer, box)
    e_omm, f_omm = group_energy_forces(context, GROUP_CONTACTS)

    entries = entries_from_contacts_section(sections["contacts"])
    assert entries, "contacts section unexpectedly empty"
    e_ref, f_ref = analytic_contacts(entries, positions)
    assert_matches(f"SBP-prot rcut={rcut}", e_omm, e_ref, f_omm, f_ref)


def test_sb_chain_legacy_fixed_cutoff_equivalence(sbp_work: str) -> None:
    """Legacy functype 1 (fixed 1.1 nm) and adjustable functype 2 at rcut=1.1
    must give the same contact energy and forces."""
    system_adj, sections, positions, box = build_sb_system(sbp_work, 1.1, protein_only=True)
    system_leg, _, _, _ = build_sb_system(
        sbp_work, None, protein_only=True, contact_type=1
    )

    adjustable = contact_forces_in_system(system_adj, adjustable=True)
    assert any(f.getNumBonds() > 0 for f in adjustable)
    assert sections["contacts"], "contacts section unexpectedly empty"
    for contacts in adjustable:
        contacts.setForceGroup(GROUP_CONTACTS)
    legacy = contact_forces_in_system(system_leg, adjustable=False)
    assert any(f.getNumBonds() > 0 for f in legacy)
    for force in legacy:
        force.setForceGroup(GROUP_CONTACTS)

    ctx_adj = make_reference_context(system_adj, positions * unit.nanometer, box)
    ctx_leg = make_reference_context(system_leg, positions * unit.nanometer, box)
    e_adj, f_adj = group_energy_forces(ctx_adj, GROUP_CONTACTS)
    e_leg, f_leg = group_energy_forces(ctx_leg, GROUP_CONTACTS)

    de = abs(e_adj - e_leg)
    df = float(np.abs(f_adj - f_leg).max())
    print(f"[SBP legacy-vs-adj@1.1] E={e_leg:.9f} dE={de:.2e} maxdF={df:.2e}")
    assert e_adj == pytest.approx(e_leg, rel=ENERGY_RTOL, abs=ENERGY_ATOL)
    np.testing.assert_allclose(f_adj, f_leg, rtol=FORCE_RTOL, atol=FORCE_ATOL)


# --------------------------------------------------------------------------- MB tests
def _check_mb_chain(work: str, method: str, rcut: float, cutoff_contacts: float) -> None:
    system, sections, positions, box = build_mb_system(work, method, rcut, cutoff_contacts)

    mbp_forces = [
        f for f in system.getForces()
        if isinstance(f, mm.CustomCVForce) and "ctg_" in f.getEnergyFunction()
    ]
    assert len(mbp_forces) == 1
    mbp_forces[0].setForceGroup(GROUP_MBP)
    shared = contact_forces_in_system(system, adjustable=True)
    for force in shared:
        force.setForceGroup(GROUP_CONTACTS)

    context = make_reference_context(system, positions * unit.nanometer, box)

    # Shared contacts (present in both states) sit in their own force.
    shared_rows = sections.get("contacts", [])
    assert (len(shared_rows) > 0) == (len(shared) > 0)
    # Non-empty guards: if contacts were silently lost, both sides of the
    # comparisons below would be zero and the test would vacuously pass.
    if cutoff_contacts > 0:
        assert shared_rows, "cutoff_contacts>0 must produce shared contacts"
    else:
        assert not shared_rows, "cutoff_contacts=0 must route all contacts to multi_contacts"
    if shared:
        e_omm, f_omm = group_energy_forces(context, GROUP_CONTACTS)
        e_ref, f_ref = analytic_contacts(entries_from_contacts_section(shared_rows), positions)
        assert_matches(
            f"MB-{method} rcut={rcut} cc={cutoff_contacts} shared({len(shared_rows)})",
            e_omm, e_ref, f_omm, f_ref,
        )

    # Per-state contacts are mixed by the MBP CustomCVForce.
    per_state: dict[int, list] = {1: [], 2: []}
    for a1, a2, _nstates, state, functype, sigma, eps, rcut_ in sections.get("multi_contacts", []):
        assert functype == "2"
        per_state[int(state)].append(
            (int(a1) - 1, int(a2) - 1, float(sigma), float(eps), float(rcut_))
        )
    assert per_state[1], "state-1 multi_contacts unexpectedly empty"
    assert per_state[2], "state-2 multi_contacts unexpectedly empty"
    e1, f1 = analytic_contacts(per_state[1], positions)
    e2, f2 = analytic_contacts(per_state[2], positions)

    mb_rows = sections["multiple_basin"]
    assert mb_rows[0][1].upper() == method
    mix = mixing_from_fields(method, mb_rows[0])
    e_mix, w1, w2 = mix(e1, e2)
    assert 0.01 < w2 < 0.99, f"mixing weights saturated: w2={w2}"

    e_omm, f_omm = group_energy_forces(context, GROUP_MBP)
    f_mix = w1 * f1 + w2 * f2
    print(f"    [MB-{method} rcut={rcut} cc={cutoff_contacts}] "
          f"multi={len(per_state[1]) + len(per_state[2])} "
          f"E1={e1:.4f} E2={e2:.4f} w1={w1:.4f} w2={w2:.4f}")
    assert_matches(
        f"MB-{method} rcut={rcut} cc={cutoff_contacts} mixed", e_omm, e_mix, f_omm, f_mix
    )


@pytest.mark.parametrize("rcut", RCUTS)
@pytest.mark.parametrize("method", ["EXP", "HAM"])
def test_mb_chain_multi_contacts(mb_work: str, method: str, rcut: float) -> None:
    """All contacts land in multi_contacts when the combination cutoff is 0."""
    _check_mb_chain(mb_work, method, rcut, cutoff_contacts=0.0)


@pytest.mark.parametrize("rcut", [1.1, 1.5])
@pytest.mark.parametrize("method", ["EXP", "HAM"])
def test_mb_chain_shared_and_multi_contacts(mb_work: str, method: str, rcut: float) -> None:
    """cutoff_contacts=0.06 moves similar-sigma pairs into the shared contacts
    section; both the shared and the per-state paths are checked."""
    _check_mb_chain(mb_work, method, rcut, cutoff_contacts=0.06)
