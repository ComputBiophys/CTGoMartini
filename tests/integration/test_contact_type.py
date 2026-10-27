"""Full-chain validation of the ``contact_type`` generator switch.

Uses the same GlnBP fixtures and OpenMM-vs-analytic methodology as
``test_contact_rcut`` (generator -> itp -> create_system -> Reference
platform single point):

- SBP with contact_type=3 (10-12 LJ): protein-only energy/forces vs the
  independent analytic 10-12 sum, at rcut 1.1 and 1.5 nm.
- MB (EXP) with contact_type=3: both the shared contacts (functype 3) and
  the per-state multi_contacts (functype 3, mixed by the MBP CustomCVForce)
  vs analytic, at rcut 1.1 and 1.5 nm.
- SBP contact_type=1 (fixed 1.1 nm) vs contact_type=2 at rcut=1.1: same
  physics through different forces - forces bitwise identical, energies
  equal to within 1 ulp on the Reference platform.
"""

import os
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
from tests.integration.test_contact_rcut import (
    FIXTURES_CONTACTS,
    FIXTURES_MB,
    GROUP_CONTACTS,
    GROUP_MBP,
    MB_PARAMETERS,
    N_PROTEIN_ATOMS,
    assert_matches,
    contact_forces_in_system,
    group_energy_forces,
    make_reference_context,
    mixing_from_fields,
    parse_itp_sections,
)

RCUTS = [1.1, 1.5]

# Measured deviations on the Reference platform: dE <= 4.6e-13 kJ/mol,
# max|dF| <= 4.1e-13 kJ/mol/nm; tolerances carry a >100x margin.
ENERGY_RTOL = 1e-12
ENERGY_ATOL = 1e-9
FORCE_ATOL = 1e-9
FORCE_RTOL = 1e-12

CONTACT_EXPR_FRAGMENT_10_12 = "corr = (C12/rcut^12 - C10/rcut^10)"


def contact_energy_force_10_12(
    r: float, c10: float, c12: float, rcut: float
) -> tuple[float, float]:
    """Shifted 10-12 contact potential U and dU/dr; both exactly 0 beyond rcut."""
    if r > rcut:
        return 0.0, 0.0
    u = (c12 / r**12 - c10 / r**10) - (c12 / rcut**12 - c10 / rcut**10)
    dudr = -12 * c12 / r**13 + 10 * c10 / r**11
    return u, dudr


def analytic_typed_contacts(
    entries: list[tuple[int, int, str, float, float, float]], positions: np.ndarray
) -> tuple[float, np.ndarray]:
    """Sum shifted contacts of mixed functypes: entries are
    (i, j, functype, sigma, epsilon, rcut) with 0-based atom indices
    (functype '1' implies rcut=1.1). Returns (energy, per-atom forces)."""
    energy = 0.0
    forces = np.zeros_like(positions)
    for i, j, functype, sigma, eps, rcut in entries:
        c12 = 4 * eps * sigma**12
        delta = positions[j] - positions[i]
        r = float(np.linalg.norm(delta))
        if functype == "3":
            u, dudr = contact_energy_force_10_12(r, 4 * eps * sigma**10, c12, rcut)
        else:
            c6 = 4 * eps * sigma**6
            u, dudr = contact_energy_force_6_12_typed(r, c6, c12, rcut)
        energy += u
        fvec = dudr * delta / r
        forces[i] += fvec
        forces[j] -= fvec
    return energy, forces


def contact_energy_force_6_12_typed(
    r: float, c6: float, c12: float, rcut: float
) -> tuple[float, float]:
    """Shifted 6-12 contact potential U and dU/dr; both exactly 0 beyond rcut."""
    if r > rcut:
        return 0.0, 0.0
    u = (c12 / r**12 - c6 / r**6) - (c12 / rcut**12 - c6 / rcut**6)
    dudr = -12 * c12 / r**13 + 6 * c6 / r**7
    return u, dudr


def entries_from_typed_contacts_section(
    rows: list[list[str]],
) -> list[tuple[int, int, str, float, float, float]]:
    entries = []
    for fields in rows:
        a1, a2, functype, sigma, eps = fields[:5]
        rcut = float(fields[5]) if len(fields) > 5 else 1.1
        entries.append((int(a1) - 1, int(a2) - 1, functype, float(sigma), float(eps), rcut))
    return entries


def contact_forces_10_12_in_system(system) -> list:
    """functype-3 contact CustomBondForces carry the C10 correction term."""
    return [
        force
        for force in system.getForces()
        if isinstance(force, mm.CustomBondForce)
        and CONTACT_EXPR_FRAGMENT_10_12 in force.getEnergyFunction()
    ]


@pytest.fixture(scope="module")
def sbp_work_ct(tmp_path_factory) -> str:
    work = tmp_path_factory.mktemp("sbp_ct") / "Contacts"
    shutil.copytree(os.path.join(FIXTURES_CONTACTS, "GlnBP_Open/Contacts"), work)
    os.rename(os.path.join(work, "system.top"), os.path.join(work, "open.top"))
    return str(work)


@pytest.fixture(scope="module")
def mb_work_ct(tmp_path_factory) -> str:
    work = str(tmp_path_factory.mktemp("mb_ct"))
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


def build_sb_system_typed(work: str, rcut: float | None, contact_type: int):
    """Protein-only SBP chain (no virtual sites) with an explicit contact_type."""
    gro = app.GromacsGroFile(os.path.join(work, "minimized.gro"))
    positions = np.array(gro.getPositions().value_in_unit(unit.nanometer), dtype=float)
    box = gro.getPeriodicBoxVectors()
    with WorkingDirectoryContext(work):
        sbmol = create_sb_topology(
            [["open.top", "open"]], "sbct", contact_rcut=rcut, contact_type=contact_type
        )
        sbmol._topology.pop("virtual_sitesn", None)
        write_itp(sbmol)
        with open("run.top", "w") as handle:
            handle.write(
                '#include "FF/martini_v3.0.0.itp"\n'
                '#include "sbct_params.itp"\n#include "sbct.itp"\n'
                "[ system ]\nTitle\n[ molecules ]\nsbct 1\n"
            )
        top = MartiniTopFile("run.top", periodicBoxVectors=box)
        system = top.create_system(nonbonded_cutoff=1.1 * unit.nanometer, epsilon_r=15.0)
        sections = parse_itp_sections("sbct.itp")
    return system, sections, positions[:N_PROTEIN_ATOMS], box


def build_mb_system_typed(
    work: str, method: str, rcut: float, cutoff_contacts: float, contact_type: int
):
    """Contacts-only MB chain with an explicit contact_type (see test_contact_rcut)."""
    pdb = app.PDBFile(os.path.join(work, "Open/Open_cg.pdb"))
    positions = np.array(
        pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer), dtype=float
    )
    assert len(positions) == N_PROTEIN_ATOMS
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
            contact_type=contact_type,
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


# --------------------------------------------------------------------------- SBP type 3 (10-12)
@pytest.mark.parametrize("rcut", RCUTS)
def test_sb_type3_10_12_protein_energy_forces(sbp_work_ct: str, rcut: float) -> None:
    system, sections, positions, box = build_sb_system_typed(sbp_work_ct, rcut, 3)
    assert system.getNumParticles() == N_PROTEIN_ATOMS
    contacts = contact_forces_10_12_in_system(system)
    assert len(contacts) == 1
    contacts[0].setForceGroup(GROUP_CONTACTS)

    context = make_reference_context(system, positions * unit.nanometer, box)
    e_omm, f_omm = group_energy_forces(context, GROUP_CONTACTS)

    rows = sections["contacts"]
    assert rows, "contacts section unexpectedly empty"
    assert {row[2] for row in rows} == {"3"}
    entries = entries_from_typed_contacts_section(rows)
    assert len(entries) == contacts[0].getNumBonds()
    e_ref, f_ref = analytic_typed_contacts(entries, positions)
    assert_matches(f"SBP-10-12 rcut={rcut}", e_omm, e_ref, f_omm, f_ref)


# ------------------------------------------------------- SBP type 1 vs type 2 @1.1
def test_sb_type1_matches_type2_at_1_1(sbp_work_ct: str) -> None:
    """Fixed-cutoff functype 1 and adjustable functype 2 at rcut=1.1 give the
    same contact physics: forces bitwise identical, energy equal to 1 ulp."""
    system_t2, sections_t2, positions, box = build_sb_system_typed(sbp_work_ct, 1.1, 2)
    system_t1, sections_t1, _, _ = build_sb_system_typed(sbp_work_ct, None, 1)

    assert {row[2] for row in sections_t2["contacts"]} == {"2"}
    assert {row[2] for row in sections_t1["contacts"]} == {"1"}
    assert [row[:2] + row[3:5] for row in sections_t1["contacts"]] == [
        row[:2] + row[3:5] for row in sections_t2["contacts"]
    ]

    for force in contact_forces_in_system(system_t2, adjustable=True):
        force.setForceGroup(GROUP_CONTACTS)
    for force in contact_forces_in_system(system_t1, adjustable=False):
        force.setForceGroup(GROUP_CONTACTS)

    ctx_t2 = make_reference_context(system_t2, positions * unit.nanometer, box)
    ctx_t1 = make_reference_context(system_t1, positions * unit.nanometer, box)
    e_t2, f_t2 = group_energy_forces(ctx_t2, GROUP_CONTACTS)
    e_t1, f_t1 = group_energy_forces(ctx_t1, GROUP_CONTACTS)

    de = abs(e_t2 - e_t1)
    print(f"[SBP type1-vs-type2@1.1] E1={e_t1:.12f} dE={de:.2e}")
    np.testing.assert_array_equal(f_t2, f_t1)
    assert de <= 1e-12


# ------------------------------------------------- MB type 3 (10-12), EXP mixing
@pytest.mark.parametrize("rcut", RCUTS)
def test_mb_type3_10_12_exp_energy_forces(mb_work_ct: str, rcut: float) -> None:
    """cutoff_contacts=0.06: shared functype-3 contacts in their own force and
    per-state functype-3 multi_contacts mixed by the MBP CustomCVForce."""
    system, sections, positions, box = build_mb_system_typed(
        mb_work_ct, "EXP", rcut, cutoff_contacts=0.06, contact_type=3
    )

    mbp_forces = [
        f
        for f in system.getForces()
        if isinstance(f, mm.CustomCVForce) and "ctg_" in f.getEnergyFunction()
    ]
    assert len(mbp_forces) == 1
    mbp_forces[0].setForceGroup(GROUP_MBP)
    shared = contact_forces_10_12_in_system(system)
    for force in shared:
        force.setForceGroup(GROUP_CONTACTS)

    context = make_reference_context(system, positions * unit.nanometer, box)

    shared_rows = sections.get("contacts", [])
    assert shared_rows, "cutoff_contacts=0.06 must produce shared contacts"
    assert {row[2] for row in shared_rows} == {"3"}
    e_omm, f_omm = group_energy_forces(context, GROUP_CONTACTS)
    e_ref, f_ref = analytic_typed_contacts(
        entries_from_typed_contacts_section(shared_rows), positions
    )
    assert_matches(
        f"MB-EXP-10-12 rcut={rcut} shared({len(shared_rows)})", e_omm, e_ref, f_omm, f_ref
    )

    per_state: dict[int, list] = {1: [], 2: []}
    for row in sections.get("multi_contacts", []):
        a1, a2, _nstates, state, functype, sigma, eps, rcut_ = row
        assert functype == "3"
        per_state[int(state)].append(
            (int(a1) - 1, int(a2) - 1, functype, float(sigma), float(eps), float(rcut_))
        )
    assert per_state[1], "state-1 multi_contacts unexpectedly empty"
    assert per_state[2], "state-2 multi_contacts unexpectedly empty"
    e1, f1 = analytic_typed_contacts(per_state[1], positions)
    e2, f2 = analytic_typed_contacts(per_state[2], positions)

    mb_rows = sections["multiple_basin"]
    mix = mixing_from_fields("EXP", mb_rows[0])
    e_mix, w1, w2 = mix(e1, e2)
    assert 0.01 < w2 < 0.99, f"mixing weights saturated: w2={w2}"

    e_omm, f_omm = group_energy_forces(context, GROUP_MBP)
    print(f"    [MB-EXP-10-12 rcut={rcut}] E1={e1:.4f} E2={e2:.4f} w1={w1:.4f} w2={w2:.4f}")
    assert_matches(
        f"MB-EXP-10-12 rcut={rcut} mixed", e_omm, e_mix, f_omm, w1 * f1 + w2 * f2
    )
