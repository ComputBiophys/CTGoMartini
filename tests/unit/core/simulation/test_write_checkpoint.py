"""Tests for XML checkpoint writing in ctgomartini.simulation.base."""

from __future__ import annotations

from pathlib import Path

import openmm as mm
import openmm.unit as unit
from openmm.app import Simulation, Topology

from ctgomartini.simulation.base import write_checkpoint
from ctgomartini.topology.interactions import EXPInteraction


def _two_atom_topology() -> Topology:
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue("MOL", chain)
    topology.addAtom("B1", None, residue)
    topology.addAtom("B2", None, residue)
    return topology


def _bond_force(k: float) -> mm.HarmonicBondForce:
    """Harmonic bond between particles 0 and 1, r0 = 1.0 nm."""
    force = mm.HarmonicBondForce()
    force.addBond(0, 1, 1.0, k)
    return force


def test_write_checkpoint_preserves_runtime_global_parameters(tmp_path: Path) -> None:
    mbp_force = EXPInteraction().addForce(
        {'1': [_bond_force(100.0)], '2': [_bond_force(200.0)]},
        5.0,
        [10.0, 20.0],
        molname='GBP',
    )
    system = mm.System()
    system.addParticle(1.0)
    system.addParticle(1.0)
    system.addForce(mbp_force)

    platform = mm.Platform.getPlatformByName("Reference")
    simulation = Simulation(
        _two_atom_topology(), system, mm.VerletIntegrator(0.001), platform
    )
    simulation.context.setPositions(
        [mm.Vec3(0, 0, 0), mm.Vec3(0.3, 0, 0)] * unit.nanometer
    )
    simulation.context.setParameter('ctg_GBP_beta', 2.0)
    simulation.context.setParameter('ctg_GBP_C1', 30.0)

    checkpoint_file = tmp_path / "state.chk"
    write_checkpoint(simulation, str(checkpoint_file))

    # A Context rebuilt from the same System starts at the topology defaults;
    # the values adjusted at runtime must be restored from the checkpoint.
    new_context = mm.Context(system, mm.VerletIntegrator(0.001), platform)
    assert new_context.getParameter('ctg_GBP_beta') == 5.0
    with open(checkpoint_file) as f:
        new_context.setState(mm.XmlSerializer.deserialize(f.read()))

    assert new_context.getParameter('ctg_GBP_beta') == 2.0
    assert new_context.getParameter('ctg_GBP_C1') == 30.0
    assert new_context.getParameter('ctg_GBP_C2') == 20.0
