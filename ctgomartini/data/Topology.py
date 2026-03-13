#!/usr/bin/env python3
"""Topology and simulation setup utilities for CTGoMartini.

This module provides functions for loading molecular structures, generating
restraints, configuring OpenMM platforms, and setting up molecular dynamics
simulations using the Martini force field.

Authors: Song Yang
Date: 20250505
"""

from __future__ import annotations

import argparse
import datetime
import os
import signal
import sys
import warnings
from typing import TYPE_CHECKING, Any

import MDAnalysis as mda
import numpy as np
import openmm as mm
import openmm.unit as u
from openmm.app import GromacsGroFile, PDBFile, Simulation

from ctgomartini.api import MartiniTopFile
from ctgomartini.utils import read_inputs

if TYPE_CHECKING:
    from openmm import System
    from openmm.app import Topology

warnings.filterwarnings("ignore")

# Type alias for input parameters object
Inputs = Any


def report_time(start_time: datetime.datetime) -> None:
    """Report elapsed time between start and current time.

    Args:
        start_time: The starting datetime to calculate elapsed time from.
    """
    end_time = datetime.datetime.now()
    elapsed = (end_time - start_time).total_seconds()

    start_time_str = start_time.__format__("%Y-%m-%d %H:%M:%S")
    end_time_str = end_time.__format__("%Y-%m-%d %H:%M:%S")
    print(f"  Start Time: {start_time_str}")
    print(f"    End Time: {end_time_str}")
    print(f"Elapsed Time: {elapsed:.2f}")


def load_structure(str_file: str) -> tuple[GromacsGroFile | PDBFile, mm.Vec3 | None]:
    """Load molecular structure from file.

    Args:
        str_file: Path to structure file (.gro or .pdb).

    Returns:
        A tuple containing:
            - The loaded structure object (GromacsGroFile or PDBFile).
            - The periodic box vectors, or None if not found.

    Raises:
        Exception: If the file format is not supported.
    """
    if str_file.split(".")[-1] == "gro":
        conf = GromacsGroFile(str_file)
        box_vectors = conf.getPeriodicBoxVectors()
    elif str_file.split(".")[-1] == "pdb":
        conf = PDBFile(str_file)
        box_vectors = conf.getTopology().getPeriodicBoxVectors()
    else:
        raise Exception("Unsupported structure file: ", str_file)
    return conf, box_vectors


def gen_restraints(
    str_file: str,
    atomname: str,
    fc: float = 1000.0,
    rest_file: str = "restraints.txt",
) -> None:
    """Generate restraint file for selected atoms.

    Args:
        str_file: Structure file (gro/pdb) containing atoms to restrain.
        atomname: Atom name to apply restraints to (e.g., 'CA').
        fc: Force constant in kJ/(mol·nm²). Defaults to 1000.
        rest_file: Output restraint file path. Defaults to "restraints.txt".
    """
    universe = mda.Universe(str_file)
    sel = universe.select_atoms(f"name {atomname}")

    newlines = ["; atomindex functype(1) fc_x fc_y fc_z\n"]
    newlines.append("; atomid start from 1\n")
    for i in sel.indices:
        newlines.append(f"{i+1:>5} 1 {fc} {fc} {fc}\n")

    with open(rest_file, "w") as g:
        g.writelines(newlines)


def add_restraints(system: System, inputs: Inputs) -> System:
    """Add positional restraints to the system.

    Args:
        system: The OpenMM system to add restraints to.
        inputs: Input parameters object containing restraint settings.
            Expected attributes: rest, rest_ref, rest_file.

    Returns:
        The system with restraints added (if enabled).
    """
    crd, _ = load_structure(inputs.rest_ref)
    if inputs.rest == "yes":
        # Create custom external force for anisotropic positional restraints
        posres_prot = mm.CustomExternalForce(
            "1/2*kx*periodicdistance(x, 0, 0, x0, 0, 0)^2 + "
            "1/2*ky*periodicdistance(0, y, 0, 0, y0, 0)^2 + "
            "1/2*kz*periodicdistance(0, 0, z, 0, 0, z0)^2;"
        )
        posres_prot.addPerParticleParameter("kx")
        posres_prot.addPerParticleParameter("ky")
        posres_prot.addPerParticleParameter("kz")
        posres_prot.addPerParticleParameter("x0")
        posres_prot.addPerParticleParameter("y0")
        posres_prot.addPerParticleParameter("z0")

        # Parse restraint file and add particles to the force
        with open(inputs.rest_file, "r") as f:
            for line in f:
                if line.find(";") >= 0:
                    line = line.split(";")[0]
                sline = line.strip()
                if sline == "":
                    continue
                segments, functype, fcx, fcy, fcz = sline.split()[:5]
                atom1 = int(segments) - 1  # Convert to 0-based index
                fcx, fcy, fcz = float(fcx), float(fcy), float(fcz)
                assert (
                    functype == "1"
                ), f"Error: Unsupport position restraint type.\n {line}"

                # Get reference positions in nanometers
                xpos = crd.positions[atom1].value_in_unit(u.nanometers)[0]
                ypos = crd.positions[atom1].value_in_unit(u.nanometers)[1]
                zpos = crd.positions[atom1].value_in_unit(u.nanometers)[2]

                # Add restraint if any force constant is positive
                if fcx >= 0 and fcy >= 0 and fcz >= 0:
                    posres_prot.addParticle(atom1, [fcx, fcy, fcz, xpos, ypos, zpos])

        system.addForce(posres_prot)
    return system


def backup_file(file: str) -> None:
    """Create backup of existing file with incremental numbering.

    Args:
        file: Path to file that needs backup.
    """
    if os.path.isfile(file):
        i = 1
        newfile = file + f".bk{i}"
        while os.path.isfile(newfile):
            i += 1
            newfile = file + f".bk{i}"
        os.rename(file, newfile)


def write_output(output_file: str, simulation: Simulation, strfile: str) -> None:
    """Write simulation output to file in specified format.

    Args:
        output_file: Path to output file.
        simulation: Simulation object containing current state.
        strfile: Reference structure file for topology.
    """
    # Get current state information
    state = simulation.context.getState(
        getPositions=True, getVelocities=True, enforcePeriodicBox=True
    )
    crd = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    velocities = state.getVelocities(asNumpy=True).value_in_unit(
        u.angstrom / u.picosecond
    )
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)[
        [0, 1, 2], [0, 1, 2]
    ]

    # Create MDAnalysis universe and update with current simulation state
    mda_u = mda.Universe(strfile)
    mda_u.atoms.positions = crd
    mda_u.trajectory[0].velocities = True
    mda_u.dimensions[:3] = box_vectors  # only for rectangular/cubic box
    mda_u.atoms.velocities = velocities
    mda_u.atoms.write(output_file)


def write_checkpoint(simulation: Simulation, input_ochk: str) -> None:
    """Save simulation checkpoint state to XML file.

    Args:
        simulation: Simulation object to save.
        input_ochk: Path to output checkpoint file.
    """
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(input_ochk, "w") as f:
        f.write(mm.XmlSerializer.serialize(state))
    print(f"\nWrite checkpoint file: {input_ochk}")


def cleanup(signum: int, simulation: Simulation, inputs: Inputs) -> None:
    """Handle signal interrupts by saving checkpoint before exiting.

    Args:
        signum: Signal number received.
        simulation: Running simulation object.
        inputs: Input parameters object with ochk attribute.
    """
    print("Received signal", signum, ". Performing cleanup...")
    write_checkpoint(simulation, inputs.ochk)
    sys.exit(0)


signal.signal(signal.SIGTERM, cleanup)


def load_platform(inputs: Inputs) -> tuple[mm.Platform, dict[str, str]]:
    """Configure OpenMM platform based on input parameters.

    Args:
        inputs: Input parameters object containing platform settings.
            Expected attributes: platform, precision, GPU_id.

    Returns:
        A tuple containing:
            - The configured OpenMM Platform.
            - Dictionary of platform properties.

    Raises:
        Exception: If the specified platform is not supported.
    """
    if inputs.platform == "CPU":
        platform = mm.Platform.getPlatformByName("CPU")
        print("\nUsing platform: CPU, Precision: default")
        platform_properties: dict[str, str] = {}
    elif inputs.platform == "Reference":
        platform = mm.Platform.getPlatformByName("Reference")
        print("\nUsing platform: Reference, Precision: double")
        platform_properties = {}
    elif inputs.platform == "CUDA":
        platform = mm.Platform.getPlatformByName("CUDA")
        platform_properties = {"CudaPrecision": f"{inputs.precision}"}
        print(f"\nUsing platform: CUDA, Precision: {inputs.precision}")
        if inputs.GPU_id:
            platform_properties["UseBlockingSync"] = "false"
            platform_properties["DeviceIndex"] = inputs.GPU_id
            print(f"Using GPU_id: {inputs.GPU_id}")
    elif inputs.platform == "OpenCL":
        platform = mm.Platform.getPlatformByName("OpenCL")
        platform_properties = {"Precision": f"{inputs.precision}"}
        print(f"\nUsing platform: OpenCL, Precision: {inputs.precision}")
        if inputs.GPU_id:
            platform_properties["DeviceIndex"] = inputs.GPU_id
            print(f"Using GPU_id: {inputs.GPU_id}")
    else:
        raise Exception(f"Error: Unsupported platform {inputs.platform}")
    return platform, platform_properties


def generate_topology(
    inpfile: str,
) -> tuple[System, Topology, mm.Vec3]:
    """Generate OpenMM system and topology from input file.

    This is the main function to set up a CTGoMartini molecular dynamics
    simulation, loading parameters, configuring the platform, creating the
    system with restraints and barostats, and preparing the simulation context.

    Args:
        inpfile: Path to input parameter file containing simulation settings.

    Returns:
        A tuple containing:
            - The configured OpenMM System.
            - The molecular Topology.
            - The periodic box vectors.

    Workflow:
        1. Load simulation parameters from input file.
        2. Configure computational platform (CPU/GPU).
        3. Load molecular structure and topology.
        4. Create OpenMM system with forces and constraints.
        5. Add restraints if requested.
        6. Add barostat if requested.
        7. Set up integrator and simulation context.
        8. Load checkpoint or generate initial velocities.
    """
    start_time = datetime.datetime.now()

    # Load parameters
    print("Loading parameters")
    inputs = read_inputs(inpfile)

    # Platform
    platform, platform_properties = load_platform(inputs)

    # Load cord and box vectors
    conf, box_vectors = load_structure(inputs.input)

    # Load topology
    defines = inputs.defines
    top = MartiniTopFile(
        inputs.topol,
        periodicBoxVectors=box_vectors,
        defines=defines,
    )

    # Create system
    system = top.create_system(
        nonbonded_cutoff=inputs.nonbonded_cutoff * u.nanometer,
        epsilon_r=inputs.epsilon_r,
    )

    # Add restraints
    if inputs.gen_rest == "yes":
        gen_restraints(
            inputs.input, inputs.atomname, inputs.fc, inputs.gen_rest_file
        )
    if inputs.rest == "yes":
        system = add_restraints(system, inputs)

    # Add plumed
    if inputs.plumed == "yes":
        raise ValueError("Error: Plumed is not supported!")

    # Add a barostat
    if inputs.pcouple == "yes":
        if inputs.p_type == "isotropic":
            barostat = mm.MonteCarloBarostat(
                inputs.p_ref * u.bar,
                inputs.temp * u.kelvin,
                inputs.p_freq,
            )
        elif inputs.p_type == "membrane":
            barostat = mm.MonteCarloMembraneBarostat(
                inputs.p_ref * u.bar,
                inputs.p_tens * u.bar * u.nanometers,
                inputs.temp * u.kelvin,
                inputs.p_XYMode,
                inputs.p_ZMode,
                inputs.p_freq,
            )
        else:
            raise Exception("Unsupported barostat type: ", inputs.p_type)

        system.addForce(barostat)

    # Integrator
    integrator = mm.LangevinIntegrator(
        inputs.temp * u.kelvin,
        inputs.fric_coeff / u.picosecond,
        inputs.dt * u.picosecond,
    )
    if inputs.const_tol:
        integrator.setConstraintTolerance(inputs.const_tol)

    simulation = Simulation(
        top.topology, system, integrator, platform, platform_properties
    )

    # Set positions
    assert len(conf.getPositions()) == top.topology.getNumAtoms(), (
        f"Error: Number of atoms in {inputs.input} is not the same as "
        f"that from {inputs.topol}!"
    )
    simulation.context.setPositions(conf.getPositions())

    # Check the charges of the system
    if top.charges != 0:
        print(f"Warning: The charges of the system are {top.charges} instead of 0.")

    if inputs.ichk:
        try:
            with open(inputs.ichk, "r") as f:
                simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
        except Exception:
            with open(inputs.ichk, "rb") as f:
                simulation.context.loadCheckpoint(f.read())
        print(f"\nLoad checkpoint file: {inputs.ichk}")

    # System Loading Finishes!
    print("\nLoading system finishes!")
    report_time(start_time)

    start_time = datetime.datetime.now()

    # Generate initial velocities
    if inputs.gen_vel == "yes":
        print(f"\nGenerate initial velocities: {inputs.gen_temp} K")
        if inputs.gen_seed:
            simulation.context.setVelocitiesToTemperature(
                inputs.gen_temp * u.kelvin, inputs.gen_seed
            )
        else:
            simulation.context.setVelocitiesToTemperature(
                inputs.gen_temp * u.kelvin
            )

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions()
    velocities = state.getVelocities()
    box_vectors = state.getPeriodicBoxVectors()

    return system, top.topology, box_vectors
