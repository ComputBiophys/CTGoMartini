#!/usr/bin/env python3
"""
Main simulation runner for CTGoMartini.

Provides functionality for running molecular dynamics simulations using
the OpenMM engine with Martini coarse-grained force fields.

Authors: Song Yang
"""

from __future__ import annotations

import argparse
import datetime
import os
import signal
import sys
from typing import Any

import MDAnalysis as mda
import openmm as mm
import openmm.unit as u
from openmm.app import (
    DCDReporter, 
    GromacsGroFile, 
    PDBFile, 
    Simulation, 
    StateDataReporter, 
    XTCReporter
)
from openmm.app import CheckpointReporter

from ctgomartini.api import MartiniTopFile
from ctgomartini.utils import read_inputs

import warnings
warnings.filterwarnings("ignore")


def report_time(start_time: datetime.datetime) -> None:
    """
    Report elapsed simulation time.
    
    Args:
        start_time: Start time of the simulation segment.
    """
    end_time = datetime.datetime.now()
    elapsed = (end_time - start_time).total_seconds()
    
    start_str = start_time.__format__("%Y-%m-%d %H:%M:%S")
    end_str = end_time.__format__("%Y-%m-%d %H:%M:%S")
    
    print(f"  Start Time: {start_str}")
    print(f"    End Time: {end_str}")
    print(f"Elapsed Time: {elapsed:.2f}")


def load_structure(str_file: str) -> tuple[Any, Any]:
    """
    Load molecular structure from GRO or PDB file.
    
    Args:
        str_file: Path to structure file (.gro or .pdb).
        
    Returns:
        Tuple of (configuration object, periodic box vectors).
        
    Raises:
        ValueError: If file format is not supported.
    """
    ext = str_file.split('.')[-1].lower()
    
    if ext == 'gro':
        conf = GromacsGroFile(str_file)
        box_vectors = conf.getPeriodicBoxVectors()
    elif ext == 'pdb':
        conf = PDBFile(str_file)
        box_vectors = conf.getTopology().getPeriodicBoxVectors()
    else:
        raise ValueError(f'Unsupported structure file format: {ext}')
        
    return conf, box_vectors


def generate_restraints(
    str_file: str, 
    atomname: str, 
    fc: float = 1000.0, 
    rest_file: str = "restraints.txt"
) -> None:
    """
    Generate position restraint file for selected atoms.
    
    Args:
        str_file: Structure file (gro/pdb) containing atoms to restrain.
        atomname: Atom name to apply restraints to (e.g., 'BB' for backbone).
        fc: Force constant in kJ/(mol·nm²). Default 1000.
        rest_file: Output restraint file path. Default "restraints.txt".
    """
    universe = mda.Universe(str_file)
    selection = universe.select_atoms(f"name {atomname}")
    
    lines = ["; atomindex functype(1) fc_x fc_y fc_z\n"]
    lines.append("; atomid start from 1\n")
    
    for idx in selection.indices:
        lines.append(f"{idx+1:>5} 1 {fc} {fc} {fc}\n")
    
    with open(rest_file, 'w') as f:
        f.writelines(lines)


def add_restraints(
    system: mm.System, 
    inputs: Any
) -> mm.System:
    """
    Add positional restraints to the system.
    
    Args:
        system: The OpenMM system to add restraints to.
        inputs: Input parameters object with restraint settings.
        
    Returns:
        The system with restraints added.
    """
    crd, _ = load_structure(inputs.rest_ref)
    
    if inputs.rest == 'yes':
        # Anisotropic positional restraints
        posres = mm.CustomExternalForce(
            '1/2*kx*periodicdistance(x, 0, 0, x0, 0, 0)^2 + '
            '1/2*ky*periodicdistance(0, y, 0, 0, y0, 0)^2 + '
            '1/2*kz*periodicdistance(0, 0, z, 0, 0, z0)^2;'
        )
        posres.addPerParticleParameter('kx')
        posres.addPerParticleParameter('ky')
        posres.addPerParticleParameter('kz')
        posres.addPerParticleParameter('x0')
        posres.addPerParticleParameter('y0')
        posres.addPerParticleParameter('z0')
        
        # Parse restraint file
        with open(inputs.rest_file, 'r') as f:
            for line in f:
                if ';' in line:
                    line = line.split(';')[0]
                sline = line.strip()
                if not sline:
                    continue
                    
                segments, functype, fcx, fcy, fcz = sline.split()[:5]
                atom_idx = int(segments) - 1
                fcx_f = float(fcx)
                fcy_f = float(fcy)
                fcz_f = float(fcz)
                
                assert functype == '1', f'Unsupported restraint type: {line}'
                
                # Get reference positions
                xpos = crd.positions[atom_idx].value_in_unit(u.nanometers)[0]
                ypos = crd.positions[atom_idx].value_in_unit(u.nanometers)[1]
                zpos = crd.positions[atom_idx].value_in_unit(u.nanometers)[2]
                
                # Add restraint if force constants are positive
                if fcx_f >= 0 and fcy_f >= 0 and fcz_f >= 0:
                    posres.addParticle(atom_idx, [fcx_f, fcy_f, fcz_f, xpos, ypos, zpos])

        system.addForce(posres)
        
    return system


def backup_file(file: str) -> None:
    """
    Create backup of existing file with incremental numbering.
    
    Args:
        file: Path to file that needs backup.
    """
    if os.path.isfile(file):
        i = 1
        newfile = f"{file}.bk{i}"
        while os.path.isfile(newfile):
            i += 1
            newfile = f"{file}.bk{i}"
        os.rename(file, newfile)


def write_output(output_file: str, simulation: Simulation, strfile: str) -> None:
    """
    Write simulation output to file.
    
    Args:
        output_file: Path to output file.
        simulation: Simulation object containing current state.
        strfile: Reference structure file for topology.
    """
    state = simulation.context.getState(
        getPositions=True, 
        getVelocities=True,
        enforcePeriodicBox=True
    )
    
    crd = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    velocities = state.getVelocities(asNumpy=True).value_in_unit(u.angstrom / u.picosecond)
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)
    box_diag = box_vectors[[0, 1, 2], [0, 1, 2]]
    
    # Write via MDAnalysis
    mda_u = mda.Universe(strfile)
    mda_u.atoms.positions = crd
    mda_u.trajectory[0].velocities = True
    mda_u.dimensions[:3] = box_diag
    mda_u.atoms.velocities = velocities
    mda_u.atoms.write(output_file)


def write_checkpoint(simulation: Simulation, output_file: str) -> None:
    """
    Save simulation checkpoint state to XML file.
    
    Args:
        simulation: Simulation object to save.
        output_file: Path to output checkpoint file.
    """
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(output_file, 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    print(f"\nWrote checkpoint file: {output_file}")


def cleanup(signum: int, simulation: Simulation, inputs: Any) -> None:
    """
    Handle signal interrupts by saving checkpoint before exiting.
    
    Args:
        signum: Signal number.
        simulation: Running simulation.
        inputs: Input parameters object.
    """
    print(f"Received signal {signum}. Performing cleanup...")
    write_checkpoint(simulation, inputs.ochk)
    sys.exit(0)


# Register signal handler
signal.signal(signal.SIGTERM, cleanup)


def mdrun(inpfile: str) -> None:
    """
    Execute CTGoMartini molecular dynamics simulation.
    
    Args:
        inpfile: Path to input parameter file.
    """
    start_time = datetime.datetime.now()

    # Load parameters
    print("Loading parameters")
    inputs = read_inputs(inpfile)

    # Configure platform
    if inputs.platform == 'CPU':
        platform = mm.Platform.getPlatformByName("CPU")
        print("\nUsing platform: CPU, Precision: default")
        platform_properties: dict[str, str] = {}
        
    elif inputs.platform == 'Reference':
        platform = mm.Platform.getPlatformByName("Reference")
        print("\nUsing platform: Reference, Precision: double")
        platform_properties = {}
        
    elif inputs.platform == 'CUDA':
        platform = mm.Platform.getPlatformByName("CUDA")
        platform_properties = {'CudaPrecision': inputs.precision}
        print(f"\nUsing platform: CUDA, Precision: {inputs.precision}")
        if inputs.GPU_id:
            platform_properties['UseBlockingSync'] = 'false'
            platform_properties["DeviceIndex"] = inputs.GPU_id
            print(f"Using GPU_id: {inputs.GPU_id}")
            
    elif inputs.platform == 'OpenCL':
        platform = mm.Platform.getPlatformByName("OpenCL")
        platform_properties = {'Precision': inputs.precision}
        print(f"\nUsing platform: OpenCL, Precision: {inputs.precision}")
        if inputs.GPU_id:
            platform_properties["DeviceIndex"] = inputs.GPU_id
            print(f"Using GPU_id: {inputs.GPU_id}")
    else:
        raise ValueError(f"Unsupported platform: {inputs.platform}")

    # Load structure and topology
    conf, box_vectors = load_structure(inputs.input)
    
    top = MartiniTopFile(
        inputs.topol,
        periodicBoxVectors=box_vectors,
        defines=inputs.defines,
    )

    # Create system
    system = top.create_system(
        nonbonded_cutoff=inputs.nonbonded_cutoff * u.nanometer,
        epsilon_r=inputs.epsilon_r
    )

    # Add restraints
    if inputs.gen_rest == 'yes':
        generate_restraints(
            inputs.input, inputs.atomname, inputs.fc, inputs.gen_rest_file
        )
    if inputs.rest == 'yes':
        system = add_restraints(system, inputs)

    # Add PLUMED if requested
    if inputs.plumed == 'yes':
        from openmmplumed import PlumedForce
        print(f"\nAdding PLUMED: {inputs.plumed_file}")
        
        with open(inputs.plumed_file, 'r') as f:
            script = f.read()
        system.addForce(PlumedForce(script))

    # Add barostat
    if inputs.pcouple == 'yes':
        if inputs.p_type == 'isotropic':
            barostat = mm.MonteCarloBarostat(
                inputs.p_ref * u.bar,
                inputs.temp * u.kelvin,
                inputs.p_freq
            )
        elif inputs.p_type == 'membrane':
            barostat = mm.MonteCarloMembraneBarostat(
                inputs.p_ref * u.bar,
                inputs.p_tens * u.bar * u.nanometers,
                inputs.temp * u.kelvin,
                inputs.p_XYMode,
                inputs.p_ZMode,
                inputs.p_freq
            )
        else:
            raise ValueError(f'Unsupported barostat type: {inputs.p_type}')
        
        system.addForce(barostat)

    # Create integrator
    integrator = mm.LangevinIntegrator(
        inputs.temp * u.kelvin,
        inputs.fric_coeff / u.picosecond,
        inputs.dt * u.picosecond
    )
    
    if inputs.const_tol:
        integrator.setConstraintTolerance(inputs.const_tol)
        print("Set Constraint Tolerance:", integrator.getConstraintTolerance())

    # Create simulation
    simulation = Simulation(
        top.topology, system, integrator, platform, platform_properties
    )

    # Set positions
    assert len(conf.getPositions()) == top.topology.getNumAtoms(), (
        f"Error: Number of atoms in {inputs.input} does not match {inputs.topol}"
    )
    simulation.context.setPositions(conf.getPositions())

    # Check system charge
    if top.charges != 0:
        print(f'Warning: System charge is {top.charges} instead of 0.')

    # Load checkpoint if provided
    if inputs.ichk:
        try:
            with open(inputs.ichk, 'r') as f:
                simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
        except Exception:
            with open(inputs.ichk, 'rb') as f:
                simulation.context.loadCheckpoint(f.read())
        print(f"\nLoaded checkpoint file: {inputs.ichk}")

    print("\nLoading system finishes!")
    report_time(start_time)

    # Calculate initial energy
    start_time = datetime.datetime.now()
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    print(f"\nInitial system energy: {energy:.3f} kJ/mol")

    # Energy minimization
    if inputs.mini_nstep > 0:
        print(f"\nEnergy minimization: {inputs.mini_nstep} steps")
        simulation.minimizeEnergy(
            tolerance=inputs.mini_Tol,
            maxIterations=inputs.mini_nstep
        )
        energy = simulation.context.getState(
            getEnergy=True
        ).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        print(f"Minimized system energy: {energy:.3f} kJ/mol")
        report_time(start_time)

    # Generate initial velocities
    if inputs.gen_vel == 'yes':
        print(f"\nGenerating initial velocities: {inputs.gen_temp} K")
        if inputs.gen_seed:
            simulation.context.setVelocitiesToTemperature(
                inputs.gen_temp * u.kelvin, inputs.gen_seed
            )
        else:
            simulation.context.setVelocitiesToTemperature(
                inputs.gen_temp * u.kelvin
            )

    # Production MD
    if inputs.nstep > 0:
        start_time = datetime.datetime.now()

        # Select trajectory format
        if inputs.odcd and not inputs.oxtc:
            traj_reporter = DCDReporter
            traj_file = inputs.odcd
        elif inputs.oxtc and not inputs.odcd:
            traj_reporter = XTCReporter
            traj_file = inputs.oxtc
        else:
            raise ValueError("Error: Please specify either odcd or oxtc, not both!")

        # Set starting point
        if inputs.append == 'no':
            begin_step = inputs.b_step
            simulation.context.setStepCount(begin_step)
            simulation.context.setTime(inputs.dt * begin_step)
        else:
            begin_step = simulation.context.getStepCount()

        end_step = inputs.nstep
        remaining_steps = end_step - begin_step

        print(f"\nMD run: begin {begin_step}, end {end_step}, total {remaining_steps}")

        # Setup reporters
        if inputs.nstdcd > 0:
            if inputs.append == 'yes':
                simulation.reporters.append(
                    traj_reporter(traj_file, inputs.nstdcd, append=True)
                )
            else:
                backup_file(traj_file)
                simulation.reporters.append(traj_reporter(traj_file, inputs.nstdcd))

            backup_file(inputs.ochk)
            simulation.reporters.append(
                CheckpointReporter(inputs.ochk, inputs.nstdcd, writeState=True)
            )

        simulation.reporters.append(
            StateDataReporter(
                sys.stdout,
                inputs.nstout,
                step=True,
                time=True,
                potentialEnergy=True,
                temperature=True,
                progress=True,
                remainingTime=True,
                speed=True,
                totalSteps=end_step,
                separator='\t'
            )
        )

        # Run simulation
        try:
            simulation.step(remaining_steps)
        except KeyboardInterrupt:
            print("\nSimulation interrupted!")
            report_time(start_time)
            cleanup(signal.SIGINT, simulation, inputs)

    # Write final output
    if inputs.output:
        write_output(inputs.output, simulation, inputs.input)
    if inputs.output_pdb:
        write_output(inputs.output_pdb, simulation, inputs.input)

    if inputs.ochk:
        write_checkpoint(simulation, inputs.ochk)

    print("\nSimulation finished!")
    report_time(start_time)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run CTGoMartini molecular dynamics simulation using OpenMM"
    )
    parser.add_argument(
        '-i',
        dest='inpfile',
        help='Input parameter file containing simulation settings',
        required=True
    )
    args = parser.parse_args()
    mdrun(args.inpfile)
