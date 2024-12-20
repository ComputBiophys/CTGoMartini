#!/usr/bin/env python3

"""
Authors: Song Yang
Last update: 11 27, 2023
"""

import os, sys
from simtk import unit as u
from simtk import openmm as mm
from openmm.app import *
from ctgomartini.api import MartiniTopFile
from ctgomartini.func import read_inputs
import MDAnalysis as mda
import argparse
import datetime
import signal

def ReportTime(start_time):
    end_time=datetime.datetime.now()
    elapsed=(end_time-start_time).total_seconds()
    
    start_time=start_time.__format__("%Y-%m-%d %H:%M:%S")
    end_time=end_time.__format__("%Y-%m-%d %H:%M:%S")
    # print(end='\n')
    print(f"  Start Time: {start_time}")
    print(f"    End Time: {end_time}")
    print(f"Elapsed Time: {elapsed:.2f}")

def LoadStructure(str_file):
    if str_file.split('.')[-1] == 'gro':
        conf = GromacsGroFile(str_file)
        box_vectors = conf.getPeriodicBoxVectors()
    elif str_file.split('.')[-1] == 'pdb':
        conf = PDBFile(str_file)
        box_vectors = conf.getTopology().getPeriodicBoxVectors()    
    else:
        raise Exception('Unsupported structure file: ', str_file)
    return conf, box_vectors

def gen_restraints(str_file, atomname, fc=1000, rest_file="restraints.txt"):
    u=mda.Universe(str_file)
    sel=u.select_atoms(f"name {atomname}")
    
    newlines=["; atomindex functype(1) fc_x fc_y fc_z\n"]
    newlines.append(f'; atomid start from 1\n')
    for i in sel.indices:
        newlines.append(f'{i+1:>5} 1 {fc} {fc} {fc}\n')
    
    with open(rest_file,'w') as g:
        g.writelines(newlines)

def restraints(system, inputs):
    crd, _ = LoadStructure(inputs.rest_ref)
    if inputs.rest == 'yes':
        # positional restraints for protein
        # posresPROT = mm.CustomExternalForce('1/2*k*periodicdistance(x, y, z, x0, y0, z0)^2;')
        posresPROT = mm.CustomExternalForce('1/2*kx*periodicdistance(x, 0, 0, x0, 0, 0)^2 + 1/2*ky*periodicdistance(0, y, 0, 0, y0, 0)^2 + 1/2*kz*periodicdistance(0, 0, z, 0, 0, z0)^2;')
        posresPROT.addPerParticleParameter('kx')
        posresPROT.addPerParticleParameter('ky')
        posresPROT.addPerParticleParameter('kz')
        posresPROT.addPerParticleParameter('x0')
        posresPROT.addPerParticleParameter('y0')
        posresPROT.addPerParticleParameter('z0')
        for line in open(inputs.rest_file, 'r'):
            if line.find(';') >= 0: line = line.split(';')[0]
            sline = line.strip()
            if sline == '': continue
            segments, functype, fcx, fcy, fcz = sline.split()[:5]
            atom1 = int(segments) - 1
            fcx, fcy, fcz = float(fcx), float(fcy), float(fcz)
            assert functype == '1', f'Error: Unsupport position restraint type.\n {line}'
            xpos  = crd.positions[atom1].value_in_unit(u.nanometers)[0]
            ypos  = crd.positions[atom1].value_in_unit(u.nanometers)[1]
            zpos  = crd.positions[atom1].value_in_unit(u.nanometers)[2]
            if fcx >= 0 and fcy >=0 and fcz >= 0:
                posresPROT.addParticle(atom1, [fcx, fcy, fcz, xpos, ypos, zpos])

        system.addForce(posresPROT)
    return system


def BackupFile(file):
    if os.path.isfile(file):
        i = 1
        newfile = file + f'.bk{i}'
        while os.path.isfile(newfile):
            i += 1
            newfile = file + f'.bk{i}'
        os.rename(file,newfile)

def WriteOutput(output_file, simulation, strfile):
    # Get crd, volocities, box_vectors
    state = simulation.context.getState(getPositions=True,getVelocities=True)
    crd = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    velocities = state.getVelocities(asNumpy=True).value_in_unit(u.angstrom/u.picosecond)
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)[[0,1,2],[0,1,2]]
    
    # Write out output_file
    mda_u = mda.Universe(strfile)
    mda_u.atoms.positions = crd
    mda_u.trajectory[0].velocities = True
    mda_u.dimensions[:3] = box_vectors # only for rectangular/cubic box
    mda_u.atoms.velocities = velocities
    mda_u.atoms.write(output_file)

def WriteCheckPoint(simulation, input_ochk):
    # Write CheckPoint file
    state = simulation.context.getState(getPositions=True, getVelocities=True )
    with open(input_ochk, 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    print(f"\nWrite checkpoint file: {input_ochk}")

def Cleanup(signum, simulation, inputs):
    print("Received signal", signum, ". Performing cleanup...")
    WriteCheckPoint(simulation, inputs.ochk)
    sys.exit(0)
signal.signal(signal.SIGTERM, Cleanup)

def mdrun(inpfile):
    """
    inpfile: str,
        Input parameter file.
    """
    start_time=datetime.datetime.now()

    # Load parameters
    print("Loading parameters")
    inputs = read_inputs(inpfile)

    # Platform
    if inputs.platform == 'CPU':
        platform = mm.Platform.getPlatformByName("CPU")
        print("\nUsing platform: CPU, Precision: default")
        platformProperties = {}
    elif inputs.platform == 'Reference':
        platform = mm.Platform.getPlatformByName("Reference")
        print("\nUsing platform: Reference, Precision: double")        
        platformProperties = {}
    elif inputs.platform == 'CUDA':
        platform = mm.Platform.getPlatformByName("CUDA")
        platformProperties = {'CudaPrecision': '{}'.format(inputs.precision)}
        print(f"\nUsing platform: CUDA, Precision: {inputs.precision}")
        if inputs.GPU_id:
            platformProperties['UseBlockingSync']='false'
            platformProperties["DeviceIndex"] = inputs.GPU_id
            print(f"Using GPU_id: {inputs.GPU_id}")
    elif inputs.platform == 'OpenCL':
        platform = mm.Platform.getPlatformByName("OpenCL")
        platformProperties = {'Precision': '{}'.format(inputs.precision)} 
        print(f"\nUsing platform: OpenCL, Precision: {inputs.precision}")
        if inputs.GPU_id:
            platformProperties["DeviceIndex"] = inputs.GPU_id
            print(f"Using GPU_id: {inputs.GPU_id}")        
    else:
        raise Exception(f"Error: Unsupported platform {inputs.platform}")


    # Load cord and box vectors
    conf, box_vectors = LoadStructure(inputs.input)

    #Load topol
    defines = inputs.defines
    top = MartiniTopFile(
        inputs.topol,
        periodicBoxVectors=box_vectors,
        defines=defines,
    )

    # Create system
    system = top.create_system(nonbonded_cutoff=inputs.nonbonded_cutoff * u.nanometer, epsilon_r=inputs.epsilon_r)

    # Add restraints
    if inputs.gen_rest == 'yes': gen_restraints(inputs.input, inputs.atomname, inputs.fc, inputs.gen_rest_file)
    if inputs.rest == 'yes':     system = restraints(system, inputs)

    # Add a barostat
    if inputs.pcouple=='yes':
        if inputs.p_type == 'isotropic': 
            barostat = mm.MonteCarloBarostat(inputs.p_ref * u.bar, 
                                            inputs.temp * u.kelvin, inputs.p_freq)
        elif inputs.p_type == 'membrane':
            barostat = mm.MonteCarloMembraneBarostat( inputs.p_ref * u.bar,
                                                    inputs.p_tens * u.bar*u.nanometers,
                                                    inputs.temp * u.kelvin, 
                                                    inputs.p_XYMode, 
                                                    inputs.p_ZMode, 
                                                    inputs.p_freq )
        else:
            raise Exception('Unsupported barostat type: ', inputs.p_type)
        
        system.addForce(barostat)

    # Intergrator
    # integrator = mm.VerletIntegrator(inputs.dt * u.picosecond)
    #integrator = mm.LangevinMiddleIntegrator(
    integrator = mm.LangevinIntegrator(
        inputs.temp * u.kelvin, inputs.fric_coeff / u.picosecond, inputs.dt * u.picosecond
    )
    if inputs.const_tol: integrator.setConstraintTolerance(inputs.const_tol) # Set the constraint tolerance
    print(integrator.getConstraintTolerance())


    simulation = Simulation(top.topology, system, integrator, platform, platformProperties)

    # Set positions
    assert len(conf.getPositions()) == top.topology.getNumAtoms(), f"Error: Number of atoms in {inputs.input} is not the same as that from {inputs.topol}!"
    simulation.context.setPositions(conf.getPositions())

    # Check the charges of the system
    if top.charges != 0:
        print(f'Warning: The charges of the system are {top.charges} instead of 0.')

    if inputs.ichk:
        try:
            with open(inputs.ichk, 'r') as f:
                simulation.context.setState(mm.XmlSerializer.deserialize(f.read()))
        except:
            with open(inputs.ichk, 'rb') as f:
                simulation.context.loadCheckpoint(f.read())
        print(f"\nLoad checkpoint file: {inputs.ichk}")

    # System Loading Finishes!
    print("\nLoading system finishes!")
    ReportTime(start_time)
    
    start_time=datetime.datetime.now()
    # Calculate initial system energy
    state = simulation.context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    print(f"\nInitial system energy: {energy:.3f} kJ/mol")

    # Energy minimization
    if inputs.mini_nstep > 0:
        print(f"\nPlan Energy minimization: {inputs.mini_nstep} steps")
        simulation.minimizeEnergy(tolerance=inputs.mini_Tol,
                                maxIterations=inputs.mini_nstep)
        energy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        print(f"Minimized system energy: {energy:.3f} kJ/mol")
        ReportTime(start_time)
        
    # Generate initial velocities
    if inputs.gen_vel == 'yes':
        print(f"\nGenerate initial velocities: {inputs.gen_temp} K")
        if inputs.gen_seed:
            simulation.context.setVelocitiesToTemperature(inputs.gen_temp * u.kelvin, inputs.gen_seed)
        else:
            simulation.context.setVelocitiesToTemperature(inputs.gen_temp * u.kelvin)
    
    # Production
    if inputs.nstep > 0:
        start_time=datetime.datetime.now()
        # Set the output format
        if inputs.odcd and not inputs.oxtc:
            TrajReporter = DCDReporter
            otraj = inputs.odcd
        elif inputs.oxtc and not inputs.odcd:
            TrajReporter = XTCReporter
            otraj = inputs.oxtc
        else:
            raise ValueError("Error: Please specify either odcd or oxtc!")

        if inputs.append == 'no': 
            b_step = inputs.b_step
            simulation.context.setStepCount(b_step)
            simulation.context.setTime(inputs.dt * b_step)
        else:
            b_step = simulation.context.getStepCount()
        e_step = inputs.nstep
        be_step=e_step-b_step
        print(f"\nMD run: begin {b_step}, end {e_step}, total {be_step}")

        if inputs.nstdcd > 0:
            if inputs.append == 'yes':
                simulation.reporters.append(TrajReporter(otraj, inputs.nstdcd, append=True,))
            else:
                BackupFile(otraj)
                simulation.reporters.append(TrajReporter(otraj, inputs.nstdcd))
            BackupFile(inputs.ochk)
            simulation.reporters.append(CheckpointReporter(inputs.ochk, inputs.nstdcd, writeState=True))
        simulation.reporters.append(
            StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, progress=True,
                            remainingTime=True, speed=True, 
                            totalSteps=e_step, separator='\t')
        )
        
        try:
            simulation.step(be_step)
            be_step = 0 
        except KeyboardInterrupt:
            print("The task has been canceled!")
            ReportTime(start_time)
            # WriteCheckPoint(simulation, inputs.ochk)
            Cleanup(signal.SIGINT, simulation, inputs)
        

    # Write output file
    if inputs.output:
        WriteOutput(inputs.output, simulation, inputs.input)
    if inputs.output_pdb:
        WriteOutput(inputs.output_pdb, simulation, inputs.input)
    if inputs.ochk:
        WriteCheckPoint(simulation, inputs.ochk)
    
    print("The task has finished!")
    ReportTime(start_time)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='inpfile', help='Input parameter file', required=True)
    args = parser.parse_args()
    mdrun(args.inpfile)
    
