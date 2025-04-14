#!/usr/bin/env python3

"""
Authors: Song Yang
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

import warnings
warnings.filterwarnings("ignore")

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
    """Generate restraint file for selected atoms.
    
    Args:
        str_file (str): Structure file (gro/pdb) containing atoms to restrain
        atomname (str): Atom name to apply restraints to (e.g. 'CA')
        fc (int, optional): Force constant in kJ/(mol·nm²). Defaults to 1000.
        rest_file (str, optional): Output restraint file. Defaults to "restraints.txt".
    """
    u=mda.Universe(str_file)
    sel=u.select_atoms(f"name {atomname}")
    
    newlines=["; atomindex functype(1) fc_x fc_y fc_z\n"]
    newlines.append(f'; atomid start from 1\n')
    for i in sel.indices:
        newlines.append(f'{i+1:>5} 1 {fc} {fc} {fc}\n')
    
    with open(rest_file,'w') as g:
        g.writelines(newlines)

def restraints(system, inputs):
    """Add positional restraints to the system.
    
    Args:
        system (openmm.System): The system to add restraints to
        inputs (object): Input parameters object containing restraint settings
        
    Returns:
        openmm.System: The system with restraints added
    """
    crd, _ = LoadStructure(inputs.rest_ref)
    if inputs.rest == 'yes':
        # positional restraints for protein
        # posresPROT = mm.CustomExternalForce('1/2*k*periodicdistance(x, y, z, x0, y0, z0)^2;')
        # Create custom external force for anisotropic positional restraints
        posresPROT = mm.CustomExternalForce('1/2*kx*periodicdistance(x, 0, 0, x0, 0, 0)^2 + 1/2*ky*periodicdistance(0, y, 0, 0, y0, 0)^2 + 1/2*kz*periodicdistance(0, 0, z, 0, 0, z0)^2;')
        posresPROT.addPerParticleParameter('kx')
        posresPROT.addPerParticleParameter('ky')
        posresPROT.addPerParticleParameter('kz')
        posresPROT.addPerParticleParameter('x0')
        posresPROT.addPerParticleParameter('y0')
        posresPROT.addPerParticleParameter('z0')
        
        # Parse restraint file and add particles to the force
        for line in open(inputs.rest_file, 'r'):
            if line.find(';') >= 0: line = line.split(';')[0]
            sline = line.strip()
            if sline == '': continue
            segments, functype, fcx, fcy, fcz = sline.split()[:5]
            atom1 = int(segments) - 1  # Convert to 0-based index
            fcx, fcy, fcz = float(fcx), float(fcy), float(fcz)
            assert functype == '1', f'Error: Unsupport position restraint type.\n {line}'
            
            # Get reference positions in nanometers
            xpos = crd.positions[atom1].value_in_unit(u.nanometers)[0]
            ypos = crd.positions[atom1].value_in_unit(u.nanometers)[1]
            zpos = crd.positions[atom1].value_in_unit(u.nanometers)[2]
            
            # Add restraint if any force constant is positive
            if fcx >= 0 and fcy >=0 and fcz >= 0:
                posresPROT.addParticle(atom1, [fcx, fcy, fcz, xpos, ypos, zpos])

        system.addForce(posresPROT)
    return system


def BackupFile(file):
    """Create backup of existing file with incremental numbering.
    
    Args:
        file (str): Path to file that needs backup
    """
    if os.path.isfile(file):
        i = 1
        newfile = file + f'.bk{i}'
        while os.path.isfile(newfile):
            i += 1
            newfile = file + f'.bk{i}'
        os.rename(file,newfile)

def WriteOutput(output_file, simulation, strfile):
    """Write simulation output to file in specified format.
    
    Args:
        output_file (str): Path to output file
        simulation (openmm.Simulation): Simulation object containing state
        strfile (str): Reference structure file for topology
    """
    # Get current state information
    state = simulation.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
    crd = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    velocities = state.getVelocities(asNumpy=True).value_in_unit(u.angstrom/u.picosecond)
    box_vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(u.angstrom)[[0,1,2],[0,1,2]]
    
    # Create MDAnalysis universe and update with current simulation state
    mda_u = mda.Universe(strfile)
    mda_u.atoms.positions = crd
    mda_u.trajectory[0].velocities = True
    mda_u.dimensions[:3] = box_vectors # only for rectangular/cubic box
    mda_u.atoms.velocities = velocities
    mda_u.atoms.write(output_file)

def WriteCheckPoint(simulation, input_ochk):
    """Save simulation checkpoint state to XML file.
    
    Args:
        simulation (openmm.Simulation): Simulation object to save
        input_ochk (str): Path to output checkpoint file
    """
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open(input_ochk, 'w') as f:
        f.write(mm.XmlSerializer.serialize(state))
    print(f"\nWrite checkpoint file: {input_ochk}")

def Cleanup(signum, simulation, inputs):
    """Handle signal interrupts by saving checkpoint before exiting.
    
    Args:
        signum (int): Signal number
        simulation (openmm.Simulation): Running simulation
        inputs (object): Input parameters
    """
    print("Received signal", signum, ". Performing cleanup...")
    WriteCheckPoint(simulation, inputs.ochk)
    sys.exit(0)
signal.signal(signal.SIGTERM, Cleanup)

def mdrun(inpfile):
    """Main function to execute CTGoMartini molecular dynamics simulation.
    
    Args:
        inpfile (str): Path to input parameter file containing simulation settings
        
    Workflow:
        1. Load simulation parameters from input file
        2. Configure computational platform (CPU/GPU)
        3. Load molecular structure and topology
        4. Create OpenMM system with forces and constraints
        5. Perform energy minimization (if requested)
        6. Generate initial velocities (if requested)
        7. Run production MD simulation
        8. Write output files
        
    Returns:
        None: Outputs are written to files specified in input parameters
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

    # Add plumed
    if inputs.plumed == 'yes':
        from openmmplumed import PlumedForce
        print(f"\nAdd plumed: {inputs.plumed_file}")
        def SetPlumed(system, plumed_file):
            with open(plumed_file, 'r') as f:
                script = f.read()
            # print(script)
            system.addForce(PlumedForce(script))
        SetPlumed(system, inputs.plumed_file)

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
    
    # Production MD simulation
        if inputs.nstep > 0:
            start_time=datetime.datetime.now()
            
            # Determine trajectory output format (DCD or XTC)
            if inputs.odcd and not inputs.oxtc:
                TrajReporter = DCDReporter  # Binary DCD format
                otraj = inputs.odcd
            elif inputs.oxtc and not inputs.odcd:
                TrajReporter = XTCReporter  # Compressed XTC format 
                otraj = inputs.oxtc
            else:
                raise ValueError("Error: Please specify either odcd or oxtc!")
    
            # Set simulation starting point (either from input or checkpoint)
            if inputs.append == 'no': 
                b_step = inputs.b_step  # Starting step from input
                simulation.context.setStepCount(b_step)
                simulation.context.setTime(inputs.dt * b_step)
            else:
                b_step = simulation.context.getStepCount()  # Continue from current step
                
            e_step = inputs.nstep  # Total number of steps to run
            be_step = e_step - b_step  # Steps remaining to reach target
            
            print(f"\nMD run: begin {b_step}, end {e_step}, total {be_step}")
    
            # Setup trajectory and checkpoint reporters
            if inputs.nstdcd > 0:
                if inputs.append == 'yes':
                    # Append to existing trajectory file
                    simulation.reporters.append(TrajReporter(otraj, inputs.nstdcd, append=True))
                else:
                    # Create new trajectory file with backup
                    BackupFile(otraj)
                    simulation.reporters.append(TrajReporter(otraj, inputs.nstdcd))
                
                # Setup checkpoint file reporter
                BackupFile(inputs.ochk)
                simulation.reporters.append(
                    CheckpointReporter(inputs.ochk, inputs.nstdcd, writeState=True)
                )
                
            # Add console output reporter
            simulation.reporters.append(
                StateDataReporter(
                    sys.stdout, 
                    inputs.nstout,  # Reporting interval
                    step=True, 
                    time=True, 
                    potentialEnergy=True, 
                    temperature=True, 
                    progress=True,
                    remainingTime=True, 
                    speed=True, 
                    totalSteps=e_step, 
                    separator='\t'
                )
            )
            
            # Run the simulation
            try:
                simulation.step(be_step)
                be_step = 0  # Reset remaining steps counter
            except KeyboardInterrupt:
                print("The task has been canceled!")
                ReportTime(start_time)
                Cleanup(signal.SIGINT, simulation, inputs)
    
    # Write final output files
    if inputs.output:
        WriteOutput(inputs.output, simulation, inputs.input)
    if inputs.output_pdb:
        WriteOutput(inputs.output_pdb, simulation, inputs.input)
    if inputs.ochk:
        WriteCheckPoint(simulation, inputs.ochk)
    
    print("The task has finished!")
    ReportTime(start_time)

if __name__ == "__main__":
    # Set up command line argument parser
    parser = argparse.ArgumentParser(
        description="Run CTGoMartini molecular dynamics simulation using OpenMM"
    )
    
    # Define required input parameter file argument
    parser.add_argument(
        '-i', 
        dest='inpfile', 
        help='Input parameter file containing simulation settings', 
        required=True
    )
    
    # Parse command line arguments
    args = parser.parse_args()
    
    # Execute the main MD simulation function with provided input file
    mdrun(args.inpfile)
    
