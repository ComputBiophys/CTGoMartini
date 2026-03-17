#!/usr/bin/env python3

"""
Authors: Song Yang
Date: 20250528
Purpose: This script is designed to run replica exchange simulation using OpenMM.
Version: 0.2.0
"""


import warnings
warnings.filterwarnings("ignore")

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


import openmmtools
import numpy as np
from openmmtools.multistate import (MultiStateReporter, MultiStateSampler,
                                    ReplicaExchangeAnalyzer,
                                    ReplicaExchangeSampler)
from openmmtools.cache import global_context_cache


from Topology import *


def run(inpfile, replica_exchange_parameters, exchange_frequency=200, output_data='output.nc', unsampled_topfiles=None):
    """Main function to execute CTGoMartini molecular dynamics simulation.
    
    Args:
        inpfile (str): Path to input parameter file containing simulation settings
        
    Returns:
        None: Outputs are written to files specified in input parameters
    """
    start_time=datetime.datetime.now()

    # Load parameters
    print("Loading parameters")
    inputs = read_inputs(inpfile)

    # Platform
    platform, platformProperties = LoadPlatform(inputs)
    global_context_cache.platform = platform
    global_context_cache.platformProperties = platformProperties

    # Load positions, velocities, and box vectors
    conf, box_vectors = LoadStructure(inputs.input)
    positions = conf.getPositions()
    velocities = None

    if inputs.ichk:
        print(f"\nLoad checkpoint file: {inputs.ichk}")
        with open(inputs.ichk, 'r') as f:
            states = mm.XmlSerializer.deserialize(f.read())
        
        positions = states.getPositions()
        velocities = states.getVelocities()
        box_vectors = states.getPeriodicBoxVectors()
        
    # system list
    system_list = []
    target_mol = replica_exchange_parameters['molname']
    for (beta, C1, C2) in zip(replica_exchange_parameters['beta'], replica_exchange_parameters['C1'], replica_exchange_parameters['C2']):
        #Load topol
        top = MartiniTopFile(
            inputs.topol,
            periodicBoxVectors=box_vectors,
            defines=inputs.defines,
        )

        # Set replica exchange parameters
        top.moleculeTypes[target_mol].multiple_basin=[
                ['True', 'exp', '2', str(beta), str(C1), str(C2)]
            ]
        
        # Create system
        system = top.create_system(nonbonded_cutoff=inputs.nonbonded_cutoff * u.nanometer, epsilon_r=inputs.epsilon_r)

        # Check the charges of the system
        if top.charges != 0:
            print(f'Warning: The charges of the system are {top.charges} instead of 0.')
        # Check the number of atoms
        assert len(positions) == top.topology.getNumAtoms(), f"Error: Number of atoms in {inputs.input} is not the same as that from {inputs.topol}!"

        # Add restraints
        if inputs.gen_rest == 'yes': gen_restraints(inputs.input, inputs.atomname, inputs.fc, inputs.gen_rest_file)
        if inputs.rest == 'yes':     system = restraints(system, inputs)

        # Add plumed
        if inputs.plumed == 'yes':
            raise ValueError("Error: Plumed is not supported!")
            # from openmmplumed import PlumedForce
            # print(f"\nAdd plumed: {inputs.plumed_file}")
            # def SetPlumed(system, plumed_file):
            #     with open(plumed_file, 'r') as f:
            #         script = f.read()
            #     # print(script)
            #     system.addForce(PlumedForce(script))
            # SetPlumed(system, inputs.plumed_file)

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

        # Add the system into the list
        system_list.append(system)

    if unsampled_topfiles is not None:
        unsampled_system_list = []
        for topfile in unsampled_topfiles:
           #Load topol
            top = MartiniTopFile(
                topfile,
                periodicBoxVectors=box_vectors,
                defines=inputs.defines,
            )
            
            # Create system
            system = top.create_system(nonbonded_cutoff=inputs.nonbonded_cutoff * u.nanometer, epsilon_r=inputs.epsilon_r)
            
            # Check the charges of the system
            if top.charges != 0:
                print(f'Warning: The charges of the system are {top.charges} instead of 0.')
            # Check the number of atoms
            assert len(positions) == top.topology.getNumAtoms(), f"Error: Number of atoms in {inputs.input} is not the same as that from {inputs.topol}!"

            # Add restraints
            if inputs.gen_rest == 'yes': gen_restraints(inputs.input, inputs.atomname, inputs.fc, inputs.gen_rest_file)
            if inputs.rest == 'yes':     system = restraints(system, inputs)

            # Add plumed
            if inputs.plumed == 'yes':
                raise ValueError("Error: Plumed is not supported!")
                # from openmmplumed import PlumedForce
                # print(f"\nAdd plumed: {inputs.plumed_file}")
                # def SetPlumed(system, plumed_file):
                #     with open(plumed_file, 'r') as f:
                #         script = f.read()
                #     # print(script)
                #     system.addForce(PlumedForce(script))
                # SetPlumed(system, inputs.plumed_file)

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

            # Add the system into the list
            unsampled_system_list.append(system)

    # System Loading Finishes!
    print("\nLoading system finishes!")
    ReportTime(start_time)
    

    # Production MD
    if inputs.nstep > 0:
        total_simulation_time = inputs.nstep * inputs.dt * u.picosecond
        simulation_time_step = inputs.dt * u.picosecond

        simulation_steps = int(np.floor(total_simulation_time / simulation_time_step))
        exchange_attempts = int(np.floor(simulation_steps / exchange_frequency))

        num_replicas = len(system_list)
        sampler_states = list()
        thermodynamic_states = list()
        

        # Define thermodynamic states.
        for system in system_list:
            thermodynamic_state = openmmtools.states.ThermodynamicState(
                system=system, temperature=inputs.temp * u.kelvin
            )
            thermodynamic_states.append(thermodynamic_state)
            sampler_states.append(
                openmmtools.states.SamplerState(positions, velocities, box_vectors=box_vectors)
            ) 
        
        if unsampled_topfiles is not None:
            unsampled_thermodynamic_states = []
            for system in unsampled_system_list:
                thermodynamic_state = openmmtools.states.ThermodynamicState(
                    system=system, temperature=inputs.temp * u.kelvin
                )
                unsampled_thermodynamic_states.append(thermodynamic_state)


        # Create and configure simulation object.
        if not inputs.const_tol:
            inputs.const_tol = 1e-5
        move = openmmtools.mcmc.LangevinDynamicsMove(
            timestep=simulation_time_step,
            collision_rate=inputs.fric_coeff / u.picosecond,
            n_steps=exchange_frequency,
            reassign_velocities=False,
            constraint_tolerance=inputs.const_tol,
        )

        simulation = ReplicaExchangeSampler(
            mcmc_moves=move,
            number_of_iterations=exchange_attempts,
            replica_mixing_scheme='swap-neighbors',
        )

        if os.path.exists(output_data):
            os.remove(output_data)

        reporter = MultiStateReporter(output_data, checkpoint_interval=5)
        if unsampled_topfiles is not None:
            simulation.create(thermodynamic_states, sampler_states, reporter, unsampled_thermodynamic_states=unsampled_thermodynamic_states)
        else:
            simulation.create(thermodynamic_states, sampler_states, reporter)

        print("Running OpenMM replica exchange simulation...")
        print(f"Time step: {simulation_time_step}")
        print(f"Iterations: {exchange_attempts}",flush=True)

        # try:
        #     simulation.run()
        # except BaseException:
        #     print("Replica exchange simulation failed, try verifying your model/simulation settings.")
        simulation.run()

    return 

replica_exchange_parameters = {
    'molname': 'GlnBP',
    'beta': [1/300] * 11,
    'C1': np.linspace(-480, -80, 11),
    'C2': [0] * 11
}

unsampled_topfiles = ['system_stateA.top', 'system_stateB.top']
run('md.inp',
    replica_exchange_parameters=replica_exchange_parameters, 
    exchange_frequency=250, 
    output_data='output.nc',
    unsampled_topfiles=unsampled_topfiles)
