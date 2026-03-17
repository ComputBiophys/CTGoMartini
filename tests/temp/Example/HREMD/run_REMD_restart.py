#!/usr/bin/env python3

"""
Authors: Song Yang
Date: 20250703
Purpose: This script is designed to run replica exchange simulation using OpenMM.
Version: 0.2.0
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


import openmmtools
import numpy as np
from openmmtools.multistate import (MultiStateReporter, MultiStateSampler,
                                    ReplicaExchangeAnalyzer,
                                    ReplicaExchangeSampler)
from openmmtools.cache import global_context_cache

import warnings
warnings.filterwarnings("ignore")


from Topology import *


def run(inpfile, replica_exchange_parameters, exchange_frequency=200, output_data='output.nc'):
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

    # Production MD
    if inputs.nstep > 0:
        total_simulation_time = inputs.nstep * inputs.dt * u.picosecond
        simulation_time_step = inputs.dt * u.picosecond

        simulation_steps = int(np.floor(total_simulation_time / simulation_time_step))
        exchange_attempts = int(np.floor(simulation_steps / exchange_frequency))


        # Load in the reporter from the original simulation:
        #reporter = MultiStateReporter(output_data, open_mode="r+")
        simulation = ReplicaExchangeSampler.from_storage(output_data)
        simulation._online_analysis_interval = 10000

        # Compute how many more iterations are needed:
        n_iter_remain = exchange_attempts-simulation.iteration

        print(f"Continuing OpenMM replica exchange simulation from iteration {simulation.iteration}")
        print(f"Time step: {simulation_time_step}")
        print(f"New iterations: {n_iter_remain}",flush=True)

        simulation.extend(n_iterations=n_iter_remain)

    return 

run('md.inp',
    replica_exchange_parameters=None, 
    exchange_frequency=250, 
    output_data='output.nc')
