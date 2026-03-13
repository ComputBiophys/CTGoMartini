#!/usr/bin/env python3
"""Run replica exchange molecular dynamics simulation using OpenMM.

Authors: Song Yang
Date: 20250703
Purpose: This script is designed to run replica exchange simulation using OpenMM.
Version: 0.2.0
"""

from __future__ import annotations

import argparse
import datetime
import os
import signal
import sys
import warnings
from typing import Any

import MDAnalysis as mda
import numpy as np
import openmm as mm
import openmm.unit as u
import openmmtools
from openmm.app import *
from openmmtools.cache import global_context_cache
from openmmtools.multistate import (
    MultiStateReporter,
    MultiStateSampler,
    ReplicaExchangeAnalyzer,
    ReplicaExchangeSampler,
)

from ctgomartini.api import MartiniTopFile
from ctgomartini.utils import read_inputs

from .Topology import *

warnings.filterwarnings("ignore")


def run(
    inpfile: str,
    replica_exchange_parameters: dict[str, Any] | None,
    exchange_frequency: int = 200,
    output_data: str = "output.nc",
) -> None:
    """Execute CTGoMartini replica exchange molecular dynamics simulation.

    This function loads simulation parameters from an input file, sets up the
    OpenMM platform, and continues a replica exchange simulation from a
    previous checkpoint.

    Args:
        inpfile: Path to input parameter file containing simulation settings.
        replica_exchange_parameters: Dictionary of replica exchange specific
            parameters. Currently reserved for future use.
        exchange_frequency: Number of steps between exchange attempts.
            Defaults to 200.
        output_data: Path to the output NetCDF file containing simulation
            data from previous run. Defaults to 'output.nc'.

    Returns:
        None: Outputs are written to files specified in input parameters.
    """
    start_time = datetime.datetime.now()

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
        # reporter = MultiStateReporter(output_data, open_mode="r+")
        simulation = ReplicaExchangeSampler.from_storage(output_data)

        # Compute how many more iterations are needed:
        n_iter_remain = exchange_attempts - simulation.iteration

        print(f"Continuing OpenMM replica exchange simulation from iteration {simulation.iteration}")
        print(f"Time step: {simulation_time_step}")
        print(f"New iterations: {n_iter_remain}", flush=True)

        simulation.extend(n_iterations=n_iter_remain)

    return


if __name__ == "__main__":
    run(
        "md.inp",
        replica_exchange_parameters=None,
        exchange_frequency=250,
        output_data="output.nc",
    )
