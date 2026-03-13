import numpy as np
import os
from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer
from openmm import unit
kB = (unit.MOLAR_GAS_CONSTANT_R).in_units_of(unit.kilojoule / (unit.kelvin * unit.mole))
# print(kB)

from pymbar import timeseries

def ReportExchangeRatio(output_file='output.nc'):
    reporter = MultiStateReporter(output_file, open_mode="r")
    n_accepted_matrix, n_proposed_matrix = reporter.read_mixing_statistics()

    # Because use the swap-neighbors method to exchange states, only half of the exchanges happend for state i and i+1
    ratio = n_accepted_matrix.sum(axis=0)/(n_accepted_matrix.shape[0]/2) 
    for i in range(ratio.shape[0]-1):
        j = i+1
        print(f'exchange ratio between replica {i} and {j} is {ratio[i,j]}')
    reporter.close()
ReportExchangeRatio('output.nc')

# os.system('cp -a output.nc output_bk.nc')
# ReportExchangeRatio('output_bk.nc')