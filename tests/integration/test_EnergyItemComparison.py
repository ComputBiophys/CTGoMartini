"""
Integration tests for energy item comparison.

Verifies that OpenMM and GROMACS produce consistent energies for:
    - Pull code simulations
    - Position restraints

These tests ensure force field implementations are equivalent.
"""

import os
from ctgomartini.topology import MartiniTopFile
from .function import *
from ctgomartini.simulation import generate_restraints, add_restraints
from tests.conftest import WorkingDirectoryContext


def add_position_restraints(simulation):
    
    # Add restraints
    class inputs:
        input = 'minimized.gro'
        atomname = 'BB'
        fc = 1000.0
        rest_file = 'restraints.txt'
        rest = 'yes'
        gen_rest = 'yes'
        rest_ref = 'ref.gro'

    generate_restraints(inputs.input, inputs.atomname, inputs.fc, inputs.rest_file)
    system = add_restraints(simulation.system, inputs)
    simulation.context.reinitialize()


def compare_openmm_gromacs(working_dir, strfile='md.gro', epsilon_r=15.0, **kwargs):
    print(working_dir)
    
    with WorkingDirectoryContext(os.path.join(working_dir, "openmm")):
        # strfile = "md.gro"
        topfile = "system.top"

        simulation = setup_openmm_simulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)

        try:
            if kwargs['rest']:
                add_position_restraints(simulation)
        except KeyError:
            pass

        clean_files()
        calculate_openmm_state(strfile, simulation, set_vsite=True)
        omm_energy=load_energy_data(clean=False)
        omm_forces=load_forces_data(clean=False)
    print(omm_energy)

    # gmx
    with WorkingDirectoryContext(os.path.join(working_dir, "gmx")):
        # clean_files()
        # setup_gromacs_simulation(strfile=strfile, topfile=topfile, reffile='ref.gro', CreateMDP=False, indexfile=None, double_precision=True)
        # setup_gromacs_simulation(strfile=strfile, topfile=topfile, reffile='ref.gro', CreateMDP=False, indexfile=None, double_precision=True)
        # run_gromacs_calculation()

        gmx_energy=load_energy_data(clean=False)
        gmx_forces=load_forces_data(clean=False)
    print(gmx_energy)
    
    # Compare
    print("########################################")
    result_energy=compare_energy_values(omm_energy[:,1:], gmx_energy[:,1:], isPrint=True)
    result_forces=compare_forces_values(omm_forces[:,1:], gmx_forces[:,1:], isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match.")
    

class TestEnergyItemComparison:
    """
    Test EnergyItemComparison
    """
    path = os.path.dirname(__file__)

    def test_pull_code(self):
        working_dir = os.path.join(self.path, "../fixtures/EnergyItemComparison/PullCode/GlnBP/")
        compare_openmm_gromacs(working_dir, strfile='md.gro', epsilon_r=15)

    def test_restraints(self):
        working_dir = os.path.join(self.path, "../fixtures/EnergyItemComparison/Restraints/GlnBP_Open/")
        compare_openmm_gromacs(working_dir, strfile='minimized.gro', epsilon_r=15, rest=True)
    
