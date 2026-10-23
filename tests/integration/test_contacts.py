"""
Integration tests for contact-based interactions.

Tests Go-type contact interactions between native residues,
comparing OpenMM and GROMACS implementations.

Systems tested:
    - KALP (transmembrane peptide)
    - Small proteins (1GB1, 1UBQ)
    - Conformational states (GlnBP, AdK, Beta2AR, TREK1)
"""

import os
import shutil
from ctgomartini.topology import MartiniTopFile
from .function import *
from tests.conftest import WorkingDirectoryContext


def compare_results(working_dir, tmp_dir, epsilon_r=15):
    # Run the OpenMM calculation in a scratch copy so the tracked reference
    # energy/forces files under tests/fixtures stay untouched.
    scratch_dir = os.path.join(tmp_dir, "Contacts")
    shutil.copytree(os.path.join(working_dir, "Contacts"), scratch_dir, symlinks=False)

    with WorkingDirectoryContext(scratch_dir):
        strfile = "minimized.gro"
        topfile = "system.top"

        clean_files() # Remove the forces.dat and energy.dat
        simulation = setup_openmm_simulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
        calculate_openmm_state(strfile, simulation, set_vsite=True)

        omm_energy=load_energy_data(clean=False)
        omm_forces=load_forces_data(clean=False)

    # gmx
    with WorkingDirectoryContext(os.path.join(working_dir, "LJ")):
        # Use the results from the last run
        # clean_files()
        # setup_gromacs_simulation(strfile=strfile, topfile=topfile, CreateMDP=True, indexfile=None, double_precision=True)
        # run_gromacs_calculation()

        gmx_energy=load_energy_data(clean=False)
        gmx_forces=load_forces_data(clean=False)

    # Compare
    print("########################################")
    result_energy=compare_energy_values(omm_energy[:,1:], gmx_energy[:,1:], isPrint=True)
    result_forces=compare_forces_values(omm_forces[:,1:], gmx_forces[:,1:], isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match.")


class TestMartiniTopology:
    """
    Test that the MartiniTopology class is instantiated
    """
    path = os.path.dirname(__file__)

    def test_KALP(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/KALP")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_1GB1(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/1GB1")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_1UBQ(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/1UBQ")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_GlnBP_Open(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/GlnBP_Open")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_GlnBP_Closed(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/GlnBP_Closed")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_AdK_Closed(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/AdK_Closed")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_Beta2AR_Active(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/Beta2AR_Active")
        compare_results(working_dir, tmp_path, epsilon_r = 15)

    def test_TREK1_Down(self, tmp_path):
        working_dir = os.path.join(self.path, "../fixtures/contacts/TREK1_Down")
        compare_results(working_dir, tmp_path, epsilon_r = 15)
