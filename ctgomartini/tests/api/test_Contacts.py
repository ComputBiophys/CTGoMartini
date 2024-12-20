from ctgomartini.api import MartiniTopFile
from function import *

def CompareResults(working_dir, epsilon_r = 15):
    os.chdir(os.path.join(working_dir, "Contacts"))
    strfile = "minimized.gro"
    topfile = "system.top"

    Clean() # Remove the forces.dat and energy.dat
    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy=Load_energy(clean=False)
    omm_forces=Load_forces(clean=False)
        
    # gmx
    os.chdir(os.path.join(working_dir, "LJ"))

    # Use the results from the last run
    # Clean()
    # GMX_set(strfile=strfile, topfile=topfile, CreateMDP=True, indexfile=None, double_precision=True)
    # GMX_run()

    gmx_energy=Load_energy(clean=False)
    gmx_forces=Load_forces(clean=False)   

    # Compare
    print("########################################")
    result_energy=Compare_energy(omm_energy[:,1:], gmx_energy[:,1:], isPrint=True)
    result_forces=Compare_forces(omm_forces[:,1:], gmx_forces[:,1:], isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match.")
    

class TestMartiniTopology:
    """
    Test that the MartiniTopology class is instantiated
    """
    path = os.path.dirname(__file__)
    
    def test_KALP(self):
        working_dir = os.path.join(self.path, "../data/Contacts/KALP")
        CompareResults(working_dir, epsilon_r = 15)    

    def test_1GB1(self):
        working_dir = os.path.join(self.path, "../data/Contacts/1GB1")
        CompareResults(working_dir, epsilon_r = 15)   

    def test_1UBQ(self):
        working_dir = os.path.join(self.path, "../data/Contacts/1UBQ")
        CompareResults(working_dir, epsilon_r = 15)   
    
    def test_GlnBP_Open(self):
        working_dir = os.path.join(self.path, "../data/Contacts/GlnBP_Open")
        CompareResults(working_dir, epsilon_r = 15)   

    def test_GlnBP_Closed(self):
        working_dir = os.path.join(self.path, "../data/Contacts/GlnBP_Closed")
        CompareResults(working_dir, epsilon_r = 15)   

    def test_AdK_Closed(self):
        working_dir = os.path.join(self.path, "../data/Contacts/AdK_Closed")
        CompareResults(working_dir, epsilon_r = 15)   

    def test_Beta2AR_Active(self):
        working_dir = os.path.join(self.path, "../data/Contacts/Beta2AR_Active")
        CompareResults(working_dir, epsilon_r = 15)   

    def test_TREK1_Down(self):
        working_dir = os.path.join(self.path, "../data/Contacts/TREK1_Down")
        CompareResults(working_dir, epsilon_r = 15)   