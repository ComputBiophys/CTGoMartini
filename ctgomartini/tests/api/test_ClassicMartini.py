from ctgomartini.api import MartiniTopFile
from function import *

def CompareResults(working_dir, epsilon_r = 15):
    os.chdir(os.path.join(working_dir, "openmm"))
    strfile = "minimized.gro"
    topfile = "system.top"

    Clean() # Remove the forces.dat and energy.dat
    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy=Load_energy(clean=False)
    omm_forces=Load_forces(clean=False)
        
    # gmx
    os.chdir(os.path.join(working_dir, "gmx"))

    # Use the results from the last run
    # Clean()
    # GMX_set(CreateMDP=False, double_precision=True)
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

    def test_pol_water(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/pol_water")
        CompareResults(working_dir, epsilon_r = 2.5)  
    
    def test_thy2_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/thy2_m3")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1k6u_en_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/1k6u_en_m3")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1j4n_en_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/1j4n_en_m3")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_MOL1(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/MOL1")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_PEG(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/PEG")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_simple_lipid(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/simple_lipid")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1k6u_go_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/1k6u_go_m3")
        CompareResults(working_dir, epsilon_r = 15)  

    def test_CLPR(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/mono/CLPR")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_3HT(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/mono/3HT")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_4NIAN(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/mono/4NIAN")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_2T(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/poly/2T")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_CNAP(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/poly/CNAP")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_TDMBI(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/poly/TDMBI")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_CAFF(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/poly/CAFF")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_NDMBI(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/poly/NDMBI")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_ANTH(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/small_mols_m3/poly/ANTH")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1ubq_en_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/1ubq_en_m3")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_popc_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/popc_m3")
        CompareResults(working_dir, epsilon_r = 15)  

    def test_aden_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/aden_m3")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1NACL(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/others_m3/1NACL")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1POPC(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/others_m3/1POPC")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1DOD(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/others_m3/1DOD")
        CompareResults(working_dir, epsilon_r = 15)  

    def test_2W(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/others_m3/2W")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_GlnBP_go_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/GlnBP_go_m3")
        CompareResults(working_dir, epsilon_r = 15)  

    def test_GlnBP_go_m3_pairs(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/GlnBP_go_m3_pairs")
        CompareResults(working_dir, epsilon_r = 15)  

    def test_beta_hairpin_go_m2(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/beta-hairpin_go_m2")
        CompareResults(working_dir, epsilon_r = 15)  

    def test_protein(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/protein")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_1ubq_go_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/1ubq_go_m3")
        CompareResults(working_dir, epsilon_r = 15)  
    
    def test_thym_m3(self):
        working_dir = os.path.join(self.path, "../data/ClassicMartini/thym_m3")
        CompareResults(working_dir, epsilon_r = 15)    

    # The following test all involve cholesterol.
    # It is not possible to have the constraints work
    # exactly the same in openmm as gromacs, so we have
    # temporarily removed these tests.    
    # def test_m2_chol2(self):
    #     working_dir = os.path.join(self.path, "../data/ClassicMartini/m2_chol2")
    #     CompareResults(working_dir, epsilon_r = 15)  
    
    # def test_complex_lipid(self):
    #     working_dir = os.path.join(self.path, "../data/ClassicMartini/complex_lipid")
    #     CompareResults(working_dir, epsilon_r = 15)  
    
    # def test_membrane_protein(self):
    #     working_dir = os.path.join(self.path, "../data/ClassicMartini/membrane_protein")
    #     CompareResults(working_dir, epsilon_r = 15)  
    
    # def test_m2_chol(self):
    #     working_dir = os.path.join(self.path, "../data/ClassicMartini/m2_chol")
    #     CompareResults(working_dir, epsilon_r = 15)  
    