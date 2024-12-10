from function import *
from ctgomartini.api import GenSBPTop, MartiniTopFile
from ctgomartini.func import WriteItp


def CompareResults(working_dir, epsilon_r = 15):
    os.chdir(os.path.join(working_dir, "openmm"))
    strfile = "minimized.gro"
    topfile = "system.top"

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy=Load_energy(clean=True)
    omm_forces=Load_forces(clean=True)
        
    # gmx
    os.chdir(os.path.join(working_dir, "gmx"))

    # Use the results from the last run
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
    # path = os.getcwd()

    def test_GlnBP_go_m3(self):
        working_dir = os.path.join(self.path, "../data/SingleBasinGoMartini/GlnBP_go_m3")

        # GenSBPTop
        os.chdir(os.path.join(working_dir, "openmm"))
        try:
            os.remove("gbp.top", "gbp_params.top")
        except:
            pass
        mols_list =[
            ['system_gen.top', 'gbp_open']
        ]
        sbmol_name = "gbp"
        sbmol =  GenSBPTop(mols_list, sbmol_name)
        WriteItp(sbmol)
        os.chdir(working_dir)

        CompareResults(working_dir, epsilon_r = 15)  
