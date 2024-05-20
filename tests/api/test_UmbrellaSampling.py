from ctgomartini.api import MartiniTopFile
from function import *

def Compare_OMM_GMX(working_dir, epsilon_r=15.0):
    print(working_dir)
    os.chdir(os.path.join(working_dir, "openmm"))
    strfile = "md.gro"
    topfile = "system.top"

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy=Load_energy(clean=True)
    omm_forces=Load_forces(clean=True)
    print(omm_energy)
    # gmx
    os.chdir(os.path.join(working_dir, "gmx"))

    # GMX_set(CreateMDP=False, double_precision=True)
    GMX_run()

    gmx_energy=Load_energy(clean=False)
    gmx_forces=Load_forces(clean=False)
    print(gmx_energy)
    # Compare
    print("########################################")
    result_energy=Compare_energy(omm_energy[:,1:], gmx_energy[:,1:], isPrint=True)
    result_forces=Compare_forces(omm_forces[:,1:], gmx_forces[:,1:], isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match.")
    

class TestUmbrellaSampling:
    """
    Test UmbrellaSampling
    """
    path = os.path.dirname(__file__)

    def test_GlnBP(self):
        working_dir = os.path.join(self.path, "../data/UmbrellaSampling/GlnBP/")
        Compare_OMM_GMX(working_dir, epsilon_r = 15)  
    
