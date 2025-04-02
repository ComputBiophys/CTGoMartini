from ctgomartini.api import MartiniTopFile
from function import *

from ctgomartini.data.run_ctgomartini import gen_restraints, restraints

def AddPositionRestraints(simulation):
    
    # Add restraints
    class inputs:
        input = 'minimized.gro'
        atomname = 'BB'
        fc = 1000.0
        rest_file = 'restraints.txt'
        rest = 'yes'
        gen_rest = 'yes'
        rest_ref = 'ref.gro'

    gen_restraints(inputs.input, inputs.atomname, inputs.fc, inputs.rest_file)
    system = restraints(simulation.system, inputs)
    simulation.context.reinitialize()


def Compare_OMM_GMX(working_dir, strfile='md.gro',  epsilon_r=15.0, **kwargs):
    print(working_dir)
    os.chdir(os.path.join(working_dir, "openmm"))
    # strfile = "md.gro"
    topfile = "system.top"

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)

    try:
        if kwargs['rest']:
            AddPositionRestraints(simulation)
    except KeyError:
        pass

    Clean()
    OMM_calStrfile(strfile, simulation, set_vsite=True)
    omm_energy=Load_energy(clean=False)
    omm_forces=Load_forces(clean=False)
    print(omm_energy)

    # gmx
    os.chdir(os.path.join(working_dir, "gmx"))
    # Clean()
    # GMX_set(strfile=strfile, topfile=topfile, reffile='ref.gro', CreateMDP=False, indexfile=None, double_precision=True)
    # GMX_set(strfile=strfile, topfile=topfile, reffile='ref.gro', CreateMDP=False, indexfile=None, double_precision=True)
    # GMX_run()

    gmx_energy=Load_energy(clean=False)
    gmx_forces=Load_forces(clean=False)
    print(gmx_energy)
    # Compare
    print("########################################")
    result_energy=Compare_energy(omm_energy[:,1:], gmx_energy[:,1:], isPrint=True)
    result_forces=Compare_forces(omm_forces[:,1:], gmx_forces[:,1:], isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match.")
    

class TestEnergyItemComparison:
    """
    Test EnergyItemComparison
    """
    path = os.path.dirname(__file__)

    def test_PullCode(self):
        working_dir = os.path.join(self.path, "../data/EnergyItemComparison/PullCode/GlnBP/")
        Compare_OMM_GMX(working_dir, strfile='md.gro', epsilon_r = 15)  

    def test_Restraints(self):
        working_dir = os.path.join(self.path, "../data/EnergyItemComparison/Restraints/GlnBP_Open/")
        Compare_OMM_GMX(working_dir, strfile='minimized.gro', epsilon_r = 15, rest=True)
    
