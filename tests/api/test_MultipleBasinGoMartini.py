from ctgomartini.api import MartiniTopFile
from function import *

def MBP_exp_energy_combine(energy1, energy2, beta, C1, C2):
    part1=np.exp(-beta*(energy1+C1))
    part2=np.exp(-beta*(energy2+C2))

    part=part1+part2
    energy=-np.log(part)/beta
    return energy

def MBP_exp_forces_combine(forces1, forces2, energy1, energy2, beta, C1, C2):
    part1=np.exp(-beta*(energy1+C1))
    part2=np.exp(-beta*(energy2+C2))

    part=part1+part2
    forces=part1/part*forces1+part2/part*forces2
    return forces

def MBP_ham_energy_combine(energy1, energy2, delta, C1, C2):
    dV=C2-C1
    part1=(energy1+energy2+dV)/2
    part2=(energy1-energy2-dV)/2

    energy=part1-np.sqrt(part2**2+delta**2)
    return energy

def MBP_ham_forces_combine(forces1, forces2, energy1, energy2, delta, C1, C2):
    dV=C2-C1
    alpha=(energy1-energy2-dV)/2
    C=1/2*(1-alpha/np.sqrt(alpha**2+delta**2))

    forces=C*forces1+(1-C)*forces2
    return forces

def CalculateEnergyForces(working_dir, strfile='ions.gro', epsilon_r=15):
    os.chdir(working_dir)
    strfile = strfile
    topfile = "system.top"

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy=Load_energy(clean=True)
    omm_forces=Load_forces(clean=True)    
    return omm_energy, omm_forces

    
def GetMBPParameter(working_dir, topfile='system.top'):
    """
    get MBP Parameters
    """
    os.chdir(working_dir)
    top = MartiniTopFile(topfile)
    MBP_paramters = []
    for molecule_name, molecule_type in top.moleculeTypes.items():
        if 'multiple_basin' in molecule_type._topology:
            MBP_paramters.append(molecule_type._topology['multiple_basin'][0])
    return MBP_paramters

def CompareResults(base_dir, strfile='ions.gro', epsilon_r=15):
    working_dir = os.path.join(base_dir, 'StateA')
    energy_stateA, forces_stateA = CalculateEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)

    working_dir = os.path.join(base_dir, 'StateB')
    energy_stateB, forces_stateB = CalculateEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)

    working_dir = os.path.join(base_dir, 'EXP')
    energy_exp, forces_exp = CalculateEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    exp_param = GetMBPParameter(working_dir, topfile='system.top')[0]

    working_dir = os.path.join(base_dir, 'HAM')
    energy_ham, forces_ham = CalculateEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    ham_param = GetMBPParameter(working_dir, topfile='system.top')[0]

    # Energy and forces comparison for Exponential mixing scheme
    print("Exponential mixing scheme for multiple baisn popential")
    print("Parameters: ",exp_param)
    energy_exp_cal = MBP_exp_energy_combine(energy_stateA[:,1], energy_stateB[:,1], beta=eval(exp_param[3]), C1=eval(exp_param[4]), C2=eval(exp_param[5]))
    forces_exp_cal = MBP_exp_forces_combine(forces_stateA[:,1:], forces_stateB[:,1:], energy_stateA[:,1], energy_stateB[:,1], beta=1/1000, C1=-50, C2=0)

    result_energy=Compare_energy(energy_exp_cal, energy_exp[:,1:], isPrint=True)
    result_forces=Compare_forces(forces_exp_cal, forces_exp[:,1:], isPrint=True)
    
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match for Exponential mixing scheme.")

    # Energy and forces comparison for Hamiltonian mixing scheme
    print("Hamiltonian mixing scheme for multiple baisn popential")
    print("Parameters: ",ham_param)
    energy_ham_cal = MBP_ham_energy_combine(energy_stateA[:,1], energy_stateB[:,1], delta=eval(ham_param[3]), C1=eval(ham_param[4]), C2=eval(ham_param[5]))
    forces_ham_cal = MBP_ham_forces_combine(forces_stateA[:,1:], forces_stateB[:,1:], energy_stateA[:,1], energy_stateB[:,1], delta=10, C1=-50, C2=0)

    result_energy=Compare_energy(energy_ham_cal, energy_ham[:,1:], isPrint=True)
    result_forces=Compare_forces(forces_ham_cal, forces_ham[:,1:], isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match for Hamiltonian mixing scheme.")

class TestMultipleBasinGoMartiniTopology:
    """
    Test the multiple-basin GoMartini topology
    """
    path = os.path.dirname(__file__)

    def test_GlnBP_open(self):
        base_dir = os.path.join(self.path, "../data/MultipleBasinGoMartini/GlnBP")
        CompareResults(base_dir, strfile='ions_open.gro', epsilon_r = 15)  
    
    def test_GlnBP_closed(self):
        base_dir = os.path.join(self.path, "../data/MultipleBasinGoMartini/GlnBP")
        CompareResults(base_dir, strfile='ions_closed.gro', epsilon_r = 15)      

    