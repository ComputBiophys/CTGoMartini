from ctgomartini.api import MartiniTopFile
from ctgomartini.func import WriteItp
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


def Compare_energy(energy1, energy2, isPrint=True):
    energy1 = float(energy1)
    energy2 = float(energy2)

    relative_energy_error = abs(energy1 - energy2) * \
        2 / (abs(energy1)+abs(energy2))
    abs_energy_error = abs(energy1 - energy2)

    if isPrint:
        print('Energy Compare')
        print(f'Absolute error: {abs_energy_error:.2e}')
        print(f'Relative error: {relative_energy_error:.2e}')

    if relative_energy_error <= 1e-5 or abs_energy_error <= 1e-3:
        if isPrint:
            print("Energies match!")
        result = True
    else:
        if isPrint:
            print("Error: Energies do not match!")
        result = False
    
    return result, abs_energy_error, relative_energy_error
    
def Compare_forces(forces1, forces2, isPrint=True):
    average = 0.5*np.linalg.norm(forces1, axis=1) + 0.5*np.linalg.norm(forces2, axis=1)

    relative_force_error = np.linalg.norm(forces1 - forces2, axis=1) / average
    relative_force_error = np.nan_to_num(relative_force_error, nan=0)
    max_relative_force_error, max_relative_force_error_index = relative_force_error.max(
    ), relative_force_error.argmax()

    abs_force_error = np.linalg.norm(forces1 - forces2, axis=1)
    max_abs_force_error, max_abs_force_error_index = abs_force_error.max(
    ), abs_force_error.argmax()

    atol = 1e-5
    rtol = 1e-4
    # allclose = abs_force_error-(atol+rtol*average)
    allclose = abs_force_error-(atol+rtol*np.linalg.norm(forces2, axis=1)) # use gmx forces as the standard
    mask = np.linalg.norm(forces2, axis=1) != 0 # exclude the force == 0
    max_allclose = allclose[mask].max()
    max_allclose_index = np.where(allclose == max_allclose)[0][0]
    max_abs_force_error = abs_force_error[max_allclose_index]
    max_relative_force_error = relative_force_error[max_allclose_index]
    if isPrint:
        print('###Forces Compare###')
        print(f'Max absolute error: {max_abs_force_error:.2e}')
        print(f'Max relative error: {max_relative_force_error:.2e}')
        print(f'      Max allclose: {max_allclose:.2e}')

    if max_allclose <= 0:
        if isPrint:
            print("Forces match!")
        result = True
    else:
        if isPrint:
            print("Error: Forces do not match!")
        result = False 

    return result, max_abs_force_error, max_relative_force_error, max_allclose
def CalculateOMMEnergyForces(working_dir, strfile='ions.gro', epsilon_r=15):
    os.chdir(working_dir)
    strfile = strfile
    topfile = "system.top"
    prefix = strfile.split('/')[-1].split('.')[0]

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, prefix=prefix, set_vsite=True)

    omm_energy=Load_energy(prefix=prefix, clean=True)
    omm_forces=Load_forces(prefix=prefix, clean=True)    
    return omm_energy, omm_forces

def CalculateGMXEnergyForces(working_dir, strfile='ions.gro', epsilon_r=15):
    os.chdir(working_dir)
    strfile = strfile
    topfile = "system.top"
    prefix = strfile.split('/')[-1].split('.')[0]
    print(prefix)

    GMX_set(strfile=strfile,trjfile=strfile,topfile=topfile, indexfile=None, prefix=prefix, CreateMDP=True, double_precision=True)
    GMX_run(prefix=prefix)

    gmx_energy=Load_energy(prefix=prefix, clean=False)
    gmx_forces=Load_forces(prefix=prefix, clean=False)    
    return gmx_energy, gmx_forces
def SetMBPParameter(working_dir, topfile='system.top', MBP_parameters = [['True', 'exp', '2', '1/300', '-300', '0']]):
    """
    get MBP Parameters
    """
    os.chdir(working_dir)
    top = MartiniTopFile(topfile)
    MBP_paramters = []
    for molecule_name, molecule_type in top.moleculeTypes.items():
        if 'multiple_basin' in molecule_type._topology:
            molecule_type._topology['multiple_basin'] = MBP_parameters
            WriteItp(molecule_type)

# # GMX Energy Initialization
# path = os.path.dirname(__file__)
# base_dir = os.path.join(path, "../data/MultipleBasinGoMartini_v2/GlnBP/")

# # Open
# working_dir = os.path.join(base_dir, 'Open')
# epsilon_r = 15
# for i in range(10):
#     strfile = os.path.join(base_dir,f'Strfiles/GlnBP_No{i}.gro')
#     energy_stateA, forces_stateA = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)

# # Closed
# working_dir = os.path.join(base_dir, 'Closed')
# epsilon_r = 15
# for i in range(10):
#     strfile = os.path.join(base_dir,f'Strfiles/GlnBP_No{i}.gro')
#     energy_stateA, forces_stateA = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)

def CompareResults_EXP(base_dir, strfile='ions.gro', epsilon_r=15, 
                   MBP_parameters = [['True', 'exp', '2', '1/300', '-300', '0']]):
    
    working_dir = os.path.join(base_dir, 'Open')
    os.chdir(working_dir)
    # energy_stateA, forces_stateA = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    prefix = strfile.split('/')[-1].split('.')[0]
    energy_stateA=Load_energy(prefix=prefix, clean=False)
    forces_stateA=Load_forces(prefix=prefix, clean=False)    
    os.chdir(base_dir)

    working_dir = os.path.join(base_dir, 'Closed')
    os.chdir(working_dir)
    # energy_stateB, forces_stateB = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    prefix = strfile.split('/')[-1].split('.')[0]
    energy_stateB=Load_energy(prefix=prefix, clean=False)
    forces_stateB=Load_forces(prefix=prefix, clean=False)    
    os.chdir(base_dir)
    
    working_dir = os.path.join(base_dir, 'EXP')
    SetMBPParameter(working_dir, topfile='system.top', MBP_parameters = MBP_parameters)
    energy_exp, forces_exp = CalculateOMMEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    exp_param = GetMBPParameter(working_dir, topfile='system.top')[0]


    # Energy and forces comparison for Exponential mixing scheme
    print("Exponential mixing scheme for multiple baisn popential")
    print("Parameters: ",exp_param)
    energy_exp_cal = MBP_exp_energy_combine(energy_stateA[:,1], energy_stateB[:,1], beta=eval(exp_param[3]), C1=eval(exp_param[4]), C2=eval(exp_param[5]))
    forces_exp_cal = MBP_exp_forces_combine(forces_stateA[:,1:], forces_stateB[:,1:], energy_stateA[:,1], energy_stateB[:,1], beta=eval(exp_param[3]), C1=eval(exp_param[4]), C2=eval(exp_param[5]))

    result_energy, abs_energy_error, relative_energy_error =Compare_energy(energy_exp[:,1:], energy_exp_cal, isPrint=True)
    result_forces, max_abs_force_error, max_relative_force_error, max_allclose =Compare_forces(forces_exp[:,1:], forces_exp_cal, isPrint=True)
    
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match for Exponential mixing scheme.")

    return [abs_energy_error, relative_energy_error, max_abs_force_error, max_relative_force_error, max_allclose]

def CompareResults_HAM(base_dir, strfile='ions.gro', epsilon_r=15, 
                   MBP_parameters = [['True', 'ham', '2', '300', '-300', '0']]):
    working_dir = os.path.join(base_dir, 'Open')
    os.chdir(working_dir)
    # energy_stateA, forces_stateA = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    prefix = strfile.split('/')[-1].split('.')[0]
    energy_stateA=Load_energy(prefix=prefix, clean=False)
    forces_stateA=Load_forces(prefix=prefix, clean=False)    
    os.chdir(base_dir)

    working_dir = os.path.join(base_dir, 'Closed')
    os.chdir(working_dir)
    # energy_stateB, forces_stateB = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    prefix = strfile.split('/')[-1].split('.')[0]
    energy_stateB=Load_energy(prefix=prefix, clean=False)
    forces_stateB=Load_forces(prefix=prefix, clean=False)    
    os.chdir(base_dir)


    working_dir = os.path.join(base_dir, 'HAM')
    SetMBPParameter(working_dir, topfile='system.top', MBP_parameters = MBP_parameters)
    energy_ham, forces_ham = CalculateOMMEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    ham_param = GetMBPParameter(working_dir, topfile='system.top')[0]


    # Energy and forces comparison for Hamiltonian mixing scheme
    print("Hamiltonian mixing scheme for multiple baisn popential")
    print("Parameters: ",ham_param)
    energy_ham_cal = MBP_ham_energy_combine(energy_stateA[:,1], energy_stateB[:,1], delta=eval(ham_param[3]), C1=eval(ham_param[4]), C2=eval(ham_param[5]))
    forces_ham_cal = MBP_ham_forces_combine(forces_stateA[:,1:], forces_stateB[:,1:], energy_stateA[:,1], energy_stateB[:,1], delta=eval(ham_param[3]), C1=eval(ham_param[4]), C2=eval(ham_param[5]))

    result_energy, abs_energy_error, relative_energy_error = Compare_energy(energy_ham[:,1:], energy_ham_cal, isPrint=True)
    result_forces, max_abs_force_error, max_relative_force_error, max_allclose =Compare_forces(forces_ham[:,1:], forces_ham_cal, isPrint=True)
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match for Hamiltonian mixing scheme.")
    
    return [abs_energy_error, relative_energy_error, max_abs_force_error, max_relative_force_error, max_allclose]
from multiprocessing import Pool
from functools import partial

def Compare(i, base_dir, Mixing_function, mixing_parameters):
    strfile = os.path.join(base_dir,f'Strfiles/GlnBP_No{i}.gro')
    print(strfile)
    Mixing_function(base_dir, strfile=strfile, epsilon_r = 15, MBP_parameters=mixing_parameters)

    

class TestMultipleBasinGoMartiniTopology:
    """
    Test the multiple-basin GoMartini topology
    """
    path = os.path.dirname(__file__)

    def test_GlnBP_EXP(self):
        base_dir = os.path.join(self.path, "../data/MultipleBasinGoMartini_v2/GlnBP/")
        EXP_parameters_list = [[['True', 'exp', '2', '1/300', '-300', '0']],
                                [['True', 'exp', '2', '1/500', '-300', '0']],
                                [['True', 'exp', '2', '1/300', '-600', '0']],
                                ]

        for j, EXP_parameters in enumerate(EXP_parameters_list):
            print(os.getcwd())
            print(EXP_parameters)
            Compare_per_frame = partial(Compare, base_dir=base_dir, Mixing_function=CompareResults_EXP, mixing_parameters=EXP_parameters)
            with Pool(10) as pool:
                pool.map(Compare_per_frame, range(10))
    
    def test_GlnBP_HAM(self):
        base_dir = os.path.join(self.path, "../data/MultipleBasinGoMartini_v2/GlnBP/")
        HAM_parameters_list = [[['True', 'ham', '2', '100', '-200', '0']],
                                [['True', 'ham', '2', '600', '-200', '0']],
                                [['True', 'ham', '2', '100', '-900', '0']],
                                ]

        for j, HAM_parameters in enumerate(HAM_parameters_list):
            print(os.getcwd())
            print(HAM_parameters)
            Compare_per_frame = partial(Compare, base_dir=base_dir, Mixing_function=CompareResults_HAM, mixing_parameters=HAM_parameters)
            with Pool(10) as pool:
                pool.map(Compare_per_frame, range(10))
