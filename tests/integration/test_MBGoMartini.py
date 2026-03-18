"""
Integration tests for Multiple Basin Potential (MBP).

Tests exponential (EXP) and Hamiltonian (HAM) mixing schemes for
multi-state Go-Martini models, verifying energy and force consistency.

Systems tested:
    - GlnBP (Glutamine Binding Protein)
    - Beta2AR (Beta-2 Adrenergic Receptor)
"""

import os
import pickle
from ctgomartini.topology import MartiniTopFile
from ctgomartini.utils import write_itp
from .function import *
from tests.conftest import WorkingDirectoryContext


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
    with WorkingDirectoryContext(working_dir):
        top = MartiniTopFile(topfile)
        MBP_paramters = []
        for molecule_name, molecule_type in top.moleculeTypes.items():
            if 'multiple_basin' in molecule_type._topology:
                MBP_paramters.append(molecule_type._topology['multiple_basin'][0])
    return MBP_paramters

def SetMBPParameter(working_dir, topfile='system.top', MBP_parameters = ['True', 'exp', '2', '1/300', '-300', '0']):
    """
    Set MBP Parameters and write ITP file
    """
    with WorkingDirectoryContext(working_dir):
        top = MartiniTopFile(topfile)
        for molecule_name, molecule_type in top.moleculeTypes.items():
            if 'multiple_basin' in molecule_type._topology:
                molecule_type._topology['multiple_basin'] = [MBP_parameters]
                write_itp(molecule_type)


def CalculateOMMEnergyForces(working_dir, strfile='ions.gro', epsilon_r=15):
    with WorkingDirectoryContext(working_dir):
        topfile = "system.top"
        prefix = strfile.split('/')[-1].split('.')[0]

        simulation = OMM_setSimulation(strfile, topfile, epsilon_r=epsilon_r, temperature=310.15, double_precision=True)
        OMM_calStrfile(strfile, simulation, prefix=prefix, set_vsite=True)

        omm_energy=Load_energy(prefix=prefix, clean=False)
        omm_forces=Load_forces(prefix=prefix, clean=False)    
    return omm_energy, omm_forces

def CalculateGMXEnergyForces(working_dir, strfile='ions.gro', epsilon_r=15):
    with WorkingDirectoryContext(working_dir):
        topfile = "system.top"
        prefix = strfile.split('/')[-1].split('.')[0]
        print(prefix)

        GMX_set(strfile=strfile,trjfile=strfile,topfile=topfile, indexfile=None, prefix=prefix, CreateMDP=True, double_precision=True)
        GMX_run(prefix=prefix)

        gmx_energy=Load_energy(prefix=prefix, clean=False)
        gmx_forces=Load_forces(prefix=prefix, clean=False)    
    return gmx_energy, gmx_forces




def LoadStateData(base_dir, state_name, strfile):
    """Load energy and forces data for a given state."""
    working_dir = os.path.join(base_dir, state_name)
    with WorkingDirectoryContext(working_dir):
        prefix = strfile.split('/')[-1].split('.')[0]
        energy = Load_energy(prefix=prefix, clean=False)
        forces = Load_forces(prefix=prefix, clean=False)
    return energy, forces

def CalculateMBPResults(energy_stateA, forces_stateA, energy_stateB, forces_stateB,
                        energy_mbp, forces_mbp, mbp_param):
    """Calculate and compare MBP results."""
    print("Multiple-basin mixing scheme for multiple basin potential")
    print("Parameters: ", mbp_param)

    if mbp_param[1] == 'exp':
        energy_mbp_cal = MBP_exp_energy_combine(energy_stateA[:,1], energy_stateB[:,1], 
                                                beta=eval(mbp_param[3]), C1=eval(mbp_param[4]), C2=eval(mbp_param[5]))
        forces_mbp_cal = MBP_exp_forces_combine(forces_stateA[:,1:], forces_stateB[:,1:], 
                                                energy_stateA[:,1], energy_stateB[:,1], 
                                                beta=eval(mbp_param[3]), C1=eval(mbp_param[4]), C2=eval(mbp_param[5]))
    elif mbp_param[1] == 'ham':
        energy_mbp_cal = MBP_ham_energy_combine(energy_stateA[:,1], energy_stateB[:,1], 
                                                delta=eval(mbp_param[3]), C1=eval(mbp_param[4]), C2=eval(mbp_param[5]))
        forces_mbp_cal = MBP_ham_forces_combine(forces_stateA[:,1:], forces_stateB[:,1:], 
                                                energy_stateA[:,1], energy_stateB[:,1], 
                                                delta=eval(mbp_param[3]), C1=eval(mbp_param[4]), C2=eval(mbp_param[5]))
    else:
        raise ValueError("MBP_parameters[1] must be 'exp' or 'ham'")

    result_energy, abs_energy_error, relative_energy_error = Compare_energy(
        energy_mbp[:,1:], energy_mbp_cal, isPrint=True, isReturn=True)
    result_forces, max_allclose, abs_force_error, relative_force_error = Compare_forces(
        forces_mbp[:,1:], forces_mbp_cal, isPrint=True, isReturn=True)
    
    if not (result_energy and result_forces):
        raise AssertionError("Energies or forces do not match.")
    
    return (result_energy, abs_energy_error, relative_energy_error, 
            result_forces, max_allclose, abs_force_error, relative_force_error)

def CompareResults(base_dir, strfile='strfile_No1.gro', epsilon_r=15, 
                   MBP_parameters=['True', 'exp', '2', '1/300', '-300', '0']):
    """Compare results between MBP calculation and OpenMM simulation."""
    # Load State A and State B data (these are pre-calculated)
    energy_stateA, forces_stateA = LoadStateData(base_dir, 'StateA', strfile)
    energy_stateB, forces_stateB = LoadStateData(base_dir, 'StateB', strfile)
    
    # Calculate MBP energy and forces using OpenMM
    working_dir = os.path.join(base_dir, MBP_parameters[1].upper())
    energy_mbp, forces_mbp = CalculateOMMEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)
    mbp_param = GetMBPParameter(working_dir, topfile='system.top')[0]
    
    # Compare results
    return CalculateMBPResults(energy_stateA, forces_stateA, energy_stateB, forces_stateB,
                               energy_mbp, forces_mbp, mbp_param)
    


def CalSBPState(base_dir, strfile='strfile_No1.gro', epsilon_r=15):
    working_dir = os.path.join(base_dir, 'StateA')
    with WorkingDirectoryContext(working_dir):
        energy_stateA, forces_stateA = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)

    working_dir = os.path.join(base_dir, 'StateB')
    with WorkingDirectoryContext(working_dir):
        energy_stateB, forces_stateB = CalculateGMXEnergyForces(working_dir, strfile, epsilon_r=epsilon_r)


from multiprocessing import Pool, cpu_count
from functools import partial

def Compare(i, base_dir, mixing_parameters, prefix_strfile='GlnBP'):
    strfile = os.path.join(base_dir,f'Strfiles/{prefix_strfile}_No{i}.gro')
    print(strfile)
    result_energy, abs_energy_error, relative_energy_error, result_forces, max_allclose, abs_force_error, relative_force_error = CompareResults(base_dir, strfile=strfile, epsilon_r = 15, MBP_parameters=mixing_parameters)
    return (mixing_parameters, result_energy, abs_energy_error, relative_energy_error, result_forces, max_allclose, abs_force_error, relative_force_error)


class TestMultipleBasinGoMartiniTopology:
    """
    Test the multiple-basin GoMartini topology
    """
    path = os.path.dirname(__file__)

    def test_GlnBP_EXP(self):
        base_dir = os.path.join(self.path, "../fixtures/MBGoMartini/GlnBP/")
        EXP_parameters_list = [['True', 'exp', '2', '1/300', '-300', '0'],
                                ['True', 'exp', '2', '1/500', '-300', '0'],
                                ['True', 'exp', '2', '1/300', '-600', '0'],
                                ]
        # Uncomment to generate single-basin state energy:
        # for i in range(10):
        #     strfile = os.path.join(base_dir,f'Strfiles/GlnBP_No{i}.gro')
        #     CalSBPState(base_dir, strfile=strfile, epsilon_r=15)

        for j, EXP_parameters in enumerate(EXP_parameters_list):
            print(os.getcwd())
            print(EXP_parameters)
            
            # Set MBP parameters once before parallel execution to avoid race conditions
            working_dir = os.path.join(base_dir, 'EXP')
            SetMBPParameter(working_dir, topfile='system.top', MBP_parameters=EXP_parameters)
            
            Compare_per_frame = partial(Compare, base_dir=base_dir, mixing_parameters=EXP_parameters, prefix_strfile='GlnBP')
            with Pool(max(cpu_count(),10)) as pool:
                results = pool.map(Compare_per_frame, range(10))
                
            with open(os.path.join(base_dir, f'./EXP/EXP_results_{j}.pickle'), 'wb') as f:
                pickle.dump(results, f)

    def test_GlnBP_HAM(self):
        base_dir = os.path.join(self.path, "../fixtures/MBGoMartini/GlnBP/")
        HAM_parameters_list = [['True', 'ham', '2', '100', '-200', '0'],
                                ['True', 'ham', '2', '600', '-200', '0'],
                                ['True', 'ham', '2', '100', '-900', '0'],
                                ]
        # Uncomment to generate single-basin state energy:
        # for i in range(10):
        #     strfile = os.path.join(base_dir,f'Strfiles/GlnBP_No{i}.gro')
        #     CalSBPState(base_dir, strfile=strfile, epsilon_r=15)

        for j, HAM_parameters in enumerate(HAM_parameters_list):
            print(os.getcwd())
            print(HAM_parameters)
            
            # Set MBP parameters once before parallel execution to avoid race conditions
            working_dir = os.path.join(base_dir, 'HAM')
            SetMBPParameter(working_dir, topfile='system.top', MBP_parameters=HAM_parameters)
            
            Compare_per_frame = partial(Compare, base_dir=base_dir, mixing_parameters=HAM_parameters, prefix_strfile='GlnBP')
            with Pool(max(cpu_count(),10)) as pool:
                results = pool.map(Compare_per_frame, range(10))

            with open(os.path.join(base_dir, f'./HAM/HAM_results_{j}.pickle'), 'wb') as f:
                pickle.dump(results, f)

    def test_Beta2AR_EXP(self):
        base_dir = os.path.join(self.path, "../fixtures/MBGoMartini/Beta2AR/")
        EXP_parameters_list = [['True', 'exp', '2', '1/300', '-300', '0'],
                                ['True', 'exp', '2', '1/500', '-300', '0'],
                                ['True', 'exp', '2', '1/300', '-600', '0'],
                                ]
        # Uncomment to generate single-basin state energy:
        # for i in range(10):
        #     strfile = os.path.join(base_dir,f'Strfiles/Beta2AR_No{i}.gro')
        #     CalSBPState(base_dir, strfile=strfile, epsilon_r=15)

        for j, EXP_parameters in enumerate(EXP_parameters_list):
            print(os.getcwd())
            print(EXP_parameters)
            
            # Set MBP parameters once before parallel execution to avoid race conditions
            working_dir = os.path.join(base_dir, 'EXP')
            SetMBPParameter(working_dir, topfile='system.top', MBP_parameters=EXP_parameters)
            
            Compare_per_frame = partial(Compare, base_dir=base_dir, mixing_parameters=EXP_parameters, prefix_strfile='Beta2AR')
            with Pool(max(cpu_count(),10)) as pool:
                results = pool.map(Compare_per_frame, range(10))

            with open(os.path.join(base_dir, f'./EXP/EXP_results_{j}.pickle'), 'wb') as f:
                pickle.dump(results, f)
            


    def test_Beta2AR_HAM(self):
        base_dir = os.path.join(self.path, "../fixtures/MBGoMartini/Beta2AR/")
        HAM_parameters_list = [['True', 'ham', '2', '100', '-200', '0'],
                                ['True', 'ham', '2', '600', '-200', '0'],
                                ['True', 'ham', '2', '100', '-900', '0'],
                                ]
        # Uncomment to generate single-basin state energy:
        # for i in range(10):
        #     strfile = os.path.join(base_dir,f'Strfiles/Beta2AR_No{i}.gro')
        #     CalSBPState(base_dir, strfile=strfile, epsilon_r=15)

        for j, HAM_parameters in enumerate(HAM_parameters_list):
            print(os.getcwd())
            print(HAM_parameters)
            
            # Set MBP parameters once before parallel execution to avoid race conditions
            working_dir = os.path.join(base_dir, 'HAM')
            SetMBPParameter(working_dir, topfile='system.top', MBP_parameters=HAM_parameters)
            
            Compare_per_frame = partial(Compare, base_dir=base_dir, mixing_parameters=HAM_parameters, prefix_strfile='Beta2AR')
            with Pool(max(cpu_count(),10)) as pool:
                results = pool.map(Compare_per_frame, range(10))

            with open(os.path.join(base_dir, f'./HAM/HAM_results_{j}.pickle'), 'wb') as f:
                pickle.dump(results, f)