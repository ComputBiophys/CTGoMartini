# Update: 20250112
# Update: For multiple-states
# Author: Song Yang
# Version: 3.0

import os
import numpy as np
import pandas as pd
import argparse


import MDAnalysis as mda
from multiprocessing import Pool
from functools import partial

import time, datetime


def ReportTime(start_time):
    end_time=datetime.datetime.now()
    elapsed=(end_time-start_time).total_seconds()
    
    start_time=start_time.__format__("%Y-%m-%d %H:%M:%S")
    end_time=end_time.__format__("%Y-%m-%d %H:%M:%S")
    print(f"Elapsed Time: {elapsed:.2f}")
    
start_time=datetime.datetime.now() # start time
# ReportTime(start_time) # end time


def Distance_pair(a1, a2):
    """
    a1,a2: two AtomGroups (single atom for each AtomGroup)
    """
    p1=a1.position
    p2=a2.position
    return np.linalg.norm(p1-p2)

def Distance(p1, p2):
    """
    p1,p2: np.linalg.norm(p1-p2, axis=1)
    """
    return np.linalg.norm(p1-p2, axis=1)

def SameSegidResid(atom1, atom2):
    return (atom1.segid == atom2.segid) and (atom1.resid == atom2.resid)


def Diff_distance_list(distance_list):
    diff_distance_list = []
    for i in range(len(distance_list)):
        for j in range(i+1, len(distance_list)):
            diff_distance_list.append(abs(distance_list[i]-distance_list[j]))
    return diff_distance_list


def Distances_ref(atoms_ref_list, params):
    minimum_distance     = params['minimum_distance']
    maximum_distance     = params['maximum_distance']
    minimum_difference   = params['minimum_difference']
    excluded_residues    = params['excluded_residues']

    # Length test
    assert len(atoms_ref_list) >=2, "Error: Reference structures should have at least two states!"
    for i, atoms_ref in enumerate(atoms_ref_list[1:]):
        if len(atoms_ref) != len(atoms_ref_list[0]):
            raise ValueError(f"Reference structures (1 and {i+2}) have no same number of atoms!")  
    length = len(atoms_ref_list[0])

    # Distance pairs
    distances_ref_list = [[] for _ in range(len(atoms_ref_list))]
    for i in range(length):
        for j in range(i+1, length):
            for atoms_ref in atoms_ref_list[1:]:
                assert SameSegidResid(atoms_ref_list[0][i], atoms_ref[i]), f"Error: Reference structures have no same atom ({i}) segid and resid!"
                assert SameSegidResid(atoms_ref_list[0][j], atoms_ref[j]), f"Error: Reference structures have no same atom ({j}) segid and resid!"


            atom_i_list = [atoms_ref[i] for atoms_ref in atoms_ref_list]
            atom_j_list = [atoms_ref[j] for atoms_ref in atoms_ref_list]


            # Excluded_residues: >=4 (default)
            condition=False
            if atom_i_list[0].segid != atom_j_list[0].segid:
                condition=True
            else:
                if atom_i_list[0].resid + excluded_residues <= atom_j_list[0].resid:
                    condition=True

            # Min/Max/Diff distance: 6, 50, 5 (default)
            if condition:
                distance_list = [Distance_pair(atom_i, atom_j) for atom_i, atom_j in zip(atom_i_list, atom_j_list)]
                diff_distance_min = np.min(Diff_distance_list(distance_list))

                if diff_distance_min >= minimum_difference:
                    for _, distance in enumerate(distance_list):
                        if distance >= minimum_distance and distance <= maximum_distance:
                            distances_ref_list[_].append([i,j,distance])                 

    distances_cal_list = []
    for distances_ref in distances_ref_list:
        distances_cal = pd.DataFrame(distances_ref, columns=['atomid1','atomid2','ref_d'])
        distances_cal_list.append(distances_cal)

    return distances_cal_list

def dRMS(distances, ref_distances):
    return np.sqrt(np.sum((distances-ref_distances)**2)/len(ref_distances))

def dRMS_per_frame(frame, atoms, distances_ref):
    # Loading frame
    ts=atoms.universe.trajectory[frame]
    
    p1=atoms[distances_ref.atomid1].positions
    p2=atoms[distances_ref.atomid2].positions    
    result=dRMS(Distance(p1,p2), distances_ref.ref_d)
    return [ts.time/1000, result] # unit:ns, A   

def dRMS_trajectory(atoms, distances_ref, n_cores=10, skip=1):
    dRMS_per_frame_partial=partial(dRMS_per_frame, atoms=atoms, distances_ref=distances_ref)
    length=atoms.universe.trajectory.n_frames
    with Pool(n_cores) as pool:
        results=pool.map(dRMS_per_frame_partial, list(range(0, length, skip)))
    return np.array(results)   


class MBAnalysis:
    def __init__(self):
        self.dRMStrj_A = None
        self.dRMStrj_B = None
        self.dRMStrj_C = None
        self.dRMStrj_list = None
    
    def dRMStrj(self, str_file, trj_file, str_ref_list, **kwargs):
        """
        Default parameters:
        minimum_distance     = 6.0
        maximum_distance     = 50.0
        minimum_difference   = 5
        excluded_residues    = 4
        selected_atom        = 'name BB'
        n_cores              = 10
        """
        # Parametrs
        default_params={
        'minimum_distance':6.0,
        'maximum_distance':50.0,
        'minimum_difference':5,
        'excluded_residues':4,
        'selected_atom':'name BB',
        'n_cores': 10,
        
        }
        default_params.update(kwargs)
        params = default_params
        
        minimum_distance     = params['minimum_distance']
        maximum_distance     = params['maximum_distance']
        minimum_difference   = params['minimum_difference']
        excluded_residues    = params['excluded_residues']
        selected_atom        = params['selected_atom']
        skip                 = params['skip']

        # Load reference structures
        u_ref_list = [mda.Universe(str_ref) for str_ref in str_ref_list]
        atoms_ref_list = [u_ref.select_atoms(selected_atom) for u_ref in u_ref_list]    
        
        # Cal reference distances
        distances_ref_list = Distances_ref(atoms_ref_list, params)
        
        # Trajectory analysis
        u=mda.Universe(str_file, trj_file)
        atoms=u.select_atoms(selected_atom)

        start_time=datetime.datetime.now() # start time
        n_cores=params['n_cores']
        dRMStrj_list = [dRMS_trajectory(atoms, distances_ref, n_cores=n_cores, skip=skip) for distances_ref in distances_ref_list]
        print("dRMStrj Calculation Finished!")
        ReportTime(start_time) # end time
        
        # Results
        self.dRMStrj_list = dRMStrj_list
        self.dRMStrj_A = dRMStrj_list[0]
        self.dRMStrj_B = dRMStrj_list[1]
        if len(dRMStrj_list) == 3:
            self.dRMStrj_C = dRMStrj_list[2]
        
    def dRMStrj_DataSave(self, prefix = 'dRMStrj'):
        '''
        Save dRMStrj data to .dat file
        '''
        # Save Data
        id = 'ABCDEFGHIJKLMN'
        for i, dRMStrj in enumerate(self.dRMStrj_list):
            dataName = prefix + '_state' + id[i] + '.dat'
            np.savetxt(dataName,dRMStrj)       
    

    
def input_args():
    parser=argparse.ArgumentParser(description=""" """)
    
    parser.add_argument('-s', '--str_file', required=True,
                      help='The structure file (gro pdb)')
    parser.add_argument('-f', '--trj_file', nargs='+', default=None, 
                      help='The trajectory file (xtc)')
    parser.add_argument('--skip', default=1, type=int,
                      help='Number of frames for skipping. Default is 1.')
    parser.add_argument('-r', '--ref_str', required=True, nargs='+',
                      help='The reference structure file (gro pdb)')
    parser.add_argument('-sel', '--selected_atom', default='name BB',
                        help='The selected atom for calculating the dRMS')
    parser.add_argument('-prefix', '--output_prefix', default='dRMStrj',
                        help='The prefix of output data file')


    args = parser.parse_args()
    #args = parser.parse_args('-s npt.pdb -f md_rottrans.dcd -r GBP_open_cg.pdb GBP_close_cg.pdb'.split())
    return args

    
if __name__=='__main__':
    # Data path

    # Structure and trjactory    
    args = input_args()
    str_file = args.str_file
    trj_file = args.trj_file
    skip = int(args.skip)
    str_ref_list = args.ref_str
    
    selected_atom=args.selected_atom
    output_prefix=args.output_prefix
    
    dataNameA=output_prefix+'_stateA.dat'
    dataNameB=output_prefix+'_stateB.dat'
    # Analysis
    mbana=MBAnalysis()
    mbana.dRMStrj(str_file, trj_file, str_ref_list, selected_atom=selected_atom, skip=skip)
    # Save Data
    mbana.dRMStrj_DataSave(prefix=output_prefix)

