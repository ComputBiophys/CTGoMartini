# Update: 20240112
# Update: extend the input options
# Author: Song Yang
# Version: 1.0
# Beijixing-speicla Argparse
# Add skip

import os
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
#plt.style.use('normal')

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
ReportTime(start_time) # end time



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

def SegidResidTest(atom1, atom2):
    result=(atom1.segid == atom2.segid) and (atom1.resid == atom2.resid)
    return result

def Test(atom_A1, atom_A2, atom_B1, atom_B2):
    if SegidResidTest(atom_A1, atom_B1) and SegidResidTest(atom_A2, atom_B2):
        pass
    else:
        raise Exception("Atom resids or segids are not the same!") 


def UpdateParams(kwargs):
    params={
        'minimum_distance':6.0,
        'maximum_distance':50.0,
        'minimum_difference':5,
        'excluded_residues':4,
        'selected_atom':'name BB',
        'n_cores': 10,
        
    }
    params.update(kwargs)
    return params
    
def Distances_ref(atoms_A, atoms_B, params):
    minimum_distance     = params['minimum_distance']
    maximum_distance     = params['maximum_distance']
    minimum_difference   = params['minimum_difference']
    excluded_residues    = params['excluded_residues']

    # Length test
    if len(atoms_A) == len(atoms_B):
        length=len(atoms_A)
    else:
        raise ValueError("Two-state atoms selected have no same number!")       
    
    # Distance pairs
    distances_cal_A=[]
    distances_cal_B=[]
    for i in range(length):
        for j in range(i+1, length):
            atom_A1, atom_A2=atoms_A[i], atoms_A[j]
            atom_B1, atom_B2=atoms_B[i], atoms_B[j]

            # Resid and Segid Test
            Test(atom_A1, atom_A2, atom_B1, atom_B2)

            # Excluded_residues: >=4 (default)
            condition=False
            if atom_A1.segid !=atom_A2.segid:
                condition=True
            else:
                if atom_A1.resid + excluded_residues <= atom_A2.resid:
                    condition=True

            # Min/Max/Diff distance: 6, 50, 5 (default)
            if condition:
                distance_A=Distance_pair(atom_A1, atom_A2)
                distance_B=Distance_pair(atom_B1, atom_B2)

                if abs(distance_A-distance_B) >= minimum_difference:
                    if distance_A >= minimum_distance and distance_A <= maximum_distance:
                        distances_cal_A.append([i,j,distance_A])
                    if distance_B >= minimum_distance and distance_B <= maximum_distance:
                        distances_cal_B.append([i,j,distance_B])                    

    distances_cal_A=pd.DataFrame(distances_cal_A, columns=['atomid1','atomid2','ref_d'])
    distances_cal_B=pd.DataFrame(distances_cal_B, columns=['atomid1','atomid2','ref_d'])
    
    return distances_cal_A, distances_cal_B

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
        self.dRMStrj_A=None
        self.dRMStrj_B=None
        pass
    
    def dRMStrj(self, str_file, trj_file, str_ref_stateA, str_ref_stateB, **kwargs):
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
        params=UpdateParams(kwargs)
        
        minimum_distance     = params['minimum_distance']
        maximum_distance     = params['maximum_distance']
        minimum_difference   = params['minimum_difference']
        excluded_residues    = params['excluded_residues']
        selected_atom        = params['selected_atom']
        skip                 = params['skip']

        # Load reference structures
        u_A=mda.Universe(str_ref_stateA)
        u_B=mda.Universe(str_ref_stateB)

        atoms_A=u_A.select_atoms(f'{selected_atom}')
        atoms_B=u_B.select_atoms(f'{selected_atom}')
        
        # Cal reference distances
        distances_ref_A, distances_ref_B=Distances_ref(atoms_A, atoms_B, params)
        
        # Trajectory analysis
        u=mda.Universe(str_file, trj_file)
        atoms=u.select_atoms(selected_atom)

        start_time=datetime.datetime.now() # start time
        n_cores=params['n_cores']
        dRMStrj_A=dRMS_trajectory(atoms, distances_ref_A, n_cores=n_cores, skip=skip)
        dRMStrj_B=dRMS_trajectory(atoms, distances_ref_B, n_cores=n_cores, skip=skip)
        print("dRMStrj Calculation Finished!")
        ReportTime(start_time) # end time
        
        # Results
        self.dRMStrj_A=dRMStrj_A
        self.dRMStrj_B=dRMStrj_B
        
    def dRMStrj_DataSave(self, dataNameA='dRMStrj_stateA.dat', dataNameB='dRMStrj_stateB.dat' ):
        '''
        dataNameA='dRMStrj_stateA.dat'
        dataNameB='dRMStrj_stateB.dat'        
        '''
        # Save Data
        np.savetxt(dataNameA,self.dRMStrj_A)
        np.savetxt(dataNameB,self.dRMStrj_B)           
    
    def PMF_dRMS(self, dRMStrj_array, Temp=310.15, bins=200, min_dRMS=None, max_dRMS=None):
        min_dRMS=min(dRMStrj_array) if min_dRMS is None else min_dRMS
        max_dRMS=max(dRMStrj_array) if max_dRMS is None else max_dRMS

        # PMF Cal
        results=np.histogram(dRMStrj_array, bins=bins, range=(min_dRMS, max_dRMS))
        dRMS=(results[1][:-1]+results[1][1:])/2
        Y=results[0]
        Density=Y/sum(Y)
        PMF=-8.314/1000*Temp*np.log(Density)
        PMF=PMF-PMF.min()
        
        return dRMS, PMF
    
    def Keq_Cal(self, dRMStrj_array, lp=None, rp=None, bins=200):
        # PMF Cal
        results=np.histogram(dRMStrj_array, bins=bins)
        dRMS=(results[1][:-1]+results[1][1:])/2
        Y=results[0]
        Density=Y/sum(Y)
        
        lp=min(dRMS) if lp is None else lp
        rp=max(dRMS) if rp is None else rp        
        search= (dRMS> lp) & (dRMS<rp)
        mid=Density[search].argmin()+search.argmax()
        Keq=Density[:mid].sum()/Density[mid:].sum()
        print(f"dRMS: {dRMS[mid]:.4f}, Keq: {Keq:.4f}")
        return Keq

        
    

def dRMStrj_Plot(dRMStrj, **kwargs):
    """
    params={
        'xlim':[0,10],
        'ylim':[0,10],
        'color':np.array([255,81,81])/255,
        'alpha':1,
        'lw':1,
        'filename':'repeat',
        'DataSave':".",
        'isSave': False
    }
    """
    params={
        'xlim':[0,10],
        'ylim':[0,10],
        'color':np.array([255,81,81])/255,
        'alpha':1,
        'lw':1,
        'filename':'repeat',
        'DataSave':".",
        'isSave': False
    }
    params.update(kwargs)

    plt.figure(figsize=(10,5))
    plt.plot(dRMStrj[:,0], dRMStrj[:,1], lw=params['lw'], color=params['color'], alpha=params['alpha'])
    plt.xlabel('Time (μs)')
    plt.ylabel('dRMS (Å)')
    plt.xlim(params['xlim'])
    plt.ylim(params['ylim'])
    filename=params['filename']+'.svg'
    DataSave=params['DataSave']
    filesave=os.path.join(DataSave, filename)
    if params['isSave']:
        plt.savefig(filesave, bbox_inches='tight') 
    plt.show()
    
def PMF_Plot(dRMS, PMF, **kwargs):
    """
    params={
        'xlim':[0,12],
        'ylim':[0,15],
        'color':'black',
        'alpha':1,
        'lw':2,
        'filename':'repeat',
        'DataSave':".",
        'isSave': False
    }
    """
    params={
        'xlim':[0,12],
        'ylim':[0,15],
        'color':'black',
        'alpha':1,
        'lw':2,
        'filename':'repeat',
        'DataSave':".",
        'isSave': False
    }
    params.update(kwargs)

    plt.figure(figsize=(10,5))
    plt.plot(dRMS, PMF, lw=params['lw'], color=params['color'], alpha=params['alpha'])
    plt.xlabel('dRMS')
    plt.ylabel('PMF (kJ/mol)')
    plt.xlim(params['xlim'])
    plt.ylim(params['ylim'])
    filename=params['filename']+'.svg'
    DataSave=params['DataSave']
    filesave=os.path.join(DataSave, filename)
    if params['isSave']:
        plt.savefig(filesave, bbox_inches='tight') 
    plt.show()

    
def input_args():
    parser=argparse.ArgumentParser(description=""" """)
    
    parser.add_argument('-s', '--str_file', required=True,
                      help='The structure file (gro pdb)')
    parser.add_argument('-f', '--trj_file', nargs='+', default=None, 
                      help='The trajectory file (xtc)')
    parser.add_argument('--skip', default=1, type=int,
                      help='Number of frames for skipping. Default is 1.')
    parser.add_argument('-rA', '--ref_str_A', required=True,
                      help='The reference structure file of state A (gro pdb)')
    parser.add_argument('-rB', '--ref_str_B', required=True,
                      help='The reference structure file of state B (gro pdb)')
    parser.add_argument('-sel', '--selected_atom', default='name BB',
                        help='The selected atom for calculating the dRMS')
    parser.add_argument('-prefix', '--output_prefix', default='dRMStrj',
                        help='The prefix of output data file')


    args = parser.parse_args()
    #args = parser.parse_args('-s npt.pdb -f md_rottrans.dcd -rA GBP_open_cg.pdb -rB GBP_close_cg.pdb'.split())
    return args

    
if __name__=='__main__':
    # Data path
    #DataPath="/home/ys/SongYang/MultipleBasin/Work/GlnBP/repeat-1000_-200-0"
    #os.chdir(DataPath)

    # Structure and trjactory    
    args=input_args()
    str_file=args.str_file
    trj_file=args.trj_file
    skip = int(args.skip)
    str_ref_stateA=args.ref_str_A
    str_ref_stateB=args.ref_str_B
    
    selected_atom=args.selected_atom
    output_prefix=args.output_prefix
    
    dataNameA=output_prefix+'_stateA.dat'
    dataNameB=output_prefix+'_stateB.dat'
    # Analysis
    mbana=MBAnalysis()
    mbana.dRMStrj(str_file, trj_file, str_ref_stateA, str_ref_stateB, selected_atom=selected_atom, skip=skip)
    # Save Data
    mbana.dRMStrj_DataSave(dataNameA=dataNameA, dataNameB=dataNameB)

