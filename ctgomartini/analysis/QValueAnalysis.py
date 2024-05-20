# Update: 20240109
# Q value for native contacts

import os
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

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



def GetNativeContacts(topfile, mol_name, state_id):
    from ctgomartini.api import MartiniTopFile
    """ Get native contacts and mol.atoms."""

    top = MartiniTopFile(topfile)
    mol = top.moleculeTypes[mol_name]

    native_contacts = [] # [[atom_index1, atom_index2], ...]
    for multi_contact in mol.multi_contacts:
        if int(multi_contact[3]) == state_id:
            native_contacts.append([int(multi_contact[0])-1, int(multi_contact[1])-1])
    return native_contacts, mol.atoms

def Check_proteins(mda_protein_atoms, mol_atoms):
    """Check the atom number, atom name, and atom resname for mda_sel.atoms and mol_aotms"""
    protein = mda_protein_atoms
    assert len(protein.atoms) == len(mol_atoms), f"Error: the length of proteins is not the same.\nmda: {len(protein.atoms)}, mol: {len(mol_atoms)}"
    for i in range(len(mol_atoms)):
        atom = protein.atoms[i]
        mol_atom = mol_atoms[i]
        assert atom.name == mol_atom[4], f"Error: atom names are different. {i}, {atom.name} {mol_atom[4]}"
        assert atom.residue.resname == mol_atom[3], f"Error: atom resnames are different. {i}, {atom.residue.resname} {mol_atom[3]}"


def Distance_cal_pre_frame(native_contacts, atom_positions):
    """Calculate the distance array of native contacts based on atom postions"""
    results = []
    for contact in native_contacts:
        p1 = atom_positions[contact[0]]
        p2 = atom_positions[contact[1]]
        distance = np.linalg.norm(p1-p2)

        results.append(distance)

    return np.array(results)

def Q_Cal(r_0, r, beta=5, lambda_constant=1.2):
    """Calculate the Q value for r_0 and r"""
    Q = np.mean(1/(1 + np.exp(beta * (r - lambda_constant * r_0))))
    return Q

def NativeContanceAnalysis(topfile, mol_name, state_id, 
                           ref_str, mda_sel_ref,
                           strfile, trjfile, mda_sel, skip=1,
                           beta=5, lambda_constant=1.2):
    """Calculate the Q values for one trajectory.

    Parameters
    ----------
    topfile: str,
        topology filename.
    mol_name: str,
        molecule name selected in the topology file.
    state_id: int,
        state id used in the molecule itp file corresponding to the topology file.
    ref_str: str,
        filename of the reference structure.
    mda_sel_ref: str,
        string for selecting the corresponding protein for the reference structure.
    strfile: str,
        filename of the structure analyzed.
    trjfile: str,
        filename of the trajectory analyzed.
    mda_sel: str,
        string for selecting the corresponding protein for the structure analyzed.
    skip: int, default=1,
        frames skipped for the trajectory.
    beta: float, default=10,
        parameter for calculating the Q value.
    lambda_constanct: float, default=1.2,
        paramter for calculating the Q value.

    Return
    ------
    time: numpy.array,
        array of the time
    Q_array: numpy.array,
        array of the Q values
    """
    # Obtain the native contacts and mol_atoms
    native_contacts, mol_atoms = GetNativeContacts(topfile, mol_name, state_id)

    # Obtain the u and u_ref
    u = mda.Universe(strfile, trjfile)
    protein = u.select_atoms(mda_sel)
    u_ref = mda.Universe(ref_str)
    protein_ref = u_ref.select_atoms(mda_sel_ref)

    # Check whether the mol and mda.selection are same 
    Check_proteins(protein_ref, mol_atoms)
    Check_proteins(protein, mol_atoms)

    # Calculate the distance array
    r_0 = Distance_cal_pre_frame(native_contacts, protein_ref.positions)

    r_list = []
    time_list = []
    for ts in u.trajectory[::skip]:
        r = Distance_cal_pre_frame(native_contacts, protein.positions)
        r_list.append(r)
        time_list.append(ts.time)
    time_array = np.array(time_list)

    # Calculate the Q value
    Q_list = []
    for i, r in enumerate(r_list):
        Q = Q_Cal(r_0, r, beta=10, lambda_constant=1.2)
        Q_list.append(Q)
    Q_array = np.array(Q_list)

    return time_array, Q_array
def input_args():
    parser=argparse.ArgumentParser(description="""Calculate the Q values for the trajectory""")

    parser.add_argument("--topfile", type=str, help="topology filename")
    parser.add_argument("--mol_name", type=str, help="molecule name selected in the topology file")
    parser.add_argument("--state_id", type=int, help="state id used in the molecule itp file corresponding to the topology file")
    parser.add_argument("--ref_str", type=str, help="filename of the reference structure")
    parser.add_argument("--mda_sel_ref", type=str, help="string for selecting the corresponding protein for the reference structure")
    parser.add_argument("--strfile", type=str, help="filename of the structure analyzed")
    parser.add_argument("--trjfile", type=str, help="filename of the trajectory analyzed")
    parser.add_argument("--mda_sel", type=str, help="string for selecting the corresponding protein for the structure analyzed")
    parser.add_argument("--skip", type=int, default=1, help="frames skipped for the trajectory. Default value is 1")
    parser.add_argument("--beta", type=float, default=5.0, help="parameter for calculating the Q value. Default value is 5.0")
    parser.add_argument("--lambda_constant", type=float, default=1.2, help="paramter for calculating the Q value. Default value is 1.2.")
    parser.add_argument("--outfile", type=str, default='QValue.npy', help="output filename for saving data. Default is QValue.npy")

    args = parser.parse_args()

    return args

    
if __name__=='__main__':
    args=input_args()
    
    time_array, Q_array = NativeContanceAnalysis(args.topfile, args.mol_name, args.state_id,
                                                args.ref_str, args.mda_sel_ref,
                                                args.strfile, args.trjfile, args.mda_sel, skip=args.skip,
                                                beta=args.beta, lambda_constant=args.lambda_constant)
    data = np.vstack([time_array, Q_array]).T
    np.save(args.outfile, data)
if __name__=='test':
    # Input parameters
    # Topology parameters
    top_path = "/home/ys/SongYang/MultipleBasin/Work/ContactOpt/Arc/MBContact8/repeat-sheet-copy"
    topfile = os.path.join(top_path, 'system.top')
    mol_name = 'Arc'
    state_id = 1

    # Reference parameters
    ref_str = os.path.join(top_path, 'Protein_sheet_cg.pdb')
    mda_sel_ref = 'protein'

    # Trajectory analyzed parameters
    data_path = "/home/ys/SongYang/MultipleBasin/Work/ContactOpt/Arc/MBContact8/Data/repeat-helix-long-repeat6-350_20_0"
    strfile = os.path.join(data_path, "npt_revised.pdb")
    trjfile = os.path.join(data_path, "md_rottrans_dt1ns.dcd")
    mda_sel = 'protein'

    # Q value parameters
    beta=5
    lambda_constant=1.2

    # Calculate 
    time_array, Q_array = NativeContanceAnalysis(topfile, mol_name, state_id, 
                                                ref_str, mda_sel_ref,
                                                strfile, trjfile, mda_sel, skip=1,
                                                beta=beta, lambda_constant=lambda_constant)
    ReportTime(start_time) # end time
    plt.plot(time_array, Q_array)