# update: 
# 20230507: add the argparse for input
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import distance_array

import argparse

parser = argparse.ArgumentParser(description=
'''
An example:
python Contact.py -f strfile.pdb -r vdw_enlargement_factor -o new_OV.map

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-f', '--strfile', required=True,
                  help='The structure file (pdb)')
parser.add_argument('-r', '--vdw_enlargement_factor', default=1.5,
                  help='The vdw_enlargement_factor')
parser.add_argument('-o', '--fileout', default="new_OV.map",
                  help='The OV contact filename')

args = parser.parse_args()


#####
# Input and params
strfile= args.strfile
outfile=args.fileout

excluded_residues    = 0 # >0
# vdw_enlargement_factor=(26/7)**(1/6)
vdw_enlargement_factor=float(args.vdw_enlargement_factor)
######


vdw_radius = {'SER': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'OG': 1.46},
              'PRO': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD': 1.88},
              'TYR': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.61,
                      'CD1': 1.76,
                      'CD2': 1.76,
                      'CE1': 1.76,
                      'CE2': 1.76,
                      'CZ': 1.61,
                      'OH': 1.46},
              'VAL': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG1': 1.88,
                      'CG2': 1.88},
              'TRP': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.61,
                      'CD1': 1.76,
                      'CD2': 1.61,
                      'CE2': 1.61,
                      'CE3': 1.76,
                      'NE1': 1.64,
                      'CZ2': 1.76,
                      'CZ3': 1.76,
                      'CH2': 1.76},
              'GLN': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD': 1.61,
                      'NE2': 1.64,
                      'OE1': 1.42},
              'HIS': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.61,
                      'CD2': 1.76,
                      'ND1': 1.64,
                      'CE1': 1.76,
                      'NE2': 1.64},
              'GLU': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD': 1.61,
                      'OE1': 1.46,
                      'OE2': 1.42},
              'LYS': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD': 1.88,
                      'CE': 1.88,
                      'NZ': 1.64},
              'THR': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG2': 1.88,
                      'OG1': 1.46},
              'MET': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'SD': 1.77,
                      'CE': 1.88},
              'ALA': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88},
              'CYS': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'SG': 1.77},
              'ASP': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.61,
                      'OD1': 1.46,
                      'OD2': 1.42},
              'ASN': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.61,
                      'ND2': 1.64,
                      'OD1': 1.42},
              'ARG': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD': 1.88,
                      'NE': 1.64,
                      'CZ': 1.61,
                      'NH1': 1.64,
                      'NH2': 1.64},
              'PHE': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD1': 1.61,
                      'CD2': 1.76,
                      'CE1': 1.76,
                      'CE2': 1.76,
                      'CZ': 1.76},
              'ILE': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG1': 1.88,
                      'CG2': 1.88,
                      'CD1': 1.88},
              'GLY': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42},
              'LEU': {'N': 1.64,
                      'CA': 1.88,
                      'C': 1.61,
                      'O': 1.42,
                      'CB': 1.88,
                      'CG': 1.88,
                      'CD1': 1.88,
                      'CD2': 1.88}}

def Projection(resname, atom_names):
    vdw_list=[]
    for atom in atom_names:
        if atom in vdw_radius[resname]:
            vdw=vdw_radius[resname][atom]
            vdw_list.append(vdw)
        else:
            raise Exception(f"atom {atom} not in vdw_radius")
    return np.array(vdw_list)
def VDW_array(array1,array2):
    shape=(len(array1),len(array2))
    vdw_array=np.zeros(shape=shape)
    for i in range(len(array1)):
        for j in range(len(array2)):
            vdw_array[i,j]=array1[i]+array2[j]
    return vdw_array
def RemoveH(PDBFile):
    u=mda.Universe(PDBFile)
    sel=u.select_atoms('protein and not name H* OXT')
    sel.atoms.write("_protein_noH.pdb")
    



RemoveH(strfile)
strfile='_protein_noH.pdb'
u=mda.Universe(strfile)                                                                                    
# Contact pairs
n_residues=u.residues.n_residues
pair_list=[]
for i in range(n_residues):
    for j in range(i+1, n_residues):
        # exclude residues
        if (u.residues[i].segid == u.residues[j].segid) and (u.residues[j].resid-u.residues[i].resid<=excluded_residues):
            continue

        atom_names1=u.residues[i].atoms.names
        atom_names2=u.residues[j].atoms.names
        resname1=u.residues[i].resname
        resname2=u.residues[j].resname
        vdw_array=VDW_array(Projection(resname1,atom_names1),Projection(resname2,atom_names2))

        position1=u.residues[i].atoms.positions
        position2=u.residues[j].atoms.positions        

        contacts=distance_array(position1,position2)
        judge=(contacts <= vdw_array * vdw_enlargement_factor) & (contacts > 0)
        if judge.sum()>0:
            pair=[i,j]
            pair_list.append(pair)

# Double replicate
new_pair_list=[]
for pair in pair_list:
    if pair not in new_pair_list:
        new_pair_list.append(pair)
    if pair[::-1] not in new_pair_list:
        new_pair_list.append(pair[::-1])
pair_list=new_pair_list







def ExtactContacts(contact_pair, u, sel='CA'):
    resindex1=contact_pair[0]
    resindex2=contact_pair[1]

    position1=u.residues[resindex1].atoms[u.residues[resindex1].atoms.names==sel].positions
    position2=u.residues[resindex2].atoms[u.residues[resindex2].atoms.names==sel].positions

    dist=distance_array(position1,position2)[0][0]
    return dist

lines=[]
lines.append('            I1  AA  C I(PDB)    I2  AA  C I(PDB)    DISTANCE       CMs    rCSU    aSurf    rSurf    nSurf\n')
lines.append('==========================================================================================================\n')
for i, pair in enumerate(pair_list):
    residue1=u.residues[pair[0]]
    residue2=u.residues[pair[1]]
    
    I1=pair[0]+1
    I2=pair[1]+1

    AA1=residue1.resname
    AA2=residue2.resname

    C1=residue1.segid
    C2=residue2.segid

    I1_pdb=residue1.resid
    I2_pdb=residue2.resid
    
    dist=ExtactContacts(pair, u, sel='CA')

    line=f'R {i+1:>6} {I1:>5} {AA1:>4} {C1} {I1_pdb:>4} {I2:>8} {AA2:>4} {C2} {I2_pdb:>4} {dist:>12.4f}     1 0 0 0     0   0.0000   0.0000   0.0000\n'
    lines.append(line)
    # print(line)




# Write Contact
with open(outfile,'w') as fp:
    fp.writelines(lines)
