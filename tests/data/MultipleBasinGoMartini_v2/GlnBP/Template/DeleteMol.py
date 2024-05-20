# Update: 20230706 Song Yang # revise bugs for nolipid systems
# Update: 20230709 Song Yang # revise bugs for nodiff systems
import MDAnalysis as mda
from optparse import OptionParser

def GetMolInfo(strfile):
    # parameters
    lipid_list=['POPC', 'POPS', 'POPG','POPE','POP2','PIP2']
    solute_list=['NA','NA+','CL', 'CL-', 'W']

    # analysis
    u=mda.Universe(strfile)
    mol_dict={} # mol_name: resid_list
    for solute in solute_list:
        resid_list=u.select_atoms(f'resname {solute}').residues.resids.tolist()
        if len(resid_list)==0: continue
        mol_dict[solute]=resid_list

    for lipid in lipid_list:
        resid_list=u.select_atoms(f'resname {lipid}').residues.resids
        if len(resid_list)==0: continue
        center=u.select_atoms('resname '+' '.join(lipid_list)).center_of_geometry()[2]

        ulipid_list=[]
        llipid_list=[]

        for resid in resid_list:
            if u.select_atoms(f'resname {lipid} and resid {resid} and name PO4').positions[:,2] > center:
                ulipid_list.append(resid)
            else:
                llipid_list.append(resid)

        mol_dict['u'+lipid]=ulipid_list    
        mol_dict['l'+lipid]=llipid_list

    return mol_dict

usage = "python DeleteMol.py input.pdb ref.pdb"
parser = OptionParser(usage)
(options, args) = parser.parse_args()

strfile=args[0]
reffile=args[1]

u=mda.Universe(strfile)
mol_dict=GetMolInfo(strfile)
ref_dict=GetMolInfo(reffile)
assert mol_dict.keys() == ref_dict.keys(), "Mol and reference dicts must have the same keys"

# Difference
diff_dict={}
print(f'Difference between {strfile} and {reffile}')
for key in mol_dict.keys():
    size_mol=len(mol_dict[key])
    size_ref=len(ref_dict[key])
    diff=size_mol-size_ref
    diff_dict[key]=diff
    print(f"{key:>6}: {diff:>5}")
print()

for key, value in diff_dict.items():
    if value < 0:
        raise Exception('Error: do not support <0', key, value)
del_dict={}
for key, value in diff_dict.items():
    if value > 0:
        del_dict[key]=mol_dict[key][:value]
sel_string_list=[]
print("Delete mols")
for key, value in del_dict.items():
    if key[0] in['u', 'l']:
        key=key[1:]
    value_string=' '.join(map(str, value))
    sel_string=f'(resname {key} and resid {value_string})'
    sel_string_list.append(sel_string)
    print(sel_string[1:-1])
print()

sel_string=f'all and not ({" or ".join(sel_string_list)})' if len(sel_string_list) != 0 else 'all'
u.select_atoms(sel_string).write('deleting_ionized.gro')

strfile='deleting_ionized.gro'

mol_dict=GetMolInfo(strfile)
ref_dict=GetMolInfo(reffile)
assert mol_dict.keys() == ref_dict.keys(), "Mol and reference dicts must have the same keys"

# Difference
diff_dict={}
print(f'Difference between {strfile} and {reffile}')
for key in mol_dict.keys():
    size_mol=len(mol_dict[key])
    size_ref=len(ref_dict[key])
    diff=size_mol-size_ref
    diff_dict[key]=diff
    print(f"{key:>6}: {diff:>5}")

for key, value in diff_dict.items():
    if value < 0:
        raise Exception('Error: do not support <0', key, value)
