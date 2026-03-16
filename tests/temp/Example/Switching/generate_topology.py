import sys
sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')
from ctgomartini.data.ctgomartinize import SwitchingGOMartinize

SwitchingGOMartinize(
    aa_strfile_list=['State1.pdb', 'State2.pdb'],
    aa_map_list=['State1.map', 'State2.map'],
    state_name_list=['Open', 'Closed'],
    mbmol_name='Protein',
    dict_cutoffs={
        'cutoff_BBB_angles': 15.0,
        'cutoff_BBBB_dihedrals': 30.0,
        'cutoff_SBBS_dihedrals': 30.0,
        'cutoff_contacts': 0.06
    },
    method='switching',
    dssp=None,
    ff='martini3001',
    other_params=''
)
print("Switching topology generated!")
