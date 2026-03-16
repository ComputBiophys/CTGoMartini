import sys
sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')
from ctgomartini.data.ctgomartinize import SBGOMartinize

SBGOMartinize(
    aa_strfile_list=['HAMP_Bicyclomycin.pdb'],
    aa_map_list=['HAMP_Bicyclomycin.map'],
    state_name_list=['Protein'],
    sbmol_name='Protein',
    method='SBP',
    dssp=None,
    ff='martini3001',
    other_params=''
)
print("SB topology generated!")
