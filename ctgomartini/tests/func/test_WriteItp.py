import os
from functools import partial

from ctgomartini.api import GenMBPTop, MartiniTopFile
from ctgomartini.util import SameListList
from ctgomartini.func import WriteItp

def Angles_Dihedrals_Sort(fields, category):
    if category in ['angles', 'multi_angles']:
        if int(fields[0]) >  int(fields[2]):
            fields = fields[2::-1] + fields[3:]
    elif category in ['dihedrals', 'multi_dihedrals']:
        if int(fields[0]) > int(fields[3]):
            fields = fields[3::-1] + fields[4:]
    return fields

def Comparison_Top(mbmol_ref, mbmol_test):
    categories_list = list(set(list(mbmol_ref._topology.keys()) + list(mbmol_test._topology.keys())))
    for category in categories_list:
        if category in ['angles', 'dihedrals', 'multi_angles', 'multi_dihedrals']:
            Angles_Dihedrals_Sort_partial = partial(Angles_Dihedrals_Sort, category=category)
            mbmol_ref._topology[category] = list(map(Angles_Dihedrals_Sort_partial, mbmol_ref._topology[category]))
            mbmol_test._topology[category] = list(map(Angles_Dihedrals_Sort_partial, mbmol_test._topology[category]))
        same = SameListList([mbmol_ref._topology[category], mbmol_test._topology[category]], sort=True)
        assert same is True, f"Error: comparison of {category} between test and ref is not the same!"

def Comparison_ITP(working_dir, molname, topfile):
    os.chdir(os.path.join(working_dir, "test"))
    mbmol_test = MartiniTopFile(topfile)._moleculeTypes[molname]
    
    os.chdir(os.path.join(working_dir, "ref"))
    mbmol_ref = MartiniTopFile(topfile)._moleculeTypes[molname]

    Comparison_Top(mbmol_ref, mbmol_test)


class TestMBMartiniTIP:
    """
    Test the ITP file of multiple baisn GoMartini
    """
    path = os.path.dirname(__file__)

    def test_GlnBP_ITP(self):
        working_dir = os.path.join(self.path, "../data/WriteItp/WriteItp")

        os.chdir(os.path.join(working_dir, "test"))
        topfileA = 'system_open.top'
        mol_nameA = 'gbp_open'
        topfileB = 'system_closed.top'
        mol_nameB = 'gbp_closed'
        mbmol_name = 'gbp'
        mols_list = [
            [topfileA, mol_nameA],
            [topfileB, mol_nameB]
        ]
        dict_cutoffs = { 
            'cutoff_BBB_angles': 15,
            'cutoff_BBBB_dihedrals': 30,
            'cutoff_SBBS_dihedrals': 30,
            'cutoff_contacts': 0.06 }
        mbmol =  GenMBPTop(mols_list, mbmol_name, dict_cutoffs)
        WriteItp(mbmol)
        
        Comparison_ITP(working_dir, 'gbp', 'system.top')


    

    

    