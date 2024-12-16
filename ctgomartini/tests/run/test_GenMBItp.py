import os
from functools import partial

import ctgomartini
from ctgomartini.api import GenMBPTop, MartiniTopFile
from ctgomartini.util import SameListList
from ctgomartini.func import WriteItp, ConvertLongShortElasticBonds
import subprocess

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
        # if category == 'atoms': continue
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

    def test_GenMBItp_EXP(self):
        working_dir = os.path.join(self.path, "../data/WriteItp/EXP")
        os.chdir(working_dir)

        os.system('rm -r test')
        os.system('cp -a template test')
        
        os.chdir(os.path.join(working_dir, 'test'))
        # Fetch ctgomarinize
        os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/ctgomartinize.py')} .")
        # Generate Itp
        subprocess.run("python ctgomartinize.py -s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp dssp -ff martini3001 -method exp",
                       shell=True)
        
        # Check Comparison
        Comparison_ITP(working_dir, 'gbp', 'system.top')

    def test_GenMBItp_HAM(self):
        working_dir = os.path.join(self.path, "../data/WriteItp/HAM")
        os.chdir(working_dir)

        os.system('rm -r test')
        os.system('cp -a template test')
        
        os.chdir(os.path.join(working_dir, 'test'))
        # Fetch ctgomarinize
        os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/ctgomartinize.py')} .")
        # Generate Itp
        subprocess.run("python ctgomartinize.py -s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp dssp -ff martini3001 -method ham",
                       shell=True)
        
        # Check Comparison
        Comparison_ITP(working_dir, 'gbp', 'system.top')   

