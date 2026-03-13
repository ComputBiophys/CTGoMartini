"""
End-to-end tests for multiple-basin ITP generation.

Tests the complete ctgomartinize.py workflow which:
    1. Runs Martinize2 to generate single-state topologies
    2. Generates Go contacts from contact maps
    3. Combines states into multiple-basin topology
    4. Validates output against reference

Methods tested:
    - EXP (exponential mixing)
    - HAM (harmonic mixing)
"""

import os
from functools import partial

import ctgomartini
from ctgomartini.api import GenMBPTop, MartiniTopFile
from ctgomartini.core import SameListList
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds
import subprocess

from tests.conftest import WorkingDirectoryContext


def Comparison_ITP(working_dir, molname, topfile):
    with WorkingDirectoryContext(os.path.join(working_dir, "test")):
        mbmol_test = MartiniTopFile(topfile)._moleculeTypes[molname]
    
    with WorkingDirectoryContext(os.path.join(working_dir, "ref")):
        mbmol_ref = MartiniTopFile(topfile)._moleculeTypes[molname]

    Comparison_Top(mbmol_ref, mbmol_test)

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
        same = SameListList([mbmol_ref._topology[category], mbmol_test._topology[category]], sort=True, precision=5)
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
        working_dir = os.path.join(self.path, "../fixtures/WriteItp/EXP")
        
        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test 2>/dev/null; cp -a template test')
            
            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Fetch ctgomarinize
                os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/ctgomartinize.py')} .")
                
                # Check which DSSP command is available (prefer mkdssp for compatibility)
                dssp_cmd = None
                if subprocess.run("which mkdssp", shell=True, capture_output=True).returncode == 0:
                    dssp_cmd = "mkdssp"
                elif subprocess.run("which dssp", shell=True, capture_output=True).returncode == 0:
                    dssp_cmd = "dssp"
                else:
                    raise RuntimeError("Neither 'dssp' nor 'mkdssp' command found in the system")
                
                # Generate Itp
                subprocess.run(f"python ctgomartinize.py -s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp {dssp_cmd} -ff martini3001 -method exp",
                               shell=True)
        
        # Check Comparison
        Comparison_ITP(working_dir, 'gbp', 'system.top')

    def test_GenMBItp_HAM(self):
        working_dir = os.path.join(self.path, "../fixtures/WriteItp/HAM")
        
        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test 2>/dev/null; cp -a template test')
            
            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Fetch ctgomarinize
                os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/ctgomartinize.py')} .")
                
                # Check which DSSP command is available (prefer mkdssp for compatibility)
                dssp_cmd = None
                if subprocess.run("which mkdssp", shell=True, capture_output=True).returncode == 0:
                    dssp_cmd = "mkdssp"
                elif subprocess.run("which dssp", shell=True, capture_output=True).returncode == 0:
                    dssp_cmd = "dssp"
                else:
                    raise RuntimeError("Neither 'dssp' nor 'mkdssp' command found in the system")
                
                # Generate Itp
                subprocess.run(f"python ctgomartinize.py -s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp {dssp_cmd} -ff martini3001 -method ham",
                               shell=True)
        
        # Check Comparison
        Comparison_ITP(working_dir, 'gbp', 'system.top')   

