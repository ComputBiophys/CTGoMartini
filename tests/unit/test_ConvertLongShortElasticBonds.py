"""
Unit tests for elastic bond conversion functionality.

Tests the convert_long_short_elastic_bonds function which modifies
elastic network bonds based on reference structure distances.
"""

import os
from functools import partial

from ctgomartini.api import GenMBPTop, MartiniTopFile
from ctgomartini.core import SameListList
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds
from tests.conftest import WorkingDirectoryContext


class TestMBMartiniTIP:
    """Context manager to save and restore working directory."""
    def __init__(self, new_dir):
        self.new_dir = new_dir
        self.original_dir = None
    
    def __enter__(self):
        self.original_dir = os.getcwd()
        os.chdir(self.new_dir)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.original_dir)
        return False


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
    Test elastic bond conversion for multiple-basin GoMartini.
    
    Attributes:
        path (str): Base path to test data directory.
    """
    
    path = os.path.dirname(__file__)

    def test_GlnBP_ITP_LongShortElasticBonds(self):
        """
        Test long/short elastic bond conversion for GlnBP.
        
        Converts elastic bonds based on open and closed reference structures,
        then generates and validates the resulting ITP file.
        
        Raises:
            AssertionError: If converted topology doesn't match reference.
        """
        working_dir = os.path.join(self.path, "../fixtures/WriteItp/LongShortElasticBondsConversion")
        
        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test')
            os.system('cp -a template test')
            
            # Process OpenItp
            with WorkingDirectoryContext(os.path.join(working_dir, 'test/OpenItp/')):
                prefix = 'gbp_open'
                ref_pdb = '../GBP_open_cg.pdb'
                convert_long_short_elastic_bonds(prefix, ref_pdb, convert_long_elastic_bonds=True, convert_short_elastic_bonds=False, lj_epsilon=12)

            # Process ClosedItp
            with WorkingDirectoryContext(os.path.join(working_dir, 'test/ClosedItp/')):
                prefix = 'gbp_closed'
                ref_pdb = '../GBP_closed_cg.pdb'
                convert_long_short_elastic_bonds(prefix, ref_pdb, convert_long_elastic_bonds=True, convert_short_elastic_bonds=False, lj_epsilon=12)

            # Write mb_itp
            with WorkingDirectoryContext(os.path.join(working_dir, 'test/')):
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
                write_itp(mbmol)

        # Check Comparison
        Comparison_ITP(working_dir, 'gbp', 'system.top')

    

