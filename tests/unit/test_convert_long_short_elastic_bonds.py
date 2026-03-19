"""
Unit tests for elastic bond conversion functionality.

Tests the convert_long_short_elastic_bonds function which modifies
elastic network bonds based on reference structure distances.
"""

import os
from functools import partial

from ctgomartini.topology import create_mb_topology, MartiniTopFile
from ctgomartini.topology.generator.combiner import SameListList
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds
from tests.conftest import WorkingDirectoryContext


def angles_dihedrals_sort(fields, category):
    """Sort angle/dihedral fields for consistent comparison."""
    if category in ['angles', 'multi_angles']:
        if int(fields[0]) > int(fields[2]):
            fields = fields[2::-1] + fields[3:]
    elif category in ['dihedrals', 'multi_dihedrals']:
        if int(fields[0]) > int(fields[3]):
            fields = fields[3::-1] + fields[4:]
    return fields


def compare_topology_sections(mbmol_ref, mbmol_test):
    """Compare topology sections between reference and test molecules."""
    categories_list = list(set(list(mbmol_ref._topology.keys()) + list(mbmol_test._topology.keys())))
    for category in categories_list:
        if category in ['angles', 'dihedrals', 'multi_angles', 'multi_dihedrals']:
            sort_partial = partial(angles_dihedrals_sort, category=category)
            mbmol_ref._topology[category] = list(map(sort_partial, mbmol_ref._topology[category]))
            mbmol_test._topology[category] = list(map(sort_partial, mbmol_test._topology[category]))
        same = SameListList([mbmol_ref._topology[category], mbmol_test._topology[category]], sort=True, precision=5)
        assert same is True, f"Error: comparison of {category} between test and ref is not the same!"


def compare_itp_topology(working_dir, molname, topfile):
    """Compare ITP topology between test and reference directories."""
    with WorkingDirectoryContext(os.path.join(working_dir, "test")):
        mbmol_test = MartiniTopFile(topfile)._moleculeTypes[molname]

    with WorkingDirectoryContext(os.path.join(working_dir, "ref")):
        mbmol_ref = MartiniTopFile(topfile)._moleculeTypes[molname]

    compare_topology_sections(mbmol_ref, mbmol_test)


class TestConvertLongShortElasticBonds:
    """
    Test elastic bond conversion for multiple-basin GoMartini.

    Attributes:
        path (str): Base path to test data directory.
    """

    path = os.path.dirname(__file__)

    def test_convert_elastic_bonds_glnbp(self):
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
                mbmol =  create_mb_topology(mols_list, mbmol_name, dict_cutoffs)
                write_itp(mbmol)

        # Check Comparison
        compare_itp_topology(working_dir, 'gbp', 'system.top')



