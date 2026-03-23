"""
End-to-end tests for multiple-basin ITP generation.

Tests the complete ctgomartinize.py workflow which:
    1. Runs Martinize2 to generate single-state topologies
    2. Generates Go contacts from contact maps
    3. Combines states into multiple-basin topology
    4. Validates output against reference

Methods tested:
    - EXP (exponential mixing)
    - HAM (Hamiltonian mixing)
    - EXP_multichains (multi-chain membrane protein)

Reference files:
    - ref_legacy: Generated with vermouth < 0.14.0 (older force field)
    - ref_current: Generated with vermouth >= 0.14.0 (proline angle fix)
"""

import os
from functools import partial

import ctgomartini
from ctgomartini.topology import create_mb_topology, MartiniTopFile
from ctgomartini.topology.generator.combiner import SameListList
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds
import subprocess

from tests.conftest import WorkingDirectoryContext


def get_vermouth_version() -> tuple[int, int, int]:
    """Get installed vermouth version as tuple (major, minor, patch)."""
    try:
        import vermouth
        version_str = vermouth.__version__
        parts = version_str.split('.')
        return tuple(int(p) for p in parts[:3])
    except Exception:
        return (0, 0, 0)


def is_legacy_vermouth() -> bool:
    """Check if installed vermouth is < 0.14.0 (legacy version)."""
    major, minor, _ = get_vermouth_version()
    return (major, minor) < (0, 14)


def get_ref_dir(base_dir: str) -> str:
    """Get appropriate reference directory based on vermouth version."""
    if is_legacy_vermouth():
        return os.path.join(base_dir, "ref_legacy")
    return os.path.join(base_dir, "ref_current")


def compare_itp_topology(working_dir: str, molname: str, topfile: str, ref_subdir: str = "ref") -> None:
    """Compare ITP topology between test and reference directories.
    
    Args:
        working_dir: Base working directory containing test/ and ref*/ subdirs
        molname: Molecule name to compare
        topfile: Topology file name (e.g., 'system.top')
        ref_subdir: Reference subdirectory name ('ref_legacy' or 'ref_current')
    """
    with WorkingDirectoryContext(os.path.join(working_dir, "test")):
        mbmol_test = MartiniTopFile(topfile)._moleculeTypes[molname]

    with WorkingDirectoryContext(os.path.join(working_dir, ref_subdir)):
        mbmol_ref = MartiniTopFile(topfile)._moleculeTypes[molname]

    compare_topology_sections(mbmol_ref, mbmol_test)


def category_sort(fields, category):
    if category in ['angles', 'multi_angles']:
        if int(fields[0]) >  int(fields[2]):
            fields = fields[2::-1] + fields[3:]
    elif category in ['dihedrals', 'multi_dihedrals']:
        # For dihedrals, we need to handle cases where the atom order is completely reversed
        # between REF and TEST (e.g., [a,b,c,d] vs [d,c,b,a])
        atoms = fields[:4]
        # Create a canonical form by sorting the tuple in both directions and picking the smaller one
        forward = tuple(atoms)
        backward = tuple(reversed(atoms))
        canonical = min(forward, backward)
        fields = list(canonical) + fields[4:]
    elif category in ['bonds', 'constraints', 'contacts', 'multi_contacts']:
        if int(fields[0]) > int(fields[1]):
            fields = fields[1::-1] + fields[2:]
    elif category in ['exclusions']:
        fields = list(map(int, fields))
        fields = [fields[0]] + sorted(fields[1:])
        fields = list(map(str, fields))
    else:
        raise ValueError(f'Unresolved category, {category}: {fields}')
    return fields


def compare_with_tolerance(ref_list, test_list, key_func, abs_tol=1e-4, category=None):
    """Compare two lists with absolute tolerance for numeric values.

    Args:
        ref_list: Reference list of field lists
        test_list: Test list of field lists
        key_func: Function to extract key from fields for matching
        abs_tol: Absolute tolerance for numeric comparison (default: 1e-4)
        category: Category name for special handling (e.g., ignoring functype differences)

    Returns:
        (is_equal, message) tuple
    """
    if len(ref_list) != len(test_list):
        return False, f"Count mismatch: {len(ref_list)} vs {len(test_list)}"

    # Build dictionaries
    ref_dict = {key_func(f): f for f in ref_list}
    test_dict = {key_func(f): f for f in test_list}

    ref_keys = set(ref_dict.keys())
    test_keys = set(test_dict.keys())

    only_in_ref = ref_keys - test_keys
    only_in_test = test_keys - ref_keys

    if only_in_ref:
        return False, f"Keys only in ref: {list(only_in_ref)[:3]}..."
    if only_in_test:
        return False, f"Keys only in test: {list(only_in_test)[:3]}..."

    # Compare values with tolerance
    for key in ref_keys:
        ref_fields = ref_dict[key]
        test_fields = test_dict[key]

        if len(ref_fields) != len(test_fields):
            return False, f"Field count mismatch at {key}: {len(ref_fields)} vs {len(test_fields)}"

        for i, (r, t) in enumerate(zip(ref_fields, test_fields)):

            # Try numeric comparison with tolerance
            try:
                r_val = float(r)
                t_val = float(t)
                if abs(r_val - t_val) > abs_tol:
                    return False, f"Value mismatch at {key}[{i}]: {r} vs {t} (diff: {abs(r_val - t_val):.2e})"
            except (ValueError, TypeError):
                # Non-numeric: exact comparison
                if r != t:
                    return False, f"String mismatch at {key}[{i}]: {r} vs {t}"

    return True, "OK"


def compare_topology_sections(mbmol_ref, mbmol_test):
    categories_list = list(set(list(mbmol_ref._topology.keys()) + list(mbmol_test._topology.keys())))
    for category in categories_list:
        ref_data = mbmol_ref._topology[category]
        test_data = mbmol_test._topology[category]

        # Sort angles/dihedrals for consistent comparison
        if category in ['bond', 'constraints', 'contacts', 'exclusions',
                        'angles', 'dihedrals',
                        'multi_angles', 'multi_dihedrals', 'multi_contacts']:
            category_sort_partial = partial(category_sort, category=category)
            ref_data = list(map(category_sort_partial, ref_data))
            test_data = list(map(category_sort_partial, test_data))

        # Define key functions for each category
        if category in ['atoms', 'exclusions']:
            key_func = lambda f: (f[0],) if f else ()
        elif category in ['bonds', 'constraints', 'contacts']:
            key_func = lambda f: tuple(sorted(f[:2])) if len(f) >= 2 else tuple(f)
        elif category in ['angles']:
            key_func = lambda f: tuple(f[:3]) if len(f) >= 3 else tuple(f)
        elif category in ['dihedrals']:
            key_func = lambda f: tuple(f[:4]) if len(f) >= 4 else tuple(f)
        elif category in ['multi_angles']:
            key_func = lambda f: tuple(f[:5]) if len(f) >= 5 else tuple(f)
        elif category in ['multi_dihedrals']:
            key_func = lambda f: tuple(f[:6]) if len(f) >= 6 else tuple(f)
        elif category in ['multi_contacts']:
            key_func = lambda f: tuple(f[:4]) if len(f) >= 4 else tuple(f)
        elif category in ['virtual_sitesn']:
            key_func = lambda f: (f[0],) if f else ()
        else:
            key_func = lambda f: tuple(f)

        # Use tolerance-based comparison for numeric values
        # The max observed difference between vermouth versions is ~8e-5 nm
        # Pass category to handle functype differences in angles/dihedrals
        same, msg = compare_with_tolerance(ref_data, test_data, key_func, abs_tol=1e-4, category=category)
        assert same, f"Error: comparison of {category} between test and ref is not the same: {msg}"


class TestGenMBItp:
    """
    Test the ITP file of multiple basin GoMartini.

    These tests compare generated topologies against reference files.
    Two reference versions are maintained:
        - ref_legacy: For vermouth < 0.14.0
        - ref_current: For vermouth >= 0.14.0 (with proline angle fix)
    """
    path = os.path.dirname(__file__)

    def test_generate_mb_itp_exp(self):
        """Test EXP method multiple-basin topology generation."""
        working_dir = os.path.join(self.path, "../fixtures/write_itp/EXP")
        ref_subdir = "ref_legacy" if is_legacy_vermouth() else "ref_current"

        with WorkingDirectoryContext(working_dir):
            os.system('rm -rf test 2>/dev/null; cp -a template test')

            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Generate ITP using ctgomartinize
                subprocess.run(
                    f"python -m ctgomartini.cli.ctgomartinize "
                    f"-s 1GGG_1_clean.pdb 1WDN_1_clean.pdb "
                    f"-m 1GGG_1_clean.map 1WDN_1_clean.map "
                    f"-mol gbp_open gbp_closed "
                    f"-mbmol gbp "
                    f"-ff martini3001 "
                    f"-method exp",
                    shell=True,
                    check=True
                )

        # Compare with appropriate reference
        compare_itp_topology(working_dir, 'gbp', 'system.top', ref_subdir)

    def test_generate_mb_itp_ham(self):
        """Test HAM method multiple-basin topology generation."""
        working_dir = os.path.join(self.path, "../fixtures/write_itp/HAM")
        ref_subdir = "ref_legacy" if is_legacy_vermouth() else "ref_current"

        with WorkingDirectoryContext(working_dir):
            os.system('rm -rf test 2>/dev/null; cp -a template test')

            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Generate ITP using ctgomartinize
                subprocess.run(
                    f"python -m ctgomartini.cli.ctgomartinize "
                    f"-s 1GGG_1_clean.pdb 1WDN_1_clean.pdb "
                    f"-m 1GGG_1_clean.map 1WDN_1_clean.map "
                    f"-mol gbp_open gbp_closed "
                    f"-mbmol gbp "
                    f"-ff martini3001 "
                    f"-method ham",
                    shell=True,
                    check=True
                )

        # Compare with appropriate reference
        compare_itp_topology(working_dir, 'gbp', 'system.top', ref_subdir)

    def test_generate_mb_itp_exp_multichains(self):
        """Test EXP method with multi-chain membrane protein (TRAAK).
        
        This test validates the multi-chain topology generation for the
        TRAAK potassium channel. 
        
        Note: Multi-chain Go-Martini support requires vermouth >= 0.14.0.
        For older versions, this test is skipped.
        """
        import pytest
        
        # Skip for vermouth < 0.14.0 (no multi-chain support)
        if is_legacy_vermouth():
            pytest.skip("Multi-chain Go-Martini requires vermouth >= 0.14.0")
        
        working_dir = os.path.join(self.path, "../fixtures/write_itp/EXP_multichains")
        ref_subdir = "ref_current"  # Only test with >= 0.14.0

        with WorkingDirectoryContext(working_dir):
            os.system('rm -rf test 2>/dev/null; cp -a template test')

            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Generate ITP using ctgomartinize
                # Note: Order is 7LJB (Up) first, then 4WFE (Down)
                subprocess.run(
                    f"python -m ctgomartini.cli.ctgomartinize "
                    f"-s 7LJB_clean.pdb 4WFE_clean.pdb "
                    f"-m 7LJB_clean.map 4WFE_clean.map "
                    f"-mol Up Down "
                    f"-mbmol TRAAK "
                    f"-dssp "
                    f"-ff martini3001 "
                    f"-method exp "
                    f"-other_params='-nt -cys 0.25'",
                    shell=True,
                    check=True
                )

        # Compare with reference
        compare_itp_topology(working_dir, 'TRAAK', 'system.top', ref_subdir)
