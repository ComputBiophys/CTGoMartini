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
"""

import os
from functools import partial

import ctgomartini
from ctgomartini.topology import create_mb_topology, MartiniTopFile
from ctgomartini.topology.generator.combiner import SameListList
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds
import subprocess

from tests.conftest import WorkingDirectoryContext


def compare_itp_topology(working_dir, molname, topfile):
    """Compare ITP topology between test and reference directories."""
    with WorkingDirectoryContext(os.path.join(working_dir, "test")):
        mbmol_test = MartiniTopFile(topfile)._moleculeTypes[molname]
    
    with WorkingDirectoryContext(os.path.join(working_dir, "ref")):
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
        # The max observed difference between vermouth 0.9.6 and 0.13.0 is ~8e-5 nm
        # Pass category to handle functype differences in angles/dihedrals
        same, msg = compare_with_tolerance(ref_data, test_data, key_func, abs_tol=1e-4, category=category)
        assert same, f"Error: comparison of {category} between test and ref is not the same: {msg}"





class TestGenMBItp:
    """
    Test the ITP file of multiple basin GoMartini.
    """
    path = os.path.dirname(__file__)

    def test_generate_mb_itp_exp(self):
        working_dir = os.path.join(self.path, "../fixtures/WriteItp/EXP")
        
        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test 2>/dev/null; cp -a template test')
            
            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Fetch ctgomarinize
                os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'cli/ctgomartinize.py')} .")
                
                # Generate Itp
                # Note: vermouth >= 0.15.0 uses MDTraj for DSSP, -dssp flag without argument
                subprocess.run(f"python ctgomartinize.py -s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp -ff martini3001 -method exp",
                               shell=True)
        
        # Check Comparison
        compare_itp_topology(working_dir, 'gbp', 'system.top')

    def test_generate_mb_itp_ham(self):
        working_dir = os.path.join(self.path, "../fixtures/WriteItp/HAM")
        
        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test 2>/dev/null; cp -a template test')
            
            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Fetch ctgomarinize
                os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'cli/ctgomartinize.py')} .")
                
                # Generate Itp
                # Note: vermouth >= 0.15.0 uses MDTraj for DSSP, -dssp flag without argument
                subprocess.run(f"python ctgomartinize.py -s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp -ff martini3001 -method ham",
                               shell=True)
        
        # Check Comparison
        compare_itp_topology(working_dir, 'gbp', 'system.top')
