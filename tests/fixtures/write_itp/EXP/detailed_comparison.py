#!/usr/bin/env python3
"""
Detailed comparison of REF vs TEST for EXP system.
Outputs every item and their differences.
"""

import sys
import os

sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')
from ctgomartini.api import MartiniTopFile

def load_topology(path, molname):
    """Load topology and return molecule."""
    cwd = os.getcwd()
    try:
        os.chdir(path)
        top = MartiniTopFile('system.top')
        return top._moleculeTypes[molname]._topology
    finally:
        os.chdir(cwd)

def make_key(fields, category):
    """Create a key for matching items between ref and test."""
    if not fields:
        return ()
    if category in ['atoms']:
        return (fields[0],)
    elif category in ['bonds', 'constraints', 'contacts', 'exclusions']:
        return tuple(sorted(fields[:2]))
    elif category in ['angles', 'multi_angles']:
        return tuple(fields[:3])
    elif category in ['dihedrals', 'multi_dihedrals']:
        return tuple(fields[:4])
    elif category in ['multi_contacts']:
        # atom1, atom2, state
        return (fields[0], fields[1], fields[3]) if len(fields) >= 4 else tuple(fields)
    elif category in ['virtual_sitesn']:
        return (fields[0],)
    else:
        return tuple(fields[:2]) if len(fields) >= 2 else tuple(fields)

def compare_values(ref_fields, test_fields, abs_tol=1e-4):
    """Compare two field lists and return differences."""
    if len(ref_fields) != len(test_fields):
        return False, f"Field count: {len(ref_fields)} vs {len(test_fields)}"
    
    diffs = []
    for i, (r, t) in enumerate(zip(ref_fields, test_fields)):
        try:
            r_val = float(r)
            t_val = float(t)
            diff = abs(r_val - t_val)
            if diff > abs_tol:
                diffs.append(f"[{i}] {r} vs {t} (diff: {diff:.2e})")
        except (ValueError, TypeError):
            if r != t:
                diffs.append(f"[{i}] '{r}' vs '{t}' (string)")
    
    if diffs:
        return False, "; ".join(diffs)
    return True, "OK"

def compare_category(category, ref_data, test_data, output_file):
    """Compare a single category and write detailed output."""
    output_file.write(f"\n{'='*80}\n")
    output_file.write(f"CATEGORY: {category}\n")
    output_file.write(f"{'='*80}\n")
    output_file.write(f"REF count: {len(ref_data)}\n")
    output_file.write(f"TEST count: {len(test_data)}\n\n")
    
    if not ref_data and not test_data:
        output_file.write("Both empty.\n")
        return
    
    # Build dictionaries
    ref_dict = {make_key(f, category): (i, f) for i, f in enumerate(ref_data)}
    test_dict = {make_key(f, category): (i, f) for i, f in enumerate(test_data)}
    
    ref_keys = set(ref_dict.keys())
    test_keys = set(test_dict.keys())
    
    only_in_ref = ref_keys - test_keys
    only_in_test = test_keys - ref_keys
    common = ref_keys & test_keys
    
    # Report summary
    output_file.write(f"Summary:\n")
    output_file.write(f"  Only in REF: {len(only_in_ref)}\n")
    output_file.write(f"  Only in TEST: {len(only_in_test)}\n")
    output_file.write(f"  Common: {len(common)}\n\n")
    
    # Only in REF
    if only_in_ref:
        output_file.write(f"{'-'*40}\n")
        output_file.write("ITEMS ONLY IN REF:\n")
        output_file.write(f"{'-'*40}\n")
        for key in sorted(only_in_ref):
            idx, fields = ref_dict[key]
            output_file.write(f"  REF[{idx}]: {fields}\n")
        output_file.write("\n")
    
    # Only in TEST
    if only_in_test:
        output_file.write(f"{'-'*40}\n")
        output_file.write("ITEMS ONLY IN TEST:\n")
        output_file.write(f"{'-'*40}\n")
        for key in sorted(only_in_test):
            idx, fields = test_dict[key]
            output_file.write(f"  TEST[{idx}]: {fields}\n")
        output_file.write("\n")
    
    # Common items with differences
    output_file.write(f"{'-'*40}\n")
    output_file.write("COMMON ITEMS:\n")
    output_file.write(f"{'-'*40}\n")
    
    diff_count = 0
    for key in sorted(common):
        ref_idx, ref_fields = ref_dict[key]
        test_idx, test_fields = test_dict[key]
        
        same, msg = compare_values(ref_fields, test_fields)
        
        if not same:
            diff_count += 1
            output_file.write(f"\nKey: {key}\n")
            output_file.write(f"  REF[{ref_idx}]:  {ref_fields}\n")
            output_file.write(f"  TEST[{test_idx}]: {test_fields}\n")
            output_file.write(f"  Diff: {msg}\n")
    
    if diff_count == 0:
        output_file.write("  All common items match within tolerance.\n")
    else:
        output_file.write(f"\n  Total differences in common items: {diff_count}\n")

def main():
    exp_dir = '/home/ys/CommonUse/Martini/CTGoMartini/tests/fixtures/write_itp/EXP'
    
    print("Loading REF topology...")
    ref_top = load_topology(os.path.join(exp_dir, 'ref'), 'gbp')
    
    print("Loading TEST topology...")
    test_top = load_topology(os.path.join(exp_dir, 'test'), 'gbp')
    
    output_path = os.path.join(exp_dir, 'comparison_output.txt')
    
    print(f"Writing detailed comparison to: {output_path}")
    
    with open(output_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("DETAILED COMPARISON: REF vs TEST for EXP System\n")
        f.write("="*80 + "\n")
        f.write(f"REF categories: {sorted(ref_top.keys())}\n")
        f.write(f"TEST categories: {sorted(test_top.keys())}\n")
        
        # Compare each category
        all_categories = sorted(set(ref_top.keys()) | set(test_top.keys()))
        
        for category in all_categories:
            ref_data = ref_top.get(category, [])
            test_data = test_top.get(category, [])
            compare_category(category, ref_data, test_data, f)
        
        f.write("\n" + "="*80 + "\n")
        f.write("COMPARISON COMPLETE\n")
        f.write("="*80 + "\n")
    
    print(f"Done! Output written to: {output_path}")

if __name__ == '__main__':
    main()
