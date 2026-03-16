#!/usr/bin/env python3
"""
Detailed comparison of REF vs TEST for EXP system (with sorting).
"""

import sys
import os

sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')
from ctgomartini.api import MartiniTopFile

def load_topology(path, molname):
    cwd = os.getcwd()
    try:
        os.chdir(path)
        top = MartiniTopFile('system.top')
        return top._moleculeTypes[molname]._topology
    finally:
        os.chdir(cwd)

def sort_dihedral(fields):
    """Sort dihedral atoms for consistent comparison."""
    if len(fields) >= 4:
        atoms = fields[:4]
        # Sort to get canonical form
        if int(atoms[0]) > int(atoms[3]):
            atoms = list(reversed(atoms))
        return tuple(atoms)
    return tuple(fields)

def main():
    exp_dir = '/home/ys/CommonUse/Martini/CTGoMartini/tests/fixtures/WriteItp/EXP'
    
    print("Loading REF topology...")
    ref_top = load_topology(os.path.join(exp_dir, 'ref'), 'gbp')
    
    print("Loading TEST topology...")
    test_top = load_topology(os.path.join(exp_dir, 'test'), 'gbp')
    
    output_path = os.path.join(exp_dir, 'comparison_sorted.txt')
    
    with open(output_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("DETAILED COMPARISON (WITH SORTING): REF vs TEST for EXP System\n")
        f.write("="*80 + "\n\n")
        
        # Focus on dihedrals which have ordering issues
        category = 'dihedrals'
        ref_data = ref_top.get(category, [])
        test_data = test_top.get(category, [])
        
        f.write(f"CATEGORY: {category}\n")
        f.write(f"REF count: {len(ref_data)}\n")
        f.write(f"TEST count: {len(test_data)}\n\n")
        
        # Sort dihedrals
        ref_sorted = {}
        for i, fields in enumerate(ref_data):
            key = sort_dihedral(fields)
            ref_sorted[key] = (i, fields)
        
        test_sorted = {}
        for i, fields in enumerate(test_data):
            key = sort_dihedral(fields)
            test_sorted[key] = (i, fields)
        
        ref_keys = set(ref_sorted.keys())
        test_keys = set(test_sorted.keys())
        
        only_in_ref = ref_keys - test_keys
        only_in_test = test_keys - ref_keys
        common = ref_keys & test_keys
        
        f.write(f"After sorting atoms:\n")
        f.write(f"  Only in REF: {len(only_in_ref)}\n")
        f.write(f"  Only in TEST: {len(only_in_test)}\n")
        f.write(f"  Common: {len(common)}\n\n")
        
        f.write("="*80 + "\n")
        f.write("PAIRWISE COMPARISON (sorted by dihedral key):\n")
        f.write("="*80 + "\n\n")
        
        diff_count = 0
        for key in sorted(common):
            ref_idx, ref_fields = ref_sorted[key]
            test_idx, test_fields = test_sorted[key]
            
            # Check for differences
            diffs = []
            for i, (r, t) in enumerate(zip(ref_fields, test_fields)):
                try:
                    r_val = float(r)
                    t_val = float(t)
                    diff = abs(r_val - t_val)
                    if diff > 1e-4:
                        diffs.append(f"[{i}] {r} vs {t} (diff: {diff:.2e})")
                except (ValueError, TypeError):
                    if r != t:
                        diffs.append(f"[{i}] '{r}' vs '{t}'")
            
            if diffs:
                diff_count += 1
                f.write(f"Key: {key} - DIFFERENCES:\n")
                f.write(f"  REF[{ref_idx}]:  {ref_fields}\n")
                f.write(f"  TEST[{test_idx}]: {test_fields}\n")
                f.write(f"  Diff: {'; '.join(diffs)}\n\n")
            else:
                f.write(f"Key: {key} - MATCH\n")
                f.write(f"  REF:  {ref_fields}\n")
                f.write(f"  TEST: {test_fields}\n\n")
        
        f.write("="*80 + "\n")
        f.write(f"SUMMARY: {diff_count} items have value differences\n")
        f.write("="*80 + "\n")
    
    print(f"Done! Output written to: {output_path}")

if __name__ == '__main__':
    main()
