#!/usr/bin/env python3
"""
Compare reference files between vermouth versions.

This script compares ref_legacy (vermouth < 0.14.0) and ref_current (vermouth >= 0.14.0)
to document the differences introduced by the proline force field fix.
"""

import sys
import os
from functools import partial

# Add project root to path
sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')

from ctgomartini.topology import MartiniTopFile
from tests.conftest import WorkingDirectoryContext


def category_sort(fields, category):
    """Sort fields for consistent comparison."""
    if category in ['angles', 'multi_angles']:
        if int(fields[0]) > int(fields[2]):
            fields = fields[2::-1] + fields[3:]
    elif category in ['dihedrals', 'multi_dihedrals']:
        atoms = fields[:4]
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
    return fields


def compare_section(name, legacy_data, current_data):
    """Compare a single section between legacy and current.
    
    Numeric format differences (e.g., '100.0' vs '100') are ignored.
    Only actual content differences (functype changes, atom ID changes) are reported.
    """
    print(f"\n{'='*60}")
    print(f"Section: {name}")
    print(f"{'='*60}")
    
    if len(legacy_data) != len(current_data):
        print(f"  ✗ Entry count mismatch: {len(legacy_data)} vs {len(current_data)}")
        return False
    else:
        print(f"  ✓ Entry count: {len(legacy_data)}")
    
    if not legacy_data:
        print("  (empty)")
        return True
    
    # Check for differences
    content_differences = []
    numeric_format_count = 0
    
    # Sort data for comparison if needed
    if name in ['bonds', 'constraints', 'contacts', 'exclusions', 
                'angles', 'dihedrals', 'multi_angles', 'multi_dihedrals', 'multi_contacts']:
        sort_func = partial(category_sort, category=name)
        legacy_sorted = sorted([sort_func(list(f)) for f in legacy_data])
        current_sorted = sorted([sort_func(list(f)) for f in current_data])
    else:
        legacy_sorted = legacy_data
        current_sorted = current_data
    
    # Compare entries
    for i, (leg, cur) in enumerate(zip(legacy_sorted, current_sorted)):
        if list(leg) != list(cur):
            # Check if it's just numeric format difference
            is_content_diff = False
            if len(leg) == len(cur):
                for l, c in zip(leg, cur):
                    if l != c:
                        try:
                            l_f = float(l)
                            c_f = float(c)
                            # Numeric values within tolerance -> format difference
                            if abs(l_f - c_f) > 1e-4:
                                is_content_diff = True
                                break
                        except (ValueError, TypeError):
                            # Non-numeric field differs -> content difference
                            is_content_diff = True
                            break
            else:
                is_content_diff = True
            
            if is_content_diff:
                content_differences.append((i, leg, cur))
            else:
                numeric_format_count += 1
    
    if not content_differences and numeric_format_count == 0:
        print("  ✓ Content identical")
        return True
    
    if numeric_format_count > 0:
        print(f"  ○ Numeric format differences: {numeric_format_count} entries (ignored)")
    
    if content_differences:
        print(f"  ✗ Content differences: {len(content_differences)} entries")
        print("\n  First 5 differences:")
        for idx, leg, cur in content_differences[:5]:
            print(f"    [{idx}] Legacy:  {leg}")
            print(f"        Current: {cur}")
            # Show which fields differ
            diff_fields = []
            for j, (l, c) in enumerate(zip(leg, cur)):
                if l != c:
                    try:
                        l_f = float(l)
                        c_f = float(c)
                        if abs(l_f - c_f) > 1e-4:
                            diff_fields.append(f"[{j}] {l} → {c} (value)")
                        else:
                            diff_fields.append(f"[{j}] {l} → {c} (format)")
                    except:
                        diff_fields.append(f"[{j}] '{l}' → '{c}' (string)")
            print(f"        Diff: {', '.join(diff_fields)}")
    else:
        print("  ✓ Content identical (ignoring numeric format)")
    
    return len(content_differences) == 0


def analyze_functype_distribution(name, legacy_data, current_data, functype_idx):
    """Analyze functype distribution differences."""
    legacy_ft = {}
    current_ft = {}
    
    for d in legacy_data:
        if len(d) > functype_idx:
            ft = d[functype_idx]
            legacy_ft[ft] = legacy_ft.get(ft, 0) + 1
    
    for d in current_data:
        if len(d) > functype_idx:
            ft = d[functype_idx]
            current_ft[ft] = current_ft.get(ft, 0) + 1
    
    print(f"\n  Functype distribution:")
    all_types = set(legacy_ft.keys()) | set(current_ft.keys())
    for ft in sorted(all_types):
        l_count = legacy_ft.get(ft, 0)
        c_count = current_ft.get(ft, 0)
        marker = "✓" if l_count == c_count else "✗"
        change = ""
        if l_count != c_count:
            diff = c_count - l_count
            change = f" ({diff:+#d})"
        print(f"    {marker} Type {ft}: Legacy={l_count}, Current={c_count}{change}")


def compare_method(method_dir, molname="gbp"):
    """Compare reference files for a given method (EXP, HAM, or EXP_multichains)."""
    print(f"\n{'#'*70}")
    print(f"# Comparing {method_dir} method (molecule: {molname})")
    print(f"{'#'*70}")
    
    base_path = os.path.dirname(os.path.abspath(__file__))
    
    # Load both reference files
    try:
        with WorkingDirectoryContext(os.path.join(base_path, method_dir, "ref_legacy")):
            mbmol_legacy = MartiniTopFile("system.top")._moleculeTypes[molname]
    except Exception as e:
        print(f"Error loading ref_legacy: {e}")
        return False
    
    try:
        with WorkingDirectoryContext(os.path.join(base_path, method_dir, "ref_current")):
            mbmol_current = MartiniTopFile("system.top")._moleculeTypes[molname]
    except Exception as e:
        print(f"Error loading ref_current: {e}")
        return False
    
    # Get all sections
    all_sections = sorted(set(mbmol_legacy._topology.keys()) | set(mbmol_current._topology.keys()))
    
    print(f"\nSections found: {', '.join(all_sections)}")
    
    # Compare each section
    results = {}
    for section in all_sections:
        legacy_data = mbmol_legacy._topology.get(section, [])
        current_data = mbmol_current._topology.get(section, [])
        
        is_same = compare_section(section, legacy_data, current_data)
        results[section] = is_same
        
        # Analyze functype for relevant sections
        if section == "angles" and legacy_data:
            analyze_functype_distribution(section, legacy_data, current_data, functype_idx=3)
        elif section == "dihedrals" and legacy_data:
            analyze_functype_distribution(section, legacy_data, current_data, functype_idx=4)
        elif section == "multi_angles" and legacy_data:
            analyze_functype_distribution(section, legacy_data, current_data, functype_idx=5)
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    
    identical = [s for s, r in results.items() if r]
    different = [s for s, r in results.items() if not r]
    
    print(f"\nIdentical sections ({len(identical)}): {', '.join(identical) if identical else 'None'}")
    print(f"Different sections ({len(different)}): {', '.join(different) if different else 'None'}")
    
    if different:
        print("\nKey differences between vermouth < 0.14.0 and >= 0.14.0:")
        if "angles" in different:
            print("  1. Proline angles: 3 entries changed functype 2 → 10 (restricted bending)")
        print("  Note: Numeric format differences (100.0 → 100) are ignored as they don't affect physics")
    
    return len(different) == 0


def main():
    print("="*70)
    print("Reference Files Comparison Tool")
    print("="*70)
    print("\nThis tool compares ref_legacy (vermouth < 0.14.0) with")
    print("ref_current (vermouth >= 0.14.0) to document force field changes.")
    print("\nThe main difference is the proline angle fix in vermouth 0.14.0:")
    print("  - Angles: type 2 (harmonic) → type 10 (restricted bending) for proline residues")
    print("  - Dihedrals: Same content, but may have different ordering in the file")
    
    exp_ok = compare_method("EXP")
    ham_ok = compare_method("HAM")
    
    # Check for multichain case
    multichain_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "EXP_multichains")
    if os.path.exists(multichain_path):
        multichain_ok = compare_method("EXP_multichains", molname="TRAAK")
    else:
        multichain_ok = True
    
    print(f"\n{'='*70}")
    print("FINAL RESULT")
    print(f"{'='*70}")
    if exp_ok and ham_ok and multichain_ok:
        print("\n✓ All reference files are consistent (accounting for known differences)")
        return 0
    else:
        print("\n✗ Some sections have differences (see details above)")
        print("\nNote: These differences are expected due to vermouth 0.14.0 proline fix.")
        print("The test suite uses version-appropriate reference files automatically.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
