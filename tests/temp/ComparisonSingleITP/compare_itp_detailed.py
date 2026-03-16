#!/usr/bin/env python3
"""
Detailed comparison of GROMACS ITP file sections

This script performs deep comparison of topology sections including
angles, constraints, dihedrals, and virtual_sitesn.

Usage:
    python compare_itp_detailed.py <file1.itp> <file2.itp> [options]
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def parse_itp_file(filepath: str) -> dict:
    """Parse ITP file and extract all sections"""
    sections = defaultdict(list)
    current_section = None
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            stripped = line.strip()
            
            if not stripped:
                continue
            
            section_match = re.match(r'^\[\s*(\w+)\s*\]', stripped)
            if section_match:
                current_section = section_match.group(1)
                continue
            
            if stripped.startswith(';') or stripped.startswith('#'):
                continue
            
            if current_section:
                sections[current_section].append((line_num, stripped))
    
    return dict(sections)


def parse_angles(lines: list) -> dict:
    """Parse angle data: i j k funct angle force_constant"""
    angles = {}
    for line_num, line in lines:
        parts = line.split()
        if len(parts) >= 5:
            i, j, k = int(parts[0]), int(parts[1]), int(parts[2])
            key = (i, j, k)  # Keep order for angles
            angles[key] = {
                'funct': parts[3],
                'angle': float(parts[4]),
                'fc': float(parts[5]) if len(parts) > 5 else None
            }
    return angles


def parse_constraints(lines: list) -> dict:
    """Parse constraint data: i j funct length"""
    constraints = {}
    for line_num, line in lines:
        parts = line.split()
        if len(parts) >= 4:
            i, j = int(parts[0]), int(parts[1])
            key = tuple(sorted([i, j]))
            constraints[key] = {
                'funct': parts[2],
                'length': float(parts[3])
            }
    return constraints


def parse_dihedrals(lines: list) -> dict:
    """Parse dihedral data: i j k l funct params"""
    dihedrals = {}
    for line_num, line in lines:
        # Remove comments
        if ';' in line:
            line = line.split(';')[0].strip()
        parts = line.split()
        if len(parts) >= 5:
            i, j, k, l = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
            key = (i, j, k, l)
            try:
                params = [float(p) for p in parts[5:] if p]
            except ValueError:
                params = parts[5:]  # Keep as strings if can't convert
            dihedrals[key] = {
                'funct': parts[4],
                'params': params
            }
    return dihedrals


def parse_virtual_sitesn(lines: list) -> dict:
    """Parse virtual_sitesn data: site_atom funct constructing_atoms..."""
    vsites = {}
    for line_num, line in lines:
        parts = line.split()
        if len(parts) >= 3:
            site_atom = int(parts[0])
            funct = parts[1]
            constructing = [int(p) for p in parts[2:]]
            key = (site_atom, tuple(sorted(constructing)))
            vsites[key] = {
                'site_atom': site_atom,
                'funct': funct,
                'constructing': constructing
            }
    return vsites


def parse_exclusions(lines: list) -> dict:
    """Parse exclusion data: atom excluded_atoms..."""
    exclusions = {}
    for line_num, line in lines:
        parts = line.split()
        if len(parts) >= 2:
            atom = int(parts[0])
            excluded = tuple(sorted([int(p) for p in parts[1:]]))
            exclusions[atom] = excluded
    return exclusions


def compare_section_detail(name: str, data1: dict, data2: dict, tolerance: float = 0.001) -> dict:
    """Generic detailed section comparison"""
    keys1 = set(data1.keys())
    keys2 = set(data2.keys())
    
    common = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    
    param_diffs = []
    for key in common:
        d1 = data1[key]
        d2 = data2[key]
        
        # Handle tuple values (like exclusions)
        if isinstance(d1, tuple) and isinstance(d2, tuple):
            if d1 != d2:
                param_diffs.append({
                    'key': key,
                    'diffs': {'exclusions': (d1, d2)}
                })
            continue
        
        # Compare all fields for dict values
        if isinstance(d1, dict) and isinstance(d2, dict):
            diff_fields = {}
            for field in d1:
                if field not in d2:
                    continue
                v1, v2 = d1[field], d2[field]
                if isinstance(v1, float) and isinstance(v2, float):
                    if abs(v1 - v2) > tolerance:
                        diff_fields[field] = (v1, v2)
                elif v1 != v2:
                    diff_fields[field] = (v1, v2)
            
            if diff_fields:
                param_diffs.append({
                    'key': key,
                    'diffs': diff_fields,
                    'file1': d1,
                    'file2': d2
                })
    
    return {
        'name': name,
        'total1': len(data1),
        'total2': len(data2),
        'common': len(common),
        'only_in_1': sorted(only_in_1),
        'only_in_2': sorted(only_in_2),
        'param_diffs': param_diffs
    }


def print_detailed_report(result: dict, show_examples: bool = True):
    """Print detailed comparison report"""
    
    print("\n" + "=" * 80)
    print(f"DETAILED COMPARISON: [{result['name'].upper()}]")
    print("=" * 80)
    
    print(f"\n  Total in File 1: {result['total1']}")
    print(f"  Total in File 2: {result['total2']}")
    print(f"  Common: {result['common']}")
    
    only_1 = result['only_in_1']
    only_2 = result['only_in_2']
    param_diffs = result['param_diffs']
    
    if only_1:
        print(f"\n  Only in File 1: {len(only_1)} entries")
        if show_examples and len(only_1) <= 10:
            for key in only_1:
                print(f"    {key}")
        elif show_examples:
            print(f"    Examples: {only_1[:5]}")
    
    if only_2:
        print(f"\n  Only in File 2: {len(only_2)} entries")
        if show_examples and len(only_2) <= 10:
            for key in only_2:
                print(f"    {key}")
        elif show_examples:
            print(f"    Examples: {only_2[:5]}")
    
    if param_diffs:
        print(f"\n  Parameter differences: {len(param_diffs)} entries")
        if show_examples:
            print("  " + "-" * 76)
            for diff in param_diffs[:5]:
                key = diff['key']
                print(f"    {result['name']} {key}:")
                for field, (v1, v2) in diff['diffs'].items():
                    print(f"      {field}: {v1} -> {v2}")
            if len(param_diffs) > 5:
                print(f"    ... and {len(param_diffs) - 5} more")
    
    if not (only_1 or only_2 or param_diffs):
        print(f"\n  STATUS: IDENTICAL")
    else:
        match_percent = (result['common'] / max(result['total1'], result['total2'])) * 100
        print(f"\n  Match: {match_percent:.1f}%")


def main():
    parser = argparse.ArgumentParser(
        description='Detailed comparison of ITP file sections',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s Old/Open.itp New/Open.itp
  %(prog)s Old/Open.itp New/Open.itp --section angles
        """
    )
    
    parser.add_argument('file1', help='First ITP file')
    parser.add_argument('file2', help='Second ITP file')
    parser.add_argument('--section', choices=['angles', 'constraints', 'dihedrals', 
                                               'virtual_sitesn', 'exclusions', 'all'],
                       default='all', help='Which section to compare (default: all)')
    parser.add_argument('--no-examples', action='store_true',
                       help='Do not show example differences')
    
    args = parser.parse_args()
    
    # Check files exist
    for f in [args.file1, args.file2]:
        if not Path(f).exists():
            print(f"Error: File not found: {f}")
            sys.exit(1)
    
    file1_path = Path(args.file1)
    file2_path = Path(args.file2)
    
    print(f"Parsing {file1_path.name}...")
    sections1 = parse_itp_file(str(file1_path))
    
    print(f"Parsing {file2_path.name}...")
    sections2 = parse_itp_file(str(file2_path))
    
    sections_to_compare = []
    if args.section == 'all':
        sections_to_compare = ['angles', 'constraints', 'dihedrals', 'virtual_sitesn', 'exclusions']
    else:
        sections_to_compare = [args.section]
    
    show_examples = not args.no_examples
    
    for section_name in sections_to_compare:
        if section_name not in sections1 or section_name not in sections2:
            print(f"\n[{section_name}] not found in both files, skipping...")
            continue
        
        # Parse section data
        if section_name == 'angles':
            data1 = parse_angles(sections1[section_name])
            data2 = parse_angles(sections2[section_name])
        elif section_name == 'constraints':
            data1 = parse_constraints(sections1[section_name])
            data2 = parse_constraints(sections2[section_name])
        elif section_name == 'dihedrals':
            data1 = parse_dihedrals(sections1[section_name])
            data2 = parse_dihedrals(sections2[section_name])
        elif section_name == 'virtual_sitesn':
            data1 = parse_virtual_sitesn(sections1[section_name])
            data2 = parse_virtual_sitesn(sections2[section_name])
        elif section_name == 'exclusions':
            data1 = parse_exclusions(sections1[section_name])
            data2 = parse_exclusions(sections2[section_name])
        else:
            continue
        
        # Compare
        result = compare_section_detail(section_name, data1, data2)
        print_detailed_report(result, show_examples)
    
    print("\n" + "=" * 80)
    print("COMPARISON COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
