#!/usr/bin/env python3
"""
Compare GROMACS ITP files by sections/categories

This script parses two ITP files and compares their contents section by section,
reporting differences in atom counts, bond parameters, and other topology elements.

Usage:
    python compare_itp_categories.py <file1.itp> <file2.itp> [options]

Example:
    python compare_itp_categories.py Old/Open.itp New/Open.itp
    python compare_itp_categories.py Old/Open.itp New/Open.itp --details
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def parse_itp_file(filepath: str) -> dict:
    """
    Parse ITP file and extract all sections
    
    Args:
        filepath: Path to ITP file
    
    Returns:
        Dictionary with section names as keys and content as values
    """
    sections = defaultdict(list)
    current_section = None
    header_lines = []
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            stripped = line.strip()
            
            # Skip empty lines but track them within sections
            if not stripped:
                if current_section:
                    sections[current_section].append((line_num, '', 'blank'))
                continue
            
            # Check for section header
            section_match = re.match(r'^\[\s*(\w+)\s*\]', stripped)
            if section_match:
                current_section = section_match.group(1)
                sections[current_section].append((line_num, stripped, 'header'))
                continue
            
            # Check for comments
            if stripped.startswith(';'):
                if current_section:
                    sections[current_section].append((line_num, stripped, 'comment'))
                else:
                    header_lines.append((line_num, stripped))
                continue
            
            # Check for preprocessor directives
            if stripped.startswith('#'):
                if current_section:
                    sections[current_section].append((line_num, stripped, 'directive'))
                else:
                    header_lines.append((line_num, stripped))
                continue
            
            # Regular data line
            if current_section:
                sections[current_section].append((line_num, stripped, 'data'))
            else:
                header_lines.append((line_num, stripped))
    
    return {
        'header': header_lines,
        'sections': dict(sections)
    }


def compare_atoms(section1: list, section2: list, file1_name: str, file2_name: str) -> dict:
    """Compare [ atoms ] section"""
    
    # Parse atom data
    atoms1 = {}
    atoms2 = {}
    
    for line_num, line, line_type in section1:
        if line_type != 'data':
            continue
        parts = line.split()
        if len(parts) >= 7:
            atom_id = int(parts[0])
            atoms1[atom_id] = {
                'type': parts[1],
                'resnr': int(parts[2]),
                'resname': parts[3],
                'atomname': parts[4],
                'cgnr': int(parts[5]),
                'charge': parts[6]
            }
    
    for line_num, line, line_type in section2:
        if line_type != 'data':
            continue
        parts = line.split()
        if len(parts) >= 7:
            atom_id = int(parts[0])
            atoms2[atom_id] = {
                'type': parts[1],
                'resnr': int(parts[2]),
                'resname': parts[3],
                'atomname': parts[4],
                'cgnr': int(parts[5]),
                'charge': parts[6]
            }
    
    # Compare
    ids1 = set(atoms1.keys())
    ids2 = set(atoms2.keys())
    
    common = ids1 & ids2
    only_in_1 = ids1 - ids2
    only_in_2 = ids2 - ids1
    
    differences = []
    for atom_id in common:
        a1 = atoms1[atom_id]
        a2 = atoms2[atom_id]
        if a1 != a2:
            differences.append({
                'id': atom_id,
                'file1': a1,
                'file2': a2
            })
    
    return {
        'total1': len(atoms1),
        'total2': len(atoms2),
        'common': len(common),
        'only_in_1': sorted(only_in_1),
        'only_in_2': sorted(only_in_2),
        'differences': differences
    }


def compare_bonds(section1: list, section2: list) -> dict:
    """Compare [ bonds ] section"""
    
    # Parse bond data
    bonds1 = {}
    bonds2 = {}
    
    for line_num, line, line_type in section1:
        if line_type != 'data':
            continue
        parts = line.split()
        if len(parts) >= 5:
            # Bond: atom_i atom_j funct length force_constant
            i, j = int(parts[0]), int(parts[1])
            key = tuple(sorted([i, j]))
            bonds1[key] = {
                'funct': parts[2],
                'length': float(parts[3]),
                'fc': float(parts[4])
            }
    
    for line_num, line, line_type in section2:
        if line_type != 'data':
            continue
        parts = line.split()
        if len(parts) >= 5:
            i, j = int(parts[0]), int(parts[1])
            key = tuple(sorted([i, j]))
            bonds2[key] = {
                'funct': parts[2],
                'length': float(parts[3]),
                'fc': float(parts[4])
            }
    
    # Compare
    keys1 = set(bonds1.keys())
    keys2 = set(bonds2.keys())
    
    common = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    
    param_diffs = []
    for key in common:
        b1 = bonds1[key]
        b2 = bonds2[key]
        if abs(b1['length'] - b2['length']) > 0.001 or abs(b1['fc'] - b2['fc']) > 0.1:
            param_diffs.append({
                'atoms': key,
                'file1': b1,
                'file2': b2
            })
    
    return {
        'total1': len(bonds1),
        'total2': len(bonds2),
        'common': len(common),
        'only_in_1': sorted(only_in_1),
        'only_in_2': sorted(only_in_2),
        'param_diffs': param_diffs
    }


def compare_position_restraints(section1: list, section2: list) -> dict:
    """Compare [ position_restraints ] section"""
    
    # Parse position restraint data
    restraints1 = {}
    restraints2 = {}
    
    for line_num, line, line_type in section1:
        if line_type != 'data':
            continue
        parts = line.split()
        if len(parts) >= 5:
            # Format: atom_id funct kx ky kz
            atom_id = int(parts[0])
            restraints1[atom_id] = {
                'funct': parts[1],
                'kx': parts[2],
                'ky': parts[3],
                'kz': parts[4]
            }
    
    for line_num, line, line_type in section2:
        if line_type != 'data':
            continue
        parts = line.split()
        if len(parts) >= 5:
            atom_id = int(parts[0])
            restraints2[atom_id] = {
                'funct': parts[1],
                'kx': parts[2],
                'ky': parts[3],
                'kz': parts[4]
            }
    
    # Check if using variable (POSRES_FC) vs hardcoded values
    uses_variable_1 = any(
        r['kx'].startswith('POSRES') or r['ky'].startswith('POSRES') or r['kz'].startswith('POSRES')
        for r in restraints1.values()
    )
    uses_variable_2 = any(
        r['kx'].startswith('POSRES') or r['ky'].startswith('POSRES') or r['kz'].startswith('POSRES')
        for r in restraints2.values()
    )
    
    ids1 = set(restraints1.keys())
    ids2 = set(restraints2.keys())
    
    return {
        'total1': len(restraints1),
        'total2': len(restraints2),
        'common': len(ids1 & ids2),
        'only_in_1': sorted(ids1 - ids2),
        'only_in_2': sorted(ids2 - ids1),
        'uses_variable_1': uses_variable_1,
        'uses_variable_2': uses_variable_2
    }


def compare_simple_section(section1: list, section2: list, section_name: str) -> dict:
    """Compare sections by line count and content"""
    
    data_lines1 = [line for _, line, line_type in section1 if line_type == 'data']
    data_lines2 = [line for _, line, line_type in section2 if line_type == 'data']
    
    return {
        'name': section_name,
        'count1': len(data_lines1),
        'count2': len(data_lines2),
        'identical': data_lines1 == data_lines2
    }


def print_comparison_report(result: dict, file1_name: str, file2_name: str, show_details: bool = False):
    """Print detailed comparison report"""
    
    print("\n" + "=" * 80)
    print("ITP FILES COMPARISON REPORT")
    print("=" * 80)
    print(f"\nFile 1: {file1_name}")
    print(f"File 2: {file2_name}")
    
    # Section presence comparison
    sections1 = set(result['sections1'].keys())
    sections2 = set(result['sections2'].keys())
    
    print("\n" + "-" * 80)
    print("SECTIONS PRESENT")
    print("-" * 80)
    print(f"{'Section':<30} {'File 1':>10} {'File 2':>10} {'Status':>15}")
    print("-" * 80)
    
    all_sections = sorted(sections1 | sections2)
    for section in all_sections:
        in1 = "Yes" if section in sections1 else "No"
        in2 = "Yes" if section in sections2 else "No"
        if section in sections1 and section in sections2:
            status = "Present in both"
        elif section in sections1:
            status = "Only in File 1"
        else:
            status = "Only in File 2"
        print(f"{section:<30} {in1:>10} {in2:>10} {status:>15}")
    
    # Detailed section comparisons
    print("\n" + "=" * 80)
    print("DETAILED SECTION COMPARISON")
    print("=" * 80)
    
    # Compare atoms
    if 'atoms' in sections1 and 'atoms' in sections2:
        print("\n[ atoms ]")
        print("-" * 40)
        atoms_result = result.get('atoms_comparison', {})
        print(f"  Atoms in File 1: {atoms_result.get('total1', 0)}")
        print(f"  Atoms in File 2: {atoms_result.get('total2', 0)}")
        print(f"  Common atoms: {atoms_result.get('common', 0)}")
        
        only_1 = atoms_result.get('only_in_1', [])
        only_2 = atoms_result.get('only_in_2', [])
        diffs = atoms_result.get('differences', [])
        
        if only_1:
            print(f"  Only in File 1: {len(only_1)} atoms")
            if show_details and len(only_1) <= 20:
                print(f"    IDs: {only_1}")
        if only_2:
            print(f"  Only in File 2: {len(only_2)} atoms")
            if show_details and len(only_2) <= 20:
                print(f"    IDs: {only_2}")
        if diffs:
            print(f"  Parameter differences: {len(diffs)} atoms")
            if show_details:
                for d in diffs[:10]:
                    print(f"    Atom {d['id']}: {d['file1']} vs {d['file2']}")
        if not (only_1 or only_2 or diffs):
            print("  Status: IDENTICAL")
    
    # Compare bonds
    if 'bonds' in sections1 and 'bonds' in sections2:
        print("\n[ bonds ]")
        print("-" * 40)
        bonds_result = result.get('bonds_comparison', {})
        print(f"  Bonds in File 1: {bonds_result.get('total1', 0)}")
        print(f"  Bonds in File 2: {bonds_result.get('total2', 0)}")
        print(f"  Common bonds: {bonds_result.get('common', 0)}")
        
        only_1 = bonds_result.get('only_in_1', [])
        only_2 = bonds_result.get('only_in_2', [])
        param_diffs = bonds_result.get('param_diffs', [])
        
        if only_1:
            print(f"  Only in File 1: {len(only_1)} bonds")
        if only_2:
            print(f"  Only in File 2: {len(only_2)} bonds")
        if param_diffs:
            print(f"  Parameter differences: {len(param_diffs)} bonds")
            if show_details:
                for d in param_diffs[:5]:
                    print(f"    {d['atoms']}: {d['file1']} vs {d['file2']}")
        if not (only_1 or only_2 or param_diffs):
            print("  Status: IDENTICAL")
    
    # Compare position restraints
    if 'position_restraints' in sections1 and 'position_restraints' in sections2:
        print("\n[ position_restraints ]")
        print("-" * 40)
        posres_result = result.get('position_restraints_comparison', {})
        print(f"  Restraints in File 1: {posres_result.get('total1', 0)}")
        print(f"  Restraints in File 2: {posres_result.get('total2', 0)}")
        print(f"  Common restraints: {posres_result.get('common', 0)}")
        
        only_1 = posres_result.get('only_in_1', [])
        only_2 = posres_result.get('only_in_2', [])
        
        if only_1:
            print(f"  Only in File 1: {len(only_1)}")
        if only_2:
            print(f"  Only in File 2: {len(only_2)}")
        
        # Report variable usage
        var1 = posres_result.get('uses_variable_1', False)
        var2 = posres_result.get('uses_variable_2', False)
        
        if var1 and not var2:
            print(f"  File 1 uses POSRES_FC variable")
            print(f"  File 2 uses hardcoded values")
        elif var2 and not var1:
            print(f"  File 1 uses hardcoded values")
            print(f"  File 2 uses POSRES_FC variable")
        elif var1 and var2:
            print(f"  Both use POSRES_FC variable")
        else:
            print(f"  Both use hardcoded values")
        
        if not (only_1 or only_2):
            print("  Status: Same atoms restrained")
    
    # Other sections
    print("\n" + "-" * 80)
    print("OTHER SECTIONS")
    print("-" * 80)
    
    common_sections = sections1 & sections2
    for section in sorted(common_sections):
        if section in ['atoms', 'bonds', 'position_restraints']:
            continue
        
        sec1 = result['sections1'][section]
        sec2 = result['sections2'][section]
        
        data_count1 = len([l for _, l, t in sec1 if t == 'data'])
        data_count2 = len([l for _, l, t in sec2 if t == 'data'])
        
        status = "Same count" if data_count1 == data_count2 else "DIFFERENT"
        print(f"  [{section}]: {data_count1} vs {data_count2} entries - {status}")
    
    print("\n" + "=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Compare two GROMACS ITP files by sections/categories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s Old/Open.itp New/Open.itp
  %(prog)s Old/Open.itp New/Open.itp --details
        """
    )
    
    parser.add_argument('file1', help='First ITP file')
    parser.add_argument('file2', help='Second ITP file')
    parser.add_argument('--details', action='store_true',
                        help='Show detailed differences')
    parser.add_argument('-o', '--output', default=None,
                        help='Save report to file')
    
    args = parser.parse_args()
    
    # Check files exist
    for f in [args.file1, args.file2]:
        if not Path(f).exists():
            print(f"Error: File not found: {f}")
            sys.exit(1)
    
    file1_path = Path(args.file1)
    file2_path = Path(args.file2)
    
    # Parse files
    print(f"Parsing {file1_path.name}...")
    parsed1 = parse_itp_file(str(file1_path))
    
    print(f"Parsing {file2_path.name}...")
    parsed2 = parse_itp_file(str(file2_path))
    
    # Perform comparisons
    result = {
        'sections1': parsed1['sections'],
        'sections2': parsed2['sections']
    }
    
    # Compare specific sections
    if 'atoms' in parsed1['sections'] and 'atoms' in parsed2['sections']:
        result['atoms_comparison'] = compare_atoms(
            parsed1['sections']['atoms'],
            parsed2['sections']['atoms'],
            file1_path.name,
            file2_path.name
        )
    
    if 'bonds' in parsed1['sections'] and 'bonds' in parsed2['sections']:
        result['bonds_comparison'] = compare_bonds(
            parsed1['sections']['bonds'],
            parsed2['sections']['bonds']
        )
    
    if 'position_restraints' in parsed1['sections'] and 'position_restraints' in parsed2['sections']:
        result['position_restraints_comparison'] = compare_position_restraints(
            parsed1['sections']['position_restraints'],
            parsed2['sections']['position_restraints']
        )
    
    # Print report
    print_comparison_report(result, file1_path.name, file2_path.name, args.details)
    
    # Save to file if requested
    if args.output:
        import json
        # Convert sets to lists for JSON serialization
        def convert_for_json(obj):
            if isinstance(obj, set):
                return list(obj)
            elif isinstance(obj, dict):
                return {k: convert_for_json(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_for_json(i) for i in obj]
            return obj
        
        with open(args.output, 'w') as f:
            json.dump(convert_for_json(result), f, indent=2)
        print(f"Detailed report saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
