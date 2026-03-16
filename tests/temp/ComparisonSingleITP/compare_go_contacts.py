#!/usr/bin/env python3
"""
Compare Go-Martini contact files (pair format vs nonbond_params format)

This script compares contact definitions between:
- Old style: go-table_VirtGoSites.itp (pairs format)
- New style: go_nbparams.itp (nonbond_params format)

Usage:
    python compare_go_contacts.py <go-table.itp> <go_nbparams.itp>
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict


def parse_go_table(filepath: str) -> dict:
    """
    Parse go-table_VirtGoSites.itp format (pairs section style)
    Format: Open_i Open_j funct sigma epsilon ; i j distance
    """
    contacts = {}
    pattern = re.compile(
        r'Open_(\d+)\s+Open_(\d+)\s+(\d)\s+([\d.]+)\s+([\d.]+)\s*;\s*(\d+)\s+(\d+)\s+([\d.]+)'
    )
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            
            match = pattern.match(line)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                funct = match.group(3)
                sigma = float(match.group(4))
                epsilon = float(match.group(5))
                comment_id1 = int(match.group(6))
                comment_id2 = int(match.group(7))
                distance = float(match.group(8))
                
                # Verify IDs match comments
                if id1 != comment_id1 or id2 != comment_id2:
                    print(f"  Warning [line {line_num}]: ID mismatch - atom {id1}-{id2} vs comment {comment_id1}-{comment_id2}")
                
                key = tuple(sorted([id1, id2]))
                contacts[key] = {
                    'funct': funct,
                    'sigma': sigma,
                    'epsilon': epsilon,
                    'distance': distance
                }
    
    return contacts


def parse_go_nbparams(filepath: str) -> dict:
    """
    Parse go_nbparams.itp format (nonbond_params section style)
    Format: Open_i Open_j funct sigma epsilon ;go bond distance
    """
    contacts = {}
    pattern = re.compile(
        r'Open_(\d+)\s+Open_(\d+)\s+(\d)\s+([\d.]+)\s+([\d.]+)\s*;\s*go bond\s+([\d.]+)'
    )
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('[') or line.startswith(';'):
                continue
            
            match = pattern.match(line)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                funct = match.group(3)
                sigma = float(match.group(4))
                epsilon = float(match.group(5))
                distance = float(match.group(6))
                
                key = tuple(sorted([id1, id2]))
                contacts[key] = {
                    'funct': funct,
                    'sigma': sigma,
                    'epsilon': epsilon,
                    'distance': distance
                }
    
    return contacts


def compare_contacts(contacts1: dict, contacts2: dict, file1_name: str, file2_name: str, tolerance: float = 0.001) -> dict:
    """Compare two sets of contacts"""
    
    keys1 = set(contacts1.keys())
    keys2 = set(contacts2.keys())
    
    common = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    
    # Check for parameter differences in common contacts
    param_diffs = []
    for key in sorted(common):
        c1 = contacts1[key]
        c2 = contacts2[key]
        
        diff_fields = {}
        if abs(c1['sigma'] - c2['sigma']) > tolerance:
            diff_fields['sigma'] = (c1['sigma'], c2['sigma'], abs(c1['sigma'] - c2['sigma']))
        if abs(c1['epsilon'] - c2['epsilon']) > tolerance:
            diff_fields['epsilon'] = (c1['epsilon'], c2['epsilon'], abs(c1['epsilon'] - c2['epsilon']))
        if abs(c1['distance'] - c2['distance']) > tolerance:
            diff_fields['distance'] = (c1['distance'], c2['distance'], abs(c1['distance'] - c2['distance']))
        
        if diff_fields:
            param_diffs.append({
                'key': key,
                'diffs': diff_fields
            })
    
    return {
        'total1': len(contacts1),
        'total2': len(contacts2),
        'common': len(common),
        'only_in_1': sorted(only_in_1),
        'only_in_2': sorted(only_in_2),
        'param_diffs': param_diffs
    }


def print_report(result: dict, file1_name: str, file2_name: str, show_all: bool = False):
    """Print detailed comparison report"""
    
    print("\n" + "=" * 80)
    print("GO-MARTINI CONTACTS COMPARISON")
    print("=" * 80)
    print(f"\nFile 1 (Old): {file1_name}")
    print(f"File 2 (New): {file2_name}")
    
    print("\n" + "-" * 80)
    print("SUMMARY")
    print("-" * 80)
    print(f"  Contacts in File 1: {result['total1']}")
    print(f"  Contacts in File 2: {result['total2']}")
    print(f"  Common contacts:    {result['common']}")
    print(f"  Only in File 1:     {len(result['only_in_1'])}")
    print(f"  Only in File 2:     {len(result['only_in_2'])}")
    print(f"  Parameter diffs:    {len(result['param_diffs'])}")
    
    match_percent = (result['common'] / max(result['total1'], result['total2'])) * 100
    print(f"  Match percentage:   {match_percent:.1f}%")
    
    # Show contacts only in File 1
    if result['only_in_1']:
        print("\n" + "-" * 80)
        print(f"CONTACTS ONLY IN FILE 1 ({len(result['only_in_1'])})")
        print("-" * 80)
        display_limit = 20 if not show_all else len(result['only_in_1'])
        for key in result['only_in_1'][:display_limit]:
            print(f"  Open_{key[0]} - Open_{key[1]}")
        if len(result['only_in_1']) > display_limit:
            print(f"  ... and {len(result['only_in_1']) - display_limit} more")
    
    # Show contacts only in File 2
    if result['only_in_2']:
        print("\n" + "-" * 80)
        print(f"CONTACTS ONLY IN FILE 2 ({len(result['only_in_2'])})")
        print("-" * 80)
        display_limit = 20 if not show_all else len(result['only_in_2'])
        for key in result['only_in_2'][:display_limit]:
            print(f"  Open_{key[0]} - Open_{key[1]}")
        if len(result['only_in_2']) > display_limit:
            print(f"  ... and {len(result['only_in_2']) - display_limit} more")
    
    # Show parameter differences
    if result['param_diffs']:
        print("\n" + "-" * 80)
        print(f"PARAMETER DIFFERENCES ({len(result['param_diffs'])})")
        print("-" * 80)
        print(f"{'Contact':<15} {'Parameter':<12} {'File 1':>15} {'File 2':>15} {'Diff':>12}")
        print("-" * 80)
        
        display_limit = 20 if not show_all else len(result['param_diffs'])
        for diff in result['param_diffs'][:display_limit]:
            key = diff['key']
            contact_str = f"Open_{key[0]}-{key[1]}"
            
            for param, (v1, v2, delta) in diff['diffs'].items():
                print(f"{contact_str:<15} {param:<12} {v1:>15.6f} {v2:>15.6f} {delta:>12.6f}")
        
        if len(result['param_diffs']) > display_limit:
            print(f"\n  ... and {len(result['param_diffs']) - display_limit} more")
    
    # Final status
    print("\n" + "=" * 80)
    if not (result['only_in_1'] or result['only_in_2'] or result['param_diffs']):
        print("STATUS: CONTACTS ARE IDENTICAL")
    else:
        print("STATUS: DIFFERENCES FOUND")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description='Compare Go-Martini contacts between old and new formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s Old/Open_go-table_VirtGoSites.itp New/go_nbparams.itp
  %(prog)s go-table.itp go_nbparams.itp --all
        """
    )
    
    parser.add_argument('file1', help='Old go-table_VirtGoSites.itp file')
    parser.add_argument('file2', help='New go_nbparams.itp file')
    parser.add_argument('--all', action='store_true',
                       help='Show all differences (not just first 20)')
    parser.add_argument('--tolerance', type=float, default=0.001,
                       help='Tolerance for parameter comparison (default: 0.001)')
    parser.add_argument('-o', '--output', default=None,
                       help='Save detailed results to JSON file')
    
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
    contacts1 = parse_go_table(str(file1_path))
    print(f"  Found {len(contacts1)} contacts")
    
    print(f"Parsing {file2_path.name}...")
    contacts2 = parse_go_nbparams(str(file2_path))
    print(f"  Found {len(contacts2)} contacts")
    
    if not contacts1 or not contacts2:
        print("Error: No contacts found in one or both files")
        sys.exit(1)
    
    # Compare
    print("\nComparing contacts...")
    result = compare_contacts(contacts1, contacts2, file1_path.name, file2_path.name, args.tolerance)
    
    # Print report
    print_report(result, file1_path.name, file2_path.name, args.all)
    
    # Save to file if requested
    if args.output:
        import json
        # Convert tuple keys to strings for JSON
        json_result = {}
        for key, value in result.items():
            if key in ['only_in_1', 'only_in_2']:
                json_result[key] = [f"{k[0]}-{k[1]}" for k in value]
            elif key == 'param_diffs':
                json_result[key] = [
                    {'key': f"{d['key'][0]}-{d['key'][1]}", 'diffs': d['diffs']}
                    for d in value
                ]
            else:
                json_result[key] = value
        
        with open(args.output, 'w') as f:
            json.dump(json_result, f, indent=2)
        print(f"\nResults saved to: {args.output}")
    
    return 0


if __name__ == "__main__":
    main()
