#!/usr/bin/env python3
"""
Compare two Go model non-bonded parameter files.
"""

from pathlib import Path
from collections import defaultdict


def normalize_key(atom1, atom2):
    """Normalize key by sorting atom names."""
    return tuple(sorted([atom1, atom2]))


def parse_go_params(filename):
    """Parse Go parameter file (handles both formats)."""
    pairs = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            # Skip section headers
            if line.startswith('['):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    atom1 = parts[0].strip()
                    atom2 = parts[1].strip()
                    func = int(parts[2])
                    dist = float(parts[3])
                    epsilon = float(parts[4])
                    
                    # Use normalized key (sorted)
                    key = normalize_key(atom1, atom2)
                    pairs[key] = {
                        'atom1': atom1,
                        'atom2': atom2,
                        'func': func,
                        'dist': dist,
                        'epsilon': epsilon,
                        'raw_line': line
                    }
                except (ValueError, IndexError):
                    pass
    return pairs


def extract_atom_number(atom_name):
    """Extract number from atom name like 'Up_5' -> 5."""
    try:
        return int(atom_name.split('_')[1])
    except (IndexError, ValueError):
        return None


def compare_files(file1, file2, tolerance=1e-4):
    """Compare two Go parameter files."""
    pairs1 = parse_go_params(file1)
    pairs2 = parse_go_params(file2)
    
    print("="*80)
    print("GO MODEL NON-BONDED PARAMETERS COMPARISON")
    print("="*80)
    print(f"File 1: {file1} ({len(pairs1)} pairs)")
    print(f"File 2: {file2} ({len(pairs2)} pairs)")
    print()
    
    # Find common and unique pairs
    keys1 = set(pairs1.keys())
    keys2 = set(pairs2.keys())
    
    common = keys1 & keys2
    only1 = keys1 - keys2
    only2 = keys2 - keys1
    
    print("-"*40)
    print("PAIR STATISTICS")
    print("-"*40)
    print(f"Common pairs:      {len(common)}")
    print(f"Only in file 1:    {len(only1)}")
    print(f"Only in file 2:    {len(only2)}")
    print()
    
    # Compare common pairs
    matching = 0
    dist_diffs = []
    epsilon_diffs = []
    
    for key in common:
        p1 = pairs1[key]
        p2 = pairs2[key]
        
        dist_match = abs(p1['dist'] - p2['dist']) <= tolerance
        eps_match = abs(p1['epsilon'] - p2['epsilon']) <= tolerance
        
        if dist_match and eps_match:
            matching += 1
        else:
            if not dist_match:
                dist_diffs.append((key, p1['dist'], p2['dist'], abs(p1['dist']-p2['dist'])))
            if not eps_match:
                epsilon_diffs.append((key, p1['epsilon'], p2['epsilon'], abs(p1['epsilon']-p2['epsilon'])))
    
    print("-"*40)
    print("PARAMETER COMPARISON (tolerance={})".format(tolerance))
    print("-"*40)
    print(f"Matching parameters: {matching}/{len(common)}")
    print(f"Distance mismatches: {len(dist_diffs)}")
    print(f"Epsilon mismatches:  {len(epsilon_diffs)}")
    print()
    
    if dist_diffs:
        print("Top 10 distance differences:")
        dist_diffs.sort(key=lambda x: x[3], reverse=True)
        for key, d1, d2, diff in dist_diffs[:10]:
            print(f"  {key[0]:10s} - {key[1]:10s}: {d1:.8f} vs {d2:.8f} (diff={diff:.8f})")
        print()
    
    if epsilon_diffs:
        print("Top 10 epsilon differences:")
        epsilon_diffs.sort(key=lambda x: x[3], reverse=True)
        for key, e1, e2, diff in epsilon_diffs[:10]:
            print(f"  {key[0]:10s} - {key[1]:10s}: {e1:.8f} vs {e2:.8f} (diff={diff:.8f})")
        print()
    
    # Analyze pairs only in file 2 (prepared)
    if only2:
        print("-"*40)
        print(f"PAIRS ONLY IN FILE 2 (showing first 20 of {len(only2)})")
        print("-"*40)
        
        # Sort by atom numbers
        def sort_key(item):
            n1 = extract_atom_number(item[0])
            n2 = extract_atom_number(item[1])
            return (n1 or 0, n2 or 0)
        
        sorted_only2 = sorted(only2, key=sort_key)
        
        print(f"{'Atom 1':<12} {'Atom 2':<12} {'Dist (nm)':<14} {'|i-j|':<10} {'Type':<15}")
        print("-"*60)
        for key in sorted_only2[:20]:
            p = pairs2[key]
            n1 = extract_atom_number(key[0])
            n2 = extract_atom_number(key[1])
            diff = abs(n1 - n2) if n1 and n2 else 0
            
            if diff < 10:
                range_type = "short-range"
            elif diff < 50:
                range_type = "medium-range"
            else:
                range_type = "long-range"
            
            print(f"{key[0]:<12} {key[1]:<12} {p['dist']:<14.8f} {diff:<10} {range_type:<15}")
        
        # Analyze range distribution
        atom_nums = []
        range_counts = {'short': 0, 'medium': 0, 'long': 0}
        for key in only2:
            n1 = extract_atom_number(key[0])
            n2 = extract_atom_number(key[1])
            if n1 and n2:
                atom_nums.extend([n1, n2])
                diff = abs(n1 - n2)
                if diff < 10:
                    range_counts['short'] += 1
                elif diff < 50:
                    range_counts['medium'] += 1
                else:
                    range_counts['long'] += 1
        
        if atom_nums:
            print()
            print(f"Atom number range: {min(atom_nums)} - {max(atom_nums)}")
            print(f"Range distribution for file 2 unique pairs:")
            print(f"  Short-range  (|i-j| < 10):  {range_counts['short']}")
            print(f"  Medium-range (10-49):       {range_counts['medium']}")
            print(f"  Long-range   (|i-j| >= 50): {range_counts['long']}")
        
        # Show some long-range examples
        long_range = [k for k in only2 if extract_atom_number(k[0]) and extract_atom_number(k[1]) 
                      and abs(extract_atom_number(k[0]) - extract_atom_number(k[1])) >= 200]
        if long_range:
            long_range.sort(key=lambda k: abs(extract_atom_number(k[0]) - extract_atom_number(k[1])), reverse=True)
            print()
            print("Top 10 longest range pairs (only in file 2):")
            print(f"{'Atom 1':<12} {'Atom 2':<12} {'Dist (nm)':<14} {'|i-j|':<10}")
            print("-"*48)
            for key in long_range[:10]:
                p = pairs2[key]
                n1 = extract_atom_number(key[0])
                n2 = extract_atom_number(key[1])
                diff = abs(n1 - n2)
                print(f"{key[0]:<12} {key[1]:<12} {p['dist']:<14.8f} {diff:<10}")
        
        if len(only2) > 20:
            print(f"\n... and {len(only2)-20} more pairs in file 2")
        print()
    
    # Analyze pairs only in file 1 (worked)
    if only1:
        print("-"*40)
        print(f"PAIRS ONLY IN FILE 1 (showing first 15 of {len(only1)})")
        print("-"*40)
        
        def sort_key(item):
            n1 = extract_atom_number(item[0])
            n2 = extract_atom_number(item[1])
            return (n1 or 0, n2 or 0)
        
        sorted_only1 = sorted(only1, key=sort_key)
        
        print(f"{'Atom 1':<12} {'Atom 2':<12} {'Dist (nm)':<14} {'|i-j|':<10}")
        print("-"*48)
        for key in sorted_only1[:15]:
            p = pairs1[key]
            n1 = extract_atom_number(key[0])
            n2 = extract_atom_number(key[1])
            diff = abs(n1 - n2) if n1 and n2 else 0
            print(f"{key[0]:<12} {key[1]:<12} {p['dist']:<14.8f} {diff:<10}")
        
        # Range analysis for only1
        range_counts = {'short': 0, 'medium': 0, 'long': 0}
        for key in only1:
            n1 = extract_atom_number(key[0])
            n2 = extract_atom_number(key[1])
            if n1 and n2:
                diff = abs(n1 - n2)
                if diff < 10:
                    range_counts['short'] += 1
                elif diff < 50:
                    range_counts['medium'] += 1
                else:
                    range_counts['long'] += 1
        
        print()
        print(f"Range distribution for file 1 unique pairs:")
        print(f"  Short-range  (|i-j| < 10):  {range_counts['short']}")
        print(f"  Medium-range (10-49):       {range_counts['medium']}")
        print(f"  Long-range   (|i-j| >= 50): {range_counts['long']}")
        print()


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Compare Go model parameter files')
    parser.add_argument('file1', help='First file (e.g., worked/Up/go_nbparams.itp)')
    parser.add_argument('file2', help='Second file (e.g., prepared/Up/Up_go-table_VirtGoSites.itp)')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-4,
                        help='Tolerance for floating point comparison (default: 1e-4)')
    
    args = parser.parse_args()
    
    compare_files(args.file1, args.file2, args.tolerance)


if __name__ == '__main__':
    main()
