#!/usr/bin/env python3
"""
Detailed comparison of two ITP files for TRAAK.
Analyzes contacts, multi_contacts, exclusions, and their patterns.
"""

from pathlib import Path
from collections import defaultdict, Counter
import argparse


def parse_section(filename: str, section_name: str):
    """Parse a section from ITP file."""
    lines = []
    in_section = False
    with open(filename) as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped == f'[ {section_name} ]':
                in_section = True
                continue
            if in_section and line_stripped.startswith('['):
                break
            if in_section and line_stripped and not line_stripped.startswith(';'):
                lines.append(line_stripped)
    return lines


def parse_contacts(lines):
    """Parse contacts lines into (i, j, func, dist, epsilon) tuples."""
    contacts = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 5:
            try:
                contacts.append((
                    int(parts[0]),   # atom i
                    int(parts[1]),   # atom j
                    int(parts[2]),   # function
                    float(parts[3]), # distance
                    float(parts[4])  # epsilon
                ))
            except (ValueError, IndexError):
                pass
    return contacts


def parse_multi_contacts(lines):
    """Parse multi_contacts lines into (i, j, state, dist, epsilon) tuples."""
    contacts = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 5:
            try:
                contacts.append((
                    int(parts[0]),   # atom i
                    int(parts[1]),   # atom j
                    int(parts[2]),   # state
                    float(parts[3]), # distance
                    float(parts[4])  # epsilon
                ))
            except (ValueError, IndexError):
                pass
    return contacts


def analyze_contacts(contacts, name):
    """Analyze contact patterns."""
    print(f"\n{'='*60}")
    print(f"ANALYSIS: {name}")
    print('='*60)
    print(f"Total contacts: {len(contacts)}")
    
    if not contacts:
        return
    
    # Atom index range
    all_atoms = set()
    for c in contacts:
        all_atoms.add(c[0])
        all_atoms.add(c[1])
    
    min_atom = min(all_atoms)
    max_atom = max(all_atoms)
    print(f"Atom index range: {min_atom} - {max_atom}")
    
    # Chain analysis (assuming chain A: 1-1727, chain B: 1728+)
    # Actually need to figure out the split point
    # For TRAAK with 2 chains of 259 residues each
    # and ~6.67 atoms per residue on average
    # Chain A: ~1-1727, Chain B: ~1728-3454
    
    chain_split = 1728  # Approximate
    
    intra_a = sum(1 for c in contacts if c[0] < chain_split and c[1] < chain_split)
    intra_b = sum(1 for c in contacts if c[0] >= chain_split and c[1] >= chain_split)
    cross_ab = sum(1 for c in contacts 
                   if (c[0] < chain_split <= c[1]) or (c[1] < chain_split <= c[0]))
    
    print(f"Intra-chain A (atoms < {chain_split}): {intra_a}")
    print(f"Intra-chain B (atoms >= {chain_split}): {intra_b}")
    print(f"Cross-chain A-B: {cross_ab}")
    
    # Distance distribution
    distances = [c[3] for c in contacts]
    print(f"Distance range: {min(distances):.4f} - {max(distances):.4f} nm")
    print(f"Mean distance: {sum(distances)/len(distances):.4f} nm")
    
    # Show first few atom pairs
    print("\nFirst 10 contacts:")
    for c in contacts[:10]:
        chain_i = 'A' if c[0] < chain_split else 'B'
        chain_j = 'A' if c[1] < chain_split else 'B'
        print(f"  {c[0]} ({chain_i}) - {c[1]} ({chain_j}): {c[3]:.4f} nm")
    
    # Show cross-chain contacts if any
    if cross_ab > 0:
        print(f"\nCross-chain contacts (first 10):")
        cross_contacts = [c for c in contacts 
                         if (c[0] < chain_split <= c[1]) or (c[1] < chain_split <= c[0])]
        for c in cross_contacts[:10]:
            chain_i = 'A' if c[0] < chain_split else 'B'
            chain_j = 'A' if c[1] < chain_split else 'B'
            print(f"  {c[0]} ({chain_i}) - {c[1]} ({chain_j}): {c[3]:.4f} nm")


def compare_atom_pairs(contacts1, contacts2, name):
    """Compare atom pairs between two contact sets."""
    pairs1 = set((c[0], c[1]) for c in contacts1)
    pairs2 = set((c[0], c[1]) for c in contacts2)
    
    common = pairs1 & pairs2
    only1 = pairs1 - pairs2
    only2 = pairs2 - pairs1
    
    print(f"\n{'='*60}")
    print(f"ATOM PAIR COMPARISON: {name}")
    print('='*60)
    print(f"File 1 unique pairs: {len(pairs1)}")
    print(f"File 2 unique pairs: {len(pairs2)}")
    print(f"Common pairs: {len(common)}")
    print(f"Only in file 1: {len(only1)}")
    print(f"Only in file 2: {len(only2)}")
    
    if only1 and len(only1) <= 20:
        print(f"\nPairs only in file 1:")
        for p in sorted(only1)[:20]:
            print(f"  {p[0]} - {p[1]}")
    
    if only2 and len(only2) <= 20:
        print(f"\nPairs only in file 2:")
        for p in sorted(only2)[:20]:
            print(f"  {p[0]} - {p[1]}")


def compare_distances(contacts1, contacts2, tolerance=1e-4):
    """Compare distances for common atom pairs."""
    dict1 = {(c[0], c[1]): c[3] for c in contacts1}
    dict2 = {(c[0], c[1]): c[3] for c in contacts2}
    
    common_pairs = set(dict1.keys()) & set(dict2.keys())
    
    matching = 0
    different = []
    
    for pair in common_pairs:
        d1 = dict1[pair]
        d2 = dict2[pair]
        if abs(d1 - d2) <= tolerance:
            matching += 1
        else:
            different.append((pair, d1, d2, abs(d1-d2)))
    
    print(f"\n{'='*60}")
    print(f"DISTANCE COMPARISON (tolerance={tolerance})")
    print('='*60)
    print(f"Common pairs: {len(common_pairs)}")
    print(f"Matching distances: {matching}")
    print(f"Different distances: {len(different)}")
    
    if different:
        # Sort by difference magnitude
        different.sort(key=lambda x: x[3], reverse=True)
        print("\nTop 10 largest distance differences:")
        for pair, d1, d2, diff in different[:10]:
            print(f"  {pair[0]} - {pair[1]}: {d1:.6f} vs {d2:.6f} (diff={diff:.6f})")


def main():
    parser = argparse.ArgumentParser(description='Detailed ITP comparison for TRAAK')
    parser.add_argument('file1', help='First ITP file (e.g., worked/TRAAK.itp)')
    parser.add_argument('file2', help='Second ITP file (e.g., prepared/TRAAK.itp)')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-4,
                        help='Tolerance for distance comparison (default: 1e-4)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("DETAILED ITP COMPARISON")
    print("="*60)
    print(f"File 1: {args.file1}")
    print(f"File 2: {args.file2}")
    
    # Parse contacts
    c1_lines = parse_section(args.file1, 'contacts')
    c2_lines = parse_section(args.file2, 'contacts')
    contacts1 = parse_contacts(c1_lines)
    contacts2 = parse_contacts(c2_lines)
    
    analyze_contacts(contacts1, f"Contacts from {args.file1}")
    analyze_contacts(contacts2, f"Contacts from {args.file2}")
    compare_atom_pairs(contacts1, contacts2, "contacts")
    compare_distances(contacts1, contacts2, args.tolerance)
    
    # Parse multi_contacts
    mc1_lines = parse_section(args.file1, 'multi_contacts')
    mc2_lines = parse_section(args.file2, 'multi_contacts')
    multi1 = parse_multi_contacts(mc1_lines)
    multi2 = parse_multi_contacts(mc2_lines)
    
    analyze_contacts(multi1, f"Multi_contacts from {args.file1}")
    analyze_contacts(multi2, f"Multi_contacts from {args.file2}")
    compare_atom_pairs(multi1, multi2, "multi_contacts")
    compare_distances(multi1, multi2, args.tolerance)
    
    # Parse exclusions
    e1 = parse_section(args.file1, 'exclusions')
    e2 = parse_section(args.file2, 'exclusions')
    print(f"\n{'='*60}")
    print("EXCLUSIONS COMPARISON")
    print('='*60)
    print(f"File 1: {len(e1)} exclusion lines")
    print(f"File 2: {len(e2)} exclusion lines")
    
    # Parse as sets of atom pairs
    excl1 = set()
    excl2 = set()
    for line in e1:
        parts = line.split()
        if len(parts) >= 2:
            try:
                excl1.add((int(parts[0]), int(parts[1])))
            except:
                pass
    for line in e2:
        parts = line.split()
        if len(parts) >= 2:
            try:
                excl2.add((int(parts[0]), int(parts[1])))
            except:
                pass
    
    print(f"File 1: {len(excl1)} exclusion pairs")
    print(f"File 2: {len(excl2)} exclusion pairs")
    print(f"Common: {len(excl1 & excl2)}")
    print(f"Only in file 1: {len(excl1 - excl2)}")
    print(f"Only in file 2: {len(excl2 - excl1)}")


if __name__ == '__main__':
    main()
