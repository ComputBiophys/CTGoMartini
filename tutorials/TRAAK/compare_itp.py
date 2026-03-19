#!/usr/bin/env python3
"""
Compare two ITP files with tolerance for floating point values.
Allows 4 decimal places (0.0001) of tolerance by default.
"""

import re
from typing import List, Tuple, Optional
from pathlib import Path


class ITPComparator:
    """Compare ITP files with floating point tolerance."""
    
    def __init__(self, tolerance: float = 1e-4):
        self.tolerance = tolerance
        self.differences: List[dict] = []
    
    def parse_line(self, line: str) -> Tuple[str, List[str]]:
        """
        Parse a line into prefix and value parts.
        Returns (prefix, values) where prefix is the non-numeric start.
        """
        line = line.strip()
        if not line or line.startswith(';'):
            return (line, [])
        
        # Split by whitespace
        parts = line.split()
        return ("", parts)
    
    def is_numeric(self, s: str) -> bool:
        """Check if a string represents a number."""
        try:
            float(s)
            return True
        except ValueError:
            return False
    
    def values_equal(self, val1: str, val2: str) -> bool:
        """Compare two values with tolerance for numeric values."""
        # Direct string match first
        if val1 == val2:
            return True
        
        # Try numeric comparison
        try:
            num1 = float(val1)
            num2 = float(val2)
            return abs(num1 - num2) <= self.tolerance
        except ValueError:
            # Not numeric, string comparison failed above
            return False
    
    def compare_lines(self, line1: str, line2: str, line_num: int) -> Optional[dict]:
        """Compare two lines and return difference info if they don't match."""
        # Handle comments and empty lines
        s1 = line1.strip()
        s2 = line2.strip()
        
        # Skip pure comments and empty lines
        if (not s1 or s1.startswith(';')) and (not s2 or s2.startswith(';')):
            return None
        
        _, parts1 = self.parse_line(line1)
        _, parts2 = self.parse_line(line2)
        
        if len(parts1) != len(parts2):
            return {
                'line': line_num,
                'type': 'column_count',
                'content1': line1.rstrip(),
                'content2': line2.rstrip(),
                'details': f"Column count differs: {len(parts1)} vs {len(parts2)}"
            }
        
        # Compare each value
        diff_cols = []
        for i, (p1, p2) in enumerate(zip(parts1, parts2)):
            if not self.values_equal(p1, p2):
                diff_cols.append({
                    'col': i,
                    'val1': p1,
                    'val2': p2,
                    'num_diff': None
                })
                # Calculate numeric difference if both are numeric
                try:
                    num_diff = abs(float(p1) - float(p2))
                    diff_cols[-1]['num_diff'] = num_diff
                except ValueError:
                    pass
        
        if diff_cols:
            return {
                'line': line_num,
                'type': 'value_mismatch',
                'content1': line1.rstrip(),
                'content2': line2.rstrip(),
                'columns': diff_cols
            }
        
        return None
    
    def parse_sections(self, lines: List[str]) -> dict:
        """Parse ITP file into sections."""
        sections = {}
        current_section = None
        section_lines = []
        
        for i, line in enumerate(lines):
            stripped = line.strip()
            if stripped.startswith('[') and stripped.endswith(']'):
                # Save previous section
                if current_section:
                    sections[current_section] = section_lines
                # Start new section
                current_section = stripped[1:-1].strip()
                section_lines = []
            else:
                if current_section:
                    section_lines.append((i+1, line))
        
        # Save last section
        if current_section:
            sections[current_section] = section_lines
        
        return sections
    
    def compare_files(self, file1: str, file2: str) -> dict:
        """Compare two ITP files."""
        path1 = Path(file1)
        path2 = Path(file2)
        
        with open(path1) as f:
            lines1 = f.readlines()
        with open(path2) as f:
            lines2 = f.readlines()
        
        # Basic info
        result = {
            'file1': str(path1),
            'file2': str(path2),
            'lines1': len(lines1),
            'lines2': len(lines2),
            'tolerance': self.tolerance,
            'line_differences': [],
            'section_comparison': {}
        }
        
        # Compare line by line (up to min length)
        min_lines = min(len(lines1), len(lines2))
        for i in range(min_lines):
            diff = self.compare_lines(lines1[i], lines2[i], i+1)
            if diff:
                result['line_differences'].append(diff)
        
        # Handle extra lines
        if len(lines1) > len(lines2):
            for i in range(len(lines2), len(lines1)):
                result['line_differences'].append({
                    'line': i+1,
                    'type': 'extra_in_file1',
                    'content1': lines1[i].rstrip(),
                    'content2': None
                })
        elif len(lines2) > len(lines1):
            for i in range(len(lines1), len(lines2)):
                result['line_differences'].append({
                    'line': i+1,
                    'type': 'extra_in_file2',
                    'content1': None,
                    'content2': lines2[i].rstrip()
                })
        
        # Compare sections
        sections1 = self.parse_sections(lines1)
        sections2 = self.parse_sections(lines2)
        
        all_sections = set(sections1.keys()) | set(sections2.keys())
        for section in sorted(all_sections):
            if section not in sections1:
                result['section_comparison'][section] = 'only_in_file2'
            elif section not in sections2:
                result['section_comparison'][section] = 'only_in_file1'
            else:
                sec_lines1 = len(sections1[section])
                sec_lines2 = len(sections2[section])
                if sec_lines1 != sec_lines2:
                    result['section_comparison'][section] = f'different_length ({sec_lines1} vs {sec_lines2})'
                else:
                    result['section_comparison'][section] = 'same_length'
        
        return result


def print_comparison(result: dict, max_diffs: int = 50):
    """Print comparison results in a readable format."""
    print("=" * 80)
    print(f"ITP FILE COMPARISON")
    print("=" * 80)
    print(f"File 1: {result['file1']} ({result['lines1']} lines)")
    print(f"File 2: {result['file2']} ({result['lines2']} lines)")
    print(f"Tolerance: {result['tolerance']}")
    print()
    
    # Section comparison
    print("-" * 40)
    print("SECTION COMPARISON")
    print("-" * 40)
    for section, status in result['section_comparison'].items():
        print(f"  [{section}]: {status}")
    print()
    
    # Line differences
    diffs = result['line_differences']
    if not diffs:
        print("-" * 40)
        print("RESULT: Files are identical (within tolerance)")
        print("-" * 40)
        return
    
    print("-" * 40)
    print(f"LINE DIFFERENCES (showing first {min(max_diffs, len(diffs))} of {len(diffs)})")
    print("-" * 40)
    
    for i, diff in enumerate(diffs[:max_diffs]):
        print(f"\nLine {diff['line']}:")
        
        if diff['type'] == 'column_count':
            print(f"  Type: Column count mismatch")
            print(f"  File 1: {diff['content1'][:100]}")
            print(f"  File 2: {diff['content2'][:100]}")
            print(f"  {diff['details']}")
            
        elif diff['type'] == 'value_mismatch':
            print(f"  Type: Value mismatch")
            print(f"  File 1: {diff['content1'][:100]}")
            print(f"  File 2: {diff['content2'][:100]}")
            for col in diff['columns']:
                num_info = f" (diff={col['num_diff']:.6f})" if col['num_diff'] else ""
                print(f"    Column {col['col']}: '{col['val1']}' vs '{col['val2']}'{num_info}")
                
        elif diff['type'] == 'extra_in_file1':
            print(f"  Type: Extra line in file 1")
            print(f"  Content: {diff['content1'][:100]}")
            
        elif diff['type'] == 'extra_in_file2':
            print(f"  Type: Extra line in file 2")
            print(f"  Content: {diff['content2'][:100]}")
    
    if len(diffs) > max_diffs:
        print(f"\n... and {len(diffs) - max_diffs} more differences")
    
    print()
    print("-" * 40)
    print(f"TOTAL DIFFERENCES: {len(diffs)}")
    print("-" * 40)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Compare two ITP files with floating point tolerance')
    parser.add_argument('file1', help='First ITP file')
    parser.add_argument('file2', help='Second ITP file')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-4,
                        help='Tolerance for floating point comparison (default: 1e-4)')
    parser.add_argument('-n', '--max-diffs', type=int, default=50,
                        help='Maximum number of differences to display (default: 50)')
    
    args = parser.parse_args()
    
    comparator = ITPComparator(tolerance=args.tolerance)
    result = comparator.compare_files(args.file1, args.file2)
    print_comparison(result, max_diffs=args.max_diffs)


if __name__ == '__main__':
    main()
