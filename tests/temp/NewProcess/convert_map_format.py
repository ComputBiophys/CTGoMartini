#!/usr/bin/env python3
"""
Convert .map file format (rCSU web-server) to contact_map.out format (OV+rCSU)

Input format (rCSU web-server .map):
  - Columns: ID, I1, AA1, C1, I(PDB)1, I2, AA2, C2, I(PDB)2, DISTANCE, CMs(4 values), rCSU, aSurf, rSurf, nSurf

Output format (OV+rCSU contact_map.out):
  - Columns: ID, I1, AA1, C1, I(PDB)1, I2, AA2, C2, I(PDB)2, DCA, CMs(4 values), rCSU, Count, Model

Key differences:
  1. Different header format
  2. Input has aSurf, rSurf, nSurf columns; output does not
  3. Output has Count and Model columns; input does not
  4. CMs format may differ slightly

Usage:
    python convert_map_format.py <input.map> [output.out] [options]

Examples:
    # Basic usage (auto-generate output filename)
    python convert_map_format.py input.map
    
    # Specify output file
    python convert_map_format.py input.map output.out
    
    # Specify PDB name for header
    python convert_map_format.py input.map output.out --pdb-name myprotein.pdb
    
    # Specify number of residues
    python convert_map_format.py input.map output.out --nresidues 300
"""

import re
import sys
import argparse
from pathlib import Path


def parse_map_file(filepath: str) -> list[dict]:
    """
    Parse rCSU web-server format .map file
    
    Args:
        filepath: Path to input .map file
    
    Returns:
        List of contact dictionaries
    """
    contacts = []
    
    # Regex pattern for data lines
    # Format: R ID I1 AA C I(PDB) I2 AA C I(PDB) DISTANCE CMs rCSU aSurf rSurf nSurf
    pattern = re.compile(
        r'^\s*R\s+(\d+)\s+'           # R ID
        r'(\d+)\s+([A-Z]{3})\s+([A-Z])\s+(\d+)\s+'  # I1 AA C I(PDB)
        r'(\d+)\s+([A-Z]{3})\s+([A-Z])\s+(\d+)\s+'  # I2 AA C I(PDB)
        r'(\d+\.?\d*)\s+'            # DISTANCE
        r'(\d)\s+(\d)\s+(\d)\s+(\d)\s+'  # CMs (4 values)
        r'(-?\d+)\s+'                # rCSU
        r'(\d+\.?\d*)\s+'            # aSurf
        r'(\d+\.?\d*)\s+'            # rSurf
        r'(\d+\.?\d*)'               # nSurf
    )
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '====' in line:
                continue
            
            match = pattern.match(line)
            if match:
                contact = {
                    'id': int(match.group(1)),
                    'i1': int(match.group(2)),
                    'aa1': match.group(3),
                    'chain1': match.group(4),
                    'ipdb1': int(match.group(5)),
                    'i2': int(match.group(6)),
                    'aa2': match.group(7),
                    'chain2': match.group(8),
                    'ipdb2': int(match.group(9)),
                    'distance': float(match.group(10)),
                    'cm_ov': int(match.group(11)),
                    'cm_csu': int(match.group(12)),
                    'cm_ocsu': int(match.group(13)),
                    'cm_rcsu': int(match.group(14)),
                    'rcsu_net': int(match.group(15)),
                    'asurf': float(match.group(16)),
                    'rsurf': float(match.group(17)),
                    'nsurf': float(match.group(18)),
                }
                contacts.append(contact)
    
    return contacts


def generate_output_header(pdb_name: str = "input.pdb", nresidues: int = 220) -> str:
    """
    Generate header for output file
    
    Args:
        pdb_name: PDB filename to display in header
        nresidues: Number of residues
    
    Returns:
        Header string
    """
    header = f"""                         CONTACT MAPS FROM PDB FILES                          
                                                                              
 This software is written by:                                                
       Rodrigo Azevedo Moreira da Silva                                      
                                                                              
 Copyright (c) 2020 - IPPT-PAN                                              
       Institute of Fundamental Techonological Research                     
       Polish Academy of Sciences                                          
 MIT LICENSE, check out LICENSE for more informations.                      
                                                                              
Reading file:    uploads/xxx/{pdb_name}
pdb natoms:      1714
pdb nresidues:   {nresidues}
Memory usage:     6.88 MB
Fibonacci grid:  610
ALPHA:           1.24
WATER_RADIUS:    2.80

Residue-Residue Contacts

ID       - atom identification
I1,I2    - serial residue id
AA       - 3-letter code of aminoacid
C        - chain
I(PDB)   - residue number in PDB file
DCA      - distance between CA
CMs      - OV , CSU , oCSU , rCSU
           (CSU does not take into account chemical properties of atoms)
rCSU     - net contact from rCSU
Count    - number of contacts between residues
MODEL    - model number

      ID    I1  AA  C I(PDB)     I2  AA  C I(PDB)        DCA       CMs    rCSU   Count Model
============================================================================================
"""
    return header


def convert_to_output_format(contacts: list[dict]) -> list[str]:
    """
    Convert parsed contacts to output format
    
    Args:
        contacts: List of contact dictionaries
    
    Returns:
        List of formatted output lines
    """
    output_lines = []
    
    for contact in contacts:
        # Calculate Count based on rCSU_net
        count = abs(contact['rcsu_net'])
        if count == 0:
            count = 0
        
        # Model is always 0
        model = 0
        
        # Format output line
        line = (
            f"R {contact['id']:>6}  "
            f"{contact['i1']:>3}  {contact['aa1']:>3}  {contact['chain1']:>1} {contact['ipdb1']:>4}   "
            f"{contact['i2']:>4}  {contact['aa2']:>3}  {contact['chain2']:>1} {contact['ipdb2']:>4}   "
            f"{contact['distance']:>10.4f}   "
            f"{contact['cm_ov']:>1} {contact['cm_csu']:>1} {contact['cm_ocsu']:>1} {contact['cm_rcsu']:>1}  "
            f"{contact['rcsu_net']:>5}   "
            f"{count:>5}    {model:>1}"
        )
        output_lines.append(line)
    
    return output_lines


def main():
    parser = argparse.ArgumentParser(
        description='Convert rCSU web-server .map format to OV+rCSU contact_map.out format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.map
  %(prog)s input.map output.out
  %(prog)s input.map output.out --pdb-name protein.pdb --nresidues 250
        """
    )
    
    parser.add_argument('input', help='Input .map file from rCSU web-server')
    parser.add_argument('output', nargs='?', default=None,
                        help='Output .out file (default: input_filename.out)')
    parser.add_argument('--pdb-name', default="input.pdb",
                        help='PDB filename to display in header (default: input.pdb)')
    parser.add_argument('--nresidues', type=int, default=None,
                        help='Number of residues (default: auto-detect from max residue ID)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite output file if it exists')
    
    args = parser.parse_args()
    
    # Check input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}")
        sys.exit(1)
    
    # Determine output filename
    if args.output:
        output_path = Path(args.output)
    else:
        # Auto-generate output filename
        output_path = input_path.with_suffix('.out')
    
    # Check if output file exists
    if output_path.exists() and not args.force:
        print(f"Error: Output file already exists: {output_path}")
        print("Use -f or --force to overwrite")
        sys.exit(1)
    
    print("=" * 80)
    print("Map File Format Converter")
    print("=" * 80)
    print()
    print(f"Input file:  {input_path.absolute()}")
    print(f"Output file: {output_path.absolute()}")
    print()
    
    # Parse input file
    print("Parsing input file...")
    contacts = parse_map_file(str(input_path))
    print(f"  Parsed {len(contacts)} contacts")
    print()
    
    if not contacts:
        print("Error: No contacts found. Please check the file format.")
        sys.exit(1)
    
    # Determine number of residues
    if args.nresidues:
        max_residue = args.nresidues
    else:
        max_residue = max(c['i2'] for c in contacts)
        print(f"Auto-detected max residue ID: {max_residue}")
    
    # Generate output
    print("Generating output file...")
    header = generate_output_header(pdb_name=args.pdb_name, nresidues=max_residue)
    output_lines = convert_to_output_format(contacts)
    
    # Write output file
    with open(output_path, 'w') as f:
        f.write(header)
        for line in output_lines:
            f.write(line + '\n')
    
    print(f"  Output file generated: {output_path.absolute()}")
    print()
    
    # Display sample output
    print("=" * 80)
    print("Sample output (first 5 lines):")
    print("=" * 80)
    print(header.split('\n')[0])  # Title
    print("...")
    print("      ID    I1  AA  C I(PDB)     I2  AA  C I(PDB)        DCA       CMs    rCSU   Count Model")
    print("============================================================================================")
    for line in output_lines[:5]:
        print(line)
    if len(output_lines) > 5:
        print("...")
    print()
    
    # Statistics
    print("=" * 80)
    print("Conversion Statistics:")
    print("=" * 80)
    print(f"  Total contacts: {len(contacts)}")
    print(f"  Max residue ID: {max_residue}")
    
    # Contact type statistics
    ov_contacts = sum(1 for c in contacts if c['cm_ov'] == 1)
    rcsu_contacts = sum(1 for c in contacts if c['cm_rcsu'] == 1)
    print(f"  OV contacts: {ov_contacts}")
    print(f"  rCSU contacts: {rcsu_contacts}")
    print()
    
    print("Conversion completed successfully!")


if __name__ == "__main__":
    main()
