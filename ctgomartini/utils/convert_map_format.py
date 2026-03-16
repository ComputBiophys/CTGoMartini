#!/usr/bin/env python3
"""
Convert .map file format (rCSU web-server) to contact_map.out format (OV+rCSU).

Input format (rCSU web-server .map):
  - Columns: ID, I1, AA1, C1, I(PDB)1, I2, AA2, C2, I(PDB)2, DISTANCE, CMs(4 values), rCSU, aSurf, rSurf, nSurf

Output format (OV+rCSU contact_map.out):
  - Columns: ID, I1, AA1, C1, I(PDB)1, I2, AA2, C2, I(PDB)2, DCA, CMs(4 values), rCSU, Count, Model

Key differences:
  1. Different header format
  2. Input has aSurf, rSurf, nSurf columns; output does not
  3. Output has Count and Model columns; input does not
  4. CMs format may differ slightly

This module provides functions to convert rCSU web-server contact map files
to the format expected by martinize2 -go option (vermouth >= 0.15.0).
"""

import re
from pathlib import Path
from typing import Any


def parse_map_file(filepath: str) -> list[dict[str, Any]]:
    """
    Parse rCSU web-server format .map file.

    Args:
        filepath: Path to input .map file.

    Returns:
        List of contact dictionaries with parsed contact information.

    Raises:
        FileNotFoundError: If input file does not exist.
        ValueError: If file format is invalid.
    """
    contacts = []

    # Regex pattern for data lines
    # Format: R ID I1 AA C I(PDB) I2 AA C I(PDB) DISTANCE CMs rCSU aSurf rSurf nSurf
    pattern = re.compile(
        r'^\s*R\s+(\d+)\s+'  # R ID
        r'(\d+)\s+([A-Z]{3})\s+([A-Z])\s+(\d+)\s+'  # I1 AA C I(PDB)
        r'(\d+)\s+([A-Z]{3})\s+([A-Z])\s+(\d+)\s+'  # I2 AA C I(PDB)
        r'(\d+\.?\d*)\s+'  # DISTANCE
        r'(\d)\s+(\d)\s+(\d)\s+(\d)\s+'  # CMs (4 values)
        r'(-?\d+)\s+'  # rCSU
        r'(\d+\.?\d*)\s+'  # aSurf
        r'(\d+\.?\d*)\s+'  # rSurf
        r'(\d+\.?\d*)'  # nSurf
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
    Generate header for output file.

    Args:
        pdb_name: PDB filename to display in header.
        nresidues: Number of residues.

    Returns:
        Header string for contact_map.out file.
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


def convert_to_output_format(contacts: list[dict[str, Any]]) -> list[str]:
    """
    Convert parsed contacts to output format.

    Args:
        contacts: List of contact dictionaries.

    Returns:
        List of formatted output lines.
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


def convert_map_format(
    input_file: str,
    output_file: str | None = None,
    pdb_name: str = "input.pdb",
    nresidues: int | None = None,
    force: bool = False,
) -> str:
    """
    Convert rCSU web-server .map format to OV+rCSU contact_map.out format.

    This is the main function that combines parsing, header generation,
    and format conversion.

    Args:
        input_file: Path to input .map file from rCSU web-server.
        output_file: Path to output .out file. If None, auto-generate from input filename.
        pdb_name: PDB filename to display in header.
        nresidues: Number of residues. If None, auto-detect from max residue ID.
        force: Overwrite output file if it exists.

    Returns:
        Path to the generated output file.

    Raises:
        FileNotFoundError: If input file does not exist.
        FileExistsError: If output file exists and force is False.
        ValueError: If no contacts found in input file.

    Example:
        >>> convert_map_format('input.map', 'output.out', pdb_name='protein.pdb')
        'output.out'
    """
    # Check input file exists
    input_path = Path(input_file)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Determine output filename
    if output_file:
        output_path = Path(output_file)
    else:
        # Auto-generate output filename
        output_path = input_path.with_suffix('.out')

    # Check if output file exists
    if output_path.exists() and not force:
        raise FileExistsError(f"Output file already exists: {output_path}. Use force=True to overwrite.")

    # Parse input file
    contacts = parse_map_file(str(input_path))

    if not contacts:
        raise ValueError("No contacts found. Please check the file format.")

    # Determine number of residues
    if nresidues is None:
        nresidues = max(c['i2'] for c in contacts)

    # Generate output
    header = generate_output_header(pdb_name=pdb_name, nresidues=nresidues)
    output_lines = convert_to_output_format(contacts)

    # Write output file
    with open(output_path, 'w') as f:
        f.write(header)
        for line in output_lines:
            f.write(line + '\n')

    return str(output_path)


def main():
    """Command-line interface for map format conversion."""
    import argparse
    import sys

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
    parser.add_argument('output', nargs='?',
                        help='Output .out file (default: input_filename.out)')
    parser.add_argument('--pdb-name', default="input.pdb",
                        help='PDB filename to display in header (default: input.pdb)')
    parser.add_argument('--nresidues', type=int,
                        help='Number of residues (default: auto-detect from max residue ID)')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Overwrite output file if it exists')

    args = parser.parse_args()

    try:
        output_path = convert_map_format(
            input_file=args.input,
            output_file=args.output,
            pdb_name=args.pdb_name,
            nresidues=args.nresidues,
            force=args.force,
        )
        print(f"Conversion completed: {output_path}")
    except (FileNotFoundError, FileExistsError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
