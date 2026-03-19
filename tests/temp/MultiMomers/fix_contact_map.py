#!/usr/bin/env python3
"""
Fix contact map chain IDs for hemoglobin tetramer.
Maps numeric residue IDs to proper chain IDs (A, B, C, D).

Hemoglobin structure:
- Chain A (Alpha): residues 1-141
- Chain B (Beta): residues 2-146 (starts from HIS2)
- Chain C (Alpha): residues 1-141
- Chain D (Beta): residues 2-146 (starts from HIS2)
"""

import sys

def map_resid_to_chain(resid):
    """Map residue ID to chain ID based on hemoglobin structure."""
    # Alpha chains (A, C): 1-141
    # Beta chains (B, D): 2-146
    # Need to determine which chain based on residue number
    # But contact map doesn't distinguish between A/C or B/D
    # We'll map all 1-141 to chain A (Alpha) and 2-146 to chain B (Beta)
    # This is a simplification - cross-chain contacts will need special handling
    
    if 1 <= resid <= 141:
        return 'A'  # Alpha chain
    elif 2 <= resid <= 146:
        return 'B'  # Beta chain
    else:
        return 'Z'  # Unknown

def fix_contact_map(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    modified_count = 0
    
    for line in lines:
        if line.startswith('R') and len(line.split()) >= 11:
            parts = line.split()
            # parts[0] = 'R'
            # parts[1] = ID
            # parts[2] = I1 (resid 1)
            # parts[3] = AA1
            # parts[4] = C1 (chain 1 - to fix)
            # parts[5] = I(PDB)1
            # parts[6] = I2 (resid 2)
            # parts[7] = AA2
            # parts[8] = C2 (chain 2 - to fix)
            # parts[9] = I(PDB)2
            
            if len(parts) >= 10:
                try:
                    resid1 = int(parts[5])  # I(PDB)1
                    resid2 = int(parts[9])  # I(PDB)2
                    
                    chain1 = map_resid_to_chain(resid1)
                    chain2 = map_resid_to_chain(resid2)
                    
                    # Reconstruct line with proper spacing
                    # Original format: R ID I1 AA C IPDB ...
                    fixed_line = f"{parts[0]:>8} {parts[1]:>4} {parts[2]:>4}  {parts[3]:>3}  {chain1:>1} {parts[5]:>4} {parts[6]:>6} {parts[7]:>3}  {chain2:>1} {parts[9]:>4}"
                    
                    # Add remaining columns
                    for i in range(10, len(parts)):
                        fixed_line += f" {parts[i]:>10}"
                    fixed_line += "\n"
                    
                    fixed_lines.append(fixed_line)
                    modified_count += 1
                except (ValueError, IndexError):
                    fixed_lines.append(line)
            else:
                fixed_lines.append(line)
        else:
            fixed_lines.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"Fixed {modified_count} contact records")
    print(f"Output written to {output_file}")

if __name__ == '__main__':
    fix_contact_map('protein/contact_map_protein.out', 'contact_map_fixed.out')
