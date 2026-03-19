#!/usr/bin/env python3
"""
Fix contact map chain IDs for hemoglobin tetramer.

Original contact map uses sequential residue numbering:
- Chain A (Alpha): residues 1-141
- Chain B (Beta): residues 142-287 (maps to 1-146)
- Chain C (Alpha): residues 288-428 (maps to 1-141)
- Chain D (Beta): residues 429-574 (maps to 2-146)

We need to convert sequential numbering to actual PDB chain IDs and residue numbers.
"""

def convert_sequential_to_chain_resid(seq_id):
    """
    Convert sequential residue ID to (chain, resid).
    
    Hemoglobin structure in contact map:
    - 1-141: Chain A (Alpha), resid = seq_id
    - 142-287: Chain B (Beta), resid = seq_id - 141
    - 288-428: Chain C (Alpha), resid = seq_id - 287
    - 429-574: Chain D (Beta), resid = seq_id - 428
    """
    if 1 <= seq_id <= 141:
        return 'A', seq_id
    elif 142 <= seq_id <= 287:
        return 'B', seq_id - 141
    elif 288 <= seq_id <= 428:
        return 'C', seq_id - 287
    elif 429 <= seq_id <= 574:
        return 'D', seq_id - 428
    else:
        return 'Z', seq_id  # Unknown

def fix_contact_map(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    contact_count = 0
    chain_distribution = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
    
    for line in lines:
        if not line.startswith('R'):
            fixed_lines.append(line)
            continue
            
        parts = line.split()
        if len(parts) < 11:
            fixed_lines.append(line)
            continue
        
        try:
            # Parse original contact map format
            # R ID I1 AA1 C1 IPDB1 I2 AA2 C2 IPDB2 DCA ...
            # Note: C1 and C2 in original are not real chain IDs
            seq_id1 = int(parts[2])  # I1 (sequential residue ID)
            seq_id2 = int(parts[6])  # I2 (sequential residue ID)
            
            # Convert to actual chain and resid
            chain1, resid1 = convert_sequential_to_chain_resid(seq_id1)
            chain2, resid2 = convert_sequential_to_chain_resid(seq_id2)
            
            chain_distribution[chain1] = chain_distribution.get(chain1, 0) + 1
            chain_distribution[chain2] = chain_distribution.get(chain2, 0) + 1
            
            # Reconstruct line with proper chain IDs
            # Format: R ID I1 AA1 C IPDB I2 AA2 C IPDB DCA ...
            fixed_line = (
                f"{parts[0]:>8} {parts[1]:>4} {parts[2]:>4}  {parts[3]:>3}  {chain1:>1} {resid1:>4}"
                f" {parts[6]:>6} {parts[7]:>3}  {chain2:>1} {resid2:>4}"
            )
            
            # Add remaining columns (DCA, CMs, etc.)
            for i in range(10, len(parts)):
                fixed_line += f" {parts[i]:>10}"
            fixed_line += "\n"
            
            fixed_lines.append(fixed_line)
            contact_count += 1
            
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse line: {line.strip()[:50]}... ({e})")
            fixed_lines.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"Fixed {contact_count} contact records")
    print(f"Chain distribution: {chain_distribution}")
    print(f"Output written to {output_file}")

if __name__ == '__main__':
    fix_contact_map('protein/contact_map_protein.out', 'contact_map_fixed.out')
