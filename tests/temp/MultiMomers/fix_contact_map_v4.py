#!/usr/bin/env python3
"""
Fix contact map chain IDs and residue IDs for Martinize2 merged hemoglobin.

Martinize2 with -merge all renumbers residues sequentially:
- Chain A: 1-141 (unchanged)
- Chain B: 142-286 (was 2-146, offset by +140)
- Chain C: 287-427 (was 1-141, offset by +286)
- Chain D: 428-572 (was 2-146, offset by +426)
"""

def convert_to_martinize_resid(chain, resid):
    """Convert original PDB resid to Martinize2 merged resid."""
    if chain == 'A':
        return resid  # 1-141
    elif chain == 'B':
        return resid + 140  # 2-146 -> 142-286
    elif chain == 'C':
        return resid + 286  # 1-141 -> 287-427
    elif chain == 'D':
        return resid + 426  # 2-146 -> 428-572
    return resid

def fix_contact_map(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    contact_count = 0
    
    for line in lines:
        stripped = line.lstrip()
        if not stripped.startswith('R '):
            fixed_lines.append(line)
            continue
            
        parts = line.split()
        if len(parts) < 11:
            fixed_lines.append(line)
            continue
        
        try:
            # Parse: R ID I1 AA C IPDB I2 AA C IPDB DCA ...
            # parts[0] = R
            # parts[1] = ID
            # parts[2] = I1
            # parts[3] = AA1
            # parts[4] = C1 (chain)
            # parts[5] = IPDB1 (resid)
            # parts[6] = I2
            # parts[7] = AA2
            # parts[8] = C2 (chain)
            # parts[9] = IPDB2 (resid)
            
            chain1 = parts[4]
            chain2 = parts[8]
            resid1 = int(parts[5])
            resid2 = int(parts[9])
            
            # Convert to Martinize2 merged resid
            new_resid1 = convert_to_martinize_resid(chain1, resid1)
            new_resid2 = convert_to_martinize_resid(chain2, resid2)
            
            # Keep original spacing
            # Find positions of columns
            pos_c1 = line.find(f" {chain1} {resid1}")
            pos_c2 = line.find(f" {chain2} {resid2}", pos_c1 + 1)
            
            # Replace residue IDs
            new_line = line[:pos_c1+2] + f"{chain1} {new_resid1:>4}" + line[pos_c1+7:pos_c2+2] + f"{chain2} {new_resid2:>4}" + line[pos_c2+7:]
            
            fixed_lines.append(new_line)
            contact_count += 1
            
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse: {line.strip()[:60]}... ({e})")
            fixed_lines.append(line)
    
    with open(output_file, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"Fixed {contact_count} contact records")
    print(f"Output written to {output_file}")

if __name__ == '__main__':
    fix_contact_map('contact_map_fixed.out', 'contact_map_martinize.out')
