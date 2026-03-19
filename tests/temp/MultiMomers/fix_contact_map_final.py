#!/usr/bin/env python3
"""
Fix contact map for Martinize2 merged hemoglobin.
Martinize2 residue numbering:
- Chain A: 1-141
- Chain B: 142-286
- Chain C: 287-427
- Chain D: 428-572
"""

def convert_resid(chain, resid):
    if chain == 'A':
        return resid
    elif chain == 'B':
        return resid + 140
    elif chain == 'C':
        return resid + 286
    elif chain == 'D':
        return resid + 426
    return resid

with open('contact_map_fixed.out') as f:
    lines = f.readlines()

fixed = []
count = 0

for line in lines:
    if not line.strip().startswith('R '):
        fixed.append(line)
        continue
    
    parts = line.split()
    if len(parts) < 11:
        fixed.append(line)
        continue
    
    try:
        # R ID I1 AA C IPDB I2 AA C IPDB DCA ...
        new_resid1 = convert_resid(parts[4], int(parts[5]))
        new_resid2 = convert_resid(parts[8], int(parts[9]))
        
        # Reconstruct with proper spacing
        new_line = f"       {parts[0]:>4} {parts[1]:>4} {parts[2]:>4}  {parts[3]:>3}  {parts[4]:>1} {new_resid1:>4} {parts[6]:>6} {parts[7]:>3}  {parts[8]:>1} {new_resid2:>4}"
        for p in parts[10:]:
            new_line += f" {p:>10}"
        new_line += "\n"
        
        fixed.append(new_line)
        count += 1
    except:
        fixed.append(line)

with open('contact_map_final.out', 'w') as f:
    f.writelines(fixed)

print(f"Fixed {count} contacts -> contact_map_final.out")
