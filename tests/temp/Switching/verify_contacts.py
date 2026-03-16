#!/usr/bin/env python3
"""
Verification script for Go contacts conversion.

Compares go_nbparams.itp (LJ nonbond_params) with ITP [contacts] section
to ensure all contacts are correctly converted.
"""

import os


def parse_go_nbparams(filepath: str) -> set:
    """Parse go_nbparams.itp and return set of contacts."""
    contacts = set()
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('['):
                continue
            parts = line.split()
            if len(parts) >= 5 and parts[2] == '1':
                # Extract atom numbers from Open_X or Closed_X
                atom1 = int(parts[0].split('_')[1])
                atom2 = int(parts[1].split('_')[1])
                sigma = float(parts[3])
                key = tuple(sorted([atom1, atom2]))
                contacts.add((key, sigma))
    return contacts


def parse_itp_contacts(filepath: str) -> set:
    """Parse ITP file [contacts] section."""
    contacts = set()
    in_contacts = False
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line == '[ contacts ]':
                in_contacts = True
                continue
            if in_contacts and line.startswith('['):
                break
            if in_contacts and line and line[0].isdigit():
                parts = line.split()
                if len(parts) >= 5 and parts[2] == '1':
                    atom1 = int(parts[0])
                    atom2 = int(parts[1])
                    sigma = float(parts[3])
                    # Convert virtual_site ID to atom ID (512 offset)
                    atom1_real = atom1 - 512 if atom1 > 512 else atom1
                    atom2_real = atom2 - 512 if atom2 > 512 else atom2
                    key = tuple(sorted([atom1_real, atom2_real]))
                    contacts.add((key, sigma))
    return contacts


def verify_state(state_name: str) -> bool:
    """Verify contacts for a single state."""
    print(f"\n【{state_name} 状态】")
    
    nb_file = f'{state_name}/go_nbparams.itp'
    itp_file = f'{state_name}.itp'
    
    if not os.path.exists(nb_file):
        print(f"  ❌ {nb_file} not found")
        return False
    
    if not os.path.exists(itp_file):
        print(f"  ❌ {itp_file} not found")
        return False
    
    # Parse both files
    nb_contacts = parse_go_nbparams(nb_file)
    itp_contacts = parse_itp_contacts(itp_file)
    
    print(f"  go_nbparams.itp 接触对: {len(nb_contacts)}")
    print(f"  {state_name}.itp contacts: {len(itp_contacts)}")
    
    # Find differences
    missing = nb_contacts - itp_contacts
    extra = itp_contacts - nb_contacts
    
    if len(missing) == 0 and len(extra) == 0:
        print(f"  ✅ 完全匹配 ({len(nb_contacts)} contacts)")
        return True
    else:
        print(f"  ❌ 不匹配:")
        print(f"     - 遗漏: {len(missing)} contacts")
        print(f"     - 多余: {len(extra)} contacts")
        if missing:
            print(f"     - 示例遗漏: {list(missing)[:3]}")
        return False


def main():
    """Run verification for all states."""
    print("=" * 70)
    print("Go Contacts Verification")
    print("=" * 70)
    print()
    print("Comparing go_nbparams.itp (LJ) with ITP [contacts] (bonds)")
    print()
    
    states = ['Open', 'Closed']
    results = []
    
    for state in states:
        results.append(verify_state(state))
    
    print()
    print("=" * 70)
    if all(results):
        print("✅ 所有状态验证通过！")
    else:
        print("❌ 部分状态验证失败")
    print("=" * 70)


if __name__ == '__main__':
    main()
