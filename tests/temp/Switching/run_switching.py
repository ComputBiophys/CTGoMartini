#!/usr/bin/env python3
"""
Test script for Switching Go-Martini mode with vermouth >= 0.15.0

This script demonstrates the new workflow:
1. Generate Go contacts using martinize2 -go option
2. Convert LJ nonbond_params to [contacts] for computational efficiency and accuracy
3. Generate separate ITP files for each state (Open and Closed)
"""

import sys
import os

# Add ctgomartini to path
sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')

from ctgomartini.data.ctgomartinize import SwitchingGOMartinize


def main():
    """Run Switching Go-Martini topology generation."""
    
    print("=" * 70)
    print("Switching Go-Martini Test (vermouth >= 0.15.0)")
    print("=" * 70)
    print()
    
    # Configuration
    config = {
        'aa_strfile_list': ['State1.pdb', 'State2.pdb'],
        'aa_map_list': ['State1.map', 'State2.map'],
        'state_name_list': ['Open', 'Closed'],
        'mbmol_name': 'Protein',
        'dict_cutoffs': {
            'cutoff_BBB_angles': 15.0,
            'cutoff_BBBB_dihedrals': 30.0,
            'cutoff_SBBS_dihedrals': 30.0,
            'cutoff_contacts': 0.06
        },
        'method': 'switching',
        'dssp': None,  # Use MDTraj (no external DSSP required)
        'ff': 'martini3001',
        'other_params': ''
    }
    
    print("Configuration:")
    print(f"  States: {config['state_name_list']}")
    print(f"  Method: {config['method']}")
    print(f"  DSSP: MDTraj (default)")
    print(f"  Force Field: {config['ff']}")
    print()
    
    # Run SwitchingGOMartinize
    SwitchingGOMartinize(**config)
    
    print()
    print("=" * 70)
    print("Switching Go-Martini generation completed!")
    print("=" * 70)
    print()
    print("Generated files:")
    print("  - Open/Open.itp, Open/Open_params.itp")
    print("  - Closed/Closed.itp, Closed/Closed_params.itp")
    print()
    print("Each ITP contains:")
    print("  - [atoms] section (coarse-grained beads)")
    print("  - [contacts] section (Go contacts from go_nbparams.itp)")
    print("  - Standard bonded interactions (bonds, angles, dihedrals)")


if __name__ == '__main__':
    main()
