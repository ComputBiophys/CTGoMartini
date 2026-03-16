"""Generate MB HAM topology for protein conformational transition."""
import os
import sys

# Add the package to the path
sys.path.insert(0, '/home/ys/CommonUse/Martini/CTGoMartini')

from ctgomartini.data.ctgomartinize import MBGOMartinize


def main():
    """Generate topology using MB HAM method."""
    print("="*70)
    print("Generating MB HAM Topology")
    print("="*70)
    
    # Configuration
    aa_strfile_list = ["State1.pdb", "State2.pdb"]
    aa_map_list = ["State1.map", "State2.map"]
    state_name_list = ["Open", "Closed"]
    mbmol_name = "Protein"
    
    # Distance cutoffs for determining shared interactions
    dict_cutoffs = {
        'cutoff_BBB_angles': 15.0,
        'cutoff_BBBB_dihedrals': 30.0,
        'cutoff_SBBS_dihedrals': 30.0,
        'cutoff_contacts': 0.06
    }
    
    # Generate topology with HAM method
    MBGOMartinize(
        aa_strfile_list=aa_strfile_list,
        aa_map_list=aa_map_list,
        state_name_list=state_name_list,
        mbmol_name=mbmol_name,
        dict_cutoffs=dict_cutoffs,
        method='ham',
        dssp=None,
        ff='martini3001',
        other_params=''
    )
    
    print("\n" + "="*70)
    print("Topology generated successfully!")
    print("Generated files:")
    print("  - Protein.itp: Main topology file")
    print("  - Protein_params.itp: Force field parameters")
    print("="*70)


if __name__ == "__main__":
    main()
