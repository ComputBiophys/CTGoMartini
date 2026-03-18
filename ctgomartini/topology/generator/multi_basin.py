"""
Multiple-basin topology generation module for CTGoMartini.

This module provides functionality for generating multiple-basin (MB) topologies
for Martini Go-models, allowing simulation of conformational transitions.
"""

from __future__ import annotations

from typing import Any

from .combiner import (
    extract_contacts_from_topology, 
    combine_atoms,
    combine_bonds_constraints,
    combine_exclusions,
    combine_contacts,
    combine_angles,
    combine_dihedrals,
    differentiate_angles, 
    differentiate_dihedrals, 
    SameListList
)
from ..builder import MartiniTopFile


def create_mb_topology(
    mols_list: list[list[str]],
    mbmol_name: str,
    dict_cutoffs: dict[str, float] | None = None
) -> Any:
    """
    Generate a multiple-basin (MB) topology for Martini Go-models.
    
    Combines multiple molecular topologies into a single multiple-basin topology,
    allowing simulation of molecules that can exist in multiple conformational states.
    
    Args:
        mols_list: List of [topfile, mol_name] pairs for each molecule to combine.
            Example: [['topfileA.itp', 'molA'], ['topfileB.itp', 'molB']]
        mbmol_name: Name for the new multiple-basin molecule.
        dict_cutoffs: Dictionary of cutoff values for different interaction types.
            Keys: 'cutoff_BBB_angles', 'cutoff_BBBB_dihedrals', 
                  'cutoff_SBBS_dihedrals', 'cutoff_contacts'.
            Default is all zeros (no cutoffs).

    Returns:
        MBMolecule: A multiple-basin molecule object with combined topology.
        The object has attributes:
            - mbmol_name: Name of the molecule
            - _topology: Dictionary containing all topology sections
    """
    # Set default cutoffs if not provided
    if dict_cutoffs is None:
        dict_cutoffs = {
            'cutoff_BBB_angles': 0.0,
            'cutoff_BBBB_dihedrals': 0.0,
            'cutoff_SBBS_dihedrals': 0.0,
            'cutoff_contacts': 0.0,
        }
    
    # Note on scFix (side-chain fix) influence on MB combination:
    # The function handles different angle types:
    # - BBB_angles (backbone-backbone-backbone)
    # - BBS_regular_angles, SBB_regular_angles, BSS_angles
    # - BBS_scFix_angles, SBB_scFix_angles (side-chain fixed)
    # - not_BBB_angles (all other angles)
    #
    # And different dihedral types:
    # - BBBB_dihedrals (backbone dihedrals)
    # - SSSS_dihedrals (side-chain dihedrals)
    # - SBBS_scFix_dihedrals (side-chain fixed dihedrals)

    # Get the number of molecules to combine
    n_mols = len(mols_list)

    # Initialize topology components for each molecule
    mols_pairs: list[list[Any]] = []
    for (topfile, mol_name) in mols_list:
        # Load topology file and extract molecule
        top = MartiniTopFile(topfile)
        mol = top.moleculeTypes[mol_name]
        
        # Ensure all required topology sections exist
        if 'contacts' not in mol._topology:
            contacts = extract_contacts_from_topology(top, mol_name)
            if contacts:
                mol._topology['contacts'] = contacts        
        if 'bonds' not in mol._topology:
            mol._topology['bonds'] = []
        if 'constraints' not in mol._topology:
            mol._topology['constraints'] = []
        if 'angles' not in mol._topology:
            mol._topology['angles'] = []
        if 'dihedrals' not in mol._topology:
            mol._topology['dihedrals'] = []
        if 'exclusions' not in mol._topology:
            mol._topology['exclusions'] = []
        mols_pairs.append([mol_name, mol])

    # Define a class for the multiple-basin molecule
    class MBMolecule:
        """
        Class representing a multiple-basin molecule with combined topology.
        
        Attributes:
            mbmol_name: Name of the multiple-basin molecule.
            _topology: Dictionary containing all topology sections.
        """
        def __init__(self, mbmol_name: str) -> None:
            self.mbmol_name: str = mbmol_name
            self._topology: dict[str, Any] = {}
    
    # Create multiple-basin molecule instance
    mbmol = MBMolecule(mbmol_name)
    
    # Add moleculetype definition
    mbmol._topology['moleculetype'] = [[mbmol_name, '1']]
    
    # Add multiple_basin section with energy basin definitions
    energy_basin_list = [f'C{i+1}' for i in range(n_mols)]
    mbmol._topology['multiple_basin'] = [['True', 'exp', str(n_mols), 'beta'] + energy_basin_list]
    
    # Combine atoms from all molecules
    mols_atoms_pairs = [[mol_name, mol._topology['atoms']] for (mol_name, mol) in mols_pairs]
    mbatoms = combine_atoms(mbmol_name, mols_atoms_pairs)
    mbmol._topology['atoms'] = mbatoms

    # Add bonds and constraints from all molecules
    mols_bonds_list = [mol._topology['bonds'] for (_, mol) in mols_pairs]
    mols_constraints_list = [mol._topology['constraints'] for (_, mol) in mols_pairs]
    mbbonds, mbconstraints = combine_bonds_constraints(
        n_mols, mols_bonds_list, mols_constraints_list
    )
    mbmol._topology['bonds'] = mbbonds
    mbmol._topology['constraints'] = mbconstraints

    # Add angles and multi-angles
    mols_BBB_angles_list: list[Any] = []
    mols_notBBB_angles_list: list[Any] = []
    for _, mol in mols_pairs:
        mol_BBB_angles, mol_notBBB_angles = differentiate_angles(
            mol._topology['angles'], mol._topology['atoms']
        )
        mols_BBB_angles_list.append(mol_BBB_angles)
        mols_notBBB_angles_list.append(mol_notBBB_angles)

    # Combine BBB angles with cutoff, keep non-BBB angles unchanged
    mb_BBB_angles, mbmulti_BBB_angles = combine_angles(
        n_mols, mols_BBB_angles_list, dict_cutoffs['cutoff_BBB_angles']
    )
    assert SameListList(mols_notBBB_angles_list, sort=True), (
        f'Error: notBBB angles from different mols should be the same! {mols_notBBB_angles_list}'
    )
    mbangles = mb_BBB_angles + mols_notBBB_angles_list[0]
    mbmol._topology['angles'] = mbangles
    mbmol._topology['multi_angles'] = mbmulti_BBB_angles

    # Add dihedrals and multi-dihedrals
    mols_BBBB_dihedrals_list: list[Any] = []
    mols_SSSS_dihedrals_list: list[Any] = []
    mols_SBBS_dihedrals_list: list[Any] = []
    for _, mol in mols_pairs:
        mols_BBBB_dihedrals, mols_SSSS_dihedrals, mols_SBBS_dihedrals = differentiate_dihedrals(
            mol._topology['dihedrals'], mol._topology['atoms']
        )
        mols_BBBB_dihedrals_list.append(mols_BBBB_dihedrals)
        mols_SSSS_dihedrals_list.append(mols_SSSS_dihedrals)
        mols_SBBS_dihedrals_list.append(mols_SBBS_dihedrals)

    # Combine BBBB dihedrals with cutoff, keep SSSS and SBBS dihedrals unchanged
    mb_BBBB_dihdedrals, mbmulti_BBBB_dihedrals = combine_dihedrals(
        n_mols, mols_BBBB_dihedrals_list, dict_cutoffs['cutoff_BBBB_dihedrals']
    )
    assert SameListList(mols_SSSS_dihedrals_list, sort=True), (
        f'Error: SSSS dihedrals from different mols should be the same! {mols_SSSS_dihedrals_list}'
    )
    assert SameListList(mols_SBBS_dihedrals_list, sort=True), (
        f'Error: Dont support the scFix for multiple-basin Go-Martini. '
        f'SBBS dihedrals from different mols should be the same! {mols_SBBS_dihedrals_list}'
    )
    mbdihedrals = mb_BBBB_dihdedrals + mols_SSSS_dihedrals_list[0] + mols_SBBS_dihedrals_list[0]
    mbmulti_dihedrals = mbmulti_BBBB_dihedrals
    mbmol._topology['dihedrals'] = mbdihedrals
    mbmol._topology['multi_dihedrals'] = mbmulti_dihedrals    

    # Add contacts and multi-contacts
    mols_contacts_list = [mol._topology['contacts'] for (_, mol) in mols_pairs]
    mbcontacts, mbmulti_contacts = combine_contacts(
        n_mols, mols_contacts_list, dict_cutoffs['cutoff_contacts']
    )
    mbmol._topology['contacts'] = mbcontacts
    mbmol._topology['multi_contacts'] = mbmulti_contacts        

    # Add exclusions
    mols_exclusions_list = [mol._topology['exclusions'] for (_, mol) in mols_pairs]
    mbexclusions = combine_exclusions(mols_exclusions_list)
    mbmol._topology['exclusions'] = mbexclusions

    # Process other topology categories
    categories_list: list[str] = []
    for (_, mol) in mols_pairs:
        categories_list += list(mol._topology.keys())
    categories_list = list(set(categories_list))
    
    # Remove already processed categories
    processed_categories = [
        'moleculetype', 'atoms', 'bonds', 'constraints', 
        'angles', 'dihedrals', 'contacts', 'exclusions'
    ]
    for processed_category in processed_categories:
        if processed_category in processed_categories:
            categories_list.remove(processed_category)
    
    # Check remaining categories are identical and add to MB topology
    for category in categories_list:
        mols_category_list = [mol._topology[category] for (_, mol) in mols_pairs]
        assert SameListList(mols_category_list, sort=True), (
            f'Error: {category} from different mols should be the same! {mols_category_list}'
        )
        mbmol._topology[category] = mols_category_list[0]

    # Remove empty categories
    mbcategories_list = list(mbmol._topology.keys())
    for category in mbcategories_list:
        if mbmol._topology[category] == []:
            mbmol._topology.pop(category)

    return mbmol
