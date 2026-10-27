"""
Single-basin topology generation module for CTGoMartini.

This module provides functionality for generating single-basin (SB) topologies
for Martini Go-models.
"""

from __future__ import annotations

import re
from typing import Any

from .combiner import extract_contacts_from_topology, _resolve_contact_format
from ..builder import MartiniTopFile


def _change_atom_type(
    atoms: list[list[str]], 
    mol_name: str, 
    sbmol_name: str
) -> list[list[str]]:
    """
    Change atom types by replacing the molecule name prefix.
    
    Args:
        atoms: List of atom fields where each field contains atom information.
        mol_name: Original molecule name to be replaced.
        sbmol_name: New molecule name to replace the original.
        
    Returns:
        List of atom fields with updated atom types.
    """
    # Lambda function to extract molecule name and residue ID from atom type
    extract = lambda atomtype: re.findall(r'^(\w+)_(\d+)$', atomtype)[0]

    newatoms: list[list[str]] = []
    for field in atoms:
        atomtype = field[1]
        try:
            # Extract molecule name and residue ID from atom type
            mol_name_extract, mol_resid_extract = extract(atomtype)
            if mol_name_extract != mol_name:
                print(f"mol_name_extract != mol_name: {atomtype} vs {mol_name}")
            # Create new atom type with the new molecule name
            newatomtype = sbmol_name + '_' + mol_resid_extract
        except (IndexError, ValueError):
            # If atom type doesn't match expected pattern, keep it unchanged
            newatomtype = atomtype
        # Reconstruct the atom field with the new atom type
        newatoms.append([field[0], newatomtype] + field[2:])
    return newatoms


def _apply_contact_format(
    contacts: list[list[str]],
    contact_type: int,
    contact_rcut: float | None = None,
) -> list[list[str]]:
    """Convert fixed-cutoff contacts to the requested contact type.

    Rewrites each contact entry from functype 1 (fixed 1.1 nm runtime cutoff)
    to functype 2 (6-12) or 3 (10-12) with a per-contact runtime cutoff.
    Type 1 keeps the legacy fixed-cutoff format unchanged.

    Args:
        contacts: Contact entries [atomid1, atomid2, '1', sigma, epsilon].
        contact_type: Contact potential type (1, 2 or 3).
        contact_rcut: Runtime cutoff in nanometers, appended to each entry
            for types 2/3 (default 1.1 when omitted); ignored for type 1.

    Returns:
        Converted entries [atomid1, atomid2, functype, sigma, epsilon, rcut],
        or the input unchanged for type 1.
    """
    functype, rcut = _resolve_contact_format(contact_type, contact_rcut)
    if functype == '1':
        return contacts
    converted: list[list[str]] = []
    for fields in contacts:
        assert fields[2] == '1', f"Error: only support functype 1 contacts: {fields}"
        converted.append(fields[:2] + [functype] + fields[3:] + [rcut])
    return converted


def _apply_contact_rcut(
    contacts: list[list[str]], contact_rcut: float
) -> list[list[str]]:
    """Convert fixed-cutoff contacts to adjustable 6-12 contacts.

    Rewrites each contact entry from functype 1 (fixed 1.1 nm runtime cutoff)
    to functype 2 with a per-contact runtime cutoff.

    Args:
        contacts: Contact entries [atomid1, atomid2, '1', sigma, epsilon].
        contact_rcut: Runtime cutoff in nanometers, appended to each entry.

    Returns:
        Converted entries [atomid1, atomid2, '2', sigma, epsilon, rcut].
    """
    return _apply_contact_format(contacts, 2, contact_rcut)
    

def create_sb_topology(
    mols_list: list[list[str]],
    sbmol_name: str,
    contact_rcut: float | None = None,
    contact_type: int = 2,
) -> Any:
    """
    Generate a single-basin (SB) Martini topology from a list of molecules.

    Creates a unified topology for single-basin Martini models.

    Args:
        mols_list: List of [topfile, mol_name] pairs.
            Example: [['topfileA.itp', 'molA'], ['topfileB.itp', 'molB']]
        sbmol_name: Name for the new single-basin molecule.
        contact_rcut: Runtime cutoff for Go contact potentials in nanometers.
            Used by the adjustable types 2/3 (default 1.1 when omitted);
            ignored by type 1.
        contact_type: Contact potential type: 1 (6-12 LJ, fixed 1.1 nm
            cutoff, legacy format without rcut field), 2 (6-12 LJ with
            per-contact adjustable cutoff) or 3 (10-12 LJ with per-contact
            adjustable cutoff).
        
    Returns:
        SBMolecule: Single-basin molecule object with topology attributes:
            - sbmol_name: Name of the molecule
            - _topology: Dictionary containing all topology sections
            
    Notes:
        - scFix influence on the SBCombination
        - Angles: BBB_angles, BBS_regular_angles, SBB_regular_angles, BSS_angles,
          SBB_scFix_angles, BBS_scFix_angles, not_BBB_angles
        - Dihedrals: BBBB_dihedrals, SSSS_dihedrals, SBBS_scFix_dihedrals
    """
    # Get the number of molecules in the list
    n_mols = len(mols_list)
    # Initialize list to store molecule name and molecule object pairs
    mols_pairs: list[list[Any]] = []
    
    # Process each molecule in the list
    for (topfile, mol_name) in mols_list:
        # Load the topology file
        top = MartiniTopFile(topfile)
        # Get the specific molecule from the topology
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
            
        # Add the processed molecule to the pairs list
        mols_pairs.append([mol_name, mol])

    # Define a class for single-basin molecules
    class SBMolecule:
        """
        Class representing a single-basin molecule with topology information.
        
        Attributes:
            sbmol_name: Name of the single-basin molecule.
            _topology: Dictionary containing all topology sections.
        """
        def __init__(self, sbmol_name: str) -> None:
            """
            Initialize a single-basin molecule.
            
            Args:
                sbmol_name: Name of the single-basin molecule.
            """
            self.sbmol_name: str = sbmol_name
            self._topology: dict[str, Any] = {}
    
    # Create a new single-basin molecule instance
    sbmol = SBMolecule(sbmol_name)
    # Add molecule type information to the topology
    sbmol._topology['moleculetype'] = [[sbmol_name, '1']]
    
    # Get the first molecule from the pairs list as the base
    mol_name, mol = mols_pairs[0]
    
    # Copy and adapt atoms from the base molecule
    sbmol._topology['atoms'] = _change_atom_type(
        mol._topology['atoms'], mol_name, sbmol_name
    )

    # Copy simple topology sections from the base molecule
    sbmol._topology['bonds'] = mol._topology['bonds']
    sbmol._topology['constraints'] = mol._topology['constraints']
    sbmol._topology['angles'] = mol._topology['angles']
    sbmol._topology['dihedrals'] = mol._topology['dihedrals']
    sbmol._topology['contacts'] = mol._topology['contacts']
    sbmol._topology['exclusions'] = mol._topology['exclusions']

    # Convert contacts to the requested type (no-op for the legacy type 1)
    sbmol._topology['contacts'] = _apply_contact_format(
        sbmol._topology['contacts'], contact_type, contact_rcut
    )

    # Process any additional topology categories
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
        if processed_category in categories_list:
            categories_list.remove(processed_category)
    
    # Copy remaining categories from the first molecule
    for category in categories_list:
        mols_category_list = [mol._topology[category] for (_, mol) in mols_pairs]
        sbmol._topology[category] = mols_category_list[0]

    # Clean up empty topology categories
    sbcategories_list = list(sbmol._topology.keys())
    for category in sbcategories_list:
        if sbmol._topology[category] == []:
            sbmol._topology.pop(category)

    return sbmol
