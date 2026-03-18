"""Module for combining molecular topologies in multiple-basin Go-Martini simulations.

This module provides functions and classes for combining molecular topologies
from different conformational states into a single multiple-basin topology.
"""

from __future__ import annotations

import os
import re
from collections import OrderedDict
from typing import Any

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..builder import MartiniTopFile


def extract_contacts_from_topology(top: "MartiniTopFile", molecule_name: str) -> list[list[str]]:
    """Extract contact interactions from a topology file for a specific molecule.

    Args:
        top: The Martini topology file object.
        molecule_name: Name of the molecule to extract contacts for.

    Returns:
        List of contact fields, where each contact is a list of strings
        containing atom indices and interaction parameters.

    Raises:
        Exception: If unsupported contacts are found or if contact atoms are missing.
    """
    atoms = top.moleculeTypes[molecule_name]._topology['atoms']
    nonbond_params = top.forcefield._parameters['nonbond_params']

    pattern = rf'{molecule_name}_\d+'
    contact_atoms: dict[str, int] = {}  # atomtype: atomid # hard restraint
    for item in atoms:
        atomtype = item[1]
        if re.fullmatch(pattern, atomtype):
            atomid = int(item[0])
            contact_atoms[atomtype] = atomid

    contacts: list[list[str]] = []
    for fields in nonbond_params:
        atomtype1 = fields[0]
        atomtype2 = fields[1]
        judge1 = bool(re.fullmatch(pattern, atomtype1))
        judge2 = bool(re.fullmatch(pattern, atomtype2))
        if judge1 ^ judge2:
            raise Exception(f"Error: Unsupport the contact bewteen {atomtype1} {atomtype2}!")
        elif judge1 & judge2:
            try:
                atomid1 = contact_atoms[atomtype1]
                atomid2 = contact_atoms[atomtype2]
            except KeyError:
                raise Exception(f"Error: not contact_atoms! {fields}")

            atomid1, atomid2 = sorted([atomid1, atomid2])
            assert fields[2] == '1', f"Error: only support functype 1: {fields}"
            newfields = [str(atomid1), str(atomid2)] + fields[2:]
            contacts.append(newfields)
        else:
            continue

    contacts = sorted(contacts, key=lambda fields: (int(fields[0]), int(fields[1])))
    return contacts


def GetAtomNames(atomidlist: list[str], atoms: list[list[str]]) -> list[str]:
    """Get atom names from a list of atom IDs.

    Args:
        atomidlist: List of atom IDs as strings.
        atoms: List of atom field lists from the topology.

    Returns:
        List of atom names corresponding to the given atom IDs.
    """
    atoms_dict = {int(fields[0]): fields for fields in atoms}

    atomnamelist: list[str] = []
    for atomid in atomidlist:
        fields = atoms_dict[int(atomid)]
        atomname = fields[4]
        atomnamelist.append(atomname)
    return atomnamelist


def GetAngleDiehdralType(atomnamelist: list[str]) -> str:
    """Get the type string for an angle or dihedral based on atom names.

    BB (backbone) atoms are marked as 'B', SC* (side-chain) atoms as 'S'.

    Args:
        atomnamelist: List of atom names (e.g., ['BB', 'SC1', 'BB']).

    Returns:
        Type string composed of 'B' and 'S' characters.
    """
    typestr = ''
    for atomname in atomnamelist:
        if atomname == 'BB':
            typestr += 'B'
        elif atomname.startswith('SC'):
            typestr += 'S'
        else:
            typestr += 'S'
            print(f"Warning: Recognize unknown atomnames. Set this atom as side chain defaultly. {atomnamelist}")
    return typestr


def differentiate_angles(angles: list[list[str]], atoms: list[list[str]]) -> tuple[list[list[str]], list[list[str]]]:
    """Separate angles into BBB (backbone) and non-BBB angles.

    Args:
        angles: List of angle field lists.
        atoms: List of atom field lists.

    Returns:
        Tuple of (BBB_angles, notBBB_angles), where each is a list of angle fields.
    """
    BBB_angles: list[list[str]] = []
    notBBB_angles: list[list[str]] = []

    for fields in angles:
        atomidlist = fields[:3]
        atomnamelist = GetAtomNames(atomidlist, atoms)
        if GetAngleDiehdralType(atomnamelist) == 'BBB':
            BBB_angles.append(fields)
        else:
            notBBB_angles.append(fields)
    return BBB_angles, notBBB_angles


def differentiate_dihedrals(dihedrals: list[list[str]], atoms: list[list[str]]) -> tuple[list[list[str]], list[list[str]], list[list[str]]]:
    """Separate dihedrals into BBBB, SSSS, and SBBS types.

    Args:
        dihedrals: List of dihedral field lists.
        atoms: List of atom field lists.

    Returns:
        Tuple of (BBBB_dihedrals, SSSS_dihedrals, SBBS_dihedrals).

    Raises:
        AssertionError: If unsupported dihedral types are found (e.g., from scfix).
    """
    BBBB_dihedrals: list[list[str]] = []
    SSSS_dihedrals: list[list[str]] = []
    SBBS_dihedrals: list[list[str]] = []
    other_dihedrals: list[list[str]] = []

    for fields in dihedrals:
        atomidlist = fields[:4]
        atomnamelist = GetAtomNames(atomidlist, atoms)
        angle_type = GetAngleDiehdralType(atomnamelist)
        if angle_type == 'BBBB':
            BBBB_dihedrals.append(fields)
        elif angle_type == 'SBBS':
            SBBS_dihedrals.append(fields)
        elif angle_type == 'SSSS':
            SSSS_dihedrals.append(fields)
        else:
            other_dihedrals.append(fields)
    assert other_dihedrals == [], f'Error: not supported dihedral type: {other_dihedrals}'
    return BBBB_dihedrals, SSSS_dihedrals, SBBS_dihedrals


def CombineDict(dict_list: list[dict[Any, Any]]) -> OrderedDict[Any, list[Any]]:
    """Combine multiple dictionaries by merging values for common keys.

    Args:
        dict_list: List of dictionaries to combine.

    Returns:
        OrderedDict with combined values as lists for each unique key.
    """
    n_dict = len(dict_list)
    key_combined_list: list[Any] = []
    for i in range(n_dict):
        key_combined_list += list(dict_list[i].keys())
    key_combined_list = list(set(key_combined_list))
    key_combined_list = sorted(key_combined_list)

    dict_combined: OrderedDict[Any, list[Any]] = OrderedDict()
    for key in key_combined_list:
        dict_combined[key] = []
        for i in range(n_dict):
            if key in dict_list[i].keys():
                dict_combined[key].append(dict_list[i][key])
    return dict_combined


def ForceItemFloat(item: Any, precision: int | None = None) -> Any:
    """Convert an item to float if possible, with optional rounding.

    Args:
        item: The item to convert.
        precision: Number of decimal places for rounding (None for no rounding).

    Returns:
        Float value if conversion succeeds, otherwise the original item.
    """
    try:
        item = float(item)
        if precision is not None:
            item = round(item, precision)
    except (ValueError, TypeError):
        pass
    return item


def ForceListFloat(itemlist: list[Any], precision: int | None = None) -> list[Any]:
    """Recursively convert items in a nested list to floats.

    Args:
        itemlist: List of items or nested lists.
        precision: Number of decimal places for rounding.

    Returns:
        List with float conversions applied recursively.
    """
    newlist: list[Any] = []
    for item in itemlist:
        if not isinstance(item, list):
            newlist.append(ForceItemFloat(item, precision))
        else:
            newlist.append(ForceListFloat(item, precision))
    return newlist


def SameListList(listlist: list[list[Any]], typeforce: bool = True, sort: bool = False, precision: int | None = None) -> bool:
    """Check if all lists in a list are the same.

    Args:
        listlist: List of lists to compare.
        typeforce: If True, convert items to float before comparison.
        sort: If True, sort lists before comparison.
        precision: Number of decimal places for rounding (None for exact comparison).

    Returns:
        True if all lists are the same, False otherwise.
    """
    issame = True
    ref_list = ForceListFloat(listlist[0], precision=precision) if typeforce else listlist[0]
    ref_list = sorted(ref_list) if sort else ref_list
    for item in listlist[1:]:
        item = ForceListFloat(item, precision=precision) if typeforce else item
        item = sorted(item) if sort else item
        if item != ref_list:
            issame = False
    return issame


def SameList(alist: list[Any], typeforce: bool = True) -> bool:
    """Check if all items in a list are the same.

    Args:
        alist: The list to check.
        typeforce: If True, convert items to float before comparison.

    Returns:
        True if all items are the same, False otherwise.
    """
    issame = True
    alistlist = [[item] for item in alist]
    issame = SameListList(alistlist, typeforce)
    return issame


def Calculate_DiffDihedral(dihedral_list: list[float]) -> tuple[float, float]:
    """Calculate the maximum difference and mean of dihedral angles.

    Handles periodic boundary conditions by choosing the direction (clockwise
    or counter-clockwise) where the maximum difference is less than 180 degrees.

    Args:
        dihedral_list: List of dihedral angle values in degrees.

    Returns:
        Tuple of (diff_max, mean_dihedral) where diff_max is the maximum
        angular difference and mean_dihedral is the mean angle in range (-180, 180].

    Raises:
        ValueError: If the dihedral angles cannot be properly resolved.
    """
    dihedral_list = sorted(dihedral_list)
    anticlock_dihedral_list: list[float] = []
    for dihedral in dihedral_list:
        if dihedral < 0:
            dihedral += 360
        anticlock_dihedral_list.append(dihedral)
    anticlock_dihedral_list = sorted(anticlock_dihedral_list)
    anticlock_diff_max = abs(anticlock_dihedral_list[-1] - anticlock_dihedral_list[0])

    clock_dihedral_list = dihedral_list.copy()
    clock_dihedral_list = sorted(clock_dihedral_list)
    clock_diff_max = abs(clock_dihedral_list[-1] - clock_dihedral_list[0])

    if anticlock_diff_max >= 180 and clock_diff_max < 180:
        dihedral_list = clock_dihedral_list
    elif anticlock_diff_max < 180 and clock_diff_max >= 180:
        dihedral_list = anticlock_dihedral_list
    elif anticlock_diff_max == 0 and clock_diff_max == 0:
        dihedral_list = anticlock_dihedral_list
    else:
        print(anticlock_diff_max, clock_diff_max)
        raise ValueError(f'Error: something wrong with the dihedrals {dihedral_list}')

    dihedral_list = sorted(dihedral_list)
    diff_max = abs(dihedral_list[-1] - dihedral_list[0])
    mean_dihedral = sum(dihedral_list) / len(dihedral_list)
    if mean_dihedral > 180:
        mean_dihedral -= 360
    if mean_dihedral <= -180:
        mean_dihedral += 360
    return diff_max, mean_dihedral


def combine_atoms(mbmolname: str, mols_atoms_pairs: list[tuple[str, list[list[str]]]]) -> list[list[str]]:
    """Combine atoms from different states of molecules.

    Args:
        mbmolname: New prefix of atom names for virtual sites.
        mols_atoms_pairs: Atoms from different states of molecules.
            Format: [(molnameA, atomtopA), (molnameB, atomtopB), ...]

    Returns:
        Combined atom topology as a list of atom field lists.

    Raises:
        AssertionError: If number of molecules is less than 2 or atoms counts differ.
        ValueError: If atoms cannot be combined according to the rules.
    """
    n_mols = len(mols_atoms_pairs)
    n_atoms = len(mols_atoms_pairs[0][1])

    # Assert that there are not less than 2 mols
    assert n_mols >= 2, "Error: The number of mols must more than or equal 2"

    # Assert that different molecules have same number of atoms
    for pair in mols_atoms_pairs[1:]:
        assert len(pair[1]) == len(mols_atoms_pairs[0][1])

    mols_atoms_dict_list: list[dict[tuple[int, ...], list[str]]] = []
    for pair in mols_atoms_pairs:
        atoms = pair[1]
        mols_atoms_dict_list.append({(int(atoms[i][0]),): atoms[i] for i in range(n_atoms)})
    mols_atoms_dict_combined = CombineDict(mols_atoms_dict_list)

    mbmol_atomtop: list[list[str]] = []

    def Extract(atomtype: str) -> tuple[str, str]:
        match = re.findall(r'^(\w+)_(\d+)$', atomtype)
        return match[0]

    for key, value in mols_atoms_dict_combined.items():
        assert len(value) == n_mols, f"Error: The number of molecules with the same atomid is not equal to the number of molecules. {key}"
        if SameListList(value):
            mbmol_atomtop.append(value[0].copy())
        else:
            try:
                mol_resid_extract_list: list[str] = []
                for i, atom in enumerate(value):
                    atomtype = atom[1]
                    mol_name_extract, mol_resid_extract = Extract(atomtype)
                    assert mol_name_extract == mols_atoms_pairs[i][0]
                    mol_resid_extract_list.append(mol_resid_extract)

                if SameList(mol_resid_extract_list):
                    newatomtype = f"{mbmolname}_{mol_resid_extract}"
                    newatom = value[0].copy()
                    newatom[1] = newatomtype
                    mbmol_atomtop.append(newatom)
                else:
                    raise ValueError

            except Exception:
                raise ValueError("Error: atoms from different states of one molecule cannot meet the combination rule!")

    assert len(mbmol_atomtop) == n_atoms
    return mbmol_atomtop


def combine_bonds_constraints(
    n_mols: int,
    mols_bonds_list: list[list[list[str]]],
    mols_constraints_list: list[list[list[str]]]
) -> tuple[list[list[str]], list[list[str]]]:
    """Combine bonds and constraints from multiple molecular states.

    Args:
        n_mols: Number of molecules being combined.
        mols_bonds_list: List of bond lists from each molecule.
        mols_constraints_list: List of constraint lists from each molecule.

    Returns:
        Tuple of (mbbonds, mbconstraints) as combined topology fields.
    """
    mols_connections_dict_list: list[dict[tuple[int, int], list[Any]]] = []
    for i, bonds in enumerate(mols_bonds_list):
        connection_dict: dict[tuple[int, int], list[Any]] = {}
        state = str(i + 1)
        for bond in bonds:
            assert bond[2] == "1", f"Error: bond type must be 1: {bond}"
            if int(bond[0]) > int(bond[1]):
                bond[:2] = [bond[1], bond[0]]
            key = (int(bond[0]), int(bond[1]))
            connection_dict[key] = [state] + bond
        mols_connections_dict_list.append(connection_dict)

    for i, constraints in enumerate(mols_constraints_list):
        connection_dict = {}
        state = str(i + 1)
        for constraint in constraints:
            assert constraint[2] == "1", f"Error: constraint type must be 1: {constraint}"
            if int(constraint[0]) > int(constraint[1]):
                constraint[:2] = [constraint[1], constraint[0]]
            key = (int(constraint[0]), int(constraint[1]))
            connection_dict[key] = [state] + constraint[:4] + [None]
        mols_connections_dict_list.append(connection_dict)

    mols_connections_dict_combined = CombineDict(mols_connections_dict_list)

    mbconnections: list[list[str]] = []
    for key, value in mols_connections_dict_combined.items():
        n_states_in_value = len({fields[0] for fields in value})
        assert n_states_in_value == len(value), f"Error: value repeats! {value}"
        assert n_states_in_value == n_mols, f"Error: {key} does not have {n_mols} values"
        dist_list = [float(fields[1:][3]) for fields in value]
        k_list = [float(fields[1:][4]) for fields in value if fields[1:][4] is not None]
        dist_mean = sum(dist_list) / len(dist_list)
        dist_mean_str = str(round(dist_mean, 3))
        if k_list:
            k_mean = sum(k_list) / len(k_list)
            k_mean_str = str(round(k_mean, 3))
        else:
            k_mean_str = ''

        mbconnections.append(value[0][1:][:3] + [dist_mean_str, k_mean_str])

    mbbonds: list[list[str]] = []
    mbconstraints: list[list[str]] = []
    for item in mbconnections:
        if item[4] != '':
            mbbonds.append(item)
        else:
            mbconstraints.append(item[:4])

    return mbbonds, mbconstraints


def combine_exclusions(mols_exclusions_list: list[list[list[str]]]) -> list[list[str]]:
    """Combine exclusions from multiple molecular states.

    Args:
        mols_exclusions_list: List of exclusion lists from each molecule.

    Returns:
        Combined exclusion topology fields.
    """
    exclusion_pair_list: list[tuple[int, int]] = []
    for exclusions in mols_exclusions_list:
        for fields in exclusions:
            item0 = fields[0]
            for item in fields[1:]:
                exclusion_pair_list.append(tuple(sorted([int(item0), int(item)])))
    exclusion_pair_list = sorted(list(set(exclusion_pair_list)))

    mbexclusion_dict: OrderedDict[tuple[int, ...], list[str]] = OrderedDict()
    for item in exclusion_pair_list:
        key = (item[0],)
        if key not in mbexclusion_dict:
            mbexclusion_dict[key] = [str(item[0])]
        if str(item[1]) not in mbexclusion_dict[key]:
            mbexclusion_dict[key].append(str(item[1]))

    mbexclusions = list(mbexclusion_dict.values())
    return mbexclusions


def combine_contacts(
    n_mols: int,
    mols_contacts_list: list[list[list[str]]],
    cutoff: float
) -> tuple[list[list[str]], list[list[str]]]:
    """Combine contacts from multiple molecular states.

    Args:
        n_mols: Number of molecules being combined.
        mols_contacts_list: List of contact lists from each molecule.
        cutoff: Distance cutoff for considering contacts as similar.

    Returns:
        Tuple of (mbcontacts, mbmulti_contacts) where mbcontacts are
        contacts with averaged parameters and mbmulti_contacts are
        state-specific contacts.
    """
    assert len(mols_contacts_list) == n_mols, "Error: The number of contacts is not equal to the number of molecules."
    mols_contacts_dict_list: list[dict[tuple[int, int], list[str]]] = []
    for i, contacts in enumerate(mols_contacts_list):
        mols_contacts_dict: dict[tuple[int, int], list[str]] = {}
        state = str(i + 1)
        for fields in contacts:
            assert fields[2] == '1', f"Error: contact type is not 1. {fields}"
            if int(fields[0]) > int(fields[1]):
                fields[:2] = [fields[1], fields[0]]
            key = (int(fields[0]), int(fields[1]))
            if key not in mols_contacts_dict:
                mols_contacts_dict[key] = [state] + fields
        mols_contacts_dict_list.append(mols_contacts_dict)
    mols_contacts_dict_combined = CombineDict(mols_contacts_dict_list)

    mbcontacts: list[list[str]] = []
    mbmulti_contacts: list[list[str]] = []
    for key, value in mols_contacts_dict_combined.items():
        n_states_in_value = len({fields[0] for fields in value})
        assert n_states_in_value == len(value), f"Error: contact value repeats! {value}"
        if n_states_in_value != n_mols:
            for fields in value:
                state = fields[0]
                fields = fields[1:]
                newfields = fields[:2] + [str(n_mols), state] + fields[2:]
                mbmulti_contacts.append(newfields)
        else:
            sigma_list = [float(fields[1:][3]) for fields in value]
            epsilon_list = [float(fields[1:][4]) for fields in value]
            diff_sigma = abs(max(sigma_list) - min(sigma_list))
            if diff_sigma <= cutoff:
                mean_sigma = round(sum(sigma_list) / len(sigma_list), 10)
                mean_epsilon = round(sum(epsilon_list) / len(epsilon_list), 10)
                mbcontacts.append([str(key[0]), str(key[1]), '1', str(mean_sigma), str(mean_epsilon)])
            else:
                for fields in value:
                    state = fields[0]
                    fields = fields[1:]
                    newfields = fields[:2] + [str(n_mols), state] + fields[2:]
                    mbmulti_contacts.append(newfields)

    return mbcontacts, mbmulti_contacts


def combine_angles(
    n_mols: int,
    mols_angles_list: list[list[list[str]]],
    cutoff: float
) -> tuple[list[list[str]], list[list[str]]]:
    """Combine angles from multiple molecular states.

    Converts angle type 2 (g96 angles) to type 10 (restricted angles)
    if the angles of the same atoms from different states have different types.

    Args:
        n_mols: Number of molecules being combined.
        mols_angles_list: List of angle lists from each molecule.
        cutoff: Angle cutoff in degrees for considering angles as similar.

    Returns:
        Tuple of (mbangles, mbmulti_angles) where mbangles are angles
        with averaged parameters and mbmulti_angles are state-specific angles.
    """
    mols_angles_dict_list: list[dict[tuple[int, int, int], list[str]]] = []
    assert len(mols_angles_list) == n_mols, 'The number of molecules is not equal to the number of angles'
    for i, angles in enumerate(mols_angles_list):
        mols_angles_dict: dict[tuple[int, int, int], list[str]] = {}
        state = str(i + 1)
        for fields in angles:
            assert fields[3] in ['2', '10'], f"Error: angle type is not 2 or 10. {fields}"
            assert 0 <= float(fields[4]) <= 180, f"Error: angles should be in 0-180. {fields}"
            if int(fields[0]) > int(fields[2]):
                fields[:3] = [fields[2], fields[1], fields[0]]
            key = (int(fields[0]), int(fields[1]), int(fields[2]))
            if key not in mols_angles_dict:
                mols_angles_dict[key] = [state] + fields
        mols_angles_dict_list.append(mols_angles_dict)
    mols_angles_dict_combined = CombineDict(mols_angles_dict_list)

    def g96Torestricted(fields: list[str]) -> list[str]:
        """Convert g96 angle type (2) to restricted angle type (10)."""
        if fields[3] == '2':
            newfields = fields[:3] + ['10', fields[4], '25.0']
        elif fields[3] == '10':
            newfields = fields.copy()
        else:
            raise ValueError(f"Error: Unsupport angle type. {fields}")
        return newfields

    mbangles: list[list[str]] = []
    mbmulti_angles: list[list[str]] = []
    for key, value in mols_angles_dict_combined.items():
        n_states_in_value = len({fields[0] for fields in value})
        assert n_states_in_value == len(value), f"Error: angle value repeats! {value}"
        if n_states_in_value != n_mols:
            for fields in value:
                state = fields[0]
                fields = g96Torestricted(fields[1:])
                newfields = fields[:3] + [str(n_mols), state] + fields[3:]
                mbmulti_angles.append(newfields)
        else:
            type_list = [fields[1:][3] for fields in value]
            angle_list = [float(fields[1:][4]) for fields in value]
            diff_angle = abs(max(angle_list) - min(angle_list))
            if not SameList(type_list) or diff_angle > cutoff:
                value = [[fields[0]] + g96Torestricted(fields[1:]) for fields in value]
                type_list = [fields[1:][3] for fields in value]
                k_list = [float(fields[1:][5]) for fields in value]
            else:
                k_list = [float(fields[1:][5]) for fields in value]

            if diff_angle <= cutoff:
                mean_angle = round(sum(angle_list) / len(angle_list), 2)
                mean_k = round(sum(k_list) / len(k_list), 2)
                mbangles.append([str(key[0]), str(key[1]), str(key[2]), type_list[0], str(mean_angle), str(mean_k)])
            else:
                for fields in value:
                    state = fields[0]
                    fields = g96Torestricted(fields[1:])
                    newfields = fields[:3] + [str(n_mols), state] + fields[3:]
                    mbmulti_angles.append(newfields)
    return mbangles, mbmulti_angles


def combine_dihedrals(
    n_mols: int,
    mols_dihedrals_list: list[list[list[str]]],
    cutoff: float
) -> tuple[list[list[str]], list[list[str]]]:
    """Combine dihedrals from multiple molecular states.

    Only supports periodic dihedrals (type 1).

    Args:
        n_mols: Number of molecules being combined.
        mols_dihedrals_list: List of dihedral lists from each molecule.
        cutoff: Dihedral cutoff in degrees for considering dihedrals as similar.

    Returns:
        Tuple of (mbdihedrals, mbmulti_dihedrals) where mbdihedrals are
        dihedrals with averaged parameters and mbmulti_dihedrals are
        state-specific dihedrals.
    """
    mols_dihedrals_dict_list: list[dict[tuple[int, int, int, int], list[str]]] = []
    assert len(mols_dihedrals_list) == n_mols, 'The number of molecules is not equal to the number of angles'
    for i, dihedrals in enumerate(mols_dihedrals_list):
        mols_dihedrals_dict: dict[tuple[int, int, int, int], list[str]] = {}
        state = str(i + 1)
        for fields in dihedrals:
            assert fields[4] in ['1'], f"Error: dihedral type is not 1. {fields}"
            assert fields[7] == '1', f"Error: dihedral n is not 1. {fields}"
            assert -180 < float(fields[5]) <= 180, f"Error: dihedrals should be in -180 -- +180. {fields}"
            if int(fields[0]) > int(fields[3]):
                fields[:4] = [fields[3], fields[2], fields[1], fields[0]]
            key = (int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3]))
            if key not in mols_dihedrals_dict:
                mols_dihedrals_dict[key] = [state] + fields
        mols_dihedrals_dict_list.append(mols_dihedrals_dict)
    mols_dihedrals_dict_combined = CombineDict(mols_dihedrals_dict_list)

    mbdihedrals: list[list[str]] = []
    mbmulti_dihedrals: list[list[str]] = []
    for key, value in mols_dihedrals_dict_combined.items():
        n_states_in_value = len({fields[0] for fields in value})
        assert n_states_in_value == len(value), f"Error: one state has more than one dihedral for same atoms! {value}"
        if n_states_in_value != n_mols:
            for fields in value:
                state = fields[0]
                fields = fields[1:]
                newfields = fields[:4] + [str(n_mols), state] + fields[4:]
                mbmulti_dihedrals.append(newfields)
        else:
            dihedral_list = [float(fields[1:][5]) for fields in value]
            diff_dihedral, mean_dihedral = Calculate_DiffDihedral(dihedral_list)
            k_list = [float(fields[1:][6]) for fields in value]

            if diff_dihedral <= cutoff:
                mean_dihedral_rounded = round(mean_dihedral, 2)
                mean_k = round(sum(k_list) / len(k_list), 2)
                mbdihedrals.append([str(key[0]), str(key[1]), str(key[2]), str(key[3]),
                                    '1', str(mean_dihedral_rounded), str(mean_k), '1'])
            else:
                for fields in value:
                    state = fields[0]
                    fields = fields[1:]
                    newfields = fields[:4] + [str(n_mols), state] + fields[4:]
                    mbmulti_dihedrals.append(newfields)
    return mbdihedrals, mbmulti_dihedrals
