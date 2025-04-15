import os
import re
from collections import OrderedDict

def Extract_contacts_from_top(top, molecule_name):
    atoms = top.moleculeTypes[molecule_name]._topology['atoms']
    nonbond_params = top.forcefield._parameters['nonbond_params']

    pattern=rf'{molecule_name}_\d+'
    contact_atoms = {}  # atomtype: atomid # hard restraint
    for item in atoms:
        atomtype = item[1]
        if re.fullmatch(pattern, atomtype):
            atomid=int(item[0])
            contact_atoms[atomtype]=atomid
    
    contacts = []
    for fields in nonbond_params:
        atomtype1 = fields[0]
        atomtype2 = fields[1]
        judge1=bool(re.fullmatch(pattern, atomtype1))
        judge2=bool(re.fullmatch(pattern, atomtype2))
        if judge1 ^ judge2 :
            raise Exception(f"Error: Unsupport the contact bewteen {key}!")
        elif judge1 & judge2:
            try:
                atomid1=contact_atoms[atomtype1]
                atomid2=contact_atoms[atomtype2]
            except:
                raise Exception(f"Error: not contact_atoms! {fields}")
            
            atomid1, atomid2 = sorted([atomid1, atomid2])
            assert fields[2] == '1', f"Error: only support functype 1: {fields}"
            newfields=[str(atomid1), str(atomid2)]+fields[2:]
            contacts.append(newfields)
        else:
            continue    
        
    contacts = sorted(contacts, key=lambda fields: (int(fields[0]), int(fields[1])))
    return contacts

def GetAtomNames(atomidlist, atoms):
    atoms_dict = {int(fields[0]): fields for fields in atoms}
    
    atomnamelist = []
    for atomid in atomidlist:
        fields = atoms_dict[int(atomid)]
        atomname = fields[4]
        atomnamelist.append(atomname)
    return atomnamelist    

def GetAngleDiehdralType(atomnamelist):
    typestr = ''
    for atomname in atomnamelist:
        if atomname == 'BB':
            typestr += 'B'
        elif atomname.startswith('SC'):
            typestr += 'S'
        else:
            raise ValueError(f"Error: Unsupport atomnames other than BB and SC*: {atomnamelist}")
    return typestr

def DifferentiateAngles(angles, atoms):
    """
    Get the BBB angles and notBBB angles

    Parameters
    #########
    angles: list
        list of angle fields
    atoms: list
        list of atom fields

    Return
    ######
    BBB_angles, notBBB_angles
    """
    BBB_angles = []
    notBBB_angles = []

    for fields in angles:
        atomidlist = fields[:3]
        atomnamelist = GetAtomNames(atomidlist, atoms)
        if GetAngleDiehdralType(atomnamelist) == 'BBB':
            BBB_angles.append(fields)
        else:
            notBBB_angles.append(fields)
    return BBB_angles, notBBB_angles

def DifferentiateDihedrals(dihedrals, atoms):
    """
    Get the BBBB dihedrals, SSSS_dihedrals and SBBS dihedrals

    Parameters
    #########
    angles: list
        list of dihedrals fields
    atoms: list
        list of atom fields

    Return
    ######
    BBBB_dihedrals, SSSS_dihedrals, SBBS_dihedrals
    """
    BBBB_dihedrals = []
    SSSS_dihedrals = []
    SBBS_dihedrals = []
    other_dihedrals = []

    for fields in dihedrals:
        atomidlist = fields[:4]
        atomnamelist = GetAtomNames(atomidlist, atoms)
        if GetAngleDiehdralType(atomnamelist) == 'BBBB':
            BBBB_dihedrals.append(fields)
        elif GetAngleDiehdralType(atomnamelist) == 'SBBS':
            SBBS_dihedrals.append(fields)
        elif GetAngleDiehdralType(atomnamelist) == 'SSSS':
            SSSS_dihedrals.append(fields)            
        else:
            other_dihedrals.append(fields)
    assert other_dihedrals == [], f'Error: not supported dihedral type: {other_dihedrals}'
    return BBBB_dihedrals, SSSS_dihedrals, SBBS_dihedrals

def CombineDict(dict_list):
    n_dict = len(dict_list)
    key_combined_list = []
    for i in range(n_dict):
        key_combined_list += list(dict_list[i].keys())
    key_combined_list = list(set(key_combined_list))
    key_combined_list = sorted(key_combined_list)
    
    dict_combined = OrderedDict()
    for key in key_combined_list:
        dict_combined[key] = []
        for i in range(n_dict):
            if key in dict_list[i].keys():
                dict_combined[key].append(dict_list[i][key])
    return dict_combined

def ForceItemFloat(item, precision=None):
    try:
        item = float(item)
        if precision is not None:
            item = round(item, precision)
    except:
        pass
    return item

def ForceListFloat(itemlist, precision=None):
    newlist = []
    for item in itemlist:
        if type(item) is not list:
            newlist.append(ForceItemFloat(item, precision))
        else:
            newlist.append(ForceListFloat(item, precision))
    return newlist

def SameListList(listlist, typeforce=True, sort=False, precision=None):
    """
    Judge whether the lists in the list are the same
    
    Parameters
    #####
    listlist: list(list) or list(list(list*))
        list of some lists
    typefoce: bool, True
    precision: int or None
        number of decimal places to consider (None for exact comparison)

    Return
    ######
    True or False
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


def SameList(alist, typeforce=True):
    """
    The items in the list are the same.
    """
    issame = True
    alistlist = [[item] for item in alist]
    issame = SameListList(alistlist, typeforce)
    return issame


def Calculate_DiffDihedral(dihedral_list):
    """
    when follow the clocklike direction or anticlocklike direction, the max difference between dihedrals shold be less than 180
    """
    
    dihedral_list = sorted(dihedral_list)
    # print(dihedral_list)
    anticlock_dihedral_list = []
    for i, dihedral in enumerate(dihedral_list):
        if dihedral < 0:
            dihedral += 360
        anticlock_dihedral_list.append(dihedral)
    anticlock_dihedral_list = sorted(anticlock_dihedral_list)
    anticlock_diff_max = abs(anticlock_dihedral_list[-1] - anticlock_dihedral_list[0])
    # print(anticlock_dihedral_list)

    clock_dihedral_list = dihedral_list.copy()
    clock_dihedral_list = sorted(clock_dihedral_list)
    clock_diff_max = abs(clock_dihedral_list[-1] - clock_dihedral_list[0])
    # print(clock_dihedral_list)

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
    mean_dihedral = sum(dihedral_list)/len(dihedral_list)
    if mean_dihedral > 180: mean_dihedral -= 360
    if mean_dihedral <= -180: mean_dihedral += 360
    return diff_max, mean_dihedral


class CombineMols:
    @staticmethod
    def combine_atoms(mbmolname: str, mols_atoms_pairs: list):
        """
        Combine atoms from different states of molecules

        Parameters
        ##########
        mbmolname: str,
            New prefix of atomnames of virtual sites

        mols_atoms_pairs: list
            Atoms from different states of molecules
            [(molnameA, atomtopA), (molnameB, atomtopB), ...]

        Return
        ######
        atomtop: list
            atomtop
        """
        n_mols = len(mols_atoms_pairs)
        n_atoms = len(mols_atoms_pairs[0][1])
        
        # Assert that there are not less than 2 mols
        assert n_mols >= 2, "Error: The number of mols must more than or equal 2"

        # Assert that differet molecules have same number of atoms
        for pair in mols_atoms_pairs[1:]:
            assert len(pair[1]) == len(mols_atoms_pairs[0][1])
        
        mols_atoms_dict_list = []
        for pair in mols_atoms_pairs:
            atoms = pair[1]
            mols_atoms_dict_list.append({(int(atoms[i][0]),): atoms[i] for i in range(n_atoms)})
        mols_atoms_dict_combined = CombineDict(mols_atoms_dict_list)

        mbmol_atomtop = []
        Extract=lambda atomtype: re.findall(r'^(\w+)_(\d+)$',atomtype)[0]
        for key, value in mols_atoms_dict_combined.items():
            assert len(value) == n_mols, f"Error: The number of molecules with the same atomid is not equal to the number of molecules. {key}"
            if SameListList(value):
                mbmol_atomtop.append(value[0].copy())
            else:
                try:
                    mol_resid_extract_list = []
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

                except:
                    raise ValueError("Error: atoms from different states of one molecule cannot meet the combination rule!", mols_atom_list)
        
        assert len(mbmol_atomtop) == n_atoms
        return mbmol_atomtop    

    @staticmethod
    def combine_bonds_constraints(n_mols, mols_bonds_list, mols_constraints_list):
        """
        """

        mols_connections_dict_list = []
        for i, bonds in enumerate(mols_bonds_list):
            connection_dict = {}
            state = str(i+1)
            for bond in bonds:
                assert bond[2] == "1", f"Error: bond type must be 1: {bond}"
                if int(bond[0]) > int(bond[1]):
                    bond[:2] = [bond[1], bond[0]]
                key = tuple([int(bond[0]), int(bond[1])])
                connection_dict[key] = [state] + bond
            mols_connections_dict_list.append(connection_dict)
        
        for i, constraints in enumerate(mols_constraints_list):
            connection_dict = {}
            state = str(i+1)
            for constraint in constraints:
                assert constraint[2] == "1", f"Error: constraint type must be 1: {constraint}"
                if int(constraint[0]) > int(constraint[1]):
                    constraint[:2] = [constraint[1], constraint[0]]
                key = tuple([int(constraint[0]), int(constraint[1])])
                connection_dict[key] = [state] + constraint[:4] + [None]
            mols_connections_dict_list.append(connection_dict)
        
        mols_connections_dict_combined = CombineDict(mols_connections_dict_list)

        mbconnections = []
        for key, value in mols_connections_dict_combined.items():
            n_states_in_value = len(set([fields[0] for fields in value]))
            assert n_states_in_value == len(value), f"Error: value repeats! {value}"
            assert n_states_in_value == n_mols, f"Error: {key} does not have {n_mols} values"
            dist_list = [float(fields[1:][3]) for fields in value]
            k_list = [float(fields[1:][4]) for fields in value if fields[1:][4] is not None]
            dist_mean = sum(dist_list) / len(dist_list)
            dist_mean = str(round(dist_mean, 3))
            if k_list != []:
                k_mean = sum(k_list) / len(k_list)
                k_mean = str(round(k_mean, 3))
            else:
                k_mean = None
            
            mbconnections.append(value[0][1:][:3] + [dist_mean, k_mean])
        
        mbbonds = []
        mbconstraints = []
        for item in mbconnections:
            if item[4] is not None:
                mbbonds.append(item)
            else:
                mbconstraints.append(item[:4])

        return mbbonds, mbconstraints



    @staticmethod
    def combine_exclusions(mols_exclusions_list):
        """
        """

        exclusion_pair_list = []
        for exclusions in mols_exclusions_list:
            for fields in exclusions:
                item0 = fields[0]
                for item in fields[1:]:
                    exclusion_pair_list.append(tuple(sorted([int(item0), int(item)])))
        exclusion_pair_list = sorted(list(set(exclusion_pair_list)))
        
        mbexclusion_dict = OrderedDict()
        for item in exclusion_pair_list:
            key=(item[0],)
            if key not in mbexclusion_dict: 
                mbexclusion_dict[key]=[str(item[0])]
            if item[1] not in mbexclusion_dict[key]:
                mbexclusion_dict[key].append(str(item[1]))

        mbexclusions = list(mbexclusion_dict.values())
        return mbexclusions

    @staticmethod
    def combine_contacts(n_mols, mols_contacts_list, cutoff):
        """
        """
        assert len(mols_contacts_list) == n_mols, f"Error: The number of contacts is not equal to the number of molecules."
        mols_contacts_dict_list = []
        for i, contacts in enumerate(mols_contacts_list):
            mols_contacts_dict = {}
            state = str(i+1)
            for fields in contacts:
                assert fields[2] == '1', f"Error: contact type is not 1. {fields}"
                if int(fields[0]) > int(fields[1]):
                    fields[:2] = [fields[1], fields[0]]
                key = tuple([int(fields[0]), int(fields[1])])
                if key not in mols_contacts_dict:
                    mols_contacts_dict[key] = [state] + fields
            mols_contacts_dict_list.append(mols_contacts_dict)
        mols_contacts_dict_combined = CombineDict(mols_contacts_dict_list)

        mbcontacts = []
        mbmulti_contacts = []
        for key, value in mols_contacts_dict_combined.items():
            n_states_in_value = len(set([fields[0] for fields in value]))
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
        

    @staticmethod
    def combine_angles(n_mols, mols_angles_list, cutoff):
        """
        Combine angles from different states of molecules
        Convert the angle type 2 (g96 angles) to type 10 (restricted angles) if the angles of the same atoms from different states have different types
        
        Parameters
        ##########
        """

        mols_angles_dict_list = []
        assert len(mols_angles_list) == n_mols, 'The number of molecules is not equal to the number of angles'
        for i, angles in enumerate(mols_angles_list):
            mols_angles_dict = {}
            state = str(i+1)
            for fields in angles:
                assert fields[3] in ['2', '10'], f"Error: angle type is not 2 or 10. {fields}"
                assert float(fields[4]) >=0 and float(fields[4]) <=180, f"Error: angles should be in 0-180. {fields}"
                if int(fields[0]) > int(fields[2]):
                    fields[:3] = [fields[2], fields[1], fields[0]]
                key = tuple([int(fields[0]), int(fields[1]), int(fields[2])])
                if key not in mols_angles_dict:
                    mols_angles_dict[key] = [state] + fields
            mols_angles_dict_list.append(mols_angles_dict)
        mols_angles_dict_combined = CombineDict(mols_angles_dict_list)

        # g96_angles only for those that the angle difference is smaller than cutoff and the inital types are g96 angles (2).
        # Others should be transformed into restricted angles.
        # type_g96_angles='2'
        # type_restricted_angles='10'
        def g96Torestricted(fields):
            if fields[3] == '2':
                newfields=fields[:3]+['10',fields[4],'25.0']
            elif fields[3] == '10':
                newfields=fields.copy()
            else:
                raise ValueError(f"Error: Unsupport angle type. {fields}")
            return newfields

            
        mbangles = []
        mbmulti_angles = []
        for key, value in mols_angles_dict_combined.items():
            n_states_in_value = len(set([fields[0] for fields in value]))
            assert n_states_in_value == len(value), f"Error: angle value repeats! {value}"
            if  n_states_in_value != n_mols:
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
        

    @staticmethod
    def combine_dihedrals(n_mols, mols_dihedrals_list, cutoff):
        """
        Combine dihedrals from different states of molecules
        assert periodic dihedrals

        Parameters
        ##########
        """

        mols_dihedrals_dict_list = []
        assert len(mols_dihedrals_list) == n_mols, 'The number of molecules is not equal to the number of angles'
        for i, dihedrals in enumerate(mols_dihedrals_list):
            mols_dihedrals_dict = {}
            state = str(i+1)
            for fields in dihedrals:
                assert fields[4] in ['1'], f"Error: dihedral type is not 1. {fields}"
                assert fields[7] == '1', f"Error: dihedral n is not 1. {fields}"
                assert float(fields[5]) > -180 and float(fields[5]) <=180, f"Error: dihedrals should be in -180 -- +180. {fields}"
                if int(fields[0]) > int(fields[3]):
                    fields[:4] = [fields[3], fields[2], fields[1], fields[0]]
                key = tuple([int(fields[0]), int(fields[1]), int(fields[2]), int(fields[3])])
                if key not in mols_dihedrals_dict:
                    mols_dihedrals_dict[key] = [state] + fields
            mols_dihedrals_dict_list.append(mols_dihedrals_dict)
        mols_dihedrals_dict_combined = CombineDict(mols_dihedrals_dict_list)

        mbdihdedrals = []
        mbmulti_dihedrals = []
        for key, value in mols_dihedrals_dict_combined.items():
            n_states_in_value = len(set([fields[0] for fields in value]))
            assert n_states_in_value == len(value), f"Error: one state has more than one dihedral for same atoms! {value}"
            if  n_states_in_value != n_mols:
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
                    mean_dihedral = round(mean_dihedral, 2)
                    mean_k = round(sum(k_list) / len(k_list), 2)
                    mbdihdedrals.append([str(key[0]), str(key[1]), str(key[2]), str(key[3]), '1', str(mean_dihedral), str(mean_k), '1'])
                else:
                    for fields in value:
                        state = fields[0]
                        fields = fields[1:]
                        newfields = fields[:4] + [str(n_mols), state] + fields[4:]
                        mbmulti_dihedrals.append(newfields)
        return mbdihdedrals, mbmulti_dihedrals


