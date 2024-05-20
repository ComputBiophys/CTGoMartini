from ..util import Extract_contacts_from_top, CombineMols, DifferentiateAngles, DifferentiateDihedrals, SameListList
from .MartiniTopology import MartiniTopFile

def GenMBPTop(mols_list, mbmol_name, dict_cutoffs={
    'cutoff_BBB_angles':0.0,   
    'cutoff_BBBB_dihedrals': 0.0,
    'cutoff_SBBS_dihedrals':0.0,
    'cutoff_contacts': 0.0,  
}):
    """
    mols_list: [[topfileA, mol_nameA], [topfileB, mol_nameB], ...]  
    mbmol_name: string
    """
    # scFix influence on the MBCombination
    # # angles:
    # BBB_angles, BBS_regular_angles, SBB_regular_angles, BSS_angles, SBB_scFix_angles, BBS_scFix_angles, 
    # BBB_angles, not_BBB_angles
    # # dihedrals
    # BBBB_dihedrals, SSSS_dihedrals, SBBS_scFix_dihedrals

    # -------------------------------------------
    n_mols = len(mols_list)

    # add contacts, bonds, constraints, angles, dihedrals, exclusions
    # for the Go-Martini model, the categories including atoms, virtual_sitesn must be comprised in the topology.
    mols_pairs = []
    for (topfile, mol_name) in mols_list:
        top = MartiniTopFile(topfile)
        mol = top.moleculeTypes[mol_name]
        if 'contacts' not in mol._topology:
            contacts = Extract_contacts_from_top(top, mol_name)
            if contacts != []:
                mol._topology['contacts'] = contacts        
        if 'bonds' not in mol._topology:
            mol._toppology['bonds'] = []
        if 'constraints' not in mol._topology:
            mol._toppology['constraints'] = []
        if 'angles' not in mol._topology:
            mol._toppology['angles'] = []
        if 'dihedrals' not in mol._topology:
            mol._toppology['dihedrals'] = []
        if 'exclusions' not in mol._topology:
            mol._toppology['exclusions'] = []
        mols_pairs.append([mol_name, mol])

    # -------------------------------------------
    class MBMolecule:
        def __init__(self, mbmol_name) -> None:
            self.mbmol_name = mbmol_name
            self._topology = {}
    
    mbmol = MBMolecule(mbmol_name)
    # add moleculetype
    mbmol._topology['moleculetype'] = [[mbmol_name, '1']]
    # add multiple_basin
    energy_basin_list = [f'C{i+1}' for i in range(n_mols)]
    mbmol._topology['multiple_basin'] = [['True', 'exp', str(n_mols), 'beta'] + energy_basin_list]
    
    # combine atoms, bonds, constraints, exclusions, contacts, BBB_angles, BBBB_dihedrals
    # add mbatoms
    mols_atoms_pairs = [[mol_name, mol._topology['atoms']] for (mol_name, mol) in mols_pairs]
    mbatoms = CombineMols.combine_atoms(mbmol_name, mols_atoms_pairs)
    mbmol._topology['atoms'] = mbatoms

    # add mbbonds and mbconstraints
    mols_bonds_list = [mol._topology['bonds'] for (_, mol) in mols_pairs]
    mols_constraints_list = [mol._topology['constraints'] for (_, mol) in mols_pairs]
    mbbonds, mbconstraints = CombineMols.combine_bonds_constraints(n_mols, mols_bonds_list, mols_constraints_list)
    mbmol._topology['bonds'] = mbbonds
    mbmol._topology['constraints'] = mbconstraints

    # add mbangles, mbmulti_angles
    mols_BBB_angles_list, mols_notBBB_angles_list = [], []
    for _, mol in mols_pairs:
        mol_BBB_angles, mol_notBBB_angles = DifferentiateAngles(mol._topology['angles'], mol._topology['atoms'])
        mols_BBB_angles_list.append(mol_BBB_angles)
        mols_notBBB_angles_list.append(mol_notBBB_angles)

    mb_BBB_angles, mbmulti_BBB_angles = CombineMols.combine_angles(n_mols, mols_BBB_angles_list, dict_cutoffs['cutoff_BBB_angles'])
    assert SameListList(mols_notBBB_angles_list, sort=True), f'Error: notBBB angles from different mols should be the same! {mols_notBBB_angles_list}'
    mbangles = mb_BBB_angles + mols_notBBB_angles_list[0]
    mbmol._topology['angles'] = mbangles
    mbmol._topology['multi_angles'] = mbmulti_BBB_angles

    # add mbdihedrals, mbmulti_dihedrals
    mols_BBBB_dihedrals_list, mols_SSSS_dihedrals_list, mols_SBBS_dihedrals_list = [], [], []
    for _, mol in mols_pairs:
        mols_BBBB_dihedrals, mols_SSSS_dihedrals, mols_SBBS_dihedrals = DifferentiateDihedrals(mol._topology['dihedrals'], mol._topology['atoms'])
        mols_BBBB_dihedrals_list.append(mols_BBBB_dihedrals)
        mols_SSSS_dihedrals_list.append(mols_SSSS_dihedrals)
        mols_SBBS_dihedrals_list.append(mols_SBBS_dihedrals)

    mb_BBBB_dihdedrals, mbmulti_BBBB_dihedrals = CombineMols.combine_dihedrals(n_mols, mols_BBBB_dihedrals_list, dict_cutoffs['cutoff_BBBB_dihedrals']) 
    assert SameListList(mols_SSSS_dihedrals_list, sort=True), f'Error: SSSS dihedrals from different mols should be the same! {mols_SSSS_dihedrals_list}'
    assert SameListList(mols_SBBS_dihedrals_list, sort=True), f'Error: Dont support the scFix for multiple-basin Go-Martini. SBBS dihedrals from different mols should be the same! {mols_SBBS_dihedrals_list}'
    mbdihedrals = mb_BBBB_dihdedrals + mols_SSSS_dihedrals_list[0] + mols_SBBS_dihedrals_list[0]
    mbmulti_dihedrals = mbmulti_BBBB_dihedrals
    mbmol._topology['dihedrals'] = mbdihedrals
    mbmol._topology['multi_dihedrals'] = mbmulti_dihedrals    

    # add mbcontacts, mbmulti_contacts
    mols_contacts_list = [mol._topology['contacts'] for (_, mol) in mols_pairs]
    mbcontacts, mbmulti_contacts = CombineMols.combine_contacts(n_mols, mols_contacts_list, dict_cutoffs['cutoff_contacts'])
    mbmol._topology['contacts'] = mbcontacts
    mbmol._topology['multi_contacts'] = mbmulti_contacts        

    # add mbexclusions
    mols_exclusions_list = [mol._topology['exclusions'] for (_, mol) in mols_pairs]
    mbexclusions = CombineMols.combine_exclusions(mols_exclusions_list)
    mbmol._topology['exclusions'] = mbexclusions

    # other categories should be the same for different mols
    categories_list = []
    for (_, mol) in mols_pairs:
        categories_list += list(mol._topology.keys())
    categories_list= list(set(categories_list))
    for processed_category in ['moleculetype', 'atoms', 'bonds', 'constraints', 'angles', 'dihedrals', 'contacts', 'exclusions']:
        categories_list.remove(processed_category)
    for category in categories_list:
        mols_category_list = [mol._topology[category] for (_, mol) in mols_pairs]
        assert SameListList(mols_category_list, sort=True), f'Error: {category} from different mols should be the same! {mols_category_list}'
        mbmol._topology[category] = mols_category_list[0]

    # remove the category including nothing
    mbcategories_list = list(mbmol._topology.keys())
    for category in mbcategories_list:
        if mbmol._topology[category] == []:
            mbmol._topology.pop(category)

    return mbmol
