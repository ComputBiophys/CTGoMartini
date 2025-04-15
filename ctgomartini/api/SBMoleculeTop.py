from ..util import Extract_contacts_from_top, CombineMols, DifferentiateAngles, DifferentiateDihedrals, SameListList
from .MartiniTopology import MartiniTopFile


def ChangeAtomType(atoms, mol_name, sbmol_name):
    import re
    Extract=lambda atomtype: re.findall(r'^(\w+)_(\d+)$',atomtype)[0]

    newatoms = []
    for field in atoms:
        atomtype = field[1]
        try:
            mol_name_extract, mol_resid_extract = Extract(atomtype)
            # print(mol_name_extract, mol_resid_extract)
            if mol_name_extract != mol_name:
                print("mol_name_extract != mol_name", atomtype, mol_name)
            newatomtype = sbmol_name + '_' + mol_resid_extract
        except:
            newatomtype = atomtype
        newatoms.append([field[0], newatomtype] + field[2:])
    return newatoms
    

def GenSBPTop(mols_list, sbmol_name):
    """
    mols_list: [[topfile, mol_name]]  
    sbmol_name: string
    """
    # scFix influence on the SBCombination
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
    class SBMolecule:
        def __init__(self, sbmol_name) -> None:
            self.sbmol_name = sbmol_name
            self._topology = {}
    
    sbmol = SBMolecule(sbmol_name)
    # add moleculetype
    sbmol._topology['moleculetype'] = [[sbmol_name, '1']]
    
    mol_name, mol = mols_pairs[0]
    # atoms, bonds, constraints, exclusions, contacts, BBB_angles, BBBB_dihedrals
    # add sbatoms
    sbmol._topology['atoms'] = ChangeAtomType(mol._topology['atoms'], mol_name, sbmol_name)

    # add sbbonds and sbconstraints
    sbmol._topology['bonds'] = mol._topology['bonds']
    sbmol._topology['constraints'] = mol._topology['constraints']

    # add sbangles
    sbmol._topology['angles'] = mol._topology['angles']

    # add sbdihedrals
    sbmol._topology['dihedrals'] = mol._topology['dihedrals']

    # add sbcontacts
    sbmol._topology['contacts'] = mol._topology['contacts']

    # add mbexclusions
    sbmol._topology['exclusions'] = mol._topology['exclusions']

    # other categories should be the same for different mols
    categories_list = []
    for (_, mol) in mols_pairs:
        categories_list += list(mol._topology.keys())
    categories_list= list(set(categories_list))
    for processed_category in ['moleculetype', 'atoms', 'bonds', 'constraints', 'angles', 'dihedrals', 'contacts', 'exclusions']:
        categories_list.remove(processed_category)
    for category in categories_list:
        mols_category_list = [mol._topology[category] for (_, mol) in mols_pairs]
        sbmol._topology[category] = mols_category_list[0]

    # remove the category including nothing
    sbcategories_list = list(sbmol._topology.keys())
    for category in sbcategories_list:
        if sbmol._topology[category] == []:
            sbmol._topology.pop(category)

    return sbmol
