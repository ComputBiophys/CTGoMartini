import re

def WriteItp(mol, Itp_file='default',Params_file='default'):
    """
    mol: molecule type.
    Itp_file: default -> f'{molecule_name}.itp'
    Params_file: default -> f'{molecule_name}_params.itp'
    """
    # Output
    molecule_name = mol._topology['moleculetype'][0][0]
    if Itp_file == 'default':
        Itp_file = f'{molecule_name}.itp'
    if Params_file == 'default':
        Params_file = f'{molecule_name}_params.itp'     
    
    ## Itp_file ##
    Itp_lines = []
    # notionlines
    Itp_lines += ['; This file was generate for the multiple basin Go-Martini.\n', '\n']
    # Categories
    for category, fields_list in mol._topology.items():
        if fields_list == []:
            continue

        Itp_lines += [f'\n[ {category} ]\n']
        if category == 'multi_angles': Itp_lines += [';atomid1, atomid2, atomid3, n_states, state, functype, angle, k\n']
        if category == 'multi_dihedrals': Itp_lines += [';atomid1, atomid2, atomid3, atomid4, n_states, state, functype, phi, k, n\n']
        if category == 'multi_contacts': Itp_lines += [';atomid1, atomid2, n_states, state, functype, sigma, epsilon\n']

        for fields in fields_list:
            line = ' '.join(fields)
            line = ' ' + line +'\n'
            Itp_lines += line 
        
    with open(Itp_file, 'w') as g:
        g.writelines(Itp_lines)

    ## Params_lines ##
    Params_lines = []
    # atomtypes
    atoms_list = mol._topology['atoms']
    new_atomtypes = [item[1] for item in atoms_list if re.fullmatch(rf'{molecule_name}_\d+', item[1])]
    if new_atomtypes != []:
        new_atomtypes_list=[[atomtype, '0.0', '0.0','A','0.0','0.0'] for atomtype in new_atomtypes]
        Params_lines+=['\n[ atomtypes ]\n']
        for fields in new_atomtypes_list:
            line = ' '.join(fields)
            line = ' ' + line +'\n'
            Params_lines += line         

    with open(Params_file, 'w') as g:
        g.writelines(Params_lines)    
            
    