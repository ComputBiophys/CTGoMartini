#!/usr/bin/env python3

import os
import argparse
import subprocess
import MDAnalysis as mda
import ctgomartini
from ctgomartini.api import GenMBPTop
from ctgomartini.func import WriteItp, ConvertLongShortElasticBonds
from ctgomartini.func import Create_goVirt_for_multimer

def Martinize2(aa_strfile, dssp, ff, state_name, other_params=''):
    output = subprocess.run(f'martinize2 -f {aa_strfile} -o system.top -x {state_name}_cg.pdb -dssp {dssp} -p backbone -ff {ff} -govs-include -govs-moltype {state_name}  -cys auto {other_params}',
                   shell=True, capture_output=True, encoding='utf-8')
    print(output.args)
    if output.returncode != 0:
        stdout = output.stdout
        stderr = output.stderr
        print(f"Error: Something wrong with {output.args}")
        print()
        print("stdout:")
        print(stdout)
        print("stderr")
        print(stderr)
        raise Exception(f'Error! {os.getcwd()}')

def GetNatoms(cg_strfile, atomname='CA'):
    u = mda.Universe(cg_strfile)
    Natoms =  u.select_atoms(f'name {atomname}')[0].index
    return Natoms

def GenGoContacts(aa_strfile, cg_strfile, aa_map, state_name, go_eps=12):
    import contextlib, io

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        with contextlib.redirect_stderr(f):
            Natoms = GetNatoms(cg_strfile, atomname='CA')
            Create_goVirt_for_multimer(f'-r {aa_strfile} -s {cg_strfile} -f {aa_map} --moltype {state_name} --go_eps {go_eps} --Natoms {Natoms}')
    output = f.getvalue()

    print(f'Create_goVirt_for_multimer -r {aa_strfile} -s {cg_strfile} -f {aa_map} --moltype {state_name} --go_eps {go_eps} --Natoms {Natoms}')
    output = output.split('\n')
    for line in output:
        if line.strip().startswith('Only symmetric OV + rCSU contacts (singly counted):'):
            print(line)

def ModifyFF(forcefield_file='martini_v3.0.0.itp'):
    output = subprocess.run(r'''
sed -i "s/\[ nonbond_params \]/\#ifdef GO_VIRT\n\#include \"BB-part-def_VirtGoSites.itp\"\n\#endif\n\n\[ nonbond_params \]/" {}
echo "\n#ifdef GO_VIRT \n#include \"go-table_VirtGoSites.itp\"\n#endif" >> {}
sed -i 's/#include "martini.itp"/#include "{}"/g' system.top
'''.format(forcefield_file, forcefield_file, forcefield_file),
    shell=True, capture_output=True, encoding='utf-8')
    # print(output.args)
    if output.returncode != 0:
        stdout = output.stdout
        stderr = output.stderr
        print(f"Error: Something wrong with {output.args}")
        print()
        print("stdout:")
        print(stdout)
        print("stderr")
        print(stderr)
        raise Exception(f'Error! {os.getcwd()}')

    
def MBGOMartinize(aa_strfile_list, aa_map_list, state_name_list, mbmol_name, dict_cutoffs, method='exp', dssp='dssp', ff='martini3001'):
    working_path = os.getcwd()
    print(f'Working path: {working_path}')
    
    for aa_strfile, aa_map, state_name in zip(aa_strfile_list, aa_map_list, state_name_list):
        os.chdir(working_path)
        if os.path.exists(f'./{state_name}'):
            raise ValueError(f'Error: Directory {state_name} exists!')
            # subprocess.run(f'rm {state_name} -r', shell=True)
            # os.mkdir(state_name)
            # pass
        else:
            os.mkdir(state_name)

        os.chdir(os.path.join(working_path, state_name))
        
        print('\n############')
        print('Subworking_dir:', os.getcwd())
        
        # Martinize AA Proteins
        print(f'\nMartinize the all-atom protein ({aa_strfile}) as the CG model with the state name ({state_name})')
        Martinize2(os.path.join('../', aa_strfile), dssp, ff, state_name)

        # Generate Go-Contacts
        print(f'\nGenerate the Go-Contacts for proteins')
        GenGoContacts(os.path.join('../', aa_strfile), f'{state_name}_cg.pdb', os.path.join('../', aa_map), state_name, go_eps=12)

        # Fetch the FF file
        print(f'\nFetch the forcefield and append the Go-Contacts to the forcefields')
        assert ff == 'martini3001', f'Error: Unsupport the forcefield: {ff}'
        os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/ForceFields/Martini300/martini_v3.0.0.itp')} .")
        ModifyFF(forcefield_file='martini_v3.0.0.itp')

        # Convert Long/Short Elastic Bonds to LJ Interactions
        print('\nConvert Long/Short Elastic Bonds to LJ Interactions')
        ConvertLongShortElasticBonds(state_name, f'{state_name}_cg.pdb', convertLongElasticBonds=True, convertShortElasticBonds=False, LJ_epsilon=12)
        
    print('############')
    os.chdir(working_path)

    # Combine multiple states into the multiple-basin potential
    print(f'\nGenerate the muliple-basin potential for {mbmol_name}')
    mols_list = []
    for i, state_name in enumerate(state_name_list):
        mols_list.append([f'{state_name}/system.top', state_name])
    mbmol =  GenMBPTop(mols_list, mbmol_name, dict_cutoffs)

    # Modify the mulitple_basin parameters according to method
    if method.lower() == "exp":
        pass
    elif method.lower() == "ham":
        assert mbmol._topology['multiple_basin'][0][2] == '2', f"Error: HAM mixing scheme only supports mulitple basins for two states.\n"
        mbmol._topology['multiple_basin'][0] = ['True', 'ham', '2', 'delta', 'mbp_energy1', 'mbp_energy2']

    # Write the mbmol.itp
    print(f'\nWrite the {mbmol_name}.itp and {mbmol_name}_params.itp')
    WriteItp(mbmol)
    print('Finish!')

def SwitchingGOMartinize(aa_strfile_list, aa_map_list, state_name_list, mbmol_name, dict_cutoffs, method='switching', dssp='dssp', ff='martini3001'):
    working_path = os.getcwd()
    print(f'Working path: {working_path}')
    
    for aa_strfile, aa_map, state_name in zip(aa_strfile_list, aa_map_list, state_name_list):
        os.chdir(working_path)
        if os.path.exists(f'./{state_name}'):
            # raise ValueError(f'Error: Directory {state_name} exists!')
            subprocess.run(f'rm {state_name} -r', shell=True)
            os.mkdir(state_name)
            pass
        else:
            os.mkdir(state_name)

        os.chdir(os.path.join(working_path, state_name))
        
        print('\n############')
        print('Subworking_dir:', os.getcwd())
        
        # Martinize AA Proteins
        print(f'\nMartinize the all-atom protein ({aa_strfile}) as the CG model with the state name ({state_name})')
        Martinize2(os.path.join('../', aa_strfile), dssp, ff, state_name, other_params='-scfix')

        # Generate Go-Contacts
        print(f'\nGenerate the Go-Contacts for proteins')
        GenGoContacts(os.path.join('../', aa_strfile), f'{state_name}_cg.pdb', os.path.join('../', aa_map), state_name, go_eps=12)

        # Fetch the FF file
        print(f'\nFetch the forcefield and append the Go-Contacts to the forcefields')
        assert ff == 'martini3001', f'Error: Unsupport the forcefield: {ff}'
        os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/ForceFields/Martini300/martini_v3.0.0.itp')} .")
        ModifyFF(forcefield_file='martini_v3.0.0.itp')
        
    print('############')
    os.chdir(working_path)
    print('Finish!')

def CTGOMartinize(aa_strfile_list, aa_map_list, state_name_list, mbmol_name, dict_cutoffs, method='switching', dssp='dssp', ff='martini3001'):
    if method.lower() == 'switching':
        SwitchingGOMartinize(aa_strfile_list, aa_map_list, state_name_list, mbmol_name, dict_cutoffs, method=method, dssp=dssp, ff=ff)
    elif method.lower() in ['exp', 'ham']:
        MBGOMartinize(aa_strfile_list, aa_map_list, state_name_list, mbmol_name, dict_cutoffs, method=method, dssp=dssp, ff=ff)
    else:
        raise ValueError(f'Error: unsupport the method named {method}!')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
Generate the topology files for the Multiple-baisn Go-Martini method or the Swithing Go-Martini method. 

An example:
# Multiple-basin Go-Martini
python ctgomartinize.py -s StateA_aa.pdb StateB_aa.pdb -m StateA_aa.map StateB_aa.map -mol StateA StateB -mbmol protein -dssp dssp -ff martini3001 -method exp

# Switching Go-Martini
python ctgomartinize.py -s StateA_aa.pdb StateB_aa.pdb -m StateA_aa.map StateB_aa.map -mol StateA StateB -dssp dssp -ff martini3001 -method switching                                     
""")
    parser.add_argument('-s', dest='strfile', required=True, nargs='+', type=str,
                        help='Input structure files')
    parser.add_argument('-m', dest='mapfile', required=True, nargs='+', type=str,
                        help='Input map files')
    parser.add_argument('-mol', dest='moltype', required=True, nargs='+', type=str,
                        help='Molecule type names')
    parser.add_argument('-mbmol', dest='mbmoltype', default='mbmol', type=str,
                        help='Mulitple-basin molecule type name (default: mbmol)')    
    parser.add_argument('-dssp', dest='dssp', default='dssp', type=str,
                        help='DSSP executable for determining structure (default: dssp)')
    parser.add_argument('-ff', dest='ff', default='martini3001', type=str,
                        help='forcefield to use (default: martini3001)\nNow only support martini3001!')
    parser.add_argument('-method', dest='method', required=True, type=str,
                        help='method to use (required: exp, ham, switching)')
    parser.add_argument('-cutoff_BBB_angles', dest='cutoff_BBB_angles', default=15.0, type=float,
                        help='Cutoff of BBB angles for generating the multiple-baisn Go-Martini topology (default: 15.0 degree)')
    parser.add_argument('-cutoff_BBBB_dihedrals', dest='cutoff_BBBB_dihedrals', default=30.0, type=float,
                        help='Cutoff of BBBB dihedrals for generating the multiple-baisn Go-Martini topology (default: 30.0 degree)')        
    parser.add_argument('-cutoff_SBBS_dihedrals', dest='cutoff_SBBS_dihedrals', default=30.0, type=float,
                        help='Cutoff of SBBS dihedrals for generating the multiple-baisn Go-Martini topology (default: 30.0 degree).\nNote that this parameter is useless now.')
    parser.add_argument('-cutoff_contacts', dest='cutoff_contacts', default=0.06, type=float,
                        help='Sigma cutoff of contacts for generating the multiple-baisn Go-Martini topology (default: 0.06 nm)')          
    
    args = parser.parse_args()
    # args = parser.parse_args('-s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp dssp -ff martini3001 -method exp'.split())


    dict_cutoffs = {
        'cutoff_BBB_angles': args.cutoff_BBB_angles,
        'cutoff_BBBB_dihedrals': args.cutoff_BBBB_dihedrals,
        'cutoff_SBBS_dihedrals': args.cutoff_SBBS_dihedrals,
        'cutoff_contacts': args.cutoff_contacts 
    }
    CTGOMartinize(args.strfile, args.mapfile, args.moltype, args.mbmoltype, dict_cutoffs, method=args.method, dssp=args.dssp, ff=args.ff)
