#!/usr/bin/env python3
# 20250819 Update: fixing the switching to contacts

from __future__ import annotations

import os
import argparse
import subprocess
from typing import Any

import MDAnalysis as mda
import ctgomartini
from ctgomartini.api import GenMBPTop, GenSBPTop
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds
from ctgomartini.utils import create_go_virt_for_multimer


def Martinize2(
    aa_strfile: str,
    dssp: str,
    ff: str,
    state_name: str,
    other_params: str = '',
) -> None:
    """Run martinize2 to convert all-atom structure to coarse-grained model.

    Args:
        aa_strfile: Path to the all-atom structure file.
        dssp: Path to the DSSP executable for secondary structure determination.
        ff: Force field name to use (e.g., 'martini3001').
        state_name: Name of the state/molecule type.
        other_params: Additional parameters to pass to martinize2.

    Raises:
        Exception: If martinize2 command fails.
    """
    output = subprocess.run(
        f'martinize2 -f {aa_strfile} -o system.top -x {state_name}_cg.pdb '
        f'-dssp {dssp} -p backbone -ff {ff} -govs-include -govs-moltype {state_name} '
        f'-cys auto {other_params}',
        shell=True, capture_output=True, encoding='utf-8'
    )
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


def GetNatoms(cg_strfile: str, atomname: str = 'CA') -> int:
    """Get the number of atoms from a coarse-grained structure file.

    Args:
        cg_strfile: Path to the coarse-grained structure file.
        atomname: Name of the atom to query (default: 'CA').

    Returns:
        The number of atoms.
    """
    u = mda.Universe(cg_strfile)
    Natoms = u.select_atoms(f'name {atomname}')[0].index
    return Natoms


def GenGoContacts(
    aa_strfile: str,
    cg_strfile: str,
    aa_map: str,
    state_name: str,
    go_eps: int = 12,
) -> None:
    """Generate Go-like contacts for a protein structure.

    Args:
        aa_strfile: Path to the all-atom structure file.
        cg_strfile: Path to the coarse-grained structure file.
        aa_map: Path to the atom mapping file.
        state_name: Name of the state/molecule type.
        go_eps: Epsilon value for Go contacts (default: 12).
    """
    import contextlib
    import io

    f = io.StringIO()
    with contextlib.redirect_stdout(f):
        with contextlib.redirect_stderr(f):
            Natoms = GetNatoms(cg_strfile, atomname='CA')
            create_go_virt_for_multimer(
                f'-r {aa_strfile} -s {cg_strfile} -f {aa_map} '
                f'--moltype {state_name} --go_eps {go_eps} --Natoms {Natoms}'
            )
    output = f.getvalue()

    print(
        f'create_go_virt_for_multimer -r {aa_strfile} -s {cg_strfile} '
        f'-f {aa_map} --moltype {state_name} --go_eps {go_eps} --Natoms {Natoms}'
    )
    output = output.split('\n')
    for line in output:
        if line.strip().startswith('Only symmetric OV + rCSU contacts (singly counted):'):
            print(line)


def ModifyFF(forcefield_file: str = 'martini_v3.0.0.itp') -> None:
    """Modify the force field file to include Go virtual site definitions.

    Args:
        forcefield_file: Path to the force field file to modify.

    Raises:
        Exception: If the modification command fails.
    """
    output = subprocess.run(
        r'''
sed -i "s/\[ nonbond_params \]/\#ifdef GO_VIRT\n\#include \"BB-part-def_VirtGoSites.itp\"\n\#endif\n\n\[ nonbond_params \]/" {}
echo "\n#ifdef GO_VIRT \n#include \"go-table_VirtGoSites.itp\"\n#endif" >> {}
sed -i 's/#include "martini.itp"/#include "{}"/g' system.top
'''.format(forcefield_file, forcefield_file, forcefield_file),
        shell=True, capture_output=True, encoding='utf-8'
    )
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


def MBGOMartinize(
    aa_strfile_list: list[str],
    aa_map_list: list[str],
    state_name_list: list[str],
    mbmol_name: str,
    dict_cutoffs: dict[str, float],
    method: str = 'exp',
    dssp: str = 'dssp',
    ff: str = 'martini3001',
    other_params: str = '',
) -> None:
    """Generate multiple-basin Go-Martini topology for multiple states.

    This function processes multiple all-atom protein structures, converts them
    to coarse-grained models, generates Go contacts, and combines them into a
    single multiple-basin potential topology.

    Args:
        aa_strfile_list: List of paths to all-atom structure files.
        aa_map_list: List of paths to atom mapping files.
        state_name_list: List of state names for each structure.
        mbmol_name: Name for the combined multiple-basin molecule.
        dict_cutoffs: Dictionary of cutoff values for angles, dihedrals, and contacts.
        method: Mixing method ('exp' or 'ham', default: 'exp').
        dssp: Path to DSSP executable (default: 'dssp').
        ff: Force field to use (default: 'martini3001').
        other_params: Additional parameters for martinize2.

    Raises:
        ValueError: If a state directory already exists.
        AssertionError: If an unsupported force field is specified.
    """
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
        print(
            f'\nMartinize the all-atom protein ({aa_strfile}) as the CG model '
            f'with the state name ({state_name})'
        )
        Martinize2(os.path.join('../', aa_strfile), dssp, ff, state_name, other_params)

        # Generate Go-Contacts
        print(f'\nGenerate the Go-Contacts for proteins')
        GenGoContacts(
            os.path.join('../', aa_strfile),
            f'{state_name}_cg.pdb',
            os.path.join('../', aa_map),
            state_name,
            go_eps=12
        )

        # Fetch the FF file
        print(f'\nFetch the forcefield and append the Go-Contacts to the forcefields')
        assert ff == 'martini3001', f'Error: Unsupport the forcefield: {ff}'
        os.system(
            f"cp {os.path.join(ctgomartini.__path__[0], 'data/ForceFields/Martini300/martini_v3.0.0.itp')} ."
        )
        ModifyFF(forcefield_file='martini_v3.0.0.itp')

        # Convert Long/Short Elastic Bonds to LJ Interactions
        print('\nConvert Long/Short Elastic Bonds to LJ Interactions')
        convert_long_short_elastic_bonds(
            state_name,
            f'{state_name}_cg.pdb',
            convert_long_elastic_bonds=True,
            convert_short_elastic_bonds=False,
            lj_epsilon=12
        )

    print('############')
    os.chdir(working_path)

    # Combine multiple states into the multiple-basin potential
    print(f'\nGenerate the muliple-basin potential for {mbmol_name}')
    mols_list: list[list[str]] = []
    for i, state_name in enumerate(state_name_list):
        mols_list.append([f'{state_name}/system.top', state_name])
    mbmol = GenMBPTop(mols_list, mbmol_name, dict_cutoffs)

    # Modify the mulitple_basin parameters according to method
    if method.lower() == "exp":
        pass
    elif method.lower() == "ham":
        assert mbmol._topology['multiple_basin'][0][2] == '2', (
            f"Error: HAM mixing scheme only supports mulitple basins for two states.\n"
        )
        mbmol._topology['multiple_basin'][0] = [
            'True', 'ham', '2', 'delta', 'mbp_energy1', 'mbp_energy2'
        ]

    # Write the mbmol.itp
    print(f'\nWrite the {mbmol_name}.itp and {mbmol_name}_params.itp')
    write_itp(mbmol)
    print('Finish!')


def SBGOMartinize(
    aa_strfile_list: list[str],
    aa_map_list: list[str],
    state_name_list: list[str],
    sbmol_name: str,
    method: str = 'SBP',
    dssp: str = 'dssp',
    ff: str = 'martini3001',
    other_params: str = '',
) -> None:
    """Generate single-basin Go-Martini topology for multiple states.

    This function processes multiple all-atom protein structures, converts them
    to coarse-grained models with side-chain fixes, generates Go contacts, and
    combines them into a single-basin potential topology.

    Args:
        aa_strfile_list: List of paths to all-atom structure files.
        aa_map_list: List of paths to atom mapping files.
        state_name_list: List of state names for each structure.
        sbmol_name: Name for the combined single-basin molecule.
        method: Method to use (must be 'SBP', default: 'SBP').
        dssp: Path to DSSP executable (default: 'dssp').
        ff: Force field to use (default: 'martini3001').
        other_params: Additional parameters for martinize2.

    Raises:
        ValueError: If a state directory already exists or unsupported method.
        AssertionError: If an unsupported force field is specified.
    """
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

        # Martinize AA Proteins and scfix
        print(
            f'\nMartinize the all-atom protein ({aa_strfile}) as the CG model '
            f'with the state name ({state_name})'
        )
        Martinize2(
            os.path.join('../', aa_strfile),
            dssp,
            ff,
            state_name,
            other_params=f'-scfix {other_params}'
        )

        # Generate Go-Contacts
        print(f'\nGenerate the Go-Contacts for proteins')
        GenGoContacts(
            os.path.join('../', aa_strfile),
            f'{state_name}_cg.pdb',
            os.path.join('../', aa_map),
            state_name,
            go_eps=12
        )

        # Fetch the FF file
        print(f'\nFetch the forcefield and append the Go-Contacts to the forcefields')
        assert ff == 'martini3001', f'Error: Unsupport the forcefield: {ff}'
        os.system(
            f"cp {os.path.join(ctgomartini.__path__[0], 'data/ForceFields/Martini300/martini_v3.0.0.itp')} ."
        )
        ModifyFF(forcefield_file='martini_v3.0.0.itp')

        # # Convert Long/Short Elastic Bonds to LJ Interactions
        # print('\nConvert Long/Short Elastic Bonds to LJ Interactions')
        # convert_long_short_elastic_bonds(state_name, f'{state_name}_cg.pdb',
        #   convert_long_elastic_bonds=True, convert_short_elastic_bonds=False, lj_epsilon=12)

    print('############')
    os.chdir(working_path)

    # Combine single states into the single-basin potential
    print(f'\nGenerate the single-basin potential for {sbmol_name}')
    mols_list: list[list[str]] = []
    for i, state_name in enumerate(state_name_list):
        mols_list.append([f'{state_name}/system.top', state_name])
    sbmol = GenSBPTop(mols_list, sbmol_name)

    # Modify the mulitple_basin parameters according to method
    if method.lower() != "sbp":
        raise ValueError(f'Error: Unsupport the method: {method}')

    # Write the sbmol.itp
    print(f'\nWrite the {sbmol_name}.itp and {sbmol_name}_params.itp')
    write_itp(sbmol)
    print('Finish!')


def SwitchingGOMartinize(
    aa_strfile_list: list[str],
    aa_map_list: list[str],
    state_name_list: list[str],
    mbmol_name: str,
    dict_cutoffs: dict[str, float],
    method: str = 'switching',
    dssp: str = 'dssp',
    ff: str = 'martini3001',
    other_params: str = '',
) -> None:
    """Generate switching Go-Martini topology for multiple states.

    This function processes multiple all-atom protein structures using the
    switching method, which combines single-basin potentials for each state.

    Args:
        aa_strfile_list: List of paths to all-atom structure files.
        aa_map_list: List of paths to atom mapping files.
        state_name_list: List of state names for each structure.
        mbmol_name: Name for the combined molecule.
        dict_cutoffs: Dictionary of cutoff values for angles, dihedrals, and contacts.
        method: Method to use (default: 'switching').
        dssp: Path to DSSP executable (default: 'dssp').
        ff: Force field to use (default: 'martini3001').
        other_params: Additional parameters for martinize2.

    Raises:
        AssertionError: If an unsupported force field is specified.
    """
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
        print(
            f'\nMartinize the all-atom protein ({aa_strfile}) as the CG model '
            f'with the state name ({state_name})'
        )
        Martinize2(
            os.path.join('../', aa_strfile),
            dssp,
            ff,
            state_name,
            other_params=f'-scfix {other_params}'
        )

        # Generate Go-Contacts
        print(f'\nGenerate the Go-Contacts for proteins')
        GenGoContacts(
            os.path.join('../', aa_strfile),
            f'{state_name}_cg.pdb',
            os.path.join('../', aa_map),
            state_name,
            go_eps=12
        )

        # Fetch the FF file
        print(f'\nFetch the forcefield and append the Go-Contacts to the forcefields')
        assert ff == 'martini3001', f'Error: Unsupport the forcefield: {ff}'
        os.system(
            f"cp {os.path.join(ctgomartini.__path__[0], 'data/ForceFields/Martini300/martini_v3.0.0.itp')} ."
        )
        ModifyFF(forcefield_file='martini_v3.0.0.itp')

        # # Convert Long/Short Elastic Bonds to LJ Interactions
        # print('\nConvert Long/Short Elastic Bonds to LJ Interactions')
        # convert_long_short_elastic_bonds(state_name, f'{state_name}_cg.pdb',
        #   convert_long_elastic_bonds=True, convert_short_elastic_bonds=False, lj_epsilon=12)
        isConvertLJ2Contacts = True
        if isConvertLJ2Contacts:
            # Combine single states into the single-basin potential
            print(f'\nGenerate the single-basin potential for {state_name}')
            mols_list: list[list[str]] = []

            mols_list.append([f'system.top', state_name])
            sbmol = GenSBPTop(mols_list, state_name)

            # Write the sbmol.itp
            os.chdir(working_path)
            print(f'\nWrite the {state_name}.itp and {state_name}_params.itp')
            write_itp(sbmol)
            print('Finish!')

    print('############')
    os.chdir(working_path)
    print('Finish!')


def CTGOMartinize(
    aa_strfile_list: list[str],
    aa_map_list: list[str],
    state_name_list: list[str],
    mbmol_name: str,
    dict_cutoffs: dict[str, float],
    method: str = 'switching',
    dssp: str = 'dssp',
    ff: str = 'martini3001',
    other_params: str = '',
) -> None:
    """Main entry point for Go-Martini topology generation.

    This function dispatches to the appropriate Go-Martini method based on
    the specified method parameter.

    Args:
        aa_strfile_list: List of paths to all-atom structure files.
        aa_map_list: List of paths to atom mapping files.
        state_name_list: List of state names for each structure.
        mbmol_name: Name for the combined molecule.
        dict_cutoffs: Dictionary of cutoff values for angles, dihedrals, and contacts.
        method: Method to use ('switching', 'exp', 'ham', or 'sbp', default: 'switching').
        dssp: Path to DSSP executable (default: 'dssp').
        ff: Force field to use (default: 'martini3001').
        other_params: Additional parameters for martinize2.

    Raises:
        ValueError: If an unsupported method is specified.
    """
    if method.lower() == 'switching':
        SwitchingGOMartinize(
            aa_strfile_list, aa_map_list, state_name_list, mbmol_name,
            dict_cutoffs, method=method, dssp=dssp, ff=ff, other_params=other_params
        )
    elif method.lower() in ['exp', 'ham']:
        MBGOMartinize(
            aa_strfile_list, aa_map_list, state_name_list, mbmol_name,
            dict_cutoffs, method=method, dssp=dssp, ff=ff, other_params=other_params
        )
    elif method.lower() == 'sbp':
        SBGOMartinize(
            aa_strfile_list, aa_map_list, state_name_list, sbmol_name=mbmol_name,
            method=method, dssp=dssp, ff=ff, other_params=other_params
        )
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
    parser.add_argument('-other_params', dest='other_params', default='', type=str,
                        help='Other parameters used for martinize')
    args = parser.parse_args()
    # args = parser.parse_args('-s 1GGG_1_clean.pdb 1WDN_1_clean.pdb -m 1GGG_1_clean.map 1WDN_1_clean.map -mol gbp_open gbp_closed -mbmol gbp -dssp dssp -ff martini3001 -method exp'.split())

    dict_cutoffs: dict[str, float] = {
        'cutoff_BBB_angles': args.cutoff_BBB_angles,
        'cutoff_BBBB_dihedrals': args.cutoff_BBBB_dihedrals,
        'cutoff_SBBS_dihedrals': args.cutoff_SBBS_dihedrals,
        'cutoff_contacts': args.cutoff_contacts
    }
    CTGOMartinize(
        args.strfile, args.mapfile, args.moltype, args.mbmoltype,
        dict_cutoffs, method=args.method, dssp=args.dssp, ff=args.ff,
        other_params=args.other_params
    )
