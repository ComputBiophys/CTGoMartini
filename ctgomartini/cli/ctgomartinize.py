#!/usr/bin/env python3
"""
CTGoMartini topology generation script.

This script generates Go-Martini topology files for coarse-grained molecular dynamics
simulations. It supports multiple methods:
- Single-basin Go-Martini (SBP): Convert one structure to Go-Contact format
- Multiple-basin Go-Martini with exponential mixing (EXP): 2+ structures
- Multiple-basin Go-Martini with Hamiltonian mixing (HAM): Exactly 2 structures
- Switching Go-Martini: Generate separate topologies for each structure

The script integrates with martinize2 (vermouth >= 0.10.0) to:
1. Convert all-atom structures to coarse-grained Martini models
2. Generate Go-like contacts from various sources (auto, rCSU .map, martinize2 .out)
3. Combine multiple states into multi-basin potential topologies
"""


from __future__ import annotations

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import Any

import ctgomartini
from ctgomartini.topology import create_mb_topology, create_sb_topology
from ctgomartini.utils import write_itp, convert_long_short_elastic_bonds, extract_all_states
from ctgomartini.utils.contacts import gen_go_contacts
from ctgomartini.utils.pdb_validation import (
    validate_pdb_compatibility,
    PDBCompatibilityValidator
)
from ctgomartini.utils.constraints_to_bonds import convert_constraints_to_bonds


def Martinize2(
    aa_strfile: Path,
    ff: str,
    state_name: str,
    go_file: Path | None = None,
    go_eps: float = 12.0,
    go_low: float = 0.3,
    go_up: float = 1.1,
    dssp: str | None = None,
    other_params: str = '',
    go_res_dist: int = 0,
) -> None:
    """Run martinize2 to convert all-atom structure to coarse-grained model."""
    cmd = [
        'martinize2',
        '-f', str(aa_strfile),
        '-o', 'system.top',
        '-x', f'{state_name}_cg.pdb',
        '-p', 'backbone',
        '-ff', ff,
        '-cys', 'auto',
    ]

    if dssp and dssp != 'auto':
        cmd.extend(['-dssp', dssp])
    else:
        cmd.append('-dssp')  # Use MDTraj

    if go_file and go_file.exists():
        cmd.extend([
            '-go', str(go_file),
            '-go-low', str(go_low),
            '-go-up', str(go_up),
            '-go-eps', str(go_eps),
            '-go-res-dist', str(go_res_dist),  # Disable martinize2's filtering
        ])

    cmd.extend(['-name', state_name])

    other_params_clean = other_params.replace('-scfix', '').replace('-scfix on', '').strip()
    if other_params_clean:
        cmd.extend(other_params_clean.split())

    result = subprocess.run(cmd, capture_output=True, encoding='utf-8')
    print(' '.join(cmd))
    
    if result.returncode != 0:
        print(f"Error: martinize2 failed")
        print("stdout:", result.stdout)
        print("stderr:", result.stderr)
        raise Exception(f'martinize2 failed in {os.getcwd()}')


def ModifyFF(forcefield_file: str = 'martini_v3.0.0.itp') -> None:
    """Modify the force field file to include Go virtual site definitions."""
    # Add go_atomtypes.itp before nonbond_params
    with open(forcefield_file, 'r') as f:
        content = f.read()
    
    content = content.replace(
        '[ nonbond_params ]',
        '#ifdef GO_VIRT\n#include "go_atomtypes.itp"\n#endif\n\n[ nonbond_params ]'
    )
    
    with open(forcefield_file, 'w') as f:
        f.write(content)
    
    # Append go_nbparams.itp
    with open(forcefield_file, 'a') as f:
        f.write('\n#ifdef GO_VIRT\n#include "go_nbparams.itp"\n#endif\n')
    
    # Update system.top
    with open('system.top', 'r') as f:
        top_content = f.read()
    
    top_content = top_content.replace(
        '#include "martini.itp"',
        f'#include "{forcefield_file}"'
    )
    
    with open('system.top', 'w') as f:
        f.write(top_content)


def _process_single_state(
    working_path: Path,
    state_dir: Path,
    aa_strfile: Path,
    map_file: str | None,
    state_name: str,
    ff: str,
    dssp: str | None,
    other_params: str,
    go_eps: float,
    go_low: float,
    go_up: float,
    convert_lj: bool = False,
    min_seq_distance: int = 4,
) -> None:
    """Process a single state: generate contacts, run martinize2, setup forcefield.
    
    This is a helper function to reduce code duplication between MB/SB/Switching modes.
    """
    print(f"\n{'─' * 60}")
    print(f"  State: {state_name}")
    print(f"  Structure: {aa_strfile.name}")
    print(f"  Contacts: {map_file if map_file else 'auto'}")
    print(f"{'─' * 60}")

    # Create state directory
    state_dir.mkdir(parents=True, exist_ok=True)

    # Resolve paths relative to working_path
    aa_strfile_abs = working_path / aa_strfile
    
    # Determine map file path
    if map_file is not None and map_file.lower() != 'auto':
        map_file_path = working_path / map_file
    else:
        map_file_path = None  # auto mode

    # Change to state_dir temporarily for all operations
    original_cwd = os.getcwd()
    try:
        os.chdir(state_dir)
        
        # Generate Go contacts (now in state_dir)
        print("  [Step 1/4] Generating Go contacts...")
        go_file = gen_go_contacts(
            str(aa_strfile_abs),
            str(map_file_path) if map_file_path else None,
            state_name,
            go_eps=go_eps,
            min_seq_distance=min_seq_distance,
        )
        print(f"    ✓ Go contacts generated: {go_file}")

        # Run martinize2
        print("  [Step 2/4] Running martinize2...")
        go_file_abs = state_dir / go_file if go_file else None
        
        Martinize2(
            aa_strfile_abs,
            ff,
            state_name,
            go_file=go_file_abs,
            go_eps=go_eps,
            go_low=go_low,
            go_up=go_up,
            dssp=dssp,
            other_params=other_params,
            go_res_dist=0,  # Disable martinize2's filtering; we do it ourselves
        )
        print(f"    ✓ Coarse-grained topology created: system.top")

        # Copy force field
        print("  [Step 3/4] Setting up force field...")
        assert ff == 'martini3001', f'Error: Unsupported forcefield: {ff}'
        ff_source = Path(ctgomartini.__path__[0]) / 'data/ForceFields/Martini300/martini_v3.0.0.itp'
        shutil.copy(ff_source, '.')
        ModifyFF(forcefield_file='martini_v3.0.0.itp')
        print(f"    ✓ Force field configured with Go virtual sites")

        # Convert long/short elastic bonds if needed
        if convert_lj:
            print("  [Step 4/4] Converting long/short elastic bonds to LJ interactions...")
            convert_long_short_elastic_bonds(
                state_name,
                f'{state_name}_cg.pdb',
                convert_long_elastic_bonds=True,
                convert_short_elastic_bonds=False,
                lj_epsilon=go_eps
            )
            print(f"    ✓ Elastic bonds converted")
        else:
            print("  [Step 4/4] Skipping elastic bond conversion")
        
        print(f"  ✓ State '{state_name}' processed successfully")
    finally:
        os.chdir(original_cwd)


def MBGOMartinize(
    aa_strfile_list: list[str],
    map_file_list: list[str | None],
    state_name_list: list[str],
    mbmol_name: str,
    dict_cutoffs: dict[str, float],
    method: str = 'exp',
    dssp: str | None = None,
    ff: str = 'martini3001',
    other_params: str = '',
    go_eps: float = 12.0,
    go_low: float = 0.3,
    go_up: float = 1.1,
    constraints2bonds: float | None = None,
    min_seq_distance: int = 4,
    extract_states: bool = False,
) -> None:
    """Generate multiple-basin Go-Martini topology for multiple states."""
    working_path = Path.cwd()
    
    print("\n" + "=" * 60)
    print(f"  CTGoMartini - Multiple-Basin Mode ({method.upper()})")
    print("=" * 60)
    print(f"  Working directory: {working_path}")
    print(f"  States: {', '.join(state_name_list)}")
    print(f"  Output molecule: {mbmol_name}")
    print("=" * 60)

    # Validate PDB compatibility before processing
    print("\n  [Pre-check] Validating PDB file compatibility...")
    try:
        validate_pdb_compatibility(aa_strfile_list, state_name_list, verbose=False)
        print("    ✓ All PDB files are compatible")
    except ValueError as e:
        print(e)
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"\n  ERROR: {e}")
        sys.exit(1)

    # For MB mode: use -noscfix
    mb_other_params = f'{other_params} -noscfix'.strip()

    for aa_strfile, map_file, state_name in zip(aa_strfile_list, map_file_list, state_name_list):
        state_dir = working_path / state_name
        if state_dir.exists():
            raise ValueError(f'Error: Directory {state_name} exists!')

        _process_single_state(
            working_path=working_path,
            state_dir=state_dir,
            aa_strfile=Path(aa_strfile),
            map_file=map_file,
            state_name=state_name,
            ff=ff,
            dssp=dssp,
            other_params=mb_other_params,
            go_eps=go_eps,
            go_low=go_low,
            go_up=go_up,
            convert_lj=True,
            min_seq_distance=min_seq_distance,
        )

    print("\n" + "─" * 60)
    print(f"  Combining {len(state_name_list)} states into multiple-basin potential...")
    
    mols_list: list[list[str]] = [
        [f'{state_name}/system.top', state_name]
        for state_name in state_name_list
    ]
    mbmol = create_mb_topology(mols_list, mbmol_name, dict_cutoffs)

    # Modify multiple_basin parameters according to method
    if method.lower() == "ham":
        assert mbmol._topology['multiple_basin'][0][2] == '2', (
            "Error: HAM mixing scheme only supports multiple basins for two states."
        )
        mbmol._topology['multiple_basin'][0] = [
            'True', 'ham', '2', 'delta', 'C1', 'C2'
        ]

    write_itp(mbmol)
    print(f"  ✓ Topology written: {mbmol_name}.itp, {mbmol_name}_params.itp")

    # Extract single-state topologies if requested
    if extract_states:
        print(f"\n  [Extracting single-state topologies...")
        try:
            # Find the system.top file (should be in current directory)
            system_top = working_path / 'system.top'
            if system_top.exists():
                state_files = extract_all_states(system_top, mbmol_name)
                for state_file in state_files:
                    print(f"    ✓ Extracted: {state_file}")
            else:
                print(f"    ! Warning: system.top not found, skipping state extraction")
        except Exception as e:
            print(f"    ! Warning: Failed to extract states: {e}")

    # Post-process: convert constraints to bonds in final topology
    if constraints2bonds is not None:
        print(f"\n  [Post-processing] Converting constraints to bonds (fc={constraints2bonds})...")
        try:
            final_itp = Path(f"{mbmol_name}.itp")
            if final_itp.exists():
                n_converted = convert_constraints_to_bonds(
                    final_itp, mbmol_name, fc=constraints2bonds
                )
                if n_converted > 0:
                    print(f"    ✓ Converted {n_converted} constraints to bonds")
                else:
                    print(f"    ℹ No constraints found in {mbmol_name}")
        except Exception as e:
            print(f"    ! Warning: Failed to convert constraints: {e}")

    print("\n" + "=" * 60)
    print("  FINISHED SUCCESSFULLY")
    print("=" * 60)
    print("  Output files:")
    for state_name in state_name_list:
        print(f"    - {state_name}/{state_name}_cg.pdb")
        print(f"    - {state_name}/{state_name}.itp")
    print(f"    - {mbmol_name}.itp")
    print(f"    - {mbmol_name}_params.itp")
    if extract_states:
        for i in range(len(state_name_list)):
            state_letter = chr(ord('A') + i)
            print(f"    - {mbmol_name}_state{state_letter}.itp")
    print("=" * 60)


def SBGOMartinize(
    aa_strfile_list: list[str],
    map_file_list: list[str | None],
    state_name_list: list[str],
    method: str = 'SBP',
    dssp: str | None = None,
    ff: str = 'martini3001',
    other_params: str = '',
    go_eps: float = 12.0,
    go_low: float = 0.3,
    go_up: float = 1.1,
    constraints2bonds: float | None = None,
    min_seq_distance: int = 4,
) -> None:
    """Generate single-basin Go-Martini topology from a single PDB.
    
    SBP mode converts martinize2's Go-LJ model to CTGoMartini Go-Contact format.
    Requires exactly one PDB file and one molecule name.
    """
    # Strict validation: SBP requires exactly one PDB
    if len(aa_strfile_list) != 1:
        raise ValueError(
            f"SBP mode requires exactly one PDB file, got {len(aa_strfile_list)}.\n"
            f"Usage: ctgomartinize -s protein.pdb -m auto -mol protein -ff martini3001 -method sbp\n"
            f"For multiple structures, use EXP, HAM, or Switching mode."
        )
    if len(state_name_list) != 1:
        raise ValueError(
            f"SBP mode requires exactly one molecule name, got {len(state_name_list)}.\n"
            f"Usage: ctgomartinize -s protein.pdb -mol protein -method sbp"
        )
    
    working_path = Path.cwd()
    state_name = state_name_list[0]
    
    print("\n" + "=" * 60)
    print(f"  CTGoMartini - Single-Basin Mode (SBP)")
    print("=" * 60)
    print(f"  Working directory: {working_path}")
    print(f"  Structure: {aa_strfile_list[0]}")
    print(f"  Molecule: {state_name}")
    print("=" * 60)

    state_dir = working_path / state_name
    if state_dir.exists():
        raise ValueError(f'Error: Directory {state_name} exists!')

    _process_single_state(
        working_path=working_path,
        state_dir=state_dir,
        aa_strfile=Path(aa_strfile_list[0]),
        map_file=map_file_list[0],
        state_name=state_name,
        ff=ff,
        dssp=dssp,
        other_params=other_params,
        go_eps=go_eps,
        go_low=go_low,
        go_up=go_up,
        convert_lj=False,
        min_seq_distance=min_seq_distance,
    )

    print("\n" + "─" * 60)
    print(f"  Generating single-basin topology...")
    
    mols_list: list[list[str]] = [[f'{state_name}/system.top', state_name]]
    sbmol = create_sb_topology(mols_list, state_name)

    if method.lower() != "sbp":
        raise ValueError(f'Error: Unsupported method: {method}')

    write_itp(sbmol)
    print(f"  ✓ Topology written: {state_name}.itp, {state_name}_params.itp")

    # Post-process: convert constraints to bonds in final topology
    if constraints2bonds is not None:
        print(f"\n  [Post-processing] Converting constraints to bonds (fc={constraints2bonds})...")
        try:
            final_itp = Path(f"{state_name}.itp")
            if final_itp.exists():
                n_converted = convert_constraints_to_bonds(
                    final_itp, state_name, fc=constraints2bonds
                )
                if n_converted > 0:
                    print(f"    ✓ Converted {n_converted} constraints to bonds")
                else:
                    print(f"    ℹ No constraints found in {state_name}")
        except Exception as e:
            print(f"    ! Warning: Failed to convert constraints: {e}")

    print("\n" + "=" * 60)
    print("  FINISHED SUCCESSFULLY")
    print("=" * 60)
    print("  Output files:")
    print(f"    - {state_name}/{state_name}_cg.pdb")
    print(f"    - {state_name}.itp")
    print(f"    - {state_name}_params.itp")
    print("=" * 60)


def SwitchingGOMartinize(
    aa_strfile_list: list[str],
    map_file_list: list[str | None],
    state_name_list: list[str],
    mbmol_name: str,
    dict_cutoffs: dict[str, float],
    method: str = 'switching',
    dssp: str | None = None,
    ff: str = 'martini3001',
    other_params: str = '',
    go_eps: float = 12.0,
    go_low: float = 0.3,
    go_up: float = 1.1,
    constraints2bonds: float | None = None,
    min_seq_distance: int = 4,
) -> None:
    """Generate switching Go-Martini topology for multiple states."""
    working_path = Path.cwd()
    
    print("\n" + "=" * 60)
    print(f"  CTGoMartini - Switching Mode")
    print("=" * 60)
    print(f"  Working directory: {working_path}")
    print(f"  States: {', '.join(state_name_list)}")
    print("=" * 60)

    # Validate PDB compatibility before processing
    print("\n  [Pre-check] Validating PDB file compatibility...")
    try:
        validate_pdb_compatibility(aa_strfile_list, state_name_list, verbose=False)
        print("    ✓ All PDB files are compatible")
    except ValueError as e:
        print(e)
        sys.exit(1)
    except FileNotFoundError as e:
        print(f"\n  ERROR: {e}")
        sys.exit(1)

    for aa_strfile, map_file, state_name in zip(aa_strfile_list, map_file_list, state_name_list):
        state_dir = working_path / state_name
        # Remove and recreate for switching mode
        if state_dir.exists():
            shutil.rmtree(state_dir)
        state_dir.mkdir(parents=True)

        _process_single_state(
            working_path=working_path,
            state_dir=state_dir,
            aa_strfile=Path(aa_strfile),
            map_file=map_file,
            state_name=state_name,
            ff=ff,
            dssp=dssp,
            other_params=other_params,
            go_eps=go_eps,
            go_low=go_low,
            go_up=go_up,
            convert_lj=False,
            min_seq_distance=min_seq_distance,
        )

        # Generate single-basin topology for switching (inside state dir)
        print(f"\n  Generating single-basin topology for {state_name}...")
        mols_list: list[list[str]] = [[str(state_dir / 'system.top'), state_name]]
        sbmol = create_sb_topology(mols_list, state_name)
        write_itp(sbmol)
        print(f"  ✓ Topology written: {state_name}.itp, {state_name}_params.itp")

        # Post-process: convert constraints to bonds for this state
        if constraints2bonds is not None:
            print(f"  [Post-processing] Converting constraints to bonds (fc={constraints2bonds})...")
            try:
                state_itp = Path(f"{state_name}.itp")
                if state_itp.exists():
                    n_converted = convert_constraints_to_bonds(
                        state_itp, state_name, fc=constraints2bonds
                    )
                    if n_converted > 0:
                        print(f"    ✓ Converted {n_converted} constraints to bonds")
                    else:
                        print(f"    ℹ No constraints found")
            except Exception as e:
                print(f"    ! Warning: Failed to convert constraints: {e}")

    print("\n" + "=" * 60)
    print("  FINISHED SUCCESSFULLY")
    print("=" * 60)
    print("  Output files:")
    for state_name in state_name_list:
        print(f"    - {state_name}/{state_name}_cg.pdb")
        print(f"    - {state_name}.itp")
        print(f"    - {state_name}_params.itp")
    print("=" * 60)


def CTGOMartinize(
    aa_strfile_list: list[str],
    map_file_list: list[str | None],
    state_name_list: list[str],
    mbmol_name: str,
    dict_cutoffs: dict[str, float],
    method: str = 'switching',
    dssp: str | None = None,
    ff: str = 'martini3001',
    other_params: str = '',
    go_eps: float = 12.0,
    go_low: float = 0.3,
    go_up: float = 1.1,
    constraints2bonds: float | None = None,
    min_seq_distance: int = 4,
    extract_states: bool = False,
) -> None:
    """Main entry point for Go-Martini topology generation."""
    if method.lower() == 'switching':
        SwitchingGOMartinize(
            aa_strfile_list, map_file_list, state_name_list, mbmol_name,
            dict_cutoffs, method=method, dssp=dssp, ff=ff, other_params=other_params,
            go_eps=go_eps, go_low=go_low, go_up=go_up,
            constraints2bonds=constraints2bonds,
            min_seq_distance=min_seq_distance,
        )
    elif method.lower() in ['exp', 'ham']:
        MBGOMartinize(
            aa_strfile_list, map_file_list, state_name_list, mbmol_name,
            dict_cutoffs, method=method, dssp=dssp, ff=ff, other_params=other_params,
            go_eps=go_eps, go_low=go_low, go_up=go_up,
            constraints2bonds=constraints2bonds,
            min_seq_distance=min_seq_distance,
            extract_states=extract_states,
        )
    elif method.lower() == 'sbp':
        SBGOMartinize(
            aa_strfile_list, map_file_list, state_name_list,
            method=method, dssp=dssp, ff=ff, other_params=other_params,
            go_eps=go_eps, go_low=go_low, go_up=go_up,
            constraints2bonds=constraints2bonds,
            min_seq_distance=min_seq_distance,
        )
    else:
        raise ValueError(f'Error: unsupported method: {method}!')


def main() -> None:
    """Main entry point for command-line interface."""
    parser = argparse.ArgumentParser(
        description="Generate topology files for Multiple-basin or Switching Go-Martini method.",
        epilog="""Contact Map Options (via -m parameter):
  1. 'auto' (default): Automatically generate contacts using OVrCSU algorithm
  2. User-provided files: 
     - rCSU web-server .map files (will be converted)
     - martinize2/contact_map.out files (used directly)

Usage examples:

    # Single-basin Go-Martini (SBP) - single structure
    ctgomartinize -s protein.pdb -m auto -mol protein -ff martini3001 -dssp -method sbp
    
    # Multiple-basin Go-Martini with exponential mixing (EXP) - 2+ structures
    ctgomartinize -s StateA.pdb StateB.pdb -m auto -mol StateA StateB -mbmol protein -ff martini3001 -dssp -method exp
    
    # Multiple-basin Go-Martini with Hamiltonian mixing (HAM) - exactly 2 structures
    ctgomartinize -s StateA.pdb StateB.pdb -m auto -mol StateA StateB -mbmol protein -ff martini3001 -dssp -method ham
    
    # Switching Go-Martini - generates separate topologies for each structure
    ctgomartinize -s StateA.pdb StateB.pdb -m auto -mol StateA StateB -ff martini3001 -dssp -method switching

    # Use provided contact map files instead of auto-generation
    ctgomartinize -s StateA.pdb StateB.pdb -m StateA.map StateB.map -mol StateA StateB -mbmol protein -ff martini3001 -dssp -method exp

    # Convert constraints to bonds (allows slight flexibility, default FC=50000)
    ctgomartinize -s protein.pdb -m auto -mol protein -ff martini3001 -dssp -method sbp -constraints2bonds

    # Convert constraints with custom force constant
    ctgomartinize -s protein.pdb -m auto -mol protein -ff martini3001 -dssp -method sbp -constraints2bonds 2000

    # Extract single-state topologies for REMD unsampled states (EXP/HAM only)
    ctgomartinize -s StateA.pdb StateB.pdb -m auto -mol StateA StateB -mbmol protein -ff martini3001 -dssp -method exp -extract-states
    # Generates: protein.itp, protein_params.itp, protein_stateA.itp, protein_stateB.itp

Notes:
    - SBP mode: Requires exactly one PDB file and one molecule name
    - EXP/HAM/Switching modes: Require 2+ PDB files with matching residue sequences
    - -dssp: Use alone for MDTraj (default), or provide path to DSSP executable
    - -mbmol: Output molecule name for multi-basin modes (EXP, HAM, Switching)
    - -other_params: For params starting with dash, use -other_params="-param" (e.g., -other_params="-nt")
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-s', dest='strfile', required=True, nargs='+', type=str,
                        help='Input structure files')
    parser.add_argument('-m', dest='mapfile', required=False, nargs='+', type=str,
                        default=None,
                        help='Contact map files: "auto" for automatic generation (default), '
                             'or paths to .map/.out files')
    parser.add_argument('-mol', dest='moltype', required=True, nargs='+', type=str,
                        help='Molecule type names')
    parser.add_argument('-mbmol', dest='mbmoltype', default='mbmol', type=str,
                        help='Multiple-basin molecule type name (default: mbmol)')
    parser.add_argument('-dssp', dest='dssp', default='auto', nargs='?', const='auto', type=str,
                        help='DSSP executable path. Use "-dssp" alone for MDTraj (default), '
                             'or provide path like "-dssp /usr/bin/mkdssp"')
    parser.add_argument('-ff', dest='ff', default='martini3001', type=str,
                        help='Force field to use (default: martini3001)')
    parser.add_argument('-method', dest='method', required=True, type=str,
                        help='Method to use (required: exp, ham, switching, sbp)')
    parser.add_argument('-cutoff_BBB_angles', dest='cutoff_BBB_angles', default=15.0, type=float,
                        help='Cutoff of BBB angles (default: 15.0 degree)')
    parser.add_argument('-cutoff_BBBB_dihedrals', dest='cutoff_BBBB_dihedrals', default=30.0, type=float,
                        help='Cutoff of BBBB dihedrals (default: 30.0 degree)')
    parser.add_argument('-cutoff_SBBS_dihedrals', dest='cutoff_SBBS_dihedrals', default=30.0, type=float,
                        help='Cutoff of SBBS dihedrals (default: 30.0 degree)')
    parser.add_argument('-cutoff_contacts', dest='cutoff_contacts', default=0.06, type=float,
                        help='Sigma cutoff of contacts (default: 0.06 nm)')
    parser.add_argument('-other_params', dest='other_params', default='', type=str,
                        help='Other parameters for martinize2. '
                             'Use -other_params="-param" for params starting with dash (e.g., -other_params="-nt")')
    parser.add_argument('-go-eps', dest='go_eps', default=12.0, type=float,
                        help='Epsilon value for Go contacts (default: 12.0)')
    parser.add_argument('-go-low', dest='go_low', default=0.3, type=float,
                        help='Lower cutoff for Go contacts in nm (default: 0.3)')
    parser.add_argument('-go-up', dest='go_up', default=1.1, type=float,
                        help='Upper cutoff for Go contacts in nm (default: 1.1)')
    parser.add_argument('-min-seq-distance', dest='min_seq_distance', default=4, type=int,
                        help='Minimum sequence distance (PDB numbering) for intra-chain contacts. '
                             'Contacts within the same chain with distance < min-seq-distance are '
                             'filtered out. Inter-chain contacts are always kept. '
                             'Set to 0 to disable filtering. '
                             'Note: This replaces martinize2\'s -go-res-dist which has issues with '
                             'disulfide bonds. (default: 4)')
    parser.add_argument('-constraints2bonds', dest='constraints2bonds',
                        nargs='?', const=50000.0, default=None, type=float, metavar='FC',
                        help='Convert constraints to bonds with force constant FC (kJ/(mol·nm²)). '
                             'Use -constraints2bonds for default FC=50000, or '
                             '-constraints2bonds 2000 for custom value.')
    parser.add_argument('-extract-states', '-es', dest='extract_states',
                        action='store_true',
                        help='Extract single-state topologies from multi-basin system. '
                             'Generates {mbmol}_stateA.itp, {mbmol}_stateB.itp, etc. '
                             'for use as unsampled states in REMD simulations. '
                             'Only applicable for EXP and HAM methods.')
    args = parser.parse_args()

    # Handle -m parameter
    map_file_list: list[str | None]
    if args.mapfile is None or len(args.mapfile) == 0:
        map_file_list = ['auto'] * len(args.strfile)
    elif len(args.mapfile) == 1 and args.mapfile[0].lower() == 'auto':
        map_file_list = ['auto'] * len(args.strfile)
    elif len(args.mapfile) != len(args.strfile):
        raise ValueError(
            f"Number of map files ({len(args.mapfile)}) must match number of "
            f"structure files ({len(args.strfile)}), or use 'auto'"
        )
    else:
        map_file_list = args.mapfile

    dict_cutoffs: dict[str, float] = {
        'cutoff_BBB_angles': args.cutoff_BBB_angles,
        'cutoff_BBBB_dihedrals': args.cutoff_BBBB_dihedrals,
        'cutoff_SBBS_dihedrals': args.cutoff_SBBS_dihedrals,
        'cutoff_contacts': args.cutoff_contacts
    }
    
    CTGOMartinize(
        args.strfile, map_file_list, args.moltype, args.mbmoltype,
        dict_cutoffs, method=args.method, dssp=args.dssp, ff=args.ff,
        other_params=args.other_params, go_eps=args.go_eps,
        go_low=args.go_low, go_up=args.go_up,
        constraints2bonds=args.constraints2bonds,
        min_seq_distance=args.min_seq_distance,
        extract_states=args.extract_states,
    )


if __name__ == "__main__":
    main()
