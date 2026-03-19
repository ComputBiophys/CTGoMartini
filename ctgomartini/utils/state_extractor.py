"""
State extraction utilities for CTGoMartini.

Provides functionality to extract single-state topologies from multiple-basin
topologies for use in REMD simulations as unsampled states.
"""

from __future__ import annotations

from typing import Any
from pathlib import Path

from ..topology import MartiniTopFile
from .write_itp import write_itp


def extract_state_topology(
    topfile: str | Path,
    molecule_name: str,
    state_id: int,
    keep_original_name: bool = True,
) -> Any:
    """
    Extract a single-state topology from a multiple-basin topology.
    
    This function takes a multiple-basin topology and extracts the interactions
    for a specific state, creating a single-state molecule that can be used
    as an unsampled state in REMD simulations.
    
    Args:
        topfile: Path to the topology file (.top)
        molecule_name: Name of the multiple-basin molecule
        state_id: State number to extract (1-indexed, e.g., 1 for stateA)
        keep_original_name: If True, keep original molecule name; 
                           if False, append _stateA/_stateB etc.
        
    Returns:
        Molecule object with single-state topology
        
    Raises:
        ValueError: If molecule not found or not a multi-basin system
        ValueError: If state_id is invalid
    """
    # Load topology
    top = MartiniTopFile(topfile)
    
    if molecule_name not in top.moleculeTypes:
        available = ', '.join(top.moleculeTypes.keys())
        raise ValueError(f"Molecule '{molecule_name}' not found. Available: {available}")
    
    mol = top.moleculeTypes[molecule_name]
    mol_top = mol._topology
    
    # Check if this is a multi-basin system
    if 'multiple_basin' not in mol_top:
        raise ValueError(f"Molecule '{molecule_name}' is not a multi-basin system")
    
    # Get number of states
    n_states = int(mol_top['multiple_basin'][0][2])
    if state_id < 1 or state_id > n_states:
        raise ValueError(f"Invalid state_id {state_id}. Must be between 1 and {n_states}")
    
    # Create a copy of the molecule for single state
    # We'll modify the topology to extract only this state's interactions
    
    # Extract all state-specific interactions
    state_str = str(state_id)
    
    # Process multi_contacts: extract entries for this state
    if 'multi_contacts' in mol_top:
        state_contacts = _extract_multi_entries(
            mol_top['multi_contacts'], state_str, state_idx=2, n_states=n_states
        )
        mol_top['contacts'] = mol_top.get('contacts', []) + state_contacts
    
    # Process multi_angles: extract entries for this state
    if 'multi_angles' in mol_top:
        state_angles = _extract_multi_entries(
            mol_top['multi_angles'], state_str, state_idx=2, n_states=n_states
        )
        mol_top['angles'] = mol_top.get('angles', []) + state_angles
    
    # Process multi_dihedrals: extract entries for this state
    if 'multi_dihedrals' in mol_top:
        state_dihedrals = _extract_multi_entries(
            mol_top['multi_dihedrals'], state_str, state_idx=2, n_states=n_states
        )
        mol_top['dihedrals'] = mol_top.get('dihedrals', []) + state_dihedrals
    
    # Remove multi-basin specific sections
    for key in ['multiple_basin', 'multi_contacts', 'multi_angles', 'multi_dihedrals']:
        if key in mol_top:
            mol_top.pop(key)
    
    # Optionally update molecule type name to indicate single state
    if not keep_original_name:
        state_letter = chr(ord('A') + state_id - 1)  # 1->A, 2->B, etc.
        original_name = mol_top['moleculetype'][0][0]
        mol_top['moleculetype'][0][0] = f"{original_name}_state{state_letter}"
        mol.name = mol_top['moleculetype'][0][0]
    
    return mol


def _extract_multi_entries(
    entries: list[list[str]], 
    state: str, 
    state_idx: int,
    n_states: int
) -> list[list[str]]:
    """
    Helper to extract interaction entries for a specific state.
    
    Multi-state entries have format:
    [atoms..., n_states, state_id, params...]
    
    For single state extraction, we remove the state info columns
    and keep only entries matching the requested state.
    """
    state_entries = []
    for entry in entries:
        # Check if this entry is for the requested state
        # Format: [atoms(0..n), n_states(n), state_id(n+1), params(n+2...)]
        if int(entry[state_idx]) == n_states and int(entry[state_idx + 1]) == int(state):
            # Remove the n_states and state_id columns, keep only atoms and params
            new_entry = entry[:state_idx] + entry[state_idx + 2:]
            state_entries.append(new_entry)
    return state_entries


def write_state_itp(
    molecule: Any,
    output_file: str | Path,
) -> None:
    """
    Write a single-state molecule to an ITP file.
    
    Args:
        molecule: Molecule object with single-state topology
        output_file: Output ITP file path
    """
    write_itp(molecule, str(output_file))


def extract_all_states(
    topfile: str | Path,
    mbmol_name: str,
    output_prefix: str | None = None,
    keep_original_name: bool = True,
) -> list[str]:
    """
    Extract all single-state topologies from a multiple-basin topology.
    
    Args:
        topfile: Path to the topology file
        mbmol_name: Name of the multiple-basin molecule
        output_prefix: Prefix for output files (default: mbmol_name)
        keep_original_name: If True, keep original molecule name
        
    Returns:
        List of output ITP filenames generated
    """
    top = MartiniTopFile(topfile)
    
    if mbmol_name not in top.moleculeTypes:
        raise ValueError(f"Molecule '{mbmol_name}' not found in topology")
    
    mol_top = top.moleculeTypes[mbmol_name]._topology
    if 'multiple_basin' not in mol_top:
        raise ValueError(f"Molecule '{mbmol_name}' is not a multi-basin system")
    
    n_states = int(mol_top['multiple_basin'][0][2])
    prefix = output_prefix or mbmol_name
    
    output_files = []
    for state_id in range(1, n_states + 1):
        state_letter = chr(ord('A') + state_id - 1)
        output_file = f"{prefix}_state{state_letter}.itp"
        
        # Extract and write this state
        # Note: We need to reload the topology each time to get a fresh copy
        state_mol = extract_state_topology(topfile, mbmol_name, state_id, keep_original_name)
        write_state_itp(state_mol, output_file)
        output_files.append(output_file)
    
    return output_files
