"""Contact map utilities for Go-Martini topology generation.

This module provides functions for handling contact map files from various sources:
- rCSU web-server .map format
- martinize2 .out format  
- Automatic generation via OVrCSU algorithm
"""

from __future__ import annotations

import os
import subprocess
from typing import Literal

import ctgomartini
from ctgomartini.utils import convert_map_format


def detect_contact_file_format(filepath: str) -> Literal['map', 'out', 'unknown']:
    """Detect contact file format based on file extension.
    
    Args:
        filepath: Path to the contact file.
        
    Returns:
        Format type: 'map' (rCSU format), 'out' (martinize2 format), or 'unknown'.
        
    Raises:
        FileNotFoundError: If file does not exist.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Contact file not found: {filepath}")
    
    ext = os.path.splitext(filepath)[1].lower()
    
    if ext == '.map':
        return 'map'
    elif ext == '.out':
        return 'out'
    else:
        return 'unknown'


def run_ovrcsu(aa_strfile: str, state_name: str) -> str:
    """Run OVrCSU to generate contacts automatically.
    
    This function runs the OV+rCSU contact detection algorithm on the input
    structure file to generate a contact map.
    
    Args:
        aa_strfile: Path to the all-atom structure file.
        state_name: Name of the state for output file naming.
        
    Returns:
        Path to the generated contact map file (OV format).
        
    Raises:
        Exception: If OVrCSU execution fails.
    """
    ovrcsu_script = os.path.join(ctgomartini.__path__[0], 'utils/contact_map.py')
    output_map = f'{state_name}_OVrCSU.map'
    
    cmd = f'python {ovrcsu_script} -f {aa_strfile} -o {output_map}'
    print(f'Running OVrCSU: {cmd}')
    
    result = subprocess.run(
        cmd,
        shell=True, capture_output=True, encoding='utf-8'
    )
    
    if result.returncode != 0:
        print(f"OVrCSU stderr: {result.stderr}")
        raise Exception(f'OVrCSU failed for {state_name}')
    
    print(f'Generated OVrCSU contacts: {output_map}')
    return output_map


def gen_go_contacts(
    aa_strfile: str,
    map_file: str | None,
    state_name: str,
    go_eps: float = 12.0,
    min_seq_distance: int = 4,
) -> str | None:
    """Generate or prepare Go-like contacts for a protein structure.
    
    This function supports three modes:
    1. Auto mode (map_file='auto' or None): Run OVrCSU to generate contacts
    2. User-provided .map file: Convert rCSU web-server format to martinize2 format
    3. User-provided .out file: Use directly (already in martinize2 format)
    
    Args:
        aa_strfile: Path to the all-atom structure file.
        map_file: Path to contact file, or 'auto' for automatic generation.
                 None is treated as 'auto'.
        state_name: Name of the state/molecule type.
        go_eps: Epsilon value for Go contacts (default: 12.0).
        min_seq_distance: Minimum sequence distance (PDB numbering) for intra-chain
                         contacts. Contacts with distance < min_seq_distance within
                         the same chain are filtered out. Inter-chain contacts are
                         always kept. Set to 0 to disable filtering. (default: 4)
        
    Returns:
        Path to the contact_map.out file for martinize2, or None if no contacts.
        
    Raises:
        FileNotFoundError: If provided contact file does not exist.
        ValueError: If contact file format is unknown.
    """
    # Auto mode: run OVrCSU to generate contacts
    if map_file is None or map_file.lower() == 'auto':
        print(f'\nAuto-generating contacts for {state_name} using OVrCSU')
        ov_map = run_ovrcsu(aa_strfile, state_name)
        output_file = f'contact_map_{state_name}.out'
        convert_map_format(
            input_file=ov_map,
            output_file=output_file,
            pdb_name=os.path.basename(aa_strfile),
            force=True,
            min_seq_distance=min_seq_distance,
        )
        print(f'Generated Go contacts: {output_file}')
        return output_file
    
    # User-provided file: detect format and convert if needed
    if not os.path.exists(map_file):
        raise FileNotFoundError(f"Contact file not found: {map_file}")
    
    file_format = detect_contact_file_format(map_file)
    output_file = f'contact_map_{state_name}.out'
    
    if file_format == 'out':
        # Already in martinize2 format, use directly
        print(f'\nUsing provided contact file (martinize2 format): {map_file}')
        # Copy to expected filename
        if os.path.abspath(map_file) != os.path.abspath(output_file):
            import shutil
            shutil.copy(map_file, output_file)
        return output_file
    
    elif file_format == 'map':
        # rCSU .map format, need conversion
        print(f'\nConverting rCSU .map format to martinize2 format: {map_file}')
        convert_map_format(
            input_file=map_file,
            output_file=output_file,
            pdb_name=os.path.basename(aa_strfile),
            force=True,
            min_seq_distance=min_seq_distance,
        )
        print(f'Generated Go contacts: {output_file}')
        return output_file
    
    else:
        raise ValueError(
            f"Unknown contact file format for {map_file}. "
            "Expected .map (rCSU format) or .out (martinize2 format)."
        )
