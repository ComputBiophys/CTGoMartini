"""PDB structure validation module for CTGoMartini."""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import mdtraj as md


class AtomRecord(NamedTuple):
    """Atom information extracted from PDB."""
    serial: int      # PDB serial number
    name: str        # Atom name
    resname: str     # Residue name
    resid: int       # Residue sequence number
    chainid: str     # Chain identifier


class PDBData(NamedTuple):
    """Parsed PDB data."""
    name: str
    filepath: Path
    atoms: list[AtomRecord]
    n_residues: int
    n_chains: int
    chain_ids: list[str]


class ValidationError(Exception):
    """Raised when PDB validation fails."""
    pass


def _is_not_hydrogen(atom_name: str) -> bool:
    """Check if atom is not hydrogen (by name)."""
    return not atom_name.startswith('H')


def _parse_pdb(filepath: Path, name: str) -> PDBData:
    """Parse PDB file and extract relevant data."""
    traj = md.load_pdb(str(filepath))
    topology = traj.topology
    
    # Extract atom records (excluding hydrogen)
    atoms: list[AtomRecord] = []
    for atom in topology.atoms:
        if _is_not_hydrogen(atom.name):
            residue = atom.residue
            chain = residue.chain
            atoms.append(AtomRecord(
                serial=atom.index + 1,  # Convert to 1-based serial
                name=atom.name,
                resname=residue.name,
                resid=residue.resSeq,
                chainid=chain.chain_id if chain else ' '
            ))
    
    # Count residues and chains
    n_residues = topology.n_residues
    n_chains = topology.n_chains
    
    # Get chain IDs in order of first appearance
    chain_ids = []
    for atom in topology.atoms:
        chain_id = atom.residue.chain.chain_id if atom.residue.chain else ' '
        if chain_id not in chain_ids:
            chain_ids.append(chain_id)
    
    return PDBData(
        name=name,
        filepath=filepath,
        atoms=atoms,
        n_residues=n_residues,
        n_chains=n_chains,
        chain_ids=chain_ids
    )


def _check_atom_count(ref: PDBData, other: PDBData) -> None:
    """Check if atom counts match."""
    if len(ref.atoms) != len(other.atoms):
        raise ValidationError(
            f"Atom count mismatch: {ref.name} has {len(ref.atoms)} atoms, "
            f"{other.name} has {len(other.atoms)} atoms"
        )


def _check_residue_count(ref: PDBData, other: PDBData) -> None:
    """Check if residue counts match."""
    if ref.n_residues != other.n_residues:
        raise ValidationError(
            f"Residue count mismatch: {ref.name} has {ref.n_residues} residues, "
            f"{other.name} has {other.n_residues} residues"
        )


def _check_chain_count(ref: PDBData, other: PDBData) -> None:
    """Check if chain counts match."""
    if ref.n_chains != other.n_chains:
        raise ValidationError(
            f"Chain count mismatch: {ref.name} has {ref.n_chains} chains, "
            f"{other.name} has {other.n_chains} chains"
        )


def _check_chain_ids(ref: PDBData, other: PDBData) -> None:
    """Check if chain ID sequences match."""
    if ref.chain_ids != other.chain_ids:
        raise ValidationError(
            f"Chain ID mismatch: {ref.name} has {ref.chain_ids}, "
            f"{other.name} has {other.chain_ids}"
        )


def _check_atoms_detail(ref: PDBData, other: PDBData) -> None:
    """Check each atom's properties match."""
    for i, (a1, a2) in enumerate(zip(ref.atoms, other.atoms)):
        if a1.resname != a2.resname:
            raise ValidationError(
                f"Atom {i}: residue name mismatch - "
                f"{ref.name} has {a1.resname}, {other.name} has {a2.resname}"
            )
        if a1.resid != a2.resid:
            raise ValidationError(
                f"Atom {i}: residue ID mismatch - "
                f"{ref.name} has {a1.resid}, {other.name} has {a2.resid}"
            )
        if a1.name != a2.name:
            raise ValidationError(
                f"Atom {i}: atom name mismatch - "
                f"{ref.name} has {a1.name}, {other.name} has {a2.name}"
            )
        if a1.chainid != a2.chainid:
            raise ValidationError(
                f"Atom {i}: chain ID mismatch - "
                f"{ref.name} has {a1.chainid}, {other.name} has {a2.chainid}"
            )


def validate_pdb_compatibility(
    pdb_files: list[str | Path],
    state_names: list[str] | None = None
) -> None:
    """Validate PDB files are compatible for multi-basin topology generation.
    
    Checks are performed in order:
    1. Atom count (excluding hydrogen atoms)
    2. Residue count
    3. Chain count
    4. Chain ID sequence
    5. Each atom's resname, resid, atomname, chainid
    
    Args:
        pdb_files: List of PDB file paths
        state_names: Optional names for each structure
        
    Raises:
        ValidationError: If any check fails
        FileNotFoundError: If PDB file does not exist
    """
    if len(pdb_files) < 2:
        return
    
    if state_names is None:
        state_names = [f"Structure{i+1}" for i in range(len(pdb_files))]
    
    # Parse all PDB files
    structures: list[PDBData] = []
    for filepath, name in zip(pdb_files, state_names):
        path = Path(filepath)
        if not path.exists():
            raise FileNotFoundError(f"PDB file not found: {path.absolute()}")
        structures.append(_parse_pdb(path, name))
    
    # Compare reference (first) against all others
    ref = structures[0]
    for other in structures[1:]:
        _check_atom_count(ref, other)
        _check_residue_count(ref, other)
        _check_chain_count(ref, other)
        _check_chain_ids(ref, other)
        _check_atoms_detail(ref, other)
