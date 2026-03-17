"""
Convert constraints to bonds in Martini topology files.

This module provides functionality to convert GROMACS constraints to harmonic
bonds, allowing slight flexibility in regions that were originally constrained.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ctgomartini.core import Molecule


def convert_constraints_to_bonds_in_molecule(
    molecule,
    fc: float = 50000.0,
) -> int:
    """
    Convert all constraints to bonds within a Molecule object (in-place).

    Args:
        molecule: The Molecule object to modify.
        fc: Force constant for the new bonds in kJ/(mol·nm²).

    Returns:
        Number of constraints converted.
    """
    if "constraints" not in molecule._topology:
        return 0

    constraints = molecule._topology["constraints"]
    if not constraints:
        # Remove empty constraints section
        del molecule._topology["constraints"]
        return 0

    converted_count = 0
    for item in constraints:
        # Constraint format: [ai, aj, functype, b0]
        # Bond format: [ai, aj, functype, b0, fc]
        item.append(str(fc))
        molecule._topology["bonds"].append(item)
        converted_count += 1

    # Remove constraints section entirely
    del molecule._topology["constraints"]

    return converted_count


def convert_constraints_to_bonds(
    topology_path: str | Path,
    mol_name: str,
    fc: float = 50000.0,
) -> int:
    """
    Convert constraints to bonds in a topology file (in-place).

    Args:
        topology_path: Path to .itp or .top file.
        mol_name: Molecule type name to modify.
        fc: Force constant for the new bonds in kJ/(mol·nm²).

    Returns:
        Number of constraints converted.

    Raises:
        FileNotFoundError: If topology file doesn't exist.
        ValueError: If molecule not found in topology.
    """
    from ctgomartini.api import MartiniTopFile
    from ctgomartini.utils import write_itp

    topology_path = Path(topology_path)
    if not topology_path.exists():
        raise FileNotFoundError(f"Topology file not found: {topology_path}")

    top = MartiniTopFile(str(topology_path))

    # Handle both naming conventions (_molecule_types vs _moleculeTypes)
    mol = getattr(top, "_molecule_types", {}).get(mol_name)
    if mol is None:
        mol = getattr(top, "_moleculeTypes", {}).get(mol_name)

    if mol is None:
        available = list(
            getattr(top, "_molecule_types", {}).keys()
            or getattr(top, "_moleculeTypes", {}).keys()
        )
        raise ValueError(
            f"Molecule '{mol_name}' not found in {topology_path}. "
            f"Available: {available}"
        )

    converted = convert_constraints_to_bonds_in_molecule(mol, fc)

    if converted > 0:
        write_itp(mol, str(topology_path))

    return converted
