"""
Convert constraints to bonds in a Martini topology.

This module provides functionality to convert constraint entries to bond entries
in a Martini topology file. Constraints are converted to bonds with a force
constant of 10000 kJ/mol/nm^2, which effectively maintains the constrained
distance while allowing the simulation to use the standard bond machinery.
"""

from __future__ import annotations

from ctgomartini.api import MartiniTopFile
from ctgomartini.core import Molecule
from ctgomartini.utils import write_itp


def convert_constraints_to_bonds(mol: Molecule, force_constant: str = '10000') -> None:
    """Convert all constraints in a molecule to bonds.

    This function iterates through all constraint entries in the molecule's
    topology, appends a force constant to each constraint, and moves them
    to the bonds section. The constraints section is then removed from the
    topology.

    Args:
        mol: The molecule object containing the topology to modify.
            Must have a _topology attribute with 'constraints' and 'bonds' keys.
        force_constant: The force constant to assign to the converted bonds.
            Default is '10000' (kJ/mol/nm^2).

    Raises:
        AttributeError: If the molecule does not have a _topology attribute.
        KeyError: If the topology does not have the expected structure.

    Example:
        >>> top = MartiniTopFile('system.top')
        >>> mol = top._moleculeTypes['TREK1']
        >>> convert_constraints_to_bonds(mol)
    """
    # Convert each constraint to a bond by appending the force constant
    for item in mol._topology['constraints']:
        item.append(force_constant)
        mol._topology['bonds'].append(item)
    
    # Remove the constraints section from the topology
    mol._topology['constraints'] = []
    del mol._topology['constraints']


def main(
    topology_file: str = 'system.top',
    molecule_name: str = 'TREK1',
    force_constant: str = '10000'
) -> None:
    """Main entry point for converting constraints to bonds.

    Loads a Martini topology file, converts all constraints to bonds for a
    specified molecule, and writes the result to ITP files.

    Args:
        topology_file: Path to the GROMACS topology file (.top) to load.
            Default is 'system.top'.
        molecule_name: Name of the molecule type to process.
            Default is 'TREK1'.
        force_constant: Force constant for the converted bonds.
            Default is '10000' (kJ/mol/nm^2).

    Raises:
        FileNotFoundError: If the topology file does not exist.
        KeyError: If the specified molecule type is not found in the topology.

    Example:
        >>> main('system.top', 'TREK1')
        # Creates TREK1.itp and TREK1_params.itp
    """
    # Load the topology file
    top = MartiniTopFile(topology_file)
    
    # Get the specified molecule type
    mol = top._moleculeTypes[molecule_name]
    
    # Convert constraints to bonds
    convert_constraints_to_bonds(mol, force_constant)
    
    # Write the modified molecule to ITP files
    write_itp(mol)


if __name__ == '__main__':
    main()
