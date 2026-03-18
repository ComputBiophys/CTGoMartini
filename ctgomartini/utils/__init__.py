"""
Utilities module for CTGoMartini.

Provides utility functions for file I/O, topology writing, contact map generation,
virtual site generation, contact map conversion, and constraints conversion.
"""

from .write_itp import write_itp
from .convert_long_short_elastic_bonds import convert_long_short_elastic_bonds, bb_distance
from .convert_map_format import convert_map_format
from .constraints_to_bonds import (
    convert_constraints_to_bonds,
    convert_constraints_to_bonds_in_molecule,
)
from .contact_map import (
    compute_atom_contacts,
    compute_residue_contacts,
)

__all__ = [
    'write_itp',
    'convert_long_short_elastic_bonds',
    'bb_distance',
    'convert_map_format',
    'convert_constraints_to_bonds',
    'convert_constraints_to_bonds_in_molecule',
    'compute_atom_contacts',
    'compute_residue_contacts',
]
