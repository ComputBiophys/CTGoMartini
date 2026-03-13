"""
Utilities module for CTGoMartini.

Provides utility functions for file I/O, topology writing, input reading,
and virtual site generation.
"""

from .WriteItp import write_itp
from .ConvertLongShortElasticBonds import convert_long_short_elastic_bonds, bb_distance
from .ReadInp import read_inputs, _OpenMMReadInputs
from .Create_goVirt_for_multimer import create_go_virt_for_multimer
