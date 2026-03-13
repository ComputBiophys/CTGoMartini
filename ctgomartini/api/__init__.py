# API module for CTGoMartini
# This module provides the main API classes and functions for creating and managing
# Martini topology files and molecular structures, including support for
# Multiple Basin (MB) and Single Basin (SB) models.

from .MartiniTopology import MartiniTopFile
from .MBMoleculeTop import GenMBPTop
from .SBMoleculeTop import GenSBPTop