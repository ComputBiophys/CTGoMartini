"""Topology module for CTGoMartini.

This module provides classes and functions for working with molecular topologies
in Go-Martini coarse-grained simulations.

Main components:
    - TopologyParser: Parse GROMACS-style topology files
    - MartiniTopFile: Build OpenMM systems from topologies
    - Molecule, ForceField: Core data models
    - Interaction classes: Define bonded and non-bonded interactions
    - Generator functions: Create single-basin and multiple-basin topologies

Example:
    >>> from ctgomartini.topology import TopologyParser, MartiniTopFile
    >>> parser = TopologyParser('system.top')
    >>> top = MartiniTopFile('system.top')
    >>> system = top.create_system()
"""

# Core models
from .models import Molecule, ForceField

# Parser
from .parser import TopologyParser

# Builder
from .builder import MartiniTopFile

# Interactions (re-exported for convenience)
from .interactions import (
    # Base
    Interaction,
    InteractionError,
    FieldValidationError,
    FuncTypeError,
    register_interaction,
    get_registered_interactions,
    BondedInteractionTypes,
    # Bonds
    HarmonicBonds,
    Constraints,
    Pairs,
    ContactsLJ,
    Contacts6_12,
    Contacts10_12,
    Exclusions,
    # Angles
    DEG_TO_RAD,
    HarmonicAngles,
    G96Angles,
    RestrictedAngles,
    # Dihedrals
    PeriodicDihedrals,
    HarmonicDihedrals,
    CombinedBendingTorsionPotentials,
    RyckaertBellemansDihedrals,
    # Virtual Sites
    VirtualSite,
    VirtualSitesNCOG,
    VirtualSitesNCOM,
    VirtualSites2,
    VirtualSites2fd,
    VirtualSites3,
    VirtualSites3fd,
    VirtualSites3out,
    # Multi-state
    MultiStateError,
    MultiStateInteraction,
    MultiAllBonds,
    register_multi_handler,
    get_multi_handlers,
    MultiG96Angle,
    MultiRestrictedAngle,
    MultiPeriodicDihedral,
    MultiContactsLJ,
    MultiContacts6_12,
    MultiContacts10_12,
    # Mixing
    MixingError,
    EXPInteraction,
    HAMInteraction,
    # Registry
    NonLocalBondedInteractionDict,
    LocalBondedInteractionDict,
    build_interaction_dicts,
    rebuild_interaction_dicts,
    # Nonbonded
    Nonbonded_interaction,
    ES_self_excl_interaction,
    ES_except_interaction,
    LJ_except_interaction,
    # VSites
    VSiteManager,
    LinearSite,
    COMLinearSite,
    NonLinearSite,
    OutOfPlane,
    NormalizedInPlaneSite,
    NormalizedInPlaneTwoParticleSite,
)

# Generator functions
from .generator import (
    create_sb_topology,
    create_mb_topology,
    combine_atoms,
    combine_bonds_constraints,
    combine_exclusions,
    combine_contacts,
    combine_angles,
    combine_dihedrals,
    extract_contacts_from_topology,
    differentiate_angles,
    differentiate_dihedrals,
)

__all__ = [
    # Core models
    'Molecule',
    'ForceField',
    # Parser
    'TopologyParser',
    # Builder
    'MartiniTopFile',
    # Generator functions
    'create_sb_topology',
    'create_mb_topology',
    'combine_atoms',
    'combine_bonds_constraints',
    'combine_exclusions',
    'combine_contacts',
    'combine_angles',
    'combine_dihedrals',
    'extract_contacts_from_topology',
    'differentiate_angles',
    'differentiate_dihedrals',
    # Base interactions
    'Interaction',
    'InteractionError',
    'FieldValidationError',
    'FuncTypeError',
    'register_interaction',
    'get_registered_interactions',
    'BondedInteractionTypes',
    # Bonds
    'HarmonicBonds',
    'Constraints',
    'Pairs',
    'ContactsLJ',
    'Contacts6_12',
    'Contacts10_12',
    'Exclusions',
    # Angles
    'DEG_TO_RAD',
    'HarmonicAngles',
    'G96Angles',
    'RestrictedAngles',
    # Dihedrals
    'PeriodicDihedrals',
    'HarmonicDihedrals',
    'CombinedBendingTorsionPotentials',
    'RyckaertBellemansDihedrals',
    # Virtual Sites
    'VirtualSite',
    'VirtualSitesNCOG',
    'VirtualSitesNCOM',
    'VirtualSites2',
    'VirtualSites2fd',
    'VirtualSites3',
    'VirtualSites3fd',
    'VirtualSites3out',
    # Multi-state
    'MultiStateError',
    'MultiStateInteraction',
    'MultiAllBonds',
    'register_multi_handler',
    'get_multi_handlers',
    'MultiG96Angle',
    'MultiRestrictedAngle',
    'MultiPeriodicDihedral',
    'MultiContactsLJ',
    'MultiContacts6_12',
    'MultiContacts10_12',
    # Mixing
    'MixingError',
    'EXPInteraction',
    'HAMInteraction',
    # Registry
    'NonLocalBondedInteractionDict',
    'LocalBondedInteractionDict',
    'build_interaction_dicts',
    'rebuild_interaction_dicts',
    # Nonbonded
    'Nonbonded_interaction',
    'ES_self_excl_interaction',
    'ES_except_interaction',
    'LJ_except_interaction',
    # VSites
    'VSiteManager',
    'LinearSite',
    'COMLinearSite',
    'NonLinearSite',
    'OutOfPlane',
    'NormalizedInPlaneSite',
    'NormalizedInPlaneTwoParticleSite',
]
