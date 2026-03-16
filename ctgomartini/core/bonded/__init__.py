"""Bonded interaction classes for molecular dynamics simulations.

This module defines various bonded interaction types used in molecular dynamics
simulations, including bonds, angles, dihedrals, constraints, virtual sites,
and specialized multi-state interactions.

Examples:
    >>> from ctgomartini.core.bonded import HarmonicBonds, InteractionError
    >>> hb = HarmonicBonds()
    >>> hb.add_interaction(['1', '2', '1', '0.5', '1000.0'])
    >>> 
    >>> # Handle errors properly
    >>> try:
    ...     hb.add_interaction(['1', '2'])  # Too few fields
    ... except InteractionError as e:
    ...     print(f"Invalid interaction: {e}")
"""

from .base import (
    Interaction,
    InteractionError,
    FieldValidationError,
    FuncTypeError,
    register_interaction,
    get_registered_interactions,
    BondedInteractionTypes,
)
from .bonds import (
    HarmonicBonds,
    Constraints,
    Pairs,
    Contacts,
    Exclusions,
)
from .angles import (
    DEG_TO_RAD,
    HarmonicAngles,
    G96Angles,
    RestrictedAngles,
)
from .dihedrals import (
    PeriodicDihedrals,
    HarmonicDihedrals,
    CombinedBendingTorsionPotentials,
    RyckaertBellemansDihedrals,
)
from .virtual_sites import (
    VirtualSite,
    VirtualSitesNCOG,
    VirtualSitesNCOM,
    VirtualSites2,
    VirtualSites2fd,
    VirtualSites3,
    VirtualSites3fd,
    VirtualSites3out,
)
from .multi_state import (
    DEG_TO_RAD as MULTI_DEG_TO_RAD,
    HandlerParam,
    MultiStateHandler,
    MultiStateError,
    MultiAllBonds,
    # Built-in handlers
    G96AngleHandler,
    RestrictedAngleHandler,
    PeriodicDihedralHandler,
    ContactHandler,
)
from .mixing import (
    MixingError,
    EXPInteraction,
    HAMInteraction,
)
from .registry import (
    NonLocalBondedInteractionDict,
    LocalBondedInteractionDict,
    build_interaction_dicts,
)

__all__ = [
    # Base classes and errors
    "Interaction",
    "InteractionError",
    "FieldValidationError",
    "FuncTypeError",
    "register_interaction",
    "get_registered_interactions",
    "BondedInteractionTypes",
    # Bonds
    "HarmonicBonds",
    "Constraints",
    "Pairs",
    "Contacts",
    "Exclusions",
    # Angles
    "DEG_TO_RAD",
    "HarmonicAngles",
    "G96Angles",
    "RestrictedAngles",
    # Dihedrals
    "PeriodicDihedrals",
    "HarmonicDihedrals",
    "CombinedBendingTorsionPotentials",
    "RyckaertBellemansDihedrals",
    # Virtual Sites
    "VirtualSite",
    "VirtualSitesNCOG",
    "VirtualSitesNCOM",
    "VirtualSites2",
    "VirtualSites2fd",
    "VirtualSites3",
    "VirtualSites3fd",
    "VirtualSites3out",
    # Multi-state core
    "MULTI_DEG_TO_RAD",
    "HandlerParam",
    "MultiStateHandler",
    "MultiStateError",
    "MultiAllBonds",
    # Multi-state built-in handlers
    "G96AngleHandler",
    "RestrictedAngleHandler",
    "PeriodicDihedralHandler",
    "ContactHandler",
    # Mixing
    "MixingError",
    "EXPInteraction",
    "HAMInteraction",
    # Registry
    "NonLocalBondedInteractionDict",
    "LocalBondedInteractionDict",
    "build_interaction_dicts",
]
