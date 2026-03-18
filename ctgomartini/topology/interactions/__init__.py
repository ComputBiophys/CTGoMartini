"""Bonded interaction classes for molecular dynamics simulations.

This module defines various bonded interaction types used in molecular dynamics
simulations, including bonds, angles, dihedrals, constraints, virtual sites,
and specialized multi-state interactions.

Examples:
    >>> from ctgomartini.topology.interactions import HarmonicBonds, InteractionError
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
    ContactsLJ,
    Contacts6_12,
    Contacts10_12,
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
    MultiStateError,
    MultiStateInteraction,
    MultiAllBonds,
    register_multi_handler,
    get_multi_handlers,
    # Angle handlers
    MultiG96Angle,
    MultiRestrictedAngle,
    # Dihedral handlers
    MultiPeriodicDihedral,
    # Contact handlers
    MultiContactsLJ,
    MultiContacts6_12,
    MultiContacts10_12,
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
    rebuild_interaction_dicts,
)
from .nonbonded import (
    Interaction as NonbondedInteractionBase,
    Nonbonded_interaction,
    ES_self_excl_interaction,
    ES_except_interaction,
    LJ_except_interaction,
)
from .vsites import (
    LinearSite,
    COMLinearSite,
    NonLinearSite,
    OutOfPlane,
    NormalizedInPlaneSite,
    NormalizedInPlaneTwoParticleSite,
    VSiteManager,
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
    "ContactsLJ",
    "Contacts6_12",
    "Contacts10_12",
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
    # Virtual Site Data Classes
    "LinearSite",
    "COMLinearSite",
    "NonLinearSite",
    "OutOfPlane",
    "NormalizedInPlaneSite",
    "NormalizedInPlaneTwoParticleSite",
    "VSiteManager",
    # Multi-state core
    "MULTI_DEG_TO_RAD",
    "MultiStateError",
    "MultiStateInteraction",
    "MultiAllBonds",
    "register_multi_handler",
    "get_multi_handlers",
    # Multi-state handlers
    "MultiG96Angle",
    "MultiRestrictedAngle",
    "MultiPeriodicDihedral",
    "MultiContactsLJ",
    "MultiContacts6_12",
    "MultiContacts10_12",
    # Mixing
    "MixingError",
    "EXPInteraction",
    "HAMInteraction",
    # Registry
    "NonLocalBondedInteractionDict",
    "LocalBondedInteractionDict",
    "build_interaction_dicts",
    "rebuild_interaction_dicts",
    # Nonbonded
    "NonbondedInteractionBase",
    "Nonbonded_interaction",
    "ES_self_excl_interaction",
    "ES_except_interaction",
    "LJ_except_interaction",
]
