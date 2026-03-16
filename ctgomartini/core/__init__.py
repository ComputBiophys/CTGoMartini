"""
Core module for CTGoMartini.

Provides core classes and functions for topology parsing, force field handling,
molecule representation, and interaction definitions.
"""

from .Topology import Topology
from .Molecule import Molecule, Molecule_Category, Moleculetype, Atoms, Multiple_basin
from .ForceField import ForceField, FF_Category, Defaults, Atomtypes, Nonbond_params
from .bonded import (
    # Base classes and errors
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
    Contacts,
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
    MULTI_DEG_TO_RAD,
    HandlerParam,
    MultiStateHandler,
    MultiStateError,
    MultiAllBonds,
    G96AngleHandler,
    RestrictedAngleHandler,
    PeriodicDihedralHandler,
    ContactHandler,
    # Mixing
    MixingError,
    EXPInteraction,
    HAMInteraction,
    # Registry
    NonLocalBondedInteractionDict,
    LocalBondedInteractionDict,
    build_interaction_dicts,
)
from .Nonbonded_interaction import (
    Nonbonded_interaction, ES_self_excl_interaction,
    ES_except_interaction, LJ_except_interaction
)
from .CombineMols import (
    CombineMols, SameListList, SameList,
    Extract_contacts_from_top, GetAtomNames, GetAngleDiehdralType,
    DifferentiateAngles, DifferentiateDihedrals,
    CombineDict, ForceItemFloat, ForceListFloat,
    Calculate_DiffDihedral
)
from .vsites import (
    LinearSite, COMLinearSite, NonLinearSite,
    OutOfPlane, NormalizedInPlaneSite, NormalizedInPlaneTwoParticleSite,
    VSiteManager
)
