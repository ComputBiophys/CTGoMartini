"""
Core module for CTGoMartini.

Provides core classes and functions for topology parsing, force field handling,
molecule representation, and interaction definitions.
"""

from .Topology import Topology
from .Molecule import Molecule, Molecule_Category, Moleculetype, Atoms, Multiple_basin
from .ForceField import ForceField, FF_Category, Defaults, Atomtypes, Nonbond_params
from .Bonded_interaction import (
    Interaction, BondedInteraction_types,
    Harmonic_bonds, Harmonic_angles, G96_angles, Restricted_angles,
    Periodic_dihedrals, Harmonic_dihedrals,
    Combined_bending_torsion_potentials, Ryckaert_Bellemans_dihedrals,
    Constraints, Pairs, Contacts, Exclusions,
    VirtualSite, Virtual_sitesn_COG, Virtual_sitesn_COM,
    Virtual_sites2, Virtual_sites2fd, Virtual_sites3, Virtual_sites3fd, Virtual_sites3out,
    Multi_allbonds, Multi_harmonic_angles, Multi_g96_angles, Multi_restricted_angles,
    Multi_periodic_dihedrals, Multi_harmonic_dihedrals, Multi_contacts,
    EXP_Interaction, HAM_Interaction,
    NonLocal_BondedInteraction_dict, Local_BondedInteraction_dict
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
