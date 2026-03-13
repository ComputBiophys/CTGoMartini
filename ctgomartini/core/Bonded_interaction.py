"""Bonded interaction classes for molecular dynamics simulations.

This module defines various bonded interaction types used in molecular dynamics
simulations, including bonds, angles, dihedrals, constraints, virtual sites,
and specialized multi-state interactions.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import openmm as mm

from .vsites import (
    COMLinearSite,
    LinearSite,
    NonLinearSite,
    NormalizedInPlaneSite,
    NormalizedInPlaneTwoParticleSite,
    OutOfPlane,
    VSiteManager,
)

if TYPE_CHECKING:
    import openmm.unit as unit


class Interaction:
    """Base class for bonded interactions.
    
    This class serves as the foundation for all interaction types in the
    molecular dynamics force field. It defines the common interface and
    attributes shared by all interaction classes.
    
    Attributes:
        name: Name of the interaction type.
        description: Description of the interaction.
        category: Category name (e.g., 'bonds', 'angles', 'dihedrals').
        mm_force: OpenMM force object or None.
        type_label: List containing the field index and type identifier(s).
        contents: List storing interaction entries.
        intermolecule_sharing: Whether this interaction can be shared between molecules.
    """

    def __init__(
        self,
        name: str,
        description: str,
        category: str,
        mm_force: mm.Force | None,
        type_label: list[int | str] | None,
    ) -> None:
        """Initialize the interaction.
        
        Args:
            name: Name of the interaction type.
            description: Description of the interaction.
            category: Category name (e.g., 'bonds', 'angles').
            mm_force: OpenMM force object or None.
            type_label: List with [field_index, type_id, ...] for validation.
        """
        self.name: str = name
        self.description: str = description
        self.category: str = category
        self.mm_force: mm.Force | None = mm_force
        self.type_label: list[int | str] | None = type_label

        # default parameters
        self.contents: list[Any] = []
        self.intermolecule_sharing: bool = True

    def __str__(self) -> str:
        """Return string representation of the interaction."""
        return f"{self.name}: {self.description}"

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one interaction to mm_force.
        
        Args:
            fields: List of string values from the topology file.
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
            
        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields.
        
        Args:
            atoms: List of atom information tuples.
            fields: List of string values from the topology file.
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
            
        Returns:
            List of exception specifications.
        """
        return []


# Register interaction types
BondedInteraction_types: list[type[Interaction]] = []


# Bonds
class Harmonic_bonds(Interaction):
    """Harmonic bond potential interaction.
    
    Implements a harmonic bond potential of the form:
    V(r) = 0.5 * k * (r - r0)^2
    
    where r is the bond length, r0 is the equilibrium length,
    and k is the force constant.
    """

    def __init__(self) -> None:
        """Initialize harmonic bonds interaction."""
        super().__init__(
            name='harmonic_bonds',
            description='Harmonic bond potential 1/2*k*(r-r0)^2: atomid1, atomid2, functype, length, and k',
            category='bonds',
            mm_force=mm.HarmonicBondForce(),
            type_label=[2, "1"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one harmonic bond interaction.
        
        Args:
            fields: List of [atomid1, atomid2, functype, length, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1]
        ), "harmonic_bonds requires 5 items and the functype is 1"
        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            float(fields[3]),
            float(fields[4]),
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        return [
            [
                base_atom_index + int(fields[0]) + offset,
                base_atom_index + int(fields[1]) + offset,
                0,
                0,
                0,
            ]
        ]


BondedInteraction_types.append(Harmonic_bonds)


# Angles
class Harmonic_angles(Interaction):
    """Harmonic angle potential interaction.
    
    Implements a harmonic angle potential of the form:
    V(theta) = 0.5 * k * (theta - theta0)^2
    
    where theta is the bond angle, theta0 is the equilibrium angle,
    and k is the force constant.
    """

    def __init__(self) -> None:
        """Initialize harmonic angles interaction."""
        super().__init__(
            name='harmonic_angles',
            description='Harmonic angle potential 1/2*k*(theta-theta0)^2: atomid1, atomid2, atomid3, functype, angle, and k',
            category='angles',
            mm_force=mm.HarmonicAngleForce(),
            type_label=[3, "1"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one harmonic angle interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, functype, angle, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 6 and fields[self.type_label[0]] == self.type_label[1]
        ), "harmonic_angles requires 6 items and the functype is 1"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            float(fields[4]) * degToRad,
            float(fields[5]),
        )


BondedInteraction_types.append(Harmonic_angles)


class G96_angles(Interaction):
    """G96 angle potential interaction.
    
    Implements the GROMOS 96 angle potential of the form:
    V(theta) = 0.5 * k * (cos(theta) - cos(theta0))^2
    
    This is a more numerically stable form for certain angle ranges.
    """

    def __init__(self) -> None:
        """Initialize G96 angles interaction."""
        super().__init__(
            name='g96_angles',
            description='g96 angles 0.5*k*(cos(theta)-cos(theta0))^2: atomid1, atomid2, atomid3, functype, angle, and k',
            category='angles',
            mm_force=mm.CustomAngleForce(
                "0.5 * k * (cos(theta) - cos(theta0))^2"
            ),
            type_label=[3, "2"],
        )
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one G96 angle interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, functype, angle, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 6 and fields[self.type_label[0]] == self.type_label[1]
        ), "g96_angles requires 6 items and the functype is 2"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            [float(fields[4]) * degToRad, float(fields[5])],
        )


BondedInteraction_types.append(G96_angles)


class Restricted_angles(Interaction):
    """Restricted angle potential interaction.
    
    Implements a restricted bending potential of the form:
    V(theta) = 0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2
    
    This prevents the angle from reaching 0 or 180 degrees.
    """

    def __init__(self) -> None:
        """Initialize restricted angles interaction."""
        super().__init__(
            name='restricted_angles',
            description='Restricted angles 0.5*k*(cos(theta)-cos(theta0))^2 /sin(theta)^2: atomid1, atomid2, atomid3, functype, angle, and k',
            category='angles',
            mm_force=mm.CustomAngleForce(
                "0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2"
            ),
            type_label=[3, "10"],
        )
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one restricted angle interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, functype, angle, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 6 and fields[self.type_label[0]] == self.type_label[1]
        ), "restricted_angles requires 6 items and the functype is 10"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            [float(fields[4]) * degToRad, float(fields[5])],
        )


BondedInteraction_types.append(Restricted_angles)


# Dihedrals
class Periodic_dihedrals(Interaction):
    """Periodic dihedral (torsion) interaction.
    
    Implements a periodic dihedral potential of the form:
    V(phi) = k * (1 + cos(n * phi - phi0))
    
    where phi is the dihedral angle, phi0 is the phase angle,
    n is the periodicity, and k is the force constant.
    """

    def __init__(self) -> None:
        """Initialize periodic dihedrals interaction."""
        super().__init__(
            name='periodic_dihedrals',
            description='Periodic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, functype, theta, k, and phase',
            category='dihedrals',
            mm_force=mm.PeriodicTorsionForce(),
            type_label=[4, "1", "4", "9"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one periodic dihedral interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, atomid4, functype, theta, k, phase].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 8 and fields[self.type_label[0]] in self.type_label[1:]
        ), "periodic_dihedrals requires 8 items and the functype is 1 or 4 or 9"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            float(fields[7]),
            float(fields[5]) * degToRad,
            float(fields[6]),
        )


BondedInteraction_types.append(Periodic_dihedrals)


class Harmonic_dihedrals(Interaction):
    """Harmonic dihedral (improper) interaction.
    
    Implements a harmonic dihedral potential for improper dihedrals
    with proper periodic boundary handling.
    """

    def __init__(self) -> None:
        """Initialize harmonic dihedrals interaction."""
        super().__init__(
            name='harmonic_dihedrals',
            description='Harmonic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, functype, theta, and k',
            category='dihedrals',
            mm_force=mm.CustomTorsionForce(
                "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                % math.pi
            ),
            type_label=[4, "2"],
        )
        self.mm_force.addPerTorsionParameter("theta0")
        self.mm_force.addPerTorsionParameter("k")

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one harmonic dihedral interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, atomid4, functype, theta, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1]
        ), "harmonic_dihedrals requires 7 items and the functype is 2"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            [float(fields[5]) * degToRad, float(fields[6])],
        )


BondedInteraction_types.append(Harmonic_dihedrals)


class Combined_bending_torsion_potentials(Interaction):
    """Combined bending-torsion potential interaction.
    
    Implements a potential coupling bond bending and torsional motions:
    V = k * sin(theta0)^3 * sin(theta1)^3 * sum(ai * cos(phi)^i)
    
    This is used for more accurate conformational modeling.
    """

    def __init__(self) -> None:
        """Initialize combined bending-torsion potential interaction."""
        super().__init__(
            name='combined_bending_torsion_potentials',
            description='Combined bending-torsion potentials: atomid1, atomid2, atomid3, atomid4, functype, k, a0, a1, a2, a3, and a4',
            category='dihedrals',
            mm_force=None,
            type_label=[4, "11"],
        )
        self.mm_force = mm.CustomCompoundBondForce(
            4,
            "k*sintheta0^3*sintheta1^3*(a0 + a1*cosphi + a2*cosphi^2 + a3*cosphi^3 + a4*cosphi^4); "
            "sintheta0 = sin(angle(p1, p2, p3));"
            "sintheta1 = sin(angle(p2, p3, p4));"
            "cosphi = cos(dihedral(p1, p2, p3, p4));",
        )
        self.mm_force.addPerBondParameter("k")
        self.mm_force.addPerBondParameter("a0")
        self.mm_force.addPerBondParameter("a1")
        self.mm_force.addPerBondParameter("a2")
        self.mm_force.addPerBondParameter("a3")
        self.mm_force.addPerBondParameter("a4")

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one combined bending-torsion interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, atomid4, functype, k, a0, a1, a2, a3, a4].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 11 and fields[self.type_label[0]] == self.type_label[1]
        ), "combined_bending_torsion_potentials requires 11 items and the functype is 11"
        k = float(fields[5])
        a0 = float(fields[6])
        a1 = float(fields[7])
        a2 = float(fields[8])
        a3 = float(fields[9])
        a4 = float(fields[10])
        self.mm_force.addBond(
            [
                base_atom_index + int(fields[0]) + offset,
                base_atom_index + int(fields[1]) + offset,
                base_atom_index + int(fields[2]) + offset,
                base_atom_index + int(fields[3]) + offset,
            ],
            [k, a0, a1, a2, a3, a4],
        )


BondedInteraction_types.append(Combined_bending_torsion_potentials)


class Ryckaert_Bellemans_dihedrals(Interaction):
    """Ryckaert-Bellemans dihedral interaction.
    
    Implements the Ryckaert-Bellemans potential for aliphatic chains:
    V(phi) = sum(Ci * cos(phi)^i) for i = 0 to 5
    
    This is commonly used for alkane chains.
    """

    def __init__(self) -> None:
        """Initialize Ryckaert-Bellemans dihedrals interaction."""
        super().__init__(
            name='Ryckaert_Bellemans_dihedrals',
            description='Ryckaert_Bellemans_dihedrals: atomid1, atomid2, atomid3, atomid4, functype, C0, C1, C2, C3, C4, and C5',
            category='dihedrals',
            mm_force=mm.RBTorsionForce(),
            type_label=[4, "3", "5"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one Ryckaert-Bellemans dihedral interaction.
        
        Args:
            fields: List of [atomid1, atomid2, atomid3, atomid4, functype, C0, C1, C2, C3, C4, C5].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 11 and fields[self.type_label[0]] in self.type_label[1:]
        ), "Ryckaert_Bellemans_dihedrals requires 11 items and the functype is 3 or 5"
        c = [float(x) for x in fields[5:11]]
        if fields[self.type_label[0]] == "5":
            # Convert Fourier coefficients to RB coefficients.
            c = [
                c[1] + 0.5 * (c[0] + c[2]),
                0.5 * (-c[0] + 3 * c[2]),
                -c[1] + 4 * c[3],
                -2 * c[2],
                -4 * c[3],
                0,
            ]
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            c[0],
            c[1],
            c[2],
            c[3],
            c[4],
            c[5],
        )


BondedInteraction_types.append(Ryckaert_Bellemans_dihedrals)


# Constraints
class Constraints(Interaction):
    """Distance constraint interaction.
    
    Adds a rigid constraint between two atoms at a fixed distance.
    This is implemented as an OpenMM constraint rather than a force.
    """

    def __init__(self, sys: mm.System | None = None) -> None:
        """Initialize constraints interaction.
        
        Args:
            sys: OpenMM system to add constraints to.
        """
        super().__init__(
            name='constraints',
            description='Constraints (sys.addConstraint): atomid1, atomid2, functype, and length',
            category='constraints',
            mm_force=None,
            type_label=[2, "1"],
        )
        self.sys: mm.System | None = sys

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one constraint.
        
        Args:
            fields: List of [atomid1, atomid2, functype, length].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        # Sometimes the fields in the constraints category have 5 items,
        # which the 5th item perhaps has no meaning according to the
        # Martini-openmm project.
        assert (
            len(fields) >= 4 and fields[self.type_label[0]] == self.type_label[1]
        ), "constraints requires at least 4 items and the functype is 1"
        self.sys.addConstraint(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            float(fields[3]),
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        return [
            [
                base_atom_index + int(fields[0]) + offset,
                base_atom_index + int(fields[1]) + offset,
                0,
                0,
                0,
            ]
        ]


BondedInteraction_types.append(Constraints)


# Pairs
class Pairs(Interaction):
    """Pair interaction for 1-4 or special pairs.
    
    Implements Lennard-Jones and electrostatic interactions between
    atom pairs that are not directly bonded but need special treatment.
    """

    def __init__(self, epsilon_r: float = 15, use_sigma_eps: bool = True) -> None:
        """Initialize pairs interaction.
        
        Args:
            epsilon_r: Relative dielectric constant.
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False).
        """
        super().__init__(
            name='pairs',
            description='Pairs: atomid1, atomid2, functype, (C6, C12)',
            category='pairs',
            mm_force=None,
            type_label=[2, "1"],
        )
        self.epsilon_r: float = epsilon_r
        self.use_sigma_eps: bool = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "LJ + ES;"
            "LJ = C12/r^12 - C6/r^6;"
            "ES = f/epsilon_r*q1*q2 * 1/r;"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")
        self.mm_force.addPerBondParameter("q1")
        self.mm_force.addPerBondParameter("q2")

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Validate pair interaction fields.
        
        Args:
            fields: List of [atomid1, atomid2, functype, C6/sigma, C12/epsilon].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1]
        ), "pairs requires 5 items and the functype is 1"

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Add pair interaction and return empty exception list."""
        if self.use_sigma_eps:
            sigma = float(fields[3])
            eps = float(fields[4])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[3])
            C12 = float(fields[4])
        q1 = float(atoms[int(fields[0]) - 1][6])
        q2 = float(atoms[int(fields[1]) - 1][6])
        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            [C12, C6, q1, q2],
        )
        return []


BondedInteraction_types.append(Pairs)


# Contacts
class Contacts(Interaction):
    """Contact (soft-core) interaction.
    
    Implements a shifted Lennard-Jones potential for contacts:
    V(r) = step(rcut - r) * (LJ(r) - LJ(rcut))
    
    This smoothly goes to zero at the cutoff distance.
    """

    def __init__(
        self,
        nonbonded_cutoff: unit.Quantity = 1.1 * mm.unit.nanometer,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize contacts interaction.
        
        Args:
            nonbonded_cutoff: Cutoff distance for the potential.
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False).
        """
        super().__init__(
            name='contacts',
            description='Contacts Lenard-Jones Potential(r - cutoff): atomid1, atomid2, functype, C6/sigma, C12/epsilon',
            category='contacts',
            mm_force=None,
            type_label=[2, "1"],
        )
        self.nonbonded_cutoff: unit.Quantity = nonbonded_cutoff
        self.use_sigma_eps: bool = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={self.nonbonded_cutoff.value_in_unit(mm.unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one contact interaction.
        
        Args:
            fields: List of [atomid1, atomid2, functype, C6/sigma, C12/epsilon].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1]
        ), "contacts requires 5 items and the functype is 1"
        if self.use_sigma_eps:
            sigma = float(fields[3])
            eps = float(fields[4])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[3])
            C12 = float(fields[4])

        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            [C12, C6],
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        q1 = float(atoms[int(fields[0]) - 1][6])
        q2 = float(atoms[int(fields[1]) - 1][6])
        return [
            [
                base_atom_index + int(fields[0]) + offset,
                base_atom_index + int(fields[1]) + offset,
                q1 * q2,
                0,
                0,
            ]
        ]


BondedInteraction_types.append(Contacts)


# virtual_sitesn
class VirtualSite(Interaction):
    """Base class for virtual site interactions.
    
    Virtual sites are particles whose positions are calculated from
    the positions of other atoms rather than being integrated.
    
    Attributes:
        vsites: Virtual site manager for registering sites.
    """

    def __init__(
        self,
        vsites: VSiteManager | None,
        name: str,
        description: str,
        category: str,
        mm_force: mm.Force | None,
        type_label: list[int | str] | None,
    ) -> None:
        """Initialize virtual site interaction.
        
        Args:
            vsites: Virtual site manager.
            name: Name of the virtual site type.
            description: Description of the virtual site.
            category: Category name.
            mm_force: OpenMM force object (None for virtual sites).
            type_label: Type label for validation.
        """
        super().__init__(
            name=name,
            description=description,
            category=category,
            mm_force=mm_force,
            type_label=type_label,
        )
        self.vsites: VSiteManager | None = vsites


class Virtual_sitesn_COG(VirtualSite):
    """Center of geometry virtual site (N-body).
    
    The virtual site position is the geometric center of N atoms.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize COG virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sitesn_COG',
            description="N-body virutal site (COG): atomid, functype, atomid1, atomid2, ..., atomidn",
            category='virtual_sitesn',
            mm_force=None,
            type_label=[1, "1"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one COG virtual site.
        
        Args:
            fields: List of [atomid, functype, atomid1, atomid2, ..., atomidn].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) >= 3 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_sitesn_COG requires at least 3 items and the functype is 1"
        index = int(fields[0]) + base_atom_index + offset
        from_atoms = [int(field) + base_atom_index + offset for field in fields[2:]]

        w = 1 / len(from_atoms)
        site_dict = {atom: w for atom in from_atoms}
        site = LinearSite(site_dict)
        self.vsites.add(index, site)

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        if len(fields[2:]) == 1:
            return [
                [
                    base_atom_index + int(fields[0]) + offset,
                    base_atom_index + int(fields[2]) + offset,
                    0,
                    0,
                    0,
                ]
            ]
        else:
            return []


BondedInteraction_types.append(Virtual_sitesn_COG)


class Virtual_sitesn_COM(VirtualSite):
    """Center of mass virtual site (N-body).
    
    The virtual site position is the mass-weighted center of N atoms.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize COM virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sitesn_COM',
            description="N-body virutal site (COM): atomid, functype, atomid1, atomid2, ..., atomidn",
            category='virtual_sitesn',
            mm_force=None,
            type_label=[1, "2"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one COM virtual site.
        
        Args:
            fields: List of [atomid, functype, atomid1, atomid2, ..., atomidn].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) >= 3 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_sitesn_COM requires at least 3 items and the functype is 2"
        index = int(fields[0]) + base_atom_index + offset
        from_atoms = [int(field) + base_atom_index + offset for field in fields[2:]]

        if len(from_atoms) == 1:
            site_dict = {from_atoms[0]: 1.0}
            site = LinearSite(site_dict)
        else:
            site = COMLinearSite(from_atoms)
            # site = site.to_linear(self.sys, offset)

        self.vsites.add(index, site)

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        if len(fields[2:]) == 1:
            return [
                [
                    base_atom_index + int(fields[0]) + offset,
                    base_atom_index + int(fields[2]) + offset,
                    0,
                    0,
                    0,
                ]
            ]
        else:
            return []


BondedInteraction_types.append(Virtual_sitesn_COM)


class Virtual_sites2(VirtualSite):
    """2-body linear virtual site.
    
    The virtual site lies on the line between two atoms at a
    specified fractional position.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 2-body virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites2',
            description="2-body virutal site: atomid, atomid1, atomid2, functype, weight",
            category='virtual_sites2',
            mm_force=None,
            type_label=[3, "1"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one 2-body virtual site.
        
        Args:
            fields: List of [atomid, atomid1, atomid2, functype, weight].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_sites2 requires 5 items and the functype is 1"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        w = float(fields[4])

        site_dict = {atom1: 1 - w, atom2: w}
        site = LinearSite(site_dict)
        self.vsites.add(index, site)


BondedInteraction_types.append(Virtual_sites2)


class Virtual_sites2fd(VirtualSite):
    """2-body fixed-distance virtual site (fd = from distance).
    
    The virtual site is positioned at a fixed distance from the
    first atom along the direction to the second atom.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 2-body fd virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites2 (fd)',
            description="2-body virutal site: atomid, atomid1, atomid2, functype, distance",
            category='virtual_sites2',
            mm_force=None,
            type_label=[3, "2"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one 2-body fd virtual site.
        
        Args:
            fields: List of [atomid, atomid1, atomid2, functype, distance].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_sites2 (fd) requires 5 items and the functype is 2"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        a = float(fields[4])

        site = NormalizedInPlaneTwoParticleSite(atom1, atom2, a)
        self.vsites.add(index, site)


BondedInteraction_types.append(Virtual_sites2fd)


class Virtual_sites3(VirtualSite):
    """3-body linear virtual site.
    
    The virtual site lies in the plane of three atoms with
    specified weights.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 3-body virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites3',
            description="3-body virutal site: atomid, atomid1, atomid2, atomid3, functype, a, b",
            category='virtual_sites3',
            mm_force=None,
            type_label=[4, "1"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one 3-body virtual site.
        
        Args:
            fields: List of [atomid, atomid1, atomid2, atomid3, functype, a, b].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_site3 requires 7 items and the functype is 1"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        atom3 = int(fields[3]) + base_atom_index + offset
        a = float(fields[5])
        b = float(fields[6])
        w1 = 1 - a - b
        w2 = a
        w3 = b
        site_dict = {
            atom1: w1,
            atom2: w2,
            atom3: w3,
        }
        site = LinearSite(site_dict)
        self.vsites.add(index, site)


BondedInteraction_types.append(Virtual_sites3)


class Virtual_sites3fd(VirtualSite):
    """3-body fixed-distance virtual site in plane (fd = from distance).
    
    The virtual site is positioned in the plane of three atoms
    at a specified position relative to the first two atoms.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 3-body fd virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites3 (fd)',
            description="3-body virutal site: atomid, atomid1, atomid2, atomid3, functype, a, d",
            category='virtual_sites3',
            mm_force=None,
            type_label=[4, "2"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one 3-body fd virtual site.
        
        Args:
            fields: List of [atomid, atomid1, atomid2, atomid3, functype, a, d].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_site3 (fd) requires 7 items and the functype is 2"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        atom3 = int(fields[3]) + base_atom_index + offset
        a = float(fields[5])
        d = float(fields[6])
        site = NormalizedInPlaneSite(atom1, atom2, atom3, a, d)
        self.vsites.add(index, site)


BondedInteraction_types.append(Virtual_sites3fd)


class Virtual_sites3out(VirtualSite):
    """3-body out-of-plane virtual site.
    
    The virtual site is positioned out of the plane defined by
    three atoms using specified displacement parameters.
    """

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 3-body out virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites3 (out)',
            description="3-body virutal site: atomid, atomid1, atomid2, atomid3, functype, a, b, c",
            category='virtual_sites3',
            mm_force=None,
            type_label=[4, "4"],
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one 3-body out virtual site.
        
        Args:
            fields: List of [atomid, atomid1, atomid2, atomid3, functype, a, b, c].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 8 and fields[self.type_label[0]] == self.type_label[1]
        ), "virtual_site3 (out) requires 8 items and the functype is 4"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        atom3 = int(fields[3]) + base_atom_index + offset
        a = float(fields[5])
        b = float(fields[6])
        c = float(fields[7])
        site = OutOfPlane(atom1, atom2, atom3, a, b, c)
        self.vsites.add(index, site)


BondedInteraction_types.append(Virtual_sites3out)


# Exclusions
class Exclusions(Interaction):
    """Exclusion list interaction.
    
    Defines atom pairs that should be excluded from nonbonded
    interactions (typically 1-2, 1-3, and sometimes 1-4 pairs).
    """

    def __init__(self) -> None:
        """Initialize exclusions interaction."""
        super().__init__(
            name='exclusions',
            description='exclusions: atomid1, atomid2, ..., atomidn',
            category='exclusions',
            mm_force=None,
            type_label=None,
        )

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Process exclusion entry.
        
        Args:
            fields: List of [atomid1, atomid2, ..., atomidn].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert len(fields) >= 2, "exclusions requires at least 2 items"

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get exclusion pairs from the fields."""
        exclusions: list[list[float]] = []
        index1 = base_atom_index + int(fields[0]) + offset
        from_atoms = fields[1:]
        for atom in from_atoms:
            index2 = base_atom_index + int(atom) + offset
            exclusions.append([index1, index2, 0, 0, 0])
        return exclusions


BondedInteraction_types.append(Exclusions)


# Multi_allbonds
class Multi_allbonds(Interaction):
    """Multi-state combined bonded interactions.
    
    Supports multiple interaction types (g96 angles, restricted angles,
    periodic dihedrals, contacts) that depend on a state variable.
    """

    def __init__(
        self,
        nonbonded_cutoff: unit.Quantity = 1.1 * mm.unit.nanometer,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize multi_allbonds interaction.
        
        Args:
            nonbonded_cutoff: Cutoff distance for contacts.
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False).
        """
        super().__init__(
            name='multi_allbonds',
            description='Multi_all_bonds: multi_g96_angles, multi_restricted_angles, multi_periodic_dihedrals, multi_contacts',
            category='multi_allbonds',
            mm_force=None,
            type_label=None,
        )
        self.intermolecule_sharing = False

        self.nonbonded_cutoff: unit.Quantity = nonbonded_cutoff
        self.use_sigma_eps: bool = use_sigma_eps

        self.mm_force = mm.CustomCompoundBondForce(
            4,
            "select(delta_g96, g96_angle, 0) + select(delta_rest, restricted_angle, 0) + select(delta_pd, periodic_dihedral, 0) + "
            "select(delta_contact, contact, 0);"
            "g96_angle = 0.5 * k_g96 * (cos(theta) - cos(theta0_g96))^2;"
            "restricted_angle = 0.5 * k_rest * (cos(theta) - cos(theta0_rest))^2 / sin(theta)^2;"
            "theta = angle(p1,p2,p3);"
            "periodic_dihedral = k_pd * ( 1 + cos( n * theta_d - theta0_pd));"
            "theta_d = dihedral(p1,p2,p3,p4);"
            "theta_d = select(0, theta_d_0, 1);"
            "theta_d_0 = dihedral(p1,p2,p3,p4);"
            "contact = step(rcut - r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            "r = distance(p1,p2);"
            f"rcut={self.nonbonded_cutoff.value_in_unit(mm.unit.nanometers)};",
        )
        self.mm_force.addPerBondParameter("delta_g96")
        self.mm_force.addPerBondParameter("delta_rest")
        self.mm_force.addPerBondParameter("delta_pd")
        self.mm_force.addPerBondParameter("delta_contact")
        # g96_angle
        self.mm_force.addPerBondParameter("theta0_g96")
        self.mm_force.addPerBondParameter("k_g96")
        # restricted_angle
        self.mm_force.addPerBondParameter("theta0_rest")
        self.mm_force.addPerBondParameter("k_rest")
        # periodic_dihedral
        self.mm_force.addPerBondParameter("n")
        self.mm_force.addPerBondParameter("theta0_pd")
        self.mm_force.addPerBondParameter("k_pd")
        # contact
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")

    def add_interaction(
        self,
        state: str,
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state interaction.
        
        Args:
            state: Current state identifier.
            category: Sub-category ('multi_angles', 'multi_dihedrals', 'multi_contacts').
            fields: Interaction fields from topology.
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
            
        Raises:
            ValueError: If the fields are not recognized.
        """
        fields_used = False
        if category == 'multi_angles':
            if len(fields) == 8:
                # multi_g96_angles
                if fields[5] == '2' and fields[4] == state:
                    degToRad = math.pi / 180
                    self.mm_force.addBond(
                        [
                            base_atom_index + int(fields[0]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                            base_atom_index + int(fields[2]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                        ],
                        # delta_g96, delta_rest, delta_pd, delta_contact,
                        [1, 0, 0, 0]
                        +
                        # theta0_g96, k_g96,
                        [float(fields[6]) * degToRad, float(fields[7])]
                        +
                        # theta0_rest, k_rest,
                        [0, 0]
                        +
                        # n, theta0_pd, k_pd,
                        [0, 0, 0]
                        +
                        # C12, C6
                        [0, 0],
                    )
                    fields_used = True

                # multi_restricted_angles
                elif fields[5] == '10' and fields[4] == state:
                    degToRad = math.pi / 180
                    self.mm_force.addBond(
                        [
                            base_atom_index + int(fields[0]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                            base_atom_index + int(fields[2]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                        ],
                        # delta_g96, delta_rest, delta_pd, delta_contact,
                        [0, 1, 0, 0]
                        +
                        # theta0_g96, k_g96,
                        [0, 0]
                        +
                        # theta0_rest, k_rest,
                        [float(fields[6]) * degToRad, float(fields[7])]
                        +
                        # n, theta0_pd, k_pd,
                        [0, 0, 0]
                        +
                        # C12, C6
                        [0, 0],
                    )
                    fields_used = True

        elif category == 'multi_dihedrals':
            if len(fields) == 10:
                # multi_periodic_dihedrals
                if fields[6] == '1' and fields[5] == state:
                    degToRad = math.pi / 180
                    self.mm_force.addBond(
                        [
                            base_atom_index + int(fields[0]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                            base_atom_index + int(fields[2]) + offset,
                            base_atom_index + int(fields[3]) + offset,
                        ],
                        # delta_g96, delta_rest, delta_pd, delta_contact,
                        [0, 0, 1, 0]
                        +
                        # theta0_g96, k_g96,
                        [0, 0]
                        +
                        # theta0_rest, k_rest,
                        [0, 0]
                        +
                        # n, theta0_pd, k_pd,
                        [float(fields[9]), float(fields[7]) * degToRad, float(fields[8])]
                        +
                        # C12, C6
                        [0, 0],
                    )
                    fields_used = True

        elif category == 'multi_contacts':
            if len(fields) == 7:
                # multi_contacts
                if fields[4] == '1' and fields[3] == state:
                    if self.use_sigma_eps:
                        sigma = float(fields[5])
                        eps = float(fields[6])
                        C6 = 4 * eps * sigma ** 6
                        C12 = 4 * eps * sigma ** 12
                    else:
                        C6 = float(fields[5])
                        C12 = float(fields[6])

                    self.mm_force.addBond(
                        [
                            base_atom_index + int(fields[0]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                            base_atom_index + int(fields[0]) + offset,
                            base_atom_index + int(fields[1]) + offset,
                        ],
                        # delta_g96, delta_rest, delta_pd, delta_contact,
                        [0, 0, 0, 1]
                        +
                        # theta0_g96, k_g96,
                        [0, 0]
                        +
                        # theta0_rest, k_rest,
                        [0, 0]
                        +
                        # n, theta0_pd, k_pd,
                        [0, 0, 0]
                        +
                        # C12, C6
                        [C12, C6],
                    )
                    fields_used = True

        if not fields_used:
            raise ValueError(f'Error: Unsupport the fields: {category}: {fields}')

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        if category == 'multi_contacts':
            q1 = float(atoms[int(fields[0]) - 1][6])
            q2 = float(atoms[int(fields[1]) - 1][6])
            return [
                [
                    base_atom_index + int(fields[0]) + offset,
                    base_atom_index + int(fields[1]) + offset,
                    q1 * q2,
                    0,
                    0,
                ]
            ]
        else:
            return []


BondedInteraction_types.append(Multi_allbonds)


# Multi_angles
class Multi_harmonic_angles(Interaction):
    """Multi-state harmonic angle interaction.
    
    Harmonic angle potential that depends on a specific state.
    """

    def __init__(self) -> None:
        """Initialize multi harmonic angles interaction."""
        super().__init__(
            name='multi_harmonic_angles',
            description='Multi_harmonic angle potential 1/2*k*(theta-theta0)^2: atomid1, atomid2, atomid3, total number of states, stateid, functype, angle, and k',
            category='multi_angles',
            mm_force=mm.HarmonicAngleForce(),
            type_label=[5, "1"],
        )
        self.intermolecule_sharing = False

    def add_interaction(
        self,
        state: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state harmonic angle.
        
        Args:
            state: Current state identifier.
            fields: List of [atomid1, atomid2, atomid3, nstates, stateid, functype, angle, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 8
            and fields[self.type_label[0]] == self.type_label[1]
            and fields[4] == state
        ), f"multi_harmonic_angles requires 8 items and the functype is 1 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            float(fields[6]) * degToRad,
            float(fields[7]),
        )


BondedInteraction_types.append(Multi_harmonic_angles)


class Multi_g96_angles(Interaction):
    """Multi-state G96 angle interaction.
    
    G96 angle potential that depends on a specific state.
    """

    def __init__(self) -> None:
        """Initialize multi G96 angles interaction."""
        super().__init__(
            name='multi_g96_angles',
            description='multi_g96 angles 0.5*k*(cos(theta)-cos(theta0))^2: atomid1, atomid2, atomid3, total number of states, stateid, functype, angle, and k',
            category='multi_angles',
            mm_force=mm.CustomAngleForce(
                "0.5 * k * (cos(theta) - cos(theta0))^2"
            ),
            type_label=[5, "2"],
        )
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')
        self.intermolecule_sharing = False

    def add_interaction(
        self,
        state: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state G96 angle.
        
        Args:
            state: Current state identifier.
            fields: List of [atomid1, atomid2, atomid3, nstates, stateid, functype, angle, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 8
            and fields[self.type_label[0]] == self.type_label[1]
            and fields[4] == state
        ), f"g96_angles requires 8 items and the functype is 2 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            [float(fields[6]) * degToRad, float(fields[7])],
        )


BondedInteraction_types.append(Multi_g96_angles)


class Multi_restricted_angles(Interaction):
    """Multi-state restricted angle interaction.
    
    Restricted angle potential that depends on a specific state.
    """

    def __init__(self) -> None:
        """Initialize multi restricted angles interaction."""
        super().__init__(
            name='multi_restricted_angles',
            description='Multi_restricted angles 0.5*k*(cos(theta)-cos(theta0))^2 /sin(theta)^2: atomid1, atomid2, atomid3, total number of states, stateid, functype, angle, and k',
            category='multi_angles',
            mm_force=mm.CustomAngleForce(
                "0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2"
            ),
            type_label=[5, "10"],
        )
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')
        self.intermolecule_sharing = False

    def add_interaction(
        self,
        state: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state restricted angle.
        
        Args:
            state: Current state identifier.
            fields: List of [atomid1, atomid2, atomid3, nstates, stateid, functype, angle, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 8
            and fields[self.type_label[0]] == self.type_label[1]
            and fields[4] == state
        ), f"restricted_angles requires 8 items and the functype is 10 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            [float(fields[6]) * degToRad, float(fields[7])],
        )


BondedInteraction_types.append(Multi_restricted_angles)


# Multi_dihedrals
class Multi_periodic_dihedrals(Interaction):
    """Multi-state periodic dihedral interaction.
    
    Periodic dihedral potential that depends on a specific state.
    """

    def __init__(self) -> None:
        """Initialize multi periodic dihedrals interaction."""
        super().__init__(
            name='multi_periodic_dihedrals',
            description='Multi_periodic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, total number of states, stateid, functype, theta, k, and phase',
            category='multi_dihedrals',
            mm_force=mm.PeriodicTorsionForce(),
            type_label=[6, "1"],
        )
        self.intermolecule_sharing = False

    def add_interaction(
        self,
        state: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state periodic dihedral.
        
        Args:
            state: Current state identifier.
            fields: List of [atomid1, atomid2, atomid3, atomid4, nstates, stateid, functype, theta, k, phase].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 10
            and fields[self.type_label[0]] == self.type_label[1]
            and fields[5] == state
        ), f"periodic_dihedrals requires 10 items and the functype is 1 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            float(fields[9]),
            float(fields[7]) * degToRad,
            float(fields[8]),
        )


BondedInteraction_types.append(Multi_periodic_dihedrals)


class Multi_harmonic_dihedrals(Interaction):
    """Multi-state harmonic dihedral interaction.
    
    Harmonic dihedral potential that depends on a specific state.
    """

    def __init__(self) -> None:
        """Initialize multi harmonic dihedrals interaction."""
        super().__init__(
            name='multi_harmonic_dihedrals',
            description='Multi_harmonic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, total number of states, stateid, functype, theta, and k',
            category='multi_dihedrals',
            mm_force=mm.CustomTorsionForce(
                "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                % math.pi
            ),
            type_label=[6, "2"],
        )
        self.mm_force.addPerTorsionParameter("theta0")
        self.mm_force.addPerTorsionParameter("k")
        self.intermolecule_sharing = False

    def add_interaction(
        self,
        state: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state harmonic dihedral.
        
        Args:
            state: Current state identifier.
            fields: List of [atomid1, atomid2, atomid3, atomid4, nstates, stateid, functype, theta, k].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 9
            and fields[self.type_label[0]] == self.type_label[1]
            and fields[5] == state
        ), f"harmonic_dihedrals requires 9 items and the functype is 2 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            [float(fields[7]) * degToRad, float(fields[8])],
        )


BondedInteraction_types.append(Multi_harmonic_dihedrals)


# Multi_contacts
class Multi_contacts(Interaction):
    """Multi-state contact interaction.
    
    Contact (shifted Lennard-Jones) potential that depends on a specific state.
    """

    def __init__(
        self,
        nonbonded_cutoff: unit.Quantity = 1.1 * mm.unit.nanometer,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize multi contacts interaction.
        
        Args:
            nonbonded_cutoff: Cutoff distance for the potential.
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False).
        """
        super().__init__(
            name='multi_contacts',
            description='Multi_contacts Lenard-Jones Potential(r - cutoff): atomid1, atomid2, total number of states, stateid, functype, C6/sigma, C12/epsilon',
            category='multi_contacts',
            mm_force=None,
            type_label=[4, "1"],
        )
        self.nonbonded_cutoff: unit.Quantity = nonbonded_cutoff
        self.use_sigma_eps: bool = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={self.nonbonded_cutoff.value_in_unit(mm.unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")

        self.intermolecule_sharing = False

    def add_interaction(
        self,
        state: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one multi-state contact.
        
        Args:
            state: Current state identifier.
            fields: List of [atomid1, atomid2, nstates, stateid, functype, C6/sigma, C12/epsilon].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        assert (
            len(fields) == 7
            and fields[self.type_label[0]] == self.type_label[1]
            and fields[3] == state
        ), f"contacts requires 7 items and the functype is 1 and state is {state}"
        if self.use_sigma_eps:
            sigma = float(fields[5])
            eps = float(fields[6])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[5])
            C12 = float(fields[6])

        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            [C12, C6],
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        q1 = float(atoms[int(fields[0]) - 1][6])
        q2 = float(atoms[int(fields[1]) - 1][6])
        return [
            [
                base_atom_index + int(fields[0]) + offset,
                base_atom_index + int(fields[1]) + offset,
                q1 * q2,
                0,
                0,
            ]
        ]


BondedInteraction_types.append(Multi_contacts)


# Summary
NonLocal_BondedInteraction_dict: dict[str, list[type[Interaction]]] = {}
Local_BondedInteraction_dict: dict[str, list[type[Interaction]]] = {}
for _Interaction in BondedInteraction_types:
    interaction = _Interaction()
    category = interaction.category
    if interaction.intermolecule_sharing:
        if category not in NonLocal_BondedInteraction_dict:
            NonLocal_BondedInteraction_dict[category] = [_Interaction]
        else:
            NonLocal_BondedInteraction_dict[category].append(_Interaction)
    else:
        if category not in Local_BondedInteraction_dict:
            Local_BondedInteraction_dict[category] = [_Interaction]
        else:
            Local_BondedInteraction_dict[category].append(_Interaction)


# Multiple basin
class EXP_Interaction:
    """Exponential mixing scheme for multiple basin potential.
    
    Implements the exponential mixing formula for combining
    multiple potential energy basins.
    """

    def __init__(self) -> None:
        """Initialize exponential mixing interaction."""
        self.name: str = "exponential mixing scheme"
        self.description: str = 'exponential mixing scheme for multiple baisn popential'

    def addForce(
        self,
        mbp_force_dict: dict[int, list[mm.Force]],
        coupling_constant: float,
        basin_energy_list: list[float],
    ) -> mm.CustomCVForce:
        """Add exponential mixing force.
        
        Args:
            mbp_force_dict: Dictionary mapping state to list of forces.
            coupling_constant: Beta parameter for exponential mixing.
            basin_energy_list: List of basin energy offsets.
            
        Returns:
            The CustomCVForce implementing exponential mixing.
        """
        beta = coupling_constant
        # print("beta", beta)
        basin_energy_list = basin_energy_list

        part1_list = []
        part2_list = []
        for state, force_set in mbp_force_dict.items():
            part1_list.append(f"exp(-beta * (energy{state} + C{state}))")
            energy_combined = ' + '.join([f'state{state}_force{j+1}' for j in range(len(force_set))])
            part2_list.append(f"energy{state} = {energy_combined};")

        part1 = ' + '.join(part1_list)
        part2 = '\n'.join(part2_list)
        energy = f"""-1/beta * log({part1});\n{part2}"""
        print(energy)

        self.mm_force = mm.CustomCVForce(energy)
        for state, force_set in mbp_force_dict.items():
            for j, force in enumerate(force_set):
                self.mm_force.addCollectiveVariable(f"state{state}_force{j+1}", force)
        self.mm_force.addGlobalParameter("beta", float(beta))
        for i, basin_energy in enumerate(basin_energy_list):
            self.mm_force.addGlobalParameter(f'C{i+1}', float(basin_energy))
        return self.mm_force


class HAM_Interaction:
    """Hamiltonian mixing scheme for multiple basin potential.
    
    Implements the Hamiltonian mixing formula for combining
    two potential energy basins.
    """

    def __init__(self) -> None:
        """Initialize Hamiltonian mixing interaction."""
        self.name: str = "hamiltonian mixing scheme"
        self.description: str = 'hamiltonian mixing scheme for multiple baisn popential'

    def addForce(
        self,
        mbp_force_dict: dict[int, list[mm.Force]],
        coupling_constant: float,
        basin_energy_list: list[float],
    ) -> mm.CustomCVForce:
        """Add Hamiltonian mixing force.
        
        Args:
            mbp_force_dict: Dictionary mapping state to list of forces.
            coupling_constant: Delta parameter for Hamiltonian mixing.
            basin_energy_list: List of basin energy offsets.
            
        Returns:
            The CustomCVForce implementing Hamiltonian mixing.
            
        Raises:
            AssertionError: If the number of states is not exactly 2.
        """
        assert len(mbp_force_dict.keys()) == 2
        delta = coupling_constant
        basin_energy_list = basin_energy_list

        part2_list = []
        for state, force_set in mbp_force_dict.items():
            energy_combined = ' + '.join([f'state{state}_force{j+1}' for j in range(len(force_set))])
            part2_list.append(f"energy{state} = {energy_combined};")

        part1 = '(energy1+energy2+deltaV)/2 - sqrt(((energy1-energy2-deltaV)/2)^2+delta^2);deltaV=mbp_energy2-mbp_energy1;'
        part2 = '\n'.join(part2_list)
        energy = f"""{part1};\n{part2}"""
        print(energy)

        self.mm_force = mm.CustomCVForce(energy)
        for state, force_set in mbp_force_dict.items():
            for j, force in enumerate(force_set):
                self.mm_force.addCollectiveVariable(f"state{state}_force{j+1}", force)
        self.mm_force.addGlobalParameter("delta", float(delta))
        for i, basin_energy in enumerate(basin_energy_list):
            self.mm_force.addGlobalParameter(f'mbp_energy{i+1}', float(basin_energy))
        return self.mm_force
