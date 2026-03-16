"""Dihedral-related interaction classes."""

from __future__ import annotations

import math

import openmm as mm

from .base import Interaction, register_interaction

# Module-level constant for degree to radian conversion
DEG_TO_RAD = math.pi / 180


@register_interaction
class PeriodicDihedrals(Interaction):
    """Periodic dihedral (torsion) interaction.
    
    Implements a periodic dihedral potential of the form:
    V(phi) = k * (1 + cos(n * phi - phi0))
    
    where phi is the dihedral angle, phi0 is the phase angle,
    n is the periodicity, and k is the force constant.
    
    Field layout: [atomid1, atomid2, atomid3, atomid4, functype, theta, k, phase]
    """
    
    class Idx:
        """Field indices for periodic dihedrals."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        ATOM4 = 3
        FUNCTYPE = 4
        THETA = 5
        K = 6
        PHASE = 7
    
    EXPECTED_FIELDS = 8

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        self.mm_force.addTorsion(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM4], base_atom_index, offset),
            float(fields[idx.PHASE]),
            float(fields[idx.THETA]) * DEG_TO_RAD,
            float(fields[idx.K]),
        )


@register_interaction
class HarmonicDihedrals(Interaction):
    """Harmonic dihedral (improper) interaction.
    
    Implements a harmonic dihedral potential for improper dihedrals
    with proper periodic boundary handling.
    
    Field layout: [atomid1, atomid2, atomid3, atomid4, functype, theta, k]
    """
    
    class Idx:
        """Field indices for harmonic dihedrals."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        ATOM4 = 3
        FUNCTYPE = 4
        THETA = 5
        K = 6
    
    EXPECTED_FIELDS = 7

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        self.mm_force.addTorsion(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM4], base_atom_index, offset),
            [float(fields[idx.THETA]) * DEG_TO_RAD, float(fields[idx.K])],
        )


@register_interaction
class CombinedBendingTorsionPotentials(Interaction):
    """Combined bending-torsion potential interaction.
    
    Implements a potential coupling bond bending and torsional motions:
    V = k * sin(theta0)^3 * sin(theta1)^3 * sum(ai * cos(phi)^i)
    
    This is used for more accurate conformational modeling.
    
    Field layout: [atomid1, atomid2, atomid3, atomid4, functype, k, a0, a1, a2, a3, a4]
    """
    
    class Idx:
        """Field indices for combined bending-torsion potentials."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        ATOM4 = 3
        FUNCTYPE = 4
        K = 5
        A0 = 6
        A1 = 7
        A2 = 8
        A3 = 9
        A4 = 10
    
    EXPECTED_FIELDS = 11

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        k = float(fields[idx.K])
        a0 = float(fields[idx.A0])
        a1 = float(fields[idx.A1])
        a2 = float(fields[idx.A2])
        a3 = float(fields[idx.A3])
        a4 = float(fields[idx.A4])
        
        self.mm_force.addBond(
            self._get_atom_indices(fields, base_atom_index, offset, 4),
            [k, a0, a1, a2, a3, a4],
        )


@register_interaction
class RyckaertBellemansDihedrals(Interaction):
    """Ryckaert-Bellemans dihedral interaction.
    
    Implements the Ryckaert-Bellemans potential for aliphatic chains:
    V(phi) = sum(Ci * cos(phi)^i) for i = 0 to 5
    
    This is commonly used for alkane chains.
    
    Field layout: [atomid1, atomid2, atomid3, atomid4, functype, C0, C1, C2, C3, C4, C5]
    """
    
    class Idx:
        """Field indices for Ryckaert-Bellemans dihedrals."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        ATOM4 = 3
        FUNCTYPE = 4
        C0 = 5
        C1 = 6
        C2 = 7
        C3 = 8
        C4 = 9
        C5 = 10
    
    EXPECTED_FIELDS = 11

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        c = [float(fields[i]) for i in range(idx.C0, idx.C5 + 1)]
        
        if fields[idx.FUNCTYPE] == "5":
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
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM4], base_atom_index, offset),
            c[0], c[1], c[2], c[3], c[4], c[5],
        )


__all__ = [
    "DEG_TO_RAD",
    "PeriodicDihedrals",
    "HarmonicDihedrals",
    "CombinedBendingTorsionPotentials",
    "RyckaertBellemansDihedrals",
]
