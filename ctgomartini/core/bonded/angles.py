"""Angle-related interaction classes."""

from __future__ import annotations

import math

import openmm as mm

from .base import Interaction, register_interaction

# Module-level constant for degree to radian conversion
DEG_TO_RAD = math.pi / 180


@register_interaction
class HarmonicAngles(Interaction):
    """Harmonic angle potential interaction.
    
    Implements a harmonic angle potential of the form:
    V(theta) = 0.5 * k * (theta - theta0)^2
    
    where theta is the bond angle, theta0 is the equilibrium angle,
    and k is the force constant.
    
    Field layout: [atomid1, atomid2, atomid3, functype, angle, k]
    """
    
    class Idx:
        """Field indices for harmonic angles."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        FUNCTYPE = 3
        ANGLE = 4
        K = 5
    
    EXPECTED_FIELDS = 6

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        self.mm_force.addAngle(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset),
            float(fields[idx.ANGLE]) * DEG_TO_RAD,
            float(fields[idx.K]),
        )


@register_interaction
class G96Angles(Interaction):
    """G96 angle potential interaction.
    
    Implements the GROMOS 96 angle potential of the form:
    V(theta) = 0.5 * k * (cos(theta) - cos(theta0))^2
    
    This is a more numerically stable form for certain angle ranges.
    
    Field layout: [atomid1, atomid2, atomid3, functype, angle, k]
    """
    
    class Idx:
        """Field indices for G96 angles."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        FUNCTYPE = 3
        ANGLE = 4
        K = 5
    
    EXPECTED_FIELDS = 6

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        self.mm_force.addAngle(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset),
            [float(fields[idx.ANGLE]) * DEG_TO_RAD, float(fields[idx.K])],
        )


@register_interaction
class RestrictedAngles(Interaction):
    """Restricted angle potential interaction.
    
    Implements a restricted bending potential of the form:
    V(theta) = 0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2
    
    This prevents the angle from reaching 0 or 180 degrees.
    
    Field layout: [atomid1, atomid2, atomid3, functype, angle, k]
    """
    
    class Idx:
        """Field indices for restricted angles."""
        ATOM1 = 0
        ATOM2 = 1
        ATOM3 = 2
        FUNCTYPE = 3
        ANGLE = 4
        K = 5
    
    EXPECTED_FIELDS = 6

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        self.mm_force.addAngle(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset),
            [float(fields[idx.ANGLE]) * DEG_TO_RAD, float(fields[idx.K])],
        )


__all__ = [
    "DEG_TO_RAD",
    "HarmonicAngles",
    "G96Angles",
    "RestrictedAngles",
]
