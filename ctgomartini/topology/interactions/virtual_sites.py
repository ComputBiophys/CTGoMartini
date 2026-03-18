"""Virtual site interaction classes."""

from __future__ import annotations

from typing import Any

import openmm as mm

from .vsites import (
    COMLinearSite,
    LinearSite,
    NormalizedInPlaneSite,
    NormalizedInPlaneTwoParticleSite,
    OutOfPlane,
    VSiteManager,
)
from .base import Interaction, register_interaction


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


@register_interaction
class VirtualSitesNCOG(VirtualSite):
    """Center of geometry virtual site (N-body).
    
    The virtual site position is the geometric center of N atoms.
    
    Field layout: [atomid, functype, atomid1, atomid2, ..., atomidn]
    """
    
    class Idx:
        """Field indices for N-body COG virtual sites."""
        ATOM_ID = 0
        FUNCTYPE = 1
        FIRST_SOURCE_ATOM = 2  # Additional source atoms follow
    
    MIN_FIELDS = 3

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize COG virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sitesn_COG',
            description="N-body virtual site (COG): atomid, functype, atomid1, atomid2, ..., atomidn",
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
        if len(fields) < self.MIN_FIELDS:
            from .base import FieldValidationError
            raise FieldValidationError(
                expected=self.MIN_FIELDS,
                actual=len(fields),
                fields=fields,
                context=self.name
            )
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        from_atoms = [
            self._get_atom_index(field, base_atom_index, offset) 
            for field in fields[idx.FIRST_SOURCE_ATOM:]
        ]

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
        idx = self.Idx
        source_atoms = fields[idx.FIRST_SOURCE_ATOM:]
        if len(source_atoms) == 1:
            return [
                [
                    self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset),
                    self._get_atom_index(source_atoms[0], base_atom_index, offset),
                    0,
                    0,
                    0,
                ]
            ]
        else:
            return []


@register_interaction
class VirtualSitesNCOM(VirtualSite):
    """Center of mass virtual site (N-body).
    
    The virtual site position is the mass-weighted center of N atoms.
    
    Field layout: [atomid, functype, atomid1, atomid2, ..., atomidn]
    """
    
    class Idx:
        """Field indices for N-body COM virtual sites."""
        ATOM_ID = 0
        FUNCTYPE = 1
        FIRST_SOURCE_ATOM = 2
    
    MIN_FIELDS = 3

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize COM virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sitesn_COM',
            description="N-body virtual site (COM): atomid, functype, atomid1, atomid2, ..., atomidn",
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
        if len(fields) < self.MIN_FIELDS:
            from .base import FieldValidationError
            raise FieldValidationError(
                expected=self.MIN_FIELDS,
                actual=len(fields),
                fields=fields,
                context=self.name
            )
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        from_atoms = [
            self._get_atom_index(field, base_atom_index, offset)
            for field in fields[idx.FIRST_SOURCE_ATOM:]
        ]

        if len(from_atoms) == 1:
            site_dict = {from_atoms[0]: 1.0}
            site = LinearSite(site_dict)
        else:
            site = COMLinearSite(from_atoms)

        self.vsites.add(index, site)

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        idx = self.Idx
        source_atoms = fields[idx.FIRST_SOURCE_ATOM:]
        if len(source_atoms) == 1:
            return [
                [
                    self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset),
                    self._get_atom_index(source_atoms[0], base_atom_index, offset),
                    0,
                    0,
                    0,
                ]
            ]
        else:
            return []


@register_interaction
class VirtualSites2(VirtualSite):
    """2-body linear virtual site.
    
    The virtual site lies on the line between two atoms at a
    specified fractional position.
    
    Field layout: [atomid, atomid1, atomid2, functype, weight]
    """
    
    class Idx:
        """Field indices for 2-body virtual sites."""
        ATOM_ID = 0
        ATOM1 = 1
        ATOM2 = 2
        FUNCTYPE = 3
        WEIGHT = 4
    
    EXPECTED_FIELDS = 5

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 2-body virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites2',
            description="2-body virtual site: atomid, atomid1, atomid2, functype, weight",
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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        atom1 = self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset)
        atom2 = self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset)
        w = float(fields[idx.WEIGHT])

        site_dict = {atom1: 1 - w, atom2: w}
        site = LinearSite(site_dict)
        self.vsites.add(index, site)


@register_interaction
class VirtualSites2fd(VirtualSite):
    """2-body fixed-distance virtual site (fd = from distance).
    
    The virtual site is positioned at a fixed distance from the
    first atom along the direction to the second atom.
    
    Field layout: [atomid, atomid1, atomid2, functype, distance]
    """
    
    class Idx:
        """Field indices for 2-body fd virtual sites."""
        ATOM_ID = 0
        ATOM1 = 1
        ATOM2 = 2
        FUNCTYPE = 3
        DISTANCE = 4
    
    EXPECTED_FIELDS = 5

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 2-body fd virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites2 (fd)',
            description="2-body virtual site: atomid, atomid1, atomid2, functype, distance",
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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        atom1 = self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset)
        atom2 = self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset)
        a = float(fields[idx.DISTANCE])

        site = NormalizedInPlaneTwoParticleSite(atom1, atom2, a)
        self.vsites.add(index, site)


@register_interaction
class VirtualSites3(VirtualSite):
    """3-body linear virtual site.
    
    The virtual site lies in the plane of three atoms with
    specified weights.
    
    Field layout: [atomid, atomid1, atomid2, atomid3, functype, a, b]
    """
    
    class Idx:
        """Field indices for 3-body virtual sites."""
        ATOM_ID = 0
        ATOM1 = 1
        ATOM2 = 2
        ATOM3 = 3
        FUNCTYPE = 4
        A = 5
        B = 6
    
    EXPECTED_FIELDS = 7

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 3-body virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites3',
            description="3-body virtual site: atomid, atomid1, atomid2, atomid3, functype, a, b",
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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        atom1 = self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset)
        atom2 = self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset)
        atom3 = self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset)
        a = float(fields[idx.A])
        b = float(fields[idx.B])
        
        w1 = 1 - a - b
        w2 = a
        w3 = b
        site_dict = {atom1: w1, atom2: w2, atom3: w3}
        site = LinearSite(site_dict)
        self.vsites.add(index, site)


@register_interaction
class VirtualSites3fd(VirtualSite):
    """3-body fixed-distance virtual site in plane (fd = from distance).
    
    The virtual site is positioned in the plane of three atoms
    at a specified position relative to the first two atoms.
    
    Field layout: [atomid, atomid1, atomid2, atomid3, functype, a, d]
    """
    
    class Idx:
        """Field indices for 3-body fd virtual sites."""
        ATOM_ID = 0
        ATOM1 = 1
        ATOM2 = 2
        ATOM3 = 3
        FUNCTYPE = 4
        A = 5
        D = 6
    
    EXPECTED_FIELDS = 7

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 3-body fd virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites3 (fd)',
            description="3-body virtual site: atomid, atomid1, atomid2, atomid3, functype, a, d",
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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        atom1 = self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset)
        atom2 = self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset)
        atom3 = self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset)
        a = float(fields[idx.A])
        d = float(fields[idx.D])
        
        site = NormalizedInPlaneSite(atom1, atom2, atom3, a, d)
        self.vsites.add(index, site)


@register_interaction
class VirtualSites3out(VirtualSite):
    """3-body out-of-plane virtual site.
    
    The virtual site is positioned out of the plane defined by
    three atoms using specified displacement parameters.
    
    Field layout: [atomid, atomid1, atomid2, atomid3, functype, a, b, c]
    """
    
    class Idx:
        """Field indices for 3-body out-of-plane virtual sites."""
        ATOM_ID = 0
        ATOM1 = 1
        ATOM2 = 2
        ATOM3 = 3
        FUNCTYPE = 4
        A = 5
        B = 6
        C = 7
    
    EXPECTED_FIELDS = 8

    def __init__(self, vsites: VSiteManager | None = None) -> None:
        """Initialize 3-body out virtual site interaction."""
        super().__init__(
            vsites=vsites,
            name='virtual_sites3 (out)',
            description="3-body virtual site: atomid, atomid1, atomid2, atomid3, functype, a, b, c",
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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        index = self._get_atom_index(fields[idx.ATOM_ID], base_atom_index, offset)
        atom1 = self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset)
        atom2 = self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset)
        atom3 = self._get_atom_index(fields[idx.ATOM3], base_atom_index, offset)
        a = float(fields[idx.A])
        b = float(fields[idx.B])
        c = float(fields[idx.C])
        
        site = OutOfPlane(atom1, atom2, atom3, a, b, c)
        self.vsites.add(index, site)


__all__ = [
    "VirtualSite",
    "VirtualSitesNCOG",
    "VirtualSitesNCOM",
    "VirtualSites2",
    "VirtualSites2fd",
    "VirtualSites3",
    "VirtualSites3fd",
    "VirtualSites3out",
]
