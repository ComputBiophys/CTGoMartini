"""Bond-related interaction classes."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import openmm as mm

from .base import Interaction, register_interaction

if TYPE_CHECKING:
    import openmm.unit as unit


@register_interaction
class HarmonicBonds(Interaction):
    """Harmonic bond potential interaction.
    
    Implements a harmonic bond potential of the form:
    V(r) = 0.5 * k * (r - r0)^2
    
    where r is the bond length, r0 is the equilibrium length,
    and k is the force constant.
    
    Field layout: [atomid1, atomid2, functype, length, k]
    """
    
    class Idx:
        """Field indices for harmonic bonds."""
        ATOM1 = 0
        ATOM2 = 1
        FUNCTYPE = 2
        LENGTH = 3
        K = 4
    
    EXPECTED_FIELDS = 5

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        self.mm_force.addBond(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            float(fields[idx.LENGTH]),
            float(fields[idx.K]),
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        idx = self.Idx
        return [
            [
                self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
                self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
                0,
                0,
                0,
            ]
        ]


@register_interaction
class Constraints(Interaction):
    """Distance constraint interaction.
    
    Adds a rigid constraint between two atoms at a fixed distance.
    This is implemented as an OpenMM constraint rather than a force.
    
    Field layout: [atomid1, atomid2, functype, length] (optional 5th field ignored)
    """
    
    class Idx:
        """Field indices for constraints."""
        ATOM1 = 0
        ATOM2 = 1
        FUNCTYPE = 2
        LENGTH = 3
    
    MIN_FIELDS = 4

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
        self.sys.addConstraint(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            float(fields[idx.LENGTH]),
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        idx = self.Idx
        return [
            [
                self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
                self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
                0,
                0,
                0,
            ]
        ]


@register_interaction
class Pairs(Interaction):
    """Pair interaction for 1-4 or special pairs.
    
    Implements Lennard-Jones and electrostatic interactions between
    atom pairs that are not directly bonded but need special treatment.
    
    Field layout: [atomid1, atomid2, functype, C6/sigma, C12/epsilon]
    """
    
    class Idx:
        """Field indices for pairs."""
        ATOM1 = 0
        ATOM2 = 1
        FUNCTYPE = 2
        C6_OR_SIGMA = 3
        C12_OR_EPS = 4
    
    EXPECTED_FIELDS = 5

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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        # Actual interaction addition happens in get_exception

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Add pair interaction and return empty exception list."""
        idx = self.Idx
        if self.use_sigma_eps:
            sigma = float(fields[idx.C6_OR_SIGMA])
            eps = float(fields[idx.C12_OR_EPS])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[idx.C6_OR_SIGMA])
            C12 = float(fields[idx.C12_OR_EPS])
        q1 = float(atoms[int(fields[idx.ATOM1]) - 1][6])
        q2 = float(atoms[int(fields[idx.ATOM2]) - 1][6])
        self.mm_force.addBond(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            [C12, C6, q1, q2],
        )
        return []


@register_interaction
class ContactsLJ(Interaction):
    """Contact (soft-core) interaction.
    
    Implements a shifted Lennard-Jones potential for contacts:
    V(r) = step(rcut - r) * (LJ(r) - LJ(rcut))
    
    This smoothly goes to zero at the cutoff distance.
    rcut is fixed at 1.1 nm.
    
    Field layout: [atomid1, atomid2, functype, C6/sigma, C12/epsilon]
    """
    
    class Idx:
        """Field indices for contacts."""
        ATOM1 = 0
        ATOM2 = 1
        FUNCTYPE = 2
        C6_OR_SIGMA = 3
        C12_OR_EPS = 4
    
    EXPECTED_FIELDS = 5
    RCUT_NM: float = 1.1  # Fixed cutoff distance in nanometers

    def __init__(
        self,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize contacts interaction.
        
        Args:
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False).
        """
        super().__init__(
            name='contacts',
            description='Contacts Lenard-Jones Potential(r - cutoff): atomid1, atomid2, functype, C6/sigma, C12/epsilon',
            category='contacts',
            mm_force=None,
            type_label=[2, "1"],
        )
        self.use_sigma_eps: bool = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={self.RCUT_NM};"
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
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        if self.use_sigma_eps:
            sigma = float(fields[idx.C6_OR_SIGMA])
            eps = float(fields[idx.C12_OR_EPS])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[idx.C6_OR_SIGMA])
            C12 = float(fields[idx.C12_OR_EPS])

        self.mm_force.addBond(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
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
        idx = self.Idx
        q1 = float(atoms[int(fields[idx.ATOM1]) - 1][6])
        q2 = float(atoms[int(fields[idx.ATOM2]) - 1][6])
        return [
            [
                self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
                self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
                q1 * q2,
                0,
                0,
            ]
        ]


@register_interaction
class Contacts6_12(Interaction):
    """Contact (6-12) interaction with adjustable rcut.
    
    Implements a shifted Lennard-Jones potential for contacts:
    V(r) = step(rcut - r) * (LJ(r) - LJ(rcut))
    
    rcut is adjustable per bond (stored as per-bond parameter).
    
    Field layout: [atomid1, atomid2, functype=2, C6/sigma, C12/epsilon, rcut]
    """
    
    class Idx:
        """Field indices for adjustable 6-12 contacts."""
        ATOM1 = 0
        ATOM2 = 1
        FUNCTYPE = 2
        C6_OR_SIGMA = 3
        C12_OR_EPS = 4
        RCUT = 5
    
    EXPECTED_FIELDS = 6
    FUNCTYPE = "2"

    def __init__(
        self,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize adjustable 6-12 contacts interaction.
        
        Args:
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False).
        """
        super().__init__(
            name='contacts_6_12_adjustable',
            description='Contacts 6-12 LJ with adjustable rcut: atomid1, atomid2, functype=2, C6/sigma, C12/epsilon, rcut',
            category='contacts',
            mm_force=None,
            type_label=[2, "2"],
        )
        self.use_sigma_eps: bool = use_sigma_eps

        # rcut is per-bond parameter, allowing each contact to have different rcut
        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")
        self.mm_force.addPerBondParameter("rcut")

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one adjustable 6-12 contact interaction.
        
        Args:
            fields: List of [atomid1, atomid2, functype, C6/sigma, C12/epsilon, rcut].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        if self.use_sigma_eps:
            sigma = float(fields[idx.C6_OR_SIGMA])
            eps = float(fields[idx.C12_OR_EPS])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[idx.C6_OR_SIGMA])
            C12 = float(fields[idx.C12_OR_EPS])
        
        rcut = float(fields[idx.RCUT])

        self.mm_force.addBond(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            [C12, C6, rcut],
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        idx = self.Idx
        q1 = float(atoms[int(fields[idx.ATOM1]) - 1][6])
        q2 = float(atoms[int(fields[idx.ATOM2]) - 1][6])
        return [
            [
                self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
                self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
                q1 * q2,
                0,
                0,
            ]
        ]


@register_interaction
class Contacts10_12(Interaction):
    """Contact (10-12) interaction with adjustable rcut.
    
    Implements a shifted 10-12 Lennard-Jones potential:
    V(r) = step(rcut - r) * (C12/r^12 - C10/r^10 - C12/rcut^12 + C10/rcut^10)
    
    This is useful for specific interactions requiring different distance dependence.
    rcut is adjustable per bond (stored as per-bond parameter).
    
    Field layout: [atomid1, atomid2, functype=3, C10/sigma10, C12/epsilon, rcut]
                 or with sigma/epsilon: [atomid1, atomid2, 3, sigma, epsilon, rcut]
    """
    
    class Idx:
        """Field indices for 10-12 contacts."""
        ATOM1 = 0
        ATOM2 = 1
        FUNCTYPE = 2
        C10_OR_SIGMA = 3  # C10 or sigma (sigma^10 -> C10)
        C12_OR_EPS = 4    # C12 or epsilon
        RCUT = 5
    
    EXPECTED_FIELDS = 6
    FUNCTYPE = "3"

    def __init__(
        self,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize 10-12 contacts interaction.
        
        Args:
            use_sigma_eps: Whether to use sigma/epsilon (True) or C10/C12 (False).
        """
        super().__init__(
            name='contacts_10_12',
            description='Contacts 10-12 LJ with adjustable rcut: atomid1, atomid2, functype=3, C10/sigma, C12/epsilon, rcut',
            category='contacts',
            mm_force=None,
            type_label=[2, "3"],
        )
        self.use_sigma_eps: bool = use_sigma_eps

        # rcut is per-bond parameter
        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C10/r^10);"
            "corr = (C12/rcut^12 - C10/rcut^10);"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C10")
        self.mm_force.addPerBondParameter("rcut")

    def add_interaction(
        self,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> None:
        """Add one 10-12 contact interaction.
        
        Args:
            fields: List of [atomid1, atomid2, functype, C10/sigma, C12/epsilon, rcut].
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
        """
        self._validate_field_count(fields, self.EXPECTED_FIELDS)
        self._validate_functype(fields)
        
        idx = self.Idx
        if self.use_sigma_eps:
            sigma = float(fields[idx.C10_OR_SIGMA])
            eps = float(fields[idx.C12_OR_EPS])
            # 10-12 LJ: C10 = 4*eps*sigma^10, C12 = 4*eps*sigma^12
            C10 = 4 * eps * sigma ** 10
            C12 = 4 * eps * sigma ** 12
        else:
            C10 = float(fields[idx.C10_OR_SIGMA])
            C12 = float(fields[idx.C12_OR_EPS])
        
        rcut = float(fields[idx.RCUT])

        self.mm_force.addBond(
            self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
            self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
            [C12, C10, rcut],
        )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get the exception from the fields."""
        idx = self.Idx
        q1 = float(atoms[int(fields[idx.ATOM1]) - 1][6])
        q2 = float(atoms[int(fields[idx.ATOM2]) - 1][6])
        return [
            [
                self._get_atom_index(fields[idx.ATOM1], base_atom_index, offset),
                self._get_atom_index(fields[idx.ATOM2], base_atom_index, offset),
                q1 * q2,
                0,
                0,
            ]
        ]


@register_interaction
class Exclusions(Interaction):
    """Exclusion list interaction.
    
    Defines atom pairs that should be excluded from nonbonded
    interactions (typically 1-2, 1-3, and sometimes 1-4 pairs).
    
    Field layout: [atomid1, atomid2, ..., atomidn]
    """
    
    class Idx:
        """Field indices for exclusions."""
        ATOM1 = 0
        # Additional atoms follow
    
    MIN_FIELDS = 2

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
        if len(fields) < self.MIN_FIELDS:
            from .base import FieldValidationError
            raise FieldValidationError(
                expected=self.MIN_FIELDS,
                actual=len(fields),
                fields=fields,
                context=self.name
            )

    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get exclusion pairs from the fields."""
        exclusions: list[list[float]] = []
        index1 = self._get_atom_index(fields[self.Idx.ATOM1], base_atom_index, offset)
        from_atoms = fields[1:]
        for atom in from_atoms:
            index2 = self._get_atom_index(atom, base_atom_index, offset)
            exclusions.append([index1, index2, 0, 0, 0])
        return exclusions


__all__ = [
    "HarmonicBonds",
    "Constraints",
    "Pairs",
    "ContactsLJ",
    "Contacts6_12",
    "Contacts10_12",
    "Exclusions",
]
