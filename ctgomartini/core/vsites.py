"""
Virtual site management module for CTGoMartini.

Provides classes for defining and managing virtual sites in molecular
topologies, including linear sites (weighted averages of atoms) and
non-linear sites (out-of-plane, in-plane positions).
"""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

import openmm.unit as unit

if TYPE_CHECKING:
    import openmm as mm


class LinearSite:
    """
    Linear virtual site defined by weighted average of reference atoms.
    
    Attributes:
        atom_weights: Dictionary mapping atom indices to their weights.
    """
    
    def __init__(self, atom_weight_dict: dict[int, float]) -> None:
        """
        Initialize linear virtual site.
        
        Args:
            atom_weight_dict: Dictionary of {atom_index: weight} pairs.
        """
        self.atom_weights: dict[int, float] = atom_weight_dict

    def __repr__(self) -> str:
        return f"LinearSite({self.atom_weights})"


class COMLinearSite:
    """
    Center-of-mass linear virtual site.
    
    Initially defined by atom list, converts to LinearSite with
    mass-based weights when system is available.
    
    Attributes:
        atoms: List of atom indices comprising the COM site.
    """
    
    def __init__(self, atoms: list[int]) -> None:
        """
        Initialize COM virtual site.
        
        Args:
            atoms: List of atom indices.
        """
        self.atoms: list[int] = atoms

    def to_linear(self, system: mm.System, offset: int) -> LinearSite:
        """
        Convert to LinearSite using actual atomic masses.
        
        Args:
            system: OpenMM system containing particle masses.
            offset: Atom index offset for this molecule.
            
        Returns:
            LinearSite with mass-weighted coefficients.
        """
        masses = [
            system.getParticleMass(i + offset).value_in_unit(unit.dalton)
            for i in self.atoms
        ]
        total = sum(masses)
        weights = [m / total for m in masses]
        site_dict = {atom: weight for atom, weight in zip(self.atoms, weights)}
        return LinearSite(site_dict)

    def __repr__(self) -> str:
        return f"COMLinearSite({self.atoms})"


class NonLinearSite:
    """Base class for non-linear virtual sites."""
    pass


class OutOfPlane(NonLinearSite):
    """
    Out-of-plane virtual site.
    
    Position defined relative to three reference atoms with
    out-of-plane displacement parameters.
    
    Attributes:
        atom1, atom2, atom3: Reference atom indices.
        a, b, c: Position parameters relative to reference atoms.
    """
    
    def __init__(
        self,
        atom1: int,
        atom2: int,
        atom3: int,
        a: float,
        b: float,
        c: float
    ) -> None:
        """
        Initialize out-of-plane virtual site.
        
        Args:
            atom1, atom2, atom3: Reference atom indices.
            a, b, c: Position parameters.
        """
        self.atom1: int = atom1
        self.atom2: int = atom2
        self.atom3: int = atom3
        self.a: float = a
        self.b: float = b
        self.c: float = c

    def __repr__(self) -> str:
        return (
            f"OutOfPlane({self.atom1}, {self.atom2}, {self.atom3}, "
            f"{self.a}, {self.b}, {self.c})"
        )


class NormalizedInPlaneSite(NonLinearSite):
    """
    Normalized in-plane virtual site with three reference atoms.
    
    Attributes:
        atom1, atom2, atom3: Reference atom indices.
        a, d: In-plane position parameters.
    """
    
    def __init__(
        self,
        atom1: int,
        atom2: int,
        atom3: int,
        a: float,
        d: float
    ) -> None:
        """
        Initialize normalized in-plane virtual site.
        
        Args:
            atom1, atom2, atom3: Reference atom indices.
            a, d: Position parameters.
        """
        self.atom1: int = atom1
        self.atom2: int = atom2
        self.atom3: int = atom3
        self.a: float = a
        self.d: float = d

    def __repr__(self) -> str:
        return (
            f"NormalizedInPlaneSite({self.atom1}, {self.atom2}, {self.atom3}, "
            f"{self.a}, {self.d})"
        )


class NormalizedInPlaneTwoParticleSite(NonLinearSite):
    """
    Normalized in-plane virtual site with two reference atoms.
    
    Attributes:
        atom1, atom2: Reference atom indices.
        a: In-plane position parameter.
    """
    
    def __init__(self, atom1: int, atom2: int, a: float) -> None:
        """
        Initialize two-particle in-plane virtual site.
        
        Args:
            atom1, atom2: Reference atom indices.
            a: Position parameter.
        """
        self.atom1: int = atom1
        self.atom2: int = atom2
        self.a: float = a

    def __repr__(self) -> str:
        return (
            f"NormalizedInPlaneTwoParticleSite({self.atom1}, {self.atom2}, {self.a})"
        )


class VSiteManager:
    """
    Manager for virtual sites in a molecule.
    
    Handles adding virtual sites, converting COM sites to linear,
    and flattening nested virtual site definitions.
    
    Attributes:
        vsites: Dictionary mapping virtual site index to site object.
    """
    
    def __init__(self) -> None:
        """Initialize empty virtual site manager."""
        self.vsites: dict[int, LinearSite | NonLinearSite] = {}

    def add(self, index: int, site: LinearSite | NonLinearSite) -> None:
        """
        Add a virtual site.
        
        Args:
            index: Virtual site particle index.
            site: Virtual site definition.
            
        Raises:
            ValueError: If a site already exists at this index.
        """
        if index in self.vsites:
            raise ValueError(f"Tried to add more than one vsite for particle {index}.")
        self.vsites[index] = site

    def convert_com_to_linear(self, system: mm.System, offset: int) -> None:
        """
        Convert all COMLinearSite instances to LinearSite.
        
        Args:
            system: OpenMM system for mass lookup.
            offset: Atom index offset.
        """
        for index, site in self.vsites.items():
            if isinstance(site, COMLinearSite):
                self.vsites[index] = site.to_linear(system, offset)

    def iter(self):
        """
        Iterate over virtual sites.
        
        Yields:
            Tuples of (index, site) where site is the flattened definition.
        """
        for index, site in self.vsites.items():
            if isinstance(site, NonLinearSite):
                yield index, site
            else:
                flattened = self.flatten_site(site)
                yield index, flattened

    def flatten_site(self, site: LinearSite) -> LinearSite:
        """
        Flatten a linear site by resolving nested virtual sites.
        
        Args:
            site: Linear site that may reference other virtual sites.
            
        Returns:
            Flattened LinearSite with only real atom references.
            
        Raises:
            RuntimeError: If site depends on non-linear or COM sites.
        """
        if isinstance(site, NonLinearSite):
            raise RuntimeError(
                "A linear vsite cannot depend on a non-linear vsite."
            )

        if isinstance(site, COMLinearSite):
            raise RuntimeError(
                "All COM sites should have been converted to LinearSites."
            )

        from_atoms: defaultdict[int, list[float]] = defaultdict(list)

        for atom, weight in site.atom_weights.items():
            # If atom is itself a vsite, recursively flatten
            if atom in self.vsites:
                flattened = self.flatten_site(self.vsites[atom])
                for f_atom, f_weight in flattened.atom_weights.items():
                    from_atoms[f_atom].append(weight * f_weight)
            else:
                # Regular atom
                from_atoms[atom].append(weight)

        # Sum coefficients for each atom
        result = {atom: sum(weights) for atom, weights in from_atoms.items()}
        return LinearSite(result)
