"""Base class for bonded interactions."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, TypeVar

import openmm as mm


class InteractionError(ValueError):
    """Raised when interaction parsing or validation fails."""
    pass


class FieldValidationError(InteractionError):
    """Raised when field validation fails."""
    
    def __init__(self, expected: int, actual: int, fields: list[str], context: str = "") -> None:
        msg = f"Expected {expected} fields, got {actual}"
        if context:
            msg = f"{context}: {msg}"
        msg = f"{msg}: {fields}"
        super().__init__(msg)


class FuncTypeError(InteractionError):
    """Raised when function type validation fails."""
    
    def __init__(self, expected: str | list[str], actual: str, context: str = "") -> None:
        expected_str = expected if isinstance(expected, str) else f"one of {expected}"
        msg = f"Expected functype {expected_str}, got {actual}"
        if context:
            msg = f"{context}: {msg}"
        super().__init__(msg)


# Registry for interaction types
_INTERACTION_REGISTRY: list[type[Interaction]] = []


I = TypeVar("I", bound="Interaction")


def register_interaction(cls: type[I]) -> type[I]:
    """Decorator to register interaction classes.
    
    Usage:
        @register_interaction
        class HarmonicBonds(Interaction):
            ...
    """
    _INTERACTION_REGISTRY.append(cls)
    return cls


def get_registered_interactions() -> list[type[Interaction]]:
    """Get all registered interaction classes."""
    return _INTERACTION_REGISTRY.copy()


class Interaction(ABC):
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

        # Default parameters
        self.contents: list[Any] = []
        self.intermolecule_sharing: bool = True

    def __str__(self) -> str:
        """Return string representation of the interaction."""
        return f"{self.name}: {self.description}"

    def _validate_field_count(self, fields: list[str], expected: int) -> None:
        """Validate the number of fields.
        
        Args:
            fields: List of field values.
            expected: Expected number of fields.
            
        Raises:
            FieldValidationError: If field count doesn't match.
        """
        if len(fields) != expected:
            raise FieldValidationError(
                expected=expected,
                actual=len(fields),
                fields=fields,
                context=self.name
            )

    def _validate_functype(self, fields: list[str]) -> None:
        """Validate the function type field.
        
        Args:
            fields: List of field values.
            
        Raises:
            FuncTypeError: If functype doesn't match expected values.
        """
        if self.type_label is None:
            return
        
        field_idx = self.type_label[0]
        expected = self.type_label[1:]
        actual = fields[field_idx]
        
        if actual not in expected:
            raise FuncTypeError(
                expected=list(expected),
                actual=actual,
                context=self.name
            )

    def _get_atom_index(self, field: str, base_atom_index: int, offset: int) -> int:
        """Calculate atom index with offset.
        
        Args:
            field: Field value containing atom index.
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
            
        Returns:
            Calculated atom index.
        """
        return base_atom_index + int(field) + offset

    def _get_atom_indices(
        self, 
        fields: list[str], 
        base_atom_index: int, 
        offset: int,
        count: int
    ) -> list[int]:
        """Calculate multiple atom indices.
        
        Args:
            fields: List of field values.
            base_atom_index: Base index for atom numbering.
            offset: Offset to apply to atom indices.
            count: Number of atom indices to extract (from start of fields).
            
        Returns:
            List of calculated atom indices.
        """
        return [
            self._get_atom_index(fields[i], base_atom_index, offset)
            for i in range(count)
        ]

    @abstractmethod
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


# For backward compatibility
BondedInteractionTypes = _INTERACTION_REGISTRY

__all__ = [
    "Interaction",
    "InteractionError",
    "FieldValidationError",
    "FuncTypeError",
    "register_interaction",
    "get_registered_interactions",
    "BondedInteractionTypes",
]
