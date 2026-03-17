"""Multi-state interaction classes with lazy unified building.

This module provides multi-state interaction handlers for Multiple Basin Potential (MBP)
simulations. The MultiAllBonds class uses a two-phase approach:

1. Collection phase: add_interaction() records all interactions without building
2. Build phase: mm_force access triggers unified Force construction

This ensures the Force only contains actually used handlers, with consistent
parameter ordering.
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Any, ClassVar

import openmm as mm

from .base import Interaction, InteractionError, register_interaction

if TYPE_CHECKING:
    import openmm.unit as unit

# Module-level constant for degree to radian conversion
DEG_TO_RAD = math.pi / 180


class MultiStateError(InteractionError):
    """Raised when multi-state interaction parsing fails."""
    
    def __init__(self, state: str, category: str, fields: list[str], reason: str = "") -> None:
        msg = f"Multi-state error for state '{state}', category '{category}': {fields}"
        if reason:
            msg = f"{msg} - {reason}"
        super().__init__(msg)


# ============================================================================
# Multi-State Handler Registry
# ============================================================================

# Registry for multi-state interaction handlers
_MULTI_STATE_HANDLERS: dict[str, type[MultiStateInteraction]] = {}


def register_multi_handler(cls: type[H]) -> type[H]:
    """Decorator to register a multi-state interaction handler."""
    _MULTI_STATE_HANDLERS[cls.name] = cls
    return cls


def get_multi_handlers() -> dict[str, type[MultiStateInteraction]]:
    """Get all registered multi-state interaction handlers."""
    return _MULTI_STATE_HANDLERS.copy()


# ============================================================================
# Abstract Base Class for Multi-State Interactions
# ============================================================================

class MultiStateInteraction(ABC):
    """Abstract base class for multi-state interaction handlers.
    
    Each subclass defines a specific interaction type with:
    - Field indices for parsing topology files
    - Energy expression for OpenMM CustomCompoundBondForce
    - Parameter extraction logic
    """
    
    # Class attributes to be defined by subclasses
    name: ClassVar[str]
    category: ClassVar[str]
    functype: ClassVar[str]
    expected_fields: ClassVar[int]
    delta_flag: ClassVar[str]
    needs_sigma_eps: ClassVar[bool] = False
    
    @classmethod
    def can_handle(cls, category: str, fields: list[str]) -> bool:
        """Check if this handler can process the given fields."""
        if category != cls.category:
            return False
        if len(fields) != cls.expected_fields:
            return False
        functype_idx = cls._get_functype_index()
        if functype_idx >= len(fields):
            return False
        return fields[functype_idx] == cls.functype
    
    @classmethod
    def check_state(cls, fields: list[str], state: str) -> bool:
        """Check if the fields match the given state."""
        stateid_idx = cls._get_stateid_index()
        if stateid_idx >= len(fields):
            return False
        return fields[stateid_idx] == state
    
    @classmethod
    def _get_functype_index(cls) -> int:
        """Get field index for function type based on category."""
        if cls.category == 'multi_angles':
            return 5
        elif cls.category == 'multi_dihedrals':
            return 6
        elif cls.category == 'multi_contacts':
            return 4
        return 5
    
    @classmethod
    def _get_stateid_index(cls) -> int:
        """Get field index for state ID based on category."""
        if cls.category == 'multi_angles':
            return 4
        elif cls.category == 'multi_dihedrals':
            return 5
        elif cls.category == 'multi_contacts':
            return 3
        return 4
    
    @property
    @abstractmethod
    def energy_expr(self) -> str:
        """Return the energy expression fragment for this interaction."""
        raise NotImplementedError
    
    @property
    @abstractmethod
    def per_bond_params(self) -> list[str]:
        """Return list of per-bond parameter names."""
        raise NotImplementedError
    
    @classmethod
    @abstractmethod
    def build_atoms(cls, fields: list[str], base_atom_index: int, offset: int) -> list[int]:
        """Build the 4-atom index list for CustomCompoundBondForce."""
        raise NotImplementedError
    
    @classmethod
    @abstractmethod
    def extract_params(cls, fields: list[str]) -> dict[str, float]:
        """Extract parameters from topology fields."""
        raise NotImplementedError


# ============================================================================
# Multi-Angle Interactions
# ============================================================================

@register_multi_handler
class MultiG96Angle(MultiStateInteraction):
    """Multi-state G96 angle potential.
    
    Field layout: [ai, aj, ak, nstates, stateid, functype, theta0, k]
                 functype = 2
    """
    
    name = "g96_angle"
    category = "multi_angles"
    functype = "2"
    expected_fields = 8
    delta_flag = "delta_g96"
    
    class Idx:
        ATOM1, ATOM2, ATOM3, NSTATES, STATEID, FUNCTYPE, THETA0, K = range(8)
    
    @property
    def energy_expr(self) -> str:
        return "0.5 * k_g96 * (cos(angle(p1,p2,p3)) - cos(theta0_g96))^2"
    
    @property
    def per_bond_params(self) -> list[str]:
        return ["delta_g96", "theta0_g96", "k_g96"]
    
    @classmethod
    def build_atoms(cls, fields: list[str], base_atom_index: int, offset: int) -> list[int]:
        idx = cls.Idx
        return [
            base_atom_index + int(fields[idx.ATOM1]) + offset,
            base_atom_index + int(fields[idx.ATOM2]) + offset,
            base_atom_index + int(fields[idx.ATOM3]) + offset,
            base_atom_index + int(fields[idx.ATOM2]) + offset,
        ]
    
    @classmethod
    def extract_params(cls, fields: list[str]) -> dict[str, float]:
        idx = cls.Idx
        return {
            "delta_g96": 1.0,
            "theta0_g96": float(fields[idx.THETA0]) * DEG_TO_RAD,
            "k_g96": float(fields[idx.K]),
        }


@register_multi_handler
class MultiRestrictedAngle(MultiStateInteraction):
    """Multi-state restricted bending angle potential.
    
    Field layout: [ai, aj, ak, nstates, stateid, functype, theta0, k]
                 functype = 10
    """
    
    name = "restricted_angle"
    category = "multi_angles"
    functype = "10"
    expected_fields = 8
    delta_flag = "delta_rest"
    
    class Idx:
        ATOM1, ATOM2, ATOM3, NSTATES, STATEID, FUNCTYPE, THETA0, K = range(8)
    
    @property
    def energy_expr(self) -> str:
        return "0.5 * k_rest * (cos(angle(p1,p2,p3)) - cos(theta0_rest))^2 / sin(angle(p1,p2,p3))^2"
    
    @property
    def per_bond_params(self) -> list[str]:
        return ["delta_rest", "theta0_rest", "k_rest"]
    
    @classmethod
    def build_atoms(cls, fields: list[str], base_atom_index: int, offset: int) -> list[int]:
        idx = cls.Idx
        return [
            base_atom_index + int(fields[idx.ATOM1]) + offset,
            base_atom_index + int(fields[idx.ATOM2]) + offset,
            base_atom_index + int(fields[idx.ATOM3]) + offset,
            base_atom_index + int(fields[idx.ATOM2]) + offset,
        ]
    
    @classmethod
    def extract_params(cls, fields: list[str]) -> dict[str, float]:
        idx = cls.Idx
        return {
            "delta_rest": 1.0,
            "theta0_rest": float(fields[idx.THETA0]) * DEG_TO_RAD,
            "k_rest": float(fields[idx.K]),
        }


# ============================================================================
# Multi-Dihedral Interactions
# ============================================================================

@register_multi_handler
class MultiPeriodicDihedral(MultiStateInteraction):
    """Multi-state periodic dihedral potential.
    
    Field layout: [ai, aj, ak, al, nstates, stateid, functype, phi0, k, n]
                 functype = 1
    """
    
    name = "periodic_dihedral"
    category = "multi_dihedrals"
    functype = "1"
    expected_fields = 10
    delta_flag = "delta_pd"
    
    class Idx:
        ATOM1, ATOM2, ATOM3, ATOM4, NSTATES, STATEID, FUNCTYPE, PHI0, K, N = range(10)
    
    @property
    def energy_expr(self) -> str:
        return "k_pd * (1 + cos(n * dihedral(p1,p2,p3,p4) - theta0_pd))"
    
    @property
    def per_bond_params(self) -> list[str]:
        return ["delta_pd", "theta0_pd", "k_pd", "n"]
    
    @classmethod
    def build_atoms(cls, fields: list[str], base_atom_index: int, offset: int) -> list[int]:
        idx = cls.Idx
        return [
            base_atom_index + int(fields[idx.ATOM1]) + offset,
            base_atom_index + int(fields[idx.ATOM2]) + offset,
            base_atom_index + int(fields[idx.ATOM3]) + offset,
            base_atom_index + int(fields[idx.ATOM4]) + offset,
        ]
    
    @classmethod
    def extract_params(cls, fields: list[str]) -> dict[str, float]:
        idx = cls.Idx
        return {
            "delta_pd": 1.0,
            "theta0_pd": float(fields[idx.PHI0]) * DEG_TO_RAD,
            "k_pd": float(fields[idx.K]),
            "n": int(fields[idx.N]),
        }


# ============================================================================
# Multi-Contact Interactions
# ============================================================================

@register_multi_handler
class MultiContactLJ(MultiStateInteraction):
    """Multi-state Lennard-Jones contact potential with cutoff.
    
    Field layout: [ai, aj, nstates, stateid, functype, C6/sigma, C12/epsilon]
                 functype = 1
    """
    
    name = "contact_lj"
    category = "multi_contacts"
    functype = "1"
    expected_fields = 7
    delta_flag = "delta_contact"
    needs_sigma_eps = True
    
    class Idx:
        ATOM1, ATOM2, NSTATES, STATEID, FUNCTYPE, C6_OR_SIGMA, C12_OR_EPS = range(7)
    
    def __init__(self, nonbonded_cutoff: float = 1.1) -> None:
        self.rcut = nonbonded_cutoff
    
    @property
    def energy_expr(self) -> str:
        rcut = self.rcut
        return (
            f"step({rcut}-distance(p1,p2)) * "
            f"((C12/distance(p1,p2)^12 - C6/distance(p1,p2)^6) - "
            f"(C12/{rcut}^12 - C6/{rcut}^6))"
        )
    
    @property
    def per_bond_params(self) -> list[str]:
        return ["delta_contact", "C12", "C6"]
    
    @classmethod
    def build_atoms(cls, fields: list[str], base_atom_index: int, offset: int) -> list[int]:
        idx = cls.Idx
        a1 = base_atom_index + int(fields[idx.ATOM1]) + offset
        a2 = base_atom_index + int(fields[idx.ATOM2]) + offset
        return [a1, a2, a1, a2]
    
    @classmethod
    def extract_params(cls, fields: list[str], use_sigma_eps: bool = True) -> dict[str, float]:
        idx = cls.Idx
        
        if use_sigma_eps:
            sigma = float(fields[idx.C6_OR_SIGMA])
            eps = float(fields[idx.C12_OR_EPS])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[idx.C6_OR_SIGMA])
            C12 = float(fields[idx.C12_OR_EPS])
        
        return {
            "delta_contact": 1.0,
            "C12": C12,
            "C6": C6,
        }


# ============================================================================
# MultiAllBonds - Lazy Unified Builder
# ============================================================================

@register_interaction
class MultiAllBonds(Interaction):
    """Multi-state combined bonded interactions with lazy unified building.
    
    This class uses a two-phase approach:
    
    1. Collection phase: add_interaction() records all interactions without
       building the Force. Only the handler names and raw data are stored.
    
    2. Build phase: Accessing mm_force triggers unified construction:
       - Instantiate only used handlers (sorted by name for consistency)
       - Build energy expression with only used handlers
       - Batch add all collected bonds
    
    This ensures:
    - Force only contains actually used handlers (efficient)
    - Parameter ordering is consistent (handlers sorted by name)
    - No interactions can be added after Force is built (enforced)
    
    Example:
        >>> multi = MultiAllBonds()
        >>> # Collection phase
        >>> multi.add_interaction("1", "multi_angles", angle_fields, base, offset)
        >>> multi.add_interaction("1", "multi_contacts", contact_fields, base, offset)
        >>> # Build phase (automatic on mm_force access)
        >>> system.addForce(multi.mm_force)  # Triggers unified build
    """
    
    def __init__(
        self,
        nonbonded_cutoff: unit.Quantity | None = None,
        use_sigma_eps: bool = True,
    ) -> None:
        """Initialize MultiAllBonds in collection mode.
        
        Args:
            nonbonded_cutoff: Cutoff distance for contact interactions
            use_sigma_eps: Whether to use sigma/epsilon (True) or C6/C12 (False)
        """
        # Initialize all attributes BEFORE calling super().__init__
        # because super().__init__ sets mm_force which triggers our setter
        
        # Phase 1: Collection
        # _collected: list of (handler_name, state, fields, base_atom_index, offset)
        self._collected: list[tuple[str, str, list[str], int, int]] = []
        self._used_handlers: set[str] = set()  # Handler names encountered
        self._handler_classes: dict[str, type[MultiStateInteraction]] = get_multi_handlers()
        
        # Phase 2: Build (lazy)
        self._finalized = False
        self._handlers: dict[str, MultiStateInteraction] = {}
        self._mm_force: mm.CustomCompoundBondForce | None = None
        self._param_indices: dict[str, int] = {}
        
        # Configuration
        self.nonbonded_cutoff = nonbonded_cutoff or (1.1 * mm.unit.nanometer)
        self.use_sigma_eps = use_sigma_eps
        
        # Now safe to call super().__init__
        super().__init__(
            name='multi_allbonds',
            description='Multi-state combined bonded interactions (lazy build)',
            category='multi_allbonds',
            mm_force=None,
            type_label=None,
        )
        
        # Must set intermolecule_sharing AFTER super().__init__
        self.intermolecule_sharing = False
    
    @property
    def mm_force(self) -> mm.CustomCompoundBondForce:
        """Get the OpenMM CustomCompoundBondForce, triggering build if necessary.
        
        This is the ONLY way to trigger Force construction. The Force is built
        on first access, incorporating all collected interactions.
        
        Returns:
            The configured CustomCompoundBondForce with only used handlers
            
        Raises:
            RuntimeError: If no interactions were collected
        """
        if not self._finalized:
            self._build()
        return self._mm_force
    
    @mm_force.setter
    def mm_force(self, value: mm.CustomCompoundBondForce | None) -> None:
        """Set the mm_force (used by base class during init).
        
        Note: This setter is only valid before finalization. After
        finalization, the force is managed internally.
        """
        if not self._finalized:
            self._mm_force = value
    
    def _find_handler_name(self, category: str, fields: list[str]) -> str | None:
        """Find handler name that can process the given fields."""
        for name, handler_cls in self._handler_classes.items():
            if handler_cls.can_handle(category, fields):
                return name
        return None
    
    def _create_handler(self, handler_name: str) -> MultiStateInteraction:
        """Create handler instance with appropriate initialization."""
        handler_cls = self._handler_classes[handler_name]
        if handler_cls.category == 'multi_contacts':
            return handler_cls(
                nonbonded_cutoff=self.nonbonded_cutoff.value_in_unit(mm.unit.nanometer)
            )
        return handler_cls()
    
    def add_interaction(
        self,
        state: str,
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> bool:
        """Collect one interaction for unified building.
        
        This method records the interaction data. The Force is built later
        when mm_force is accessed.
        
        Args:
            state: State identifier to match
            category: Topology category (multi_angles, etc.)
            fields: Field values from topology file
            base_atom_index: Base index for atom numbering
            offset: Offset to apply to atom indices
            
        Returns:
            True if handler found AND state matches (interaction recorded),
            False if no handler found OR state doesn't match
        """
        # Find handler name
        handler_name = self._find_handler_name(category, fields)
        if handler_name is None:
            return False
        
        # Check state match (using class method to avoid instantiation)
        handler_cls = self._handler_classes[handler_name]
        if not handler_cls.check_state(fields, state):
            return False
        
        # State matches - record this interaction
        self._used_handlers.add(handler_name)
        self._collected.append((handler_name, state, fields, base_atom_index, offset))
        
        return True
    
    def _rebuild(self) -> None:
        """Rebuild Force to include new handlers.
        
        This is called when a new handler type is encountered after
        the Force was already built. All bonds are re-added with
        correct parameter ordering for the new handler set.
        """
        # Reset state
        self._finalized = False
        self._handlers = {}
        self._mm_force = None
        self._param_indices = {}
        
        # Rebuild with all handlers (including new ones)
        # _build() will re-add all bonds from _collected with correct params
        self._build()
    
    def _build(self) -> None:
        """Unified build phase: construct Force with all collected interactions.
        
        This is called automatically on first access to mm_force.
        After this method completes, no more interactions can be added.
        """
        if self._finalized:
            return
        
        # Define fixed handler order for backward compatibility
        # Order: angles -> dihedrals -> contacts (same as original implementation)
        _handler_order = ['g96_angle', 'restricted_angle', 'periodic_dihedral', 'contact_lj']
        
        # Instantiate used handlers in fixed order
        for name in _handler_order:
            if name in self._used_handlers:
                self._handlers[name] = self._create_handler(name)
        
        if not self._handlers:
            # No interactions collected - create empty force
            self._mm_force = mm.CustomCompoundBondForce(4, "0;")
            self._finalized = True
            return
        
        # Build energy expression with used handlers in fixed order
        select_terms = []
        all_params = []
        
        for handler in self._handlers.values():
            expr = handler.energy_expr
            select_terms.append(f"select({handler.delta_flag}, {expr}, 0)")
            all_params.extend(handler.per_bond_params)
        
        energy_expr = " + ".join(select_terms) + ";"
        self._mm_force = mm.CustomCompoundBondForce(4, energy_expr)
        
        # Add per-bond parameters
        self._param_indices = {}
        for param_name in all_params:
            self._param_indices[param_name] = self._mm_force.getNumPerBondParameters()
            self._mm_force.addPerBondParameter(param_name)
        
        # Batch add all collected bonds
        for handler_name, state, fields, base_atom_index, offset in self._collected:
            handler = self._handlers[handler_name]
            atoms = handler.build_atoms(fields, base_atom_index, offset)
            params = self._build_params(handler, fields)
            self._mm_force.addBond(atoms, params)
        
        self._finalized = True
    
    def _build_params(
        self,
        handler: MultiStateInteraction,
        fields: list[str],
    ) -> list[float]:
        """Build parameter list for a bond.
        
        Args:
            handler: The handler for this interaction
            fields: Field values from topology
            
        Returns:
            List of parameter values in correct order for all handlers
        """
        # Extract parameters from handler
        if handler.needs_sigma_eps:
            extracted = handler.extract_params(fields, self.use_sigma_eps)
        else:
            extracted = handler.extract_params(fields)
        
        # Build parameter list for all handlers (fixed order matches _build)
        params = []
        for h in self._handlers.values():
            if h.name == handler.name:
                # This handler is active
                for param_name in h.per_bond_params:
                    if param_name == h.delta_flag:
                        params.append(1.0)
                    elif param_name in extracted:
                        params.append(extracted[param_name])
                    else:
                        params.append(0.0)
            else:
                # Other handlers are inactive
                for _ in h.per_bond_params:
                    params.append(0.0)
        
        return params
    
    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
        """Get exception for multi_contacts."""
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
        return []


# ============================================================================
# Exports
# ============================================================================

__all__ = [
    # Constants and errors
    "DEG_TO_RAD",
    "MultiStateError",
    # Core classes
    "MultiStateInteraction",
    "MultiAllBonds",
    # Registry functions
    "register_multi_handler",
    "get_multi_handlers",
    # Angle handlers
    "MultiG96Angle",
    "MultiRestrictedAngle",
    # Dihedral handlers
    "MultiPeriodicDihedral",
    # Contact handlers
    "MultiContactLJ",
]
