"""Multi-state interaction classes with extensible handler architecture."""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
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


@dataclass
class HandlerParam:
    """Definition of a handler parameter."""
    name: str
    field_index: int
    converter: str = "float"  # "float", "deg_to_rad", "int"
    default: float = 0.0


class MultiStateHandler(ABC):
    """Abstract base class for multi-state interaction handlers."""
    
    # Class attributes to be defined by subclasses
    name: ClassVar[str]
    category: ClassVar[str]
    delta_flag: ClassVar[str]
    
    expected_fields: ClassVar[int]
    functypes: ClassVar[list[str]]
    atoms_template: ClassVar[list[int]]
    
    params: ClassVar[list[HandlerParam]]
    
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
        return fields[functype_idx] in cls.functypes
    
    @classmethod
    def _get_functype_index(cls) -> int:
        """Get the field index for functype based on category."""
        if cls.category == 'multi_angles':
            return 5
        elif cls.category == 'multi_dihedrals':
            return 6
        elif cls.category == 'multi_contacts':
            return 4
        return 5
    
    @classmethod
    def build_atoms(cls, fields: list[str], base_atom_index: int, offset: int) -> list[int]:
        """Build the 4-atom index list for CustomCompoundBondForce."""
        atoms = []
        for idx in cls.atoms_template:
            atoms.append(base_atom_index + int(fields[idx]) + offset)
        return atoms
    
    @classmethod
    def get_params(cls, fields: list[str]) -> list[float]:
        """Extract and convert parameter values from fields."""
        values = []
        for param in cls.params:
            raw_value = fields[param.field_index]
            if param.converter == "float":
                values.append(float(raw_value))
            elif param.converter == "deg_to_rad":
                values.append(float(raw_value) * DEG_TO_RAD)
            elif param.converter == "int":
                values.append(int(raw_value))
            else:
                values.append(float(raw_value))
        return values
    
    @classmethod
    def check_state(cls, fields: list[str], state: str) -> bool:
        """Check if the fields match the given state."""
        stateid_idx = cls._get_stateid_index()
        if stateid_idx >= len(fields):
            return False
        return fields[stateid_idx] == state
    
    @classmethod
    def _get_stateid_index(cls) -> int:
        """Get the field index for stateid based on category."""
        if cls.category == 'multi_angles':
            return 4
        elif cls.category == 'multi_dihedrals':
            return 5
        elif cls.category == 'multi_contacts':
            return 3
        return 4


# ==================== Built-in Handlers ====================

class G96AngleHandler(MultiStateHandler):
    """Handler for G96 angle interactions."""
    name = "g96_angle"
    category = "multi_angles"
    delta_flag = "delta_g96"
    expected_fields = 8
    functypes = ["2"]
    atoms_template = [0, 1, 2, 1]
    
    params = [
        HandlerParam("theta0_g96", 6, "deg_to_rad"),
        HandlerParam("k_g96", 7, "float"),
    ]


class RestrictedAngleHandler(MultiStateHandler):
    """Handler for restricted angle interactions."""
    name = "restricted_angle"
    category = "multi_angles"
    delta_flag = "delta_rest"
    expected_fields = 8
    functypes = ["10"]
    atoms_template = [0, 1, 2, 1]
    
    params = [
        HandlerParam("theta0_rest", 6, "deg_to_rad"),
        HandlerParam("k_rest", 7, "float"),
    ]


class PeriodicDihedralHandler(MultiStateHandler):
    """Handler for periodic dihedral interactions."""
    name = "periodic_dihedral"
    category = "multi_dihedrals"
    delta_flag = "delta_pd"
    expected_fields = 10
    functypes = ["1"]
    atoms_template = [0, 1, 2, 3]
    
    params = [
        HandlerParam("theta0_pd", 7, "deg_to_rad"),
        HandlerParam("k_pd", 8, "float"),
        HandlerParam("n", 9, "int"),
    ]


class ContactHandler(MultiStateHandler):
    """Handler for contact (shifted LJ) interactions."""
    name = "contact"
    category = "multi_contacts"
    delta_flag = "delta_contact"
    expected_fields = 7
    functypes = ["1"]
    atoms_template = [0, 1, 0, 1]
    
    # Field indices: 5 = C6_or_sigma, 6 = C12_or_eps
    params = [
        HandlerParam("C12", 6, "float"),
        HandlerParam("C6", 5, "float"),
    ]
    
    @classmethod
    def get_params(cls, fields: list[str], use_sigma_eps: bool = True) -> list[float]:
        """Extract and convert contact parameters from fields.
        
        Args:
            fields: The field values from topology
            use_sigma_eps: If True, interpret fields[5] as sigma and fields[6] as epsilon
                          and convert to C6/C12 using LJ formulas.
                          If False, treat fields[5] as C6 and fields[6] as C12 directly.
        """
        if use_sigma_eps:
            sigma = float(fields[5])
            eps = float(fields[6])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[5])
            C12 = float(fields[6])
        return [C12, C6]


# ==================== MultiAllBonds Class ====================

@register_interaction
class MultiAllBonds(Interaction):
    """Multi-state combined bonded interactions with extensible handlers."""
    
    _handler_registry: ClassVar[dict[str, type[MultiStateHandler]]] = {}
    _builtin_handlers: ClassVar[list[type[MultiStateHandler]]] = [
        G96AngleHandler,
        RestrictedAngleHandler,
        PeriodicDihedralHandler,
        ContactHandler,
    ]
    _registry_initialized: ClassVar[bool] = False
    
    @classmethod
    def _ensure_registry_initialized(cls) -> None:
        if not cls._registry_initialized or not cls._handler_registry:
            cls._registry_initialized = False
            cls._handler_registry.clear()
            for handler in cls._builtin_handlers:
                cls.register_handler(handler)
            cls._registry_initialized = True
    
    @classmethod
    def register_handler(cls, handler: type[MultiStateHandler]) -> None:
        if handler.name in cls._handler_registry:
            raise ValueError(f"Handler '{handler.name}' is already registered")
        cls._handler_registry[handler.name] = handler
    
    @classmethod
    def get_registered_handlers(cls) -> dict[str, type[MultiStateHandler]]:
        cls._ensure_registry_initialized()
        if not cls._registry_initialized or not cls._handler_registry:
            cls._registry_initialized = False
            cls._handler_registry.clear()
            for handler in cls._builtin_handlers:
                cls.register_handler(handler)
            cls._registry_initialized = True
        return cls._handler_registry.copy()
    
    def __init__(
        self,
        enabled_handlers: list[str] = None,
        nonbonded_cutoff: unit.Quantity = None,
        use_sigma_eps: bool = True,
    ) -> None:
        self.__class__._ensure_registry_initialized()
        
        super().__init__(
            name='multi_allbonds',
            description='Multi-state combined bonded interactions with extensible handlers',
            category='multi_allbonds',
            mm_force=None,
            type_label=None,
        )
        self.intermolecule_sharing = False
        
        if enabled_handlers is None:
            self._enabled_handlers = list(self._handler_registry.keys())
        else:
            invalid = set(enabled_handlers) - set(self._handler_registry.keys())
            if invalid:
                raise ValueError(f"Unknown handlers: {invalid}")
            self._enabled_handlers = enabled_handlers
        
        self.nonbonded_cutoff = nonbonded_cutoff or (1.1 * mm.unit.nanometer)
        self.use_sigma_eps = use_sigma_eps
        
        self._force_built = False
        self._collected_interactions: list[tuple[str, str, list[str]]] = []
        
        self._build_force()
    
    def _build_force(self) -> None:
        if self._force_built:
            return
        
        active_handlers = [
            self._handler_registry[name]
            for name in self._enabled_handlers
            if name in self._handler_registry
        ]
        
        if not active_handlers:
            raise ValueError("No handlers enabled for MultiAllBonds")
        
        # Build energy expression using inline format (no sub-expressions)
        # Format: select(delta_flag, expression, 0) for each handler
        # Variables like theta, theta_d are computed inline
        
        rcut_nm = self.nonbonded_cutoff.value_in_unit(mm.unit.nanometers)
        
        select_terms = []
        per_bond_params = []
        
        # Build parameter list and select terms
        for handler in active_handlers:
            if handler.name == "g96_angle":
                # theta is computed inline
                expr = f"0.5 * k_g96 * (cos(angle(p1,p2,p3)) - cos(theta0_g96))^2"
                select_terms.append(f"select({handler.delta_flag}, {expr}, 0)")
                per_bond_params.extend([handler.delta_flag, "theta0_g96", "k_g96"])
                
            elif handler.name == "restricted_angle":
                expr = f"0.5 * k_rest * (cos(angle(p1,p2,p3)) - cos(theta0_rest))^2 / sin(angle(p1,p2,p3))^2"
                select_terms.append(f"select({handler.delta_flag}, {expr}, 0)")
                per_bond_params.extend([handler.delta_flag, "theta0_rest", "k_rest"])
                
            elif handler.name == "periodic_dihedral":
                expr = f"k_pd * (1 + cos(n * dihedral(p1,p2,p3,p4) - theta0_pd))"
                select_terms.append(f"select({handler.delta_flag}, {expr}, 0)")
                per_bond_params.extend([handler.delta_flag, "theta0_pd", "k_pd", "n"])
                
            elif handler.name == "contact":
                # Contact with inline expressions
                expr = (f"step({rcut_nm}-distance(p1,p2)) * ((C12/distance(p1,p2)^12 - C6/distance(p1,p2)^6) - "
                        f"(C12/{rcut_nm}^12 - C6/{rcut_nm}^6))")
                select_terms.append(f"select({handler.delta_flag}, {expr}, 0)")
                per_bond_params.extend([handler.delta_flag, "C12", "C6"])
        
        energy_expr = " + ".join(select_terms) + ";"
        
        self.mm_force = mm.CustomCompoundBondForce(4, energy_expr)
        
        # Add per-bond parameters
        self._param_indices = {}
        for param_name in per_bond_params:
            self._param_indices[param_name] = self.mm_force.getNumPerBondParameters()
            self.mm_force.addPerBondParameter(param_name)
        
        self._active_handlers = active_handlers
        self._force_built = True
    
    def add_interaction(
        self,
        state: str,
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> bool:
        if not self._force_built:
            self._build_force()
        
        return self._do_add_interaction(state, category, fields, base_atom_index, offset)
    
    def _can_any_handler_process(self, category: str, fields: list[str]) -> bool:
        """Check if any handler can process these fields (regardless of state)."""
        for handler in self._active_handlers:
            if handler.can_handle(category, fields):
                return True
        return False
    
    def _do_add_interaction(
        self,
        state: str,
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> bool:
        # First check if any handler can process this category/type
        matching_handler = None
        for handler in self._active_handlers:
            if handler.can_handle(category, fields):
                matching_handler = handler
                break
        
        if matching_handler is None:
            # No handler can process this category/type
            return False
        
        # Check if state matches
        if not matching_handler.check_state(fields, state):
            # This field belongs to a different state, skip it
            return True
        
        # State matches, add the bond
        atoms = matching_handler.build_atoms(fields, base_atom_index, offset)
        
        # Build params list matching the order in _build_force
        params = []
        
        # Special handling for ContactHandler which uses sigma/epsilon conversion
        # Note: matching_handler is a class, not an instance, so use issubclass
        if issubclass(matching_handler, ContactHandler):
            # Get C12 and C6 values using ContactHandler's get_params method
            contact_params = ContactHandler.get_params(fields, use_sigma_eps=self.use_sigma_eps)
            # contact_params = [C12, C6]
            for param_name in self._param_indices.keys():
                if param_name == matching_handler.delta_flag:
                    params.append(1.0)  # Active
                elif param_name.startswith("delta_"):
                    params.append(0.0)  # Inactive
                elif param_name == "C12":
                    params.append(contact_params[0])  # C12 from get_params
                elif param_name == "C6":
                    params.append(contact_params[1])  # C6 from get_params
                else:
                    params.append(0.0)  # Default for other handlers' params
        else:
            # Standard handling for other handlers
            for param_name in self._param_indices.keys():
                # Find which handler this param belongs to
                if param_name == matching_handler.delta_flag:
                    params.append(1.0)  # Active
                elif param_name.startswith("delta_"):
                    params.append(0.0)  # Inactive
                else:
                    # Get param value from fields
                    for hp in matching_handler.params:
                        if hp.name == param_name:
                            val = fields[hp.field_index]
                            if hp.converter == "float":
                                params.append(float(val))
                            elif hp.converter == "deg_to_rad":
                                params.append(float(val) * DEG_TO_RAD)
                            elif hp.converter == "int":
                                params.append(int(val))
                            else:
                                params.append(float(val))
                            break
                    else:
                        params.append(0.0)  # Default for other handlers' params
        
        self.mm_force.addBond(atoms, params)
        return True
    
    def get_exception(
        self,
        atoms: list[tuple[Any, ...]],
        category: str,
        fields: list[str],
        base_atom_index: int = 0,
        offset: int = -1,
    ) -> list[list[float]]:
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



__all__ = [
    "DEG_TO_RAD",
    "HandlerParam",
    "MultiStateHandler",
    "MultiStateError",
    "MultiAllBonds",
    "G96AngleHandler",
    "RestrictedAngleHandler",
    "PeriodicDihedralHandler",
    "ContactHandler",
]
