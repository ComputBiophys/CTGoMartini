"""Registry for bonded interaction types."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .base import Interaction


def build_interaction_dicts() -> tuple[
    dict[str, list[type[Interaction]]],
    dict[str, list[type[Interaction]]]
]:
    """Build interaction dictionaries from registered interactions.
    
    Separates interactions into non-local (shared between molecules) and
    local (molecule-specific) dictionaries based on their
    `intermolecule_sharing` attribute.
    
    Only includes classes that inherit from Interaction and have the
    required attributes (category, intermolecule_sharing).
    
    This function is called on each access to ensure correct behavior
    in multiprocessing environments.
    
    Returns:
        Tuple of (non_local_dict, local_dict) where each dictionary
        maps category names to lists of interaction classes.
    """
    # Import here to avoid circular imports and ensure all handlers are registered
    from .base import Interaction, get_registered_interactions
    
    non_local: dict[str, list[type[Interaction]]] = {}
    local: dict[str, list[type[Interaction]]] = {}
    
    for interaction_class in get_registered_interactions():
        # Skip non-Interaction classes (e.g., EXPInteraction, HAMInteraction)
        if not issubclass(interaction_class, Interaction):
            continue
        
        # Instantiate to get category and sharing info
        try:
            interaction = interaction_class()
        except TypeError:
            # Some interactions require arguments (e.g., Constraints, Pairs)
            continue
        
        # Skip if missing required attributes
        if not hasattr(interaction, 'category') or not hasattr(interaction, 'intermolecule_sharing'):
            continue
            
        category = interaction.category
        target_dict = non_local if interaction.intermolecule_sharing else local
        
        if category not in target_dict:
            target_dict[category] = []
        target_dict[category].append(interaction_class)
    
    return non_local, local


# Module-level cache
_non_local_cache: dict[str, list[type[Interaction]]] | None = None
_local_cache: dict[str, list[type[Interaction]]] | None = None


def _get_non_local_dict() -> dict[str, list[type[Interaction]]]:
    """Get non-local interaction dictionary."""
    global _non_local_cache
    if _non_local_cache is None:
        _non_local_cache, _ = build_interaction_dicts()
    return _non_local_cache


def _get_local_dict() -> dict[str, list[type[Interaction]]]:
    """Get local interaction dictionary."""
    global _local_cache
    if _local_cache is None:
        _, _local_cache = build_interaction_dicts()
    return _local_cache


# Public API - using a class with property to allow lazy evaluation
# but still appear as module-level variables
class _RegistryDict:
    """Wrapper class that provides dict-like access with lazy loading."""
    
    def __init__(self, getter):
        self._getter = getter
    
    def __getitem__(self, key):
        return self._getter()[key]
    
    def __contains__(self, key):
        return key in self._getter()
    
    def __iter__(self):
        return iter(self._getter())
    
    def keys(self):
        return self._getter().keys()
    
    def values(self):
        return self._getter().values()
    
    def items(self):
        return self._getter().items()
    
    def get(self, key, default=None):
        return self._getter().get(key, default)
    
    def __repr__(self):
        return repr(self._getter())
    
    def __len__(self):
        return len(self._getter())


# Create wrapper instances
NonLocalBondedInteractionDict = _RegistryDict(_get_non_local_dict)
LocalBondedInteractionDict = _RegistryDict(_get_local_dict)


def rebuild_interaction_dicts() -> tuple[
    dict[str, list[type[Interaction]]],
    dict[str, list[type[Interaction]]]
]:
    """Rebuild interaction dictionaries.
    
    Useful for forcing re-initialization, e.g., after registering new handlers.
    
    Returns:
        Tuple of (non_local_dict, local_dict).
    """
    global _non_local_cache, _local_cache
    _non_local_cache = None
    _local_cache = None
    return build_interaction_dicts()


__all__ = [
    "NonLocalBondedInteractionDict",
    "LocalBondedInteractionDict",
    "build_interaction_dicts",
    "rebuild_interaction_dicts",
]
