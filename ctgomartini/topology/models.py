"""Topology models for CTGoMartini.

This module provides Molecule and ForceField classes for representing
molecular topology data and managing force field parameters.
"""

from __future__ import annotations

from typing import Any

from .interactions import BondedInteractionTypes


class Molecule:
    """
    Represents a molecule with its topology information.
    
    Attributes:
        name: Name of the molecule.
        _topology: Dictionary mapping category names to their contents.
        atoms: List of atom definitions (populated after initialize).
    """
    
    def __init__(self, name: str | None = None) -> None:
        """
        Initialize a molecule.
        
        Args:
            name: Name of the molecule.
        """
        self.name: str | None = name
        self._topology: dict[str, list[list[str]]] = {}

    def readLine(self, line: str, currentCategory: str) -> None:
        """
        Parse and classify a line into the appropriate category.
        
        Args:
            line: The line to parse.
            currentCategory: The current topology category being parsed.
        """
        line_fields = line.strip().split()
        if currentCategory not in self._topology:
            self._topology[currentCategory] = []
        self._topology[currentCategory].append(line_fields)

    def initialize(self, ff_parameters: dict[str, Any]) -> None:
        """
        Initialize molecule categories with force field parameters.
        
        Args:
            ff_parameters: Dictionary of force field parameters.
        """
        for category_name, contents in self._topology.items():
            if category_name in Molecule_Categories:
                category = Molecule_Categories[category_name](contents)
            elif category_name in Bonded_Categories:
                category = Bonded_Categories[category_name]
                category.contents = contents
            else:
                raise ValueError(f"Unknown category {category_name}")

            if category.parameter_missing:
                category.complement(ff_parameters)
                
            self.__dict__[category_name] = category.contents
            self.__dict__[f'_{category_name}'] = category


class Molecule_Category:
    """
    Base class for molecule topology categories.
    
    Attributes:
        name: Name of the category.
        description: Description of the category.
        category: Category identifier.
        contents: List of data entries in this category.
        parameter_missing: Whether parameters need to be complemented.
    """
    
    def __init__(
        self,
        name: str,
        description: str,
        category: str,
        contents: list[list[str]]
    ) -> None:
        """
        Initialize a molecule category.
        
        Args:
            name: Name of the category.
            description: Description of the category.
            category: Category identifier.
            contents: Initial data for this category.
        """
        self.name: str = name
        self.description: str = description
        self.category: str = category
        self.contents: list[list[str]] = contents
        self.parameter_missing: bool = False

    def __str__(self) -> str:
        """Return string representation."""
        return f"{self.name}: {self.description}"

    def initialize(self) -> None:
        """Initialize the category. Override in subclasses."""
        pass

    def complement(self, ff_parameters: dict[str, Any]) -> None:
        """
        Complement missing parameters from force field.
        
        Args:
            ff_parameters: Dictionary of force field parameters.
            
        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError(
            f"Do not implement complement for category ({self.name})"
        )


# Registry of molecule categories
Molecule_Categories: dict[str, type[Molecule_Category]] = {}


class Moleculetype(Molecule_Category):
    """Molecule type category defining molecule name and exclusion rules."""
    
    def __init__(self, contents: list[list[str]] | None = None) -> None:
        """
        Initialize moleculetype category.
        
        Args:
            contents: List of moleculetype definitions.
        """
        if contents is None:
            contents = []
        super().__init__(
            name="moleculetype",
            description=(
                "Defines the molecule name and nrexcl "
                "(number of bonded atoms to exclude from non-bonded interactions)"
            ),
            category="moleculetype",
            contents=contents
        )


Molecule_Categories['moleculetype'] = Moleculetype


class Atoms(Molecule_Category):
    """
    Atom definitions category.
    
    Each atom entry contains:
    atom number, atom type, residue number, residue name,
    atom name, charge group number, charge, mass.
    """
    
    def __init__(self, contents: list[list[str]] | None = None) -> None:
        """
        Initialize atoms category.
        
        Args:
            contents: List of atom definitions.
        """
        if contents is None:
            contents = []
        super().__init__(
            name="atoms",
            description=(
                "atom number; atom type; residue number; residue name; "
                "atom name; charge group number; charge (e); mass (u)"
            ),
            category="atoms",
            contents=contents
        )
        self.parameter_missing = True

    def complement(self, ff_parameters: dict[str, Any]) -> None:
        """
        Complement missing charge and mass from atomtypes.
        
        Args:
            ff_parameters: Dictionary containing 'atomtypes' parameter.
        """
        atomtypes = ff_parameters['atomtypes']
        dict_atomtypes: dict[str, list[str]] = {
            fields[0]: fields for fields in atomtypes
        }
        
        for fields in self.contents:
            if len(fields) < 5:
                raise ValueError(f"Too few fields in [ atoms ] line: {fields}")
            elif len(fields) > 8:
                raise ValueError(f"Too many fields in [ atoms ] line: {fields}")
            else:
                fields += [None] * (8 - len(fields))

            # Fill missing charge
            if fields[6] is None:
                fields[6] = dict_atomtypes[fields[1]][4]
            # Fill missing mass
            if fields[7] is None:
                fields[7] = dict_atomtypes[fields[1]][3]
            if None in fields:
                raise ValueError(f'None in the fields {fields}, {self.category}')


Molecule_Categories['atoms'] = Atoms


class Multiple_basin(Molecule_Category):
    """Multiple basin potential parameters category."""
    
    def __init__(self, contents: list[list[str]] | None = None) -> None:
        """
        Initialize multiple_basin category.
        
        Args:
            contents: List of multiple basin definitions.
            Format: [bool, method, n_states, coupling_constant, energy1, energy2, ...]
        """
        if contents is None:
            contents = []
        super().__init__(
            name="multiple_basin",
            description="bool, method, n_states, coupling_constant, energy1, energy2, ...",
            category="multiple_basin",
            contents=contents
        )


Molecule_Categories['multiple_basin'] = Multiple_basin


# Bonded interaction categories registry
Bonded_Categories: dict[str, Molecule_Category] = {}

# Initialize bonded categories from BondedInteraction_types
for _Interaction in BondedInteractionTypes:
    interaction = _Interaction()
    category_name = interaction.category

    if category_name not in Bonded_Categories:
        molecule_category = Molecule_Category(
            name=category_name,
            description=f"{category_name}: automated generation",
            category=category_name,
            contents=[]
        )
        Bonded_Categories[category_name] = molecule_category

# Register multi-state categories for MBP (handled by MultiAllBonds)
# These categories replace the legacy individual multi-state interaction classes
for _category_name in ['multi_angles', 'multi_dihedrals', 'multi_contacts']:
    if _category_name not in Bonded_Categories:
        Bonded_Categories[_category_name] = Molecule_Category(
            name=_category_name,
            description=f"{_category_name}: automated generation",
            category=_category_name,
            contents=[]
        )


class ForceField:
    """
    Force field parameter container.
    
    Attributes:
        name: Name of the force field.
        _parameters: Dictionary mapping category names to raw contents.
        atomtypes: Dictionary of atom type definitions.
        defaults: List of default parameter settings.
        nonbond_params: Dictionary of non-bonded interaction parameters.
    """
    
    def __init__(self, name: str | None = None) -> None:
        """
        Initialize force field.
        
        Args:
            name: Name of the force field.
        """
        self.name: str | None = name
        self._parameters: dict[str, list[list[str]]] = {}

    def readLine(self, line: str, currentCategory: str) -> None:
        """
        Parse a line into the appropriate force field category.
        
        Args:
            line: The line to parse.
            currentCategory: The current parameter category.
        """
        line_fields = line.strip().split()
        if currentCategory not in self._parameters:
            self._parameters[currentCategory] = []
        self._parameters[currentCategory].append(line_fields)

    def initialize(self) -> None:
        """Initialize all force field categories from parsed parameters."""
        for category_name, contents in self._parameters.items():
            if category_name in FF_Categories:
                category = FF_Categories[category_name](contents)
            else:
                try:
                    category = FF_Category(
                        name=category_name,
                        description="Automated generation",
                        category=category_name,
                        contents=contents
                    )
                except Exception as exc:
                    raise ValueError(
                        f"[ {category_name} ] is not currently supported."
                    ) from exc

            category.initialize()
            self.__dict__[category_name] = category.dict_contents
            self.__dict__[f'_{category_name}'] = category


class FF_Category:
    """
    Base class for force field parameter categories.
    
    Attributes:
        name: Name of the category.
        description: Description of the category.
        category: Category identifier.
        contents: List of raw data entries.
        dict_contents: Processed dictionary of entries.
    """
    
    def __init__(
        self,
        name: str,
        description: str,
        category: str,
        contents: list[list[str]]
    ) -> None:
        """
        Initialize a force field category.
        
        Args:
            name: Name of the category.
            description: Description of the category.
            category: Category identifier.
            contents: List of parameter entries.
        """
        self.name: str = name
        self.description: str = description
        self.category: str = category
        self.contents: list[list[str]] = contents
        self.dict_contents: dict[str, Any] | list[list[str]] = {}

    def __str__(self) -> str:
        """Return string representation."""
        return f"{self.name}: {self.description}"

    def initialize(self) -> None:
        """Initialize the category. Override in subclasses."""
        pass


# Registry of force field categories
FF_Categories: dict[str, type[FF_Category]] = {}


class Defaults(FF_Category):
    """
    Default force field parameters category.
    
    Contains GROMACS defaults (nbfunc, comb-rule, gen-pairs, fudgeLJ, fudgeQQ).
    """
    
    def __init__(self, contents: list[list[str]] | None = None) -> None:
        """
        Initialize defaults category.
        
        Args:
            contents: List of default parameter definitions.
        """
        if contents is None:
            contents = []
        super().__init__(
            name="defaults",
            description=(
                "GROMACS defaults: nbfunc, comb-rule, gen-pairs, "
                "fudgeLJ, fudgeQQ"
            ),
            category="defaults",
            contents=contents
        )
        # Force dict_contents to be contents for defaults
        self.dict_contents: list[list[str]] = self.contents


FF_Categories["defaults"] = Defaults


class Atomtypes(FF_Category):
    """
    Atom type definitions category.
    
    Each entry contains:
    atom type, (bonded type, atomic number,) mass, charge, particle type, V, W
    Bonded type and atomic number are optional.
    """
    
    def __init__(self, contents: list[list[str]] | None = None) -> None:
        """
        Initialize atomtypes category.
        
        Args:
            contents: List of atom type definitions.
        """
        if contents is None:
            contents = []
        super().__init__(
            name="atomtypes",
            description=(
                "atom type; (bonded type; atomic number;) mass (u); "
                "charge (e); particle type; V; W"
            ),
            category="atomtypes",
            contents=contents
        )

    def initialize(self) -> None:
        """Process atom type entries and fill missing optional fields."""
        for fields in self.contents:
            if len(fields) < 6:
                raise ValueError(f"Too few fields in [ atomtypes ] line: {fields}")
            
            # Check if bonded type and atomic number are missing
            if len(fields[3]) == 1:
                # Both bonded type and atomic number missing
                fields.insert(1, None)
                fields.insert(1, None)
            elif len(fields[4]) == 1 and fields[4].isalpha():
                if fields[1][0].isalpha():
                    # Atomic number is missing
                    fields.insert(2, None)
                else:
                    # Bonded type is missing
                    fields.insert(1, None)

            if len(fields) != 8:
                raise ValueError(
                    f"Wrong number of fields in [ atomtypes ] line: {fields}"
                )

            if fields[1] is not None:
                raise ValueError(f"Unsupported bonded types: {fields}")
                
            self.dict_contents[fields[0]] = fields


FF_Categories["atomtypes"] = Atomtypes


class Nonbond_params(FF_Category):
    """Non-bonded interaction parameters between atom types."""
    
    def __init__(self, contents: list[list[str]] | None = None) -> None:
        """
        Initialize nonbond_params category.
        
        Args:
            contents: List of non-bonded parameter entries.
        """
        if contents is None:
            contents = []
        super().__init__(
            name="nonbond_params",
            description="Non-bonded parameters between atom types (V; W)",
            category="nonbond_params",
            contents=contents
        )

    def initialize(self) -> None:
        """Process non-bonded parameter entries."""
        for fields in self.contents:
            if len(fields) < 5:
                raise ValueError(
                    f"Too few fields in [ nonbond_params ] line: {fields}"
                )
            if fields[2] != "1":
                raise ValueError(
                    f"Unsupported function type in [ nonbond_params ] line: {fields}"
                )
            # Store with sorted atom type pair as key
            key: tuple[str, str] = tuple(sorted(fields[:2]))
            self.dict_contents[key] = fields


FF_Categories["nonbond_params"] = Nonbond_params


__all__ = [
    "Molecule",
    "Molecule_Category",
    "Moleculetype",
    "Atoms",
    "Multiple_basin",
    "ForceField",
    "FF_Category",
    "Defaults",
    "Atomtypes",
    "Nonbond_params",
]
