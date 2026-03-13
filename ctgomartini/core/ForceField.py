"""
Force field parameter management module for CTGoMartini.

Provides classes for parsing and managing GROMACS-style force field parameters
including atom types, non-bonded parameters, and defaults.
"""

from __future__ import annotations

from typing import Any


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
