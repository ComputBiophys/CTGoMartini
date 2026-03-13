"""
GROMACS topology file parser for Martini force fields.

Parses GROMACS-style .top files with preprocessor directives and constructs
a topology tree for conversion to OpenMM systems.

Authors: Song Yang
"""

from __future__ import annotations

import os
import re
from collections import OrderedDict
from collections.abc import Mapping
from typing import Any

try:
    from shutil import which
except ImportError:
    from distutils.spawn import find_executable as which

from .Molecule import Molecule
from .ForceField import ForceField


class Topology:
    """
    Parse a GROMACS Martini topology file.
    
    Constructs a topology tree that can be converted to OpenMM MartiniTopFile.
    Supports preprocessor directives (#include, #define, #ifdef, etc.).
    """

    def __init__(
        self,
        file: str,
        includeDir: str | None = None,
        defines: dict[str, Any] | None = None
    ) -> None:
        """
        Load a topology file.

        Args:
            file: Path to the topology file to load.
            includeDir: Directory for #include files. If None, attempts to
                locate GROMACS installation.
            defines: Preprocessor definitions. Example:
                {"FLEXIBLE": None, "PositionRestraint": "1000"}
        """
        if includeDir is None:
            includeDir = _get_default_gromacs_include_dir()
            
        self._includeDirs: tuple[str, ...] = (os.path.dirname(file), includeDir)
        self._defines: OrderedDict[str, Any] = OrderedDict()
        self._genpairs: bool = True

        if defines is not None:
            for define, value in defines.items():
                self._defines[define] = value

        self._currentCategory: str | None = None
        self._currentMoleculeType: Molecule | None = None
        self._ifStack: list[bool] = []
        self._elseStack: list[bool] = []

        # Force field and molecule data
        self.forcefield: ForceField = ForceField()
        self._moleculeTypes: dict[str, Molecule] = {}
        self.molecules: list[tuple[str, int]] = []
        self.moleculeTypes: dict[str, Molecule] = {}

        # Process the topology file
        self._process_file(file)
        
        # Extract molecule types used in the system
        self.moleculeTypes = {
            molecule_name: self._moleculeTypes[molecule_name]
            for (molecule_name, _) in self.molecules
        }

        self.forcefield.initialize()
        for molecule_name, molecule_type in self.moleculeTypes.items():
            molecule_type.initialize(self.forcefield._parameters)

    def _process_file(self, file: str) -> None:
        """
        Process a topology file line by line.
        
        Args:
            file: Path to the file to process.
        """
        append = ""
        with open(file) as lines:
            for line in lines:
                if line.strip().endswith("\\"):
                    append = f"{append} {line[: line.rfind('\\')]}"
                else:
                    self._process_line(append + " " + line, file)
                    append = ""

    def _process_line(self, line: str, file: str) -> None:
        """
        Process a single line from the topology file.
        
        Args:
            line: The line to process.
            file: Current file path (for #include resolution).
        """
        # Remove comments
        if ";" in line:
            line = line[: line.index(";")]
        stripped = line.strip()
        ignore = not all(self._ifStack)
        
        if stripped.startswith("*") or len(stripped) == 0:
            # Comment or empty line
            return

        elif stripped.startswith("[") and not ignore:
            # Start of a category
            if not stripped.endswith("]"):
                raise ValueError(f"Illegal line in .top file: {line}")
            self._currentCategory = stripped[1:-1].strip()

        elif stripped.startswith("#"):
            # Preprocessor command
            self._process_preprocessor(stripped, line, file)

        elif not ignore:
            # Replace defines and process data line
            line = _replace_defines(line, self._defines)

            if self._currentCategory is None:
                raise ValueError(f"Unexpected line in .top file: {line}")
                
            self._process_category_line(line)

    def _process_preprocessor(self, stripped: str, line: str, file: str) -> None:
        """
        Process a preprocessor directive.
        
        Args:
            stripped: Stripped line content.
            line: Original line content.
            file: Current file path.
        """
        fields = stripped.split()
        command = fields[0]
        
        if len(self._ifStack) != len(self._elseStack):
            raise RuntimeError("#if/#else stack out of sync")

        if command == "#include":
            self._process_include(stripped, file)
        elif command == "#define":
            self._process_define(stripped, line)
        elif command == "#ifdef":
            self._process_ifdef(fields, True)
        elif command == "#ifndef":
            self._process_ifdef(fields, False)
        elif command == "#undef":
            self._process_undef(fields)
        elif command == "#endif":
            self._process_endif()
        elif command == "#else":
            self._process_else(stripped)

    def _process_include(self, stripped: str, file: str) -> None:
        """Process #include directive."""
        command = "#include"
        name = stripped[len(command):].strip(' \t"<>')
        searchDirs = self._includeDirs + (os.path.dirname(file),)
        
        for dir in searchDirs:
            include_file = os.path.join(dir, name)
            if os.path.isfile(include_file):
                self._process_file(include_file)
                break
        else:
            raise ValueError(f"Could not locate #include file: {name}")

    def _process_define(self, stripped: str, line: str) -> None:
        """Process #define directive."""
        fields = stripped.split()
        if len(fields) < 2:
            raise ValueError(f"Illegal line in .top file: {line}")
        name = fields[1]
        value_start = stripped.find(name, len("#define")) + len(name) + 1
        value = line[value_start:].strip()
        value = value or "1"
        self._defines[name] = value

    def _process_ifdef(self, fields: list[str], ifdef: bool) -> None:
        """Process #ifdef or #ifndef directive."""
        if len(fields) < 2:
            raise ValueError(f"Illegal line in .top file")
        name = fields[1]
        condition = (name in self._defines) if ifdef else (name not in self._defines)
        self._ifStack.append(condition)
        self._elseStack.append(False)

    def _process_undef(self, fields: list[str]) -> None:
        """Process #undef directive."""
        if len(fields) < 2:
            raise ValueError("Illegal line in .top file")
        if fields[1] in self._defines:
            self._defines.pop(fields[1])

    def _process_endif(self) -> None:
        """Process #endif directive."""
        if len(self._ifStack) == 0:
            raise ValueError("Unexpected #endif in .top file")
        del self._ifStack[-1]
        del self._elseStack[-1]

    def _process_else(self, stripped: str) -> None:
        """Process #else directive."""
        if len(self._ifStack) == 0:
            raise ValueError("Unexpected #else in .top file")
        if self._elseStack[-1]:
            raise ValueError(f"#else has already been used: {stripped}")
        self._ifStack[-1] = not self._ifStack[-1]
        self._elseStack[-1] = True

    def _process_category_line(self, line: str) -> None:
        """Process a line within a category."""
        if self._currentCategory == "molecules":
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Too few fields in [ molecules ] line: {line}")
            self.molecules.append((fields[0], int(fields[1])))
            
        elif self._currentCategory == "moleculetype":
            fields = line.split()
            if len(fields) < 1:
                raise ValueError(f"Too few fields in [ moleculetype ] line: {line}")
            mol_type = Molecule(name=fields[0])
            self._moleculeTypes[fields[0]] = mol_type
            self._currentMoleculeType = mol_type
            self._currentMoleculeType.readLine(line, self._currentCategory)
            
        elif self._currentCategory in [
            "defaults", "atomtypes", "bondtypes", "pairtypes",
            "angletypes", "dihedraltypes", "constrainttypes", "nonbond_params"
        ]:
            self.forcefield.readLine(line, self._currentCategory)
            
        elif self._currentCategory == "system":
            pass  # Ignore system description
            
        elif self._currentCategory == "position_restraints":
            raise ValueError(
                "[ position_restraints ] is not currently supported.\n"
                "Use CustomExternalForce in OpenMM instead."
            )
        else:
            if self._currentMoleculeType is not None:
                self._currentMoleculeType.readLine(line, self._currentCategory)
            else:
                raise ValueError(f"[ {self._currentCategory} ] is not currently supported.")


def _find_all_instances_in_string(string: str, substr: str) -> list[int]:
    """
    Find all indices of a substring in a string.
    
    Args:
        string: The string to search in.
        substr: The substring to find.
        
    Returns:
        List of indices where substr occurs.
    """
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx + 1)
    return indices


def _replace_defines(line: str, defines: Mapping[str, Any]) -> str:
    """
    Replace defined tokens in a line.
    
    Args:
        line: The line to process.
        defines: Dictionary of defines to replace.
        
    Returns:
        Line with defines replaced.
    """
    novarcharre = re.compile(r"\W")
    if not defines:
        return line
        
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices:
            continue
            
        inside = ""
        idx = 0
        n_to_skip = 0
        new_line = []
        
        for i, char in enumerate(line):
            if n_to_skip:
                n_to_skip -= 1
                continue
            if char in ("'\""):
                if not inside:
                    inside = char
                elif inside == char:
                    inside = ""
            if idx < len(indices) and i == indices[idx]:
                if inside:
                    new_line.append(char)
                    idx += 1
                    continue
                if i == 0 or novarcharre.match(line[i - 1]):
                    endidx = indices[idx] + len(define)
                    if endidx >= len(line) or novarcharre.match(line[endidx]):
                        new_line.extend(list(value))
                        n_to_skip = len(define) - 1
                        idx += 1
                        continue
                idx += 1
            new_line.append(char)
        line = "".join(new_line)

    return line


def _get_default_gromacs_include_dir() -> str:
    """
    Find the GROMACS include directory.
    
    Searches for GROMACS in order:
    1. GMXDATA environment variable
    2. GMXBIN environment variable
    3. pdb2gmx executable in PATH
    4. gmx executable in PATH
    5. Default location /usr/local/gromacs/share/gromacs/top
    
    Returns:
        Path to the GROMACS topology directory.
    """
    if "GMXDATA" in os.environ:
        return os.path.join(os.environ["GMXDATA"], "top")
        
    if "GMXBIN" in os.environ:
        return os.path.abspath(
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        )

    pdb2gmx_path = which("pdb2gmx")
    if pdb2gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(pdb2gmx_path), "..", "share", "gromacs", "top")
        )
    else:
        gmx_path = which("gmx")
        if gmx_path is not None:
            return os.path.abspath(
                os.path.join(os.path.dirname(gmx_path), "..", "share", "gromacs", "top")
            )

    return "/usr/local/gromacs/share/gromacs/top"
