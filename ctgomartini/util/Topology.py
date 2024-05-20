"""
Used for loading martini top files with Gromacs styles.

Authors: Song Yang
Last update: 20231001
"""
__author__ = "Song Yang"
__version__ = "1.0"

import os, re
from collections import OrderedDict, defaultdict
import distutils.spawn
from .Molecule import Molecule
from .ForceField import ForceField

class Topology():
    """Parse a Gromacs Martini top file and constructs a Topology tree which will be conveniently converted into Opemm MartiniTopFile"""

    def __init__(
            self,
            file,
            includeDir=None,
            defines=None
    ):
        """Load a top file.

        Parameters
        ----------
        file : str
            the name of the file to load
        includeDir : string=None
            A directory in which to look for other files included from the
            top file. If not specified, we will attempt to locate a gromacs
            installation on your system. When gromacs is installed in
            /usr/local, this will resolve to /usr/local/gromacs/share/gromacs/top
        defines : dict={}
            preprocessor definitions that should be predefined when parsing the file. e.g. {"FLEXIBLE": None, "PositionRestraint":"1000"}
        """

        if includeDir is None:
            includeDir = _get_default_gromacs_include_dir()
        self._includeDirs = (os.path.dirname(file), includeDir)
        self._defines = OrderedDict()
        self._genpairs = True

        if defines is not None:
            for define, value in defines.items():
                self._defines[define] = value

        self._currentCategory = None
        self._currentMoleculeType = None
        self._ifStack = []
        self._elseStack = []

        # forcefield.parameters, molecule.parameters
        self.forcefield = ForceField()
        self._moleculeTypes = {} # all moluecule types
        self.molecules = []
        self.moleculeTypes = {} # molecule types used in the system

        # Process the topology file
        self._process_file(file)
        # Extract the moleculetypes actually used in the system
        self.moleculeTypes = {molecule_name: self._moleculeTypes[molecule_name] for (molecule_name, count) in self.molecules}

        self.forcefield.initialize()
        for molecule_name, molecule_type in self.moleculeTypes.items():
            molecule_type.initialize(self.forcefield._parameters)
        

    

    def _process_file(self, file):
        append = ""
        with open(file) as lines:
            for line in lines:
                if line.strip().endswith("\\"):
                    append = "%s %s" % (append, line[: line.rfind("\\")])
                else:
                    self._process_line(append + " " + line, file)
                    append = ""

    def _process_line(self, line, file):
        """Process one line from a file."""
        if ";" in line:
            line = line[: line.index(";")]
        stripped = line.strip()
        ignore = not all(self._ifStack)
        if stripped.startswith("*") or len(stripped) == 0:
            # A comment or empty line.
            return

        elif stripped.startswith("[") and not ignore:
            # The start of a category.
            if not stripped.endswith("]"):
                raise ValueError("Illegal line in .top file: " + line)
            self._currentCategory = stripped[1:-1].strip()

        elif stripped.startswith("#"):
            # A preprocessor command.
            fields = stripped.split()
            command = fields[0]
            if len(self._ifStack) != len(self._elseStack):
                raise RuntimeError("#if/#else stack out of sync")

            if command == "#include" and not ignore:
                # Locate the file to include
                name = stripped[len(command) :].strip(' \t"<>')
                searchDirs = self._includeDirs + (os.path.dirname(file),)
                for dir in searchDirs:
                    file = os.path.join(dir, name)
                    if os.path.isfile(file):
                        # We found the file, so process it.
                        self._process_file(file)
                        break
                else:
                    raise ValueError("Could not locate #include file: " + name)
            elif command == "#define" and not ignore:
                # Add a value to our list of defines.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                valueStart = stripped.find(name, len(command)) + len(name) + 1
                value = line[valueStart:].strip()
                value = value or "1"  # Default define is 1
                self._defines[name] = value
            elif command == "#ifdef":
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                self._ifStack.append(name in self._defines)
                self._elseStack.append(False)
            elif command == "#undef":
                # Un-define a variable
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                if fields[1] in self._defines:
                    self._defines.pop(fields[1])
            elif command == "#ifndef":
                # See whether this block should be ignored.
                if len(fields) < 2:
                    raise ValueError("Illegal line in .top file: " + line)
                name = fields[1]
                self._ifStack.append(name not in self._defines)
                self._elseStack.append(False)
            elif command == "#endif":
                # Pop an entry off the if stack.
                if len(self._ifStack) == 0:
                    raise ValueError("Unexpected line in .top file: " + line)
                del self._ifStack[-1]
                del self._elseStack[-1]
            elif command == "#else":
                # Reverse the last entry on the if stack
                if len(self._ifStack) == 0:
                    raise ValueError("Unexpected line in .top file: " + line)
                if self._elseStack[-1]:
                    raise ValueError(
                        "Unexpected line in .top file: "
                        "#else has already been used " + line
                    )
                self._ifStack[-1] = not self._ifStack[-1]
                self._elseStack[-1] = True

        elif not ignore:
            # Gromacs occasionally uses #define's to introduce specific
            # parameters for individual terms (for instance, this is how
            # ff99SB-ILDN is implemented). So make sure we do the appropriate
            # pre-processor replacements necessary
            line = _replace_defines(line, self._defines)

            # A line of data for the current category
            if self._currentCategory is None:
                raise ValueError("Unexpected line in .top file: " + line)
            if self._currentCategory == "molecules":
                """Process a line in the [ molecules ] category."""
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError("Too few fields in [ molecules ] line: " + line)
                self.molecules.append((fields[0], int(fields[1])))
            elif self._currentCategory == "moleculetype":
                """Process a line in the [ moleculetype ] category."""
                fields = line.split()
                if len(fields) < 1:
                    raise ValueError("Too few fields in [ moleculetypes ] line: " + line)
                type = Molecule(name = fields[0])
                self._moleculeTypes[fields[0]] = type
                self._currentMoleculeType = type
                self._currentMoleculeType.readLine(line, self._currentCategory)
            elif self._currentCategory in ["defaults", "atomtypes", "bondtypes", "pairtypes", "angletypes", "dihedraltypes", "constrainttypes", "nonbond_params"]:
                """Process a line for FF parameters."""
                self.forcefield.readLine(line, self._currentCategory)
            elif self._currentCategory == "system":  # ignore the system description
                pass
            elif self._currentCategory == "position_restraints":
                raise ValueError(
                    "[ position_restraints ] is not currently supported.\n"
                    "The same effect can be achieved using a CustomExternalForce in OpenMM.\n"
                    "Please see github.com/maccallumlab/martini_openmm for details."
                )
            else:
                if self._currentMoleculeType is not None:
                    self._currentMoleculeType.readLine(line, self._currentCategory)
                else:
                    raise ValueError(
                        f"[ {self._currentCategory} ] is not currently supported.\n"
                    )



def _find_all_instances_in_string(string, substr):
    """Find indices of all instances of substr in string"""
    indices = []
    idx = string.find(substr, 0)
    while idx > -1:
        indices.append(idx)
        idx = string.find(substr, idx + 1)
    return indices


def _replace_defines(line, defines):
    """Replaces defined tokens in a given line"""
    novarcharre = re.compile(r"\W")
    if not defines:
        return line
    for define in reversed(defines):
        value = defines[define]
        indices = _find_all_instances_in_string(line, define)
        if not indices:
            continue
        # Check to see if it's inside of quotes
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
                else:
                    if inside == char:
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

def _get_default_gromacs_include_dir():
    """Find the location where gromacs #include files are referenced from, by
    searching for (1) gromacs environment variables, (2) for the gromacs binary
    'pdb2gmx' or 'gmx' in the PATH, or (3) just using the default gromacs
    install location, /usr/local/gromacs/share/gromacs/top"""
    if "GMXDATA" in os.environ:
        return os.path.join(os.environ["GMXDATA"], "top")
    if "GMXBIN" in os.environ:
        return os.path.abspath(
            os.path.join(os.environ["GMXBIN"], "..", "share", "gromacs", "top")
        )

    pdb2gmx_path = distutils.spawn.find_executable("pdb2gmx")
    if pdb2gmx_path is not None:
        return os.path.abspath(
            os.path.join(os.path.dirname(pdb2gmx_path), "..", "share", "gromacs", "top")
        )
    else:
        gmx_path = distutils.spawn.find_executable("gmx")
        if gmx_path is not None:
            return os.path.abspath(
                os.path.join(os.path.dirname(gmx_path), "..", "share", "gromacs", "top")
            )

    return "/usr/local/gromacs/share/gromacs/top"