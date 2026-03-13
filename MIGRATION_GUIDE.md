# CTGoMartini Migration Guide

This document describes the breaking changes in the recent refactoring and provides a mapping from old to new APIs.

## Module Structure Changes

### Directory Renaming

| Old Path | New Path | Description |
|----------|----------|-------------|
| `ctgomartini/func/` | `ctgomartini/core/` | Core functionality modules |
| `ctgomartini/util/` | `ctgomartini/utils/` | Utility functions |

### Import Statement Updates

```python
# Old imports (NO LONGER WORK)
from ctgomartini.func import Topology, Molecule
from ctgomartini.func import WriteItp, ConvertLongShortElasticBonds
from ctgomartini.util import read_inputs
from ctgomartini.util import SameListList

# New imports
from ctgomartini.core import Topology, Molecule
from ctgomartini.utils import WriteItp, ConvertLongShortElasticBonds
from ctgomartini.utils import read_inputs
from ctgomartini.core import SameListList  # Note: SameListList is in core
```

## Function Name Changes

### `data/run_ctgomartini.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `ReportTime` | `report_time` | Report elapsed simulation time |
| `LoadStructure` | `load_structure` | Load GRO or PDB structure file |
| `gen_restraints` | `generate_restraints` | Generate position restraint file |
| `restraints` | `add_restraints` | Add positional restraints to system |
| `BackupFile` | `backup_file` | Create backup of existing file |
| `WriteOutput` | `write_output` | Write simulation output |
| `WriteCheckPoint` | `write_checkpoint` | Save checkpoint to XML file |
| `Cleanup` | `cleanup` | Handle signal interrupts |
| `LoadPlatform` | `load_platform` | Configure OpenMM platform |
| `GenerateTopology` | `generate_topology` | Generate topology from input |

### `utils/WriteItp.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `WriteItp` | `write_itp` | Write GROMACS ITP files |

### `utils/ConvertLongShortElasticBonds.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `BB_Distance` | `bb_distance` | Calculate backbone distance |
| `ConvertLongShortElasticBonds` | `convert_long_short_elastic_bonds` | Convert elastic bonds |

### `utils/Create_goVirt_for_multimer.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `Create_goVirt_for_multimer` | `create_go_virt_for_multimer` | Create virtual sites for multimer |

## Class and Function Location Reference

### `ctgomartini.core` (Core Functionality)

| Name | Type | Description |
|------|------|-------------|
| `Topology` | Class | GROMACS topology parser |
| `Molecule` | Class | Molecule representation |
| `ForceField` | Class | Force field parameter handling |
| `VSiteManager` | Class | Virtual site management |
| `LinearSite` | Class | Linear virtual site definition |
| `COMLinearSite` | Class | Center-of-mass virtual site |
| `OutOfPlane` | Class | Out-of-plane virtual site |
| `NormalizedInPlaneSite` | Class | Normalized in-plane virtual site |
| `SameListList` | Function | Compare lists with tolerance |
| `CombineMols` | Class | Molecule combination utilities |
| `BondedInteraction_types` | List | Registered bonded interaction types |

### `ctgomartini.utils` (Utilities)

| Name | Type | Description |
|------|------|-------------|
| `read_inputs` | Function | Read simulation input parameters |
| `write_itp` | Function | Write GROMACS ITP files |
| `ConvertLongShortElasticBonds` | Class | Bond conversion utility |
| `Create_goVirt_for_multimer` | Function | Multimer virtual site generation |

### `ctgomartini.api` (High-level API)

| Name | Type | Description |
|------|------|-------------|
| `MartiniTopFile` | Class | Main topology file handler |
| `GenMBPTop` | Function | Generate multiple-basin topology |
| `GenSBPTop` | Function | Generate single-basin topology |

### `ctgomartini.data` (Executable Scripts)

| Name | Type | Description |
|------|------|-------------|
| `mdrun` | Function | Main MD simulation runner |
| `generate_topology` | Function | Generate OpenMM topology |
| `report_time` | Function | Report simulation time |
| `load_structure` | Function | Load molecular structure |
| `generate_restraints` | Function | Generate restraint files |
| `add_restraints` | Function | Add restraints to system |
| `backup_file` | Function | Backup existing files |
| `write_output` | Function | Write simulation output |
| `write_checkpoint` | Function | Write checkpoint files |
| `cleanup` | Function | Cleanup on interrupt |
| `load_platform` | Function | Load OpenMM platform |

## Type Annotations

All public functions now have full type annotations using Python 3.10+ syntax:

```python
def load_structure(str_file: str) -> tuple[GromacsGroFile | PDBFile, mm.Vec3 | None]:
    ...

def generate_restraints(
    str_file: str, 
    atomname: str, 
    fc: float = 1000.0, 
    rest_file: str = "restraints.txt"
) -> None:
    ...
```

## Example Migration

### Before
```python
from ctgomartini.func import Topology, Molecule, WriteItp
from ctgomartini.util import read_inputs, SameListList
from ctgomartini.data.run_ctgomartini import mdrun, LoadStructure, gen_restraints

def main():
    inputs = read_inputs('npt.inp')
    conf, box = LoadStructure(inputs.input)
    gen_restraints(inputs.input, 'BB', 1000, 'restraints.txt')
```

### After
```python
from ctgomartini.core import Topology, Molecule
from ctgomartini.utils import read_inputs, write_itp
from ctgomartini.core import SameListList
from ctgomartini.data.run_ctgomartini import mdrun, load_structure, generate_restraints

def main():
    inputs = read_inputs('npt.inp')
    conf, box = load_structure(inputs.input)
    generate_restraints(inputs.input, 'BB', 1000, 'restraints.txt')
```

## Development Tools

The package now uses modern Python packaging:
- `pyproject.toml` for package configuration
- `setup.py` is minimal (backward compatibility only)
- Full type annotations throughout
- Google-style docstrings
- MIT License

## Version

This migration guide applies to CTGoMartini >= 0.6.0
