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

## Migrating from 0.6.0 to 1.0.0

### Vermouth Upgrade (>= 0.15.0)

CTGoMartini 1.0.0 requires vermouth >= 0.15.0, which has significant changes from 0.9.6:

#### 1. DSSP Handling

| Version | Command | Behavior |
|---------|---------|----------|
| 0.9.6 | `-dssp /path/to/dssp` | Requires external DSSP executable |
| 0.15.0 | `-dssp` (no argument) | Uses MDTraj (no external dependency) |

**Migration**: DSSP parameter is now optional in `ctgomartinize.py`:
```bash
# Old (0.6.0)
python ctgomartinize.py -s protein.pdb -m contact.map -mol State1 -dssp /usr/bin/dssp

# New (1.0.0) - uses MDTraj by default
python ctgomartinize.py -s protein.pdb -m contact.map -mol State1

# Or explicitly provide DSSP path
python ctgomartinize.py -s protein.pdb -m contact.map -mol State1 -dssp /usr/bin/dssp
```

#### 2. Side-Chain Fix (`-scfix`)

| Version | Default | To Enable | To Disable |
|---------|---------|-----------|------------|
| 0.9.6 | Disabled | Add `-scfix` | (default) |
| 0.15.0 | **Enabled** | (default) | Add `-noscfix` |

**Note**: CTGoMartini 1.0.0 removes explicit `-scfix` from martinize2 calls since it's now the default.

#### 3. Go Contacts Workflow

The Go contacts generation has been completely redesigned:

**Old Workflow (0.6.0)**:
```
1. martinize2 -f protein.pdb -dssp <cmd> -govs-include ...
2. create_go_virt_for_multimer -r protein.pdb -s protein_cg.pdb -f contact.map
3. Generates: BB-part-def_VirtGoSites.itp, go-table_VirtGoSites.itp
```

**New Workflow (1.0.0)**:
```
1. convert_map_format input.map output.out  # Convert rCSU .map to martinize2 format
2. martinize2 -f protein.pdb -dssp -go contact_map.out -go-low 0.3 -go-up 1.1 -go-eps 12.0
3. Generates: go_atomtypes.itp, go_nbparams.itp
```

**Key Differences**:
- Input format: rCSU `.map` → `contact_map.out` (OV+rCSU format)
- Go parameters now specified in martinize2 command
- Output files renamed and reformatted

### New Utility Module

Added `convert_map_format` for converting rCSU web-server contact maps:

```python
from ctgomartini.utils import convert_map_format

# Convert .map to contact_map.out
convert_map_format(
    input_file='contact.map',
    output_file='contact_map.out',
    pdb_name='protein.pdb',
    force=True
)
```

### API Changes in ctgomartinize.py

| Function | Change | Details |
|----------|--------|---------|
| `Martinize2()` | DSSP parameter now optional | `dssp: str \| None = None` |
| `Martinize2()` | Removed `-scfix` | Now default in vermouth 0.15.0 |
| `GenGoContacts()` | Complete rewrite | Now converts .map format instead of calling `create_go_virt_for_multimer` |
| `ModifyFF()` | Updated includes | Now includes `go_atomtypes.itp` and `go_nbparams.itp` |

### File Output Changes

| Old File (0.6.0) | New File (1.0.0) | Purpose |
|------------------|------------------|---------|
| `BB-part-def_VirtGoSites.itp` | `go_atomtypes.itp` | Go virtual site atom types |
| `go-table_VirtGoSites.itp` | `go_nbparams.itp` | Go interaction parameters |

### Position Restraints

Removed `position_restraints` section from final ITP in SB mode (not needed for production).

## Version

This migration guide applies to CTGoMartini >= 1.0.0
