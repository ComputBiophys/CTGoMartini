# CTGoMartini Migration Guide

This document describes the breaking changes and provides a mapping from old to new APIs.

## Version 1.0.0 Major Changes

### Module Structure Overhaul

The biggest change in v1.0.0 is the consolidation of `api/` and `core/` into a unified `topology/` module:

| Old Path (0.6.x) | New Path (1.0.0) | Status |
|------------------|------------------|--------|
| `ctgomartini/api/` | `ctgomartini/topology/` | ✅ Reorganized |
| `ctgomartini/core/` | `ctgomartini/topology/` | ✅ Merged into topology |
| `ctgomartini/data/*.py` | `ctgomartini/cli/` | ✅ Moved to CLI module |
| `ctgomartini/data/run_*.py` | Deleted | ✅ Replaced by CLI commands |

### Import Statement Updates (1.0.0)

```python
# Old imports (0.6.x - NO LONGER WORK)
from ctgomartini.api import MartiniTopFile, GenMBPTop, GenSBPTop
from ctgomartini.core import Topology, Molecule, ForceField
from ctgomartini.util import read_inputs

# New imports (1.0.0)
from ctgomartini.topology import (
    TopologyParser,      # was Topology
    MartiniTopFile,
    Molecule,
    ForceField,
    create_sb_topology,  # was GenSBPTop
    create_mb_topology,  # was GenMBPTop
)
from ctgomartini.simulation import load_config  # was read_inputs
```

### New Modules in 1.0.0

| Module | Purpose | Key Components |
|--------|---------|----------------|
| `ctgomartini.simulation` | Simulation execution | `SimulationRunner`, `MDRunner`, `REMDRunner`, `load_config` |
| `ctgomartini.cli` | Command-line interface | `run_ctgomartini`, `ctgomartinize` |
| `ctgomartini.topology.generator` | Topology generation | `create_sb_topology`, `create_mb_topology`, `combine_*` |
| `ctgomartini.topology.interactions` | Force interactions | All bonded/non-bonded interaction classes |

### Function Name Changes (1.0.0)

| Old Name (0.6.x) | New Name (1.0.0) | Location |
|------------------|------------------|----------|
| `Topology` | `TopologyParser` | `topology.parser` |
| `GenMBPTop` | `create_mb_topology` | `topology.generator` |
| `GenSBPTop` | `create_sb_topology` | `topology.generator` |
| `read_inputs` | `load_config` | `simulation.config` |
| `mdrun` | `run_ctgomartini` (CLI) | `cli` |

### CLI Changes (1.0.0)

Old way (0.6.x):
```bash
python ctgomartini/data/run_ctgomartini.py -i npt.inp
python ctgomartini/data/ctgomartinize.py -s protein.pdb ...
```

New way (1.0.0):
```bash
# Installed as entry points
run_ctgomartini -i npt.inp
ctgomartinize -s protein.pdb ...

# Or as module
python -m ctgomartini.cli.run_ctgomartini -i npt.inp
```

---

## Historical Changes (0.5.x → 0.6.x)

### Directory Renaming

| Old Path | New Path | Description |
|----------|----------|-------------|
| `ctgomartini/func/` | `ctgomartini/core/` | Core functionality modules |
| `ctgomartini/util/` | `ctgomartini/utils/` | Utility functions |

### Import Statement Updates (0.6.x)

```python
# Old imports (NO LONGER WORK)
from ctgomartini.func import Topology, Molecule
from ctgomartini.func import WriteItp, ConvertLongShortElasticBonds
from ctgomartini.util import read_inputs
from ctgomartini.util import SameListList

# New imports (0.6.x)
from ctgomartini.core import Topology, Molecule
from ctgomartini.utils import WriteItp, ConvertLongShortElasticBonds
from ctgomartini.utils import read_inputs
from ctgomartini.core import SameListList  # Note: SameListList is in core
```

### Function Name Changes (0.6.x)

#### `data/run_ctgomartini.py`

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

#### `utils/WriteItp.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `WriteItp` | `write_itp` | Write GROMACS ITP files |

#### `utils/ConvertLongShortElasticBonds.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `BB_Distance` | `bb_distance` | Calculate backbone distance |
| `ConvertLongShortElasticBonds` | `convert_long_short_elastic_bonds` | Convert elastic bonds |

#### `utils/Create_goVirt_for_multimer.py`

| Old Name | New Name | Description |
|----------|----------|-------------|
| `Create_goVirt_for_multimer` | `create_go_virt_for_multimer` | Create virtual sites for multimer |

### Class and Function Location Reference (0.6.x)

#### `ctgomartini.core` (Core Functionality)

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

#### `ctgomartini.utils` (Utilities)

| Name | Type | Description |
|------|------|-------------|
| `read_inputs` | Function | Read simulation input parameters |
| `write_itp` | Function | Write GROMACS ITP files |
| `ConvertLongShortElasticBonds` | Class | Bond conversion utility |
| `Create_goVirt_for_multimer` | Function | Multimer virtual site generation |

#### `ctgomartini.api` (High-level API - 0.6.x)

| Name | Type | Description |
|------|------|-------------|
| `MartiniTopFile` | Class | Main topology file handler |
| `GenMBPTop` | Function | Generate multiple-basin topology |
| `GenSBPTop` | Function | Generate single-basin topology |

---

## Migrating from 0.6.x to 1.0.0

### Step 1: Update Topology Imports

```python
# Before (0.6.x)
from ctgomartini.api import MartiniTopFile, GenMBPTop, GenSBPTop
from ctgomartini.core import Topology, Molecule, ForceField

# After (1.0.0)
from ctgomartini.topology import (
    TopologyParser,      # Topology was renamed
    MartiniTopFile,
    Molecule,
    ForceField,
    create_sb_topology,  # GenSBPTop renamed
    create_mb_topology,  # GenMBPTop renamed
)
```

### Step 2: Update Simulation Imports

```python
# Before (0.6.x)
from ctgomartini.utils import read_inputs

# After (1.0.0)
from ctgomartini.simulation import load_config
```

### Step 3: Update Function Calls

```python
# Before (0.6.x)
from ctgomartini.api import GenMBPTop
mbmol = GenMBPTop(mols_list, mbmol_name, dict_cutoffs)

# After (1.0.0)
from ctgomartini.topology import create_mb_topology
mbmol = create_mb_topology(mols_list, mbmol_name, dict_cutoffs)
```

### Step 4: Update Scripts to CLI Commands

```bash
# Before (0.6.x)
python ctgomartini/data/run_ctgomartini.py -i npt.inp

# After (1.0.0)
run_ctgomartini -i npt.inp
```

---

## Complete Example Migration

### Before (0.6.x)
```python
from ctgomartini.api import MartiniTopFile, GenMBPTop
from ctgomartini.core import Topology, Molecule
from ctgomartini.utils import read_inputs, write_itp
from ctgomartini.core import SameListList

def main():
    inputs = read_inputs('npt.inp')
    top = MartiniTopFile(inputs.topol)
    mbmol = GenMBPTop(mols_list, 'MBMol', dict_cutoffs)
    write_itp(mbmol)
```

### After (1.0.0)
```python
from ctgomartini.topology import (
    MartiniTopFile, 
    create_mb_topology,
    write_itp
)
from ctgomartini.simulation import load_config

def main():
    inputs = load_config('npt.inp')
    top = MartiniTopFile(inputs.topol)
    mbmol = create_mb_topology(mols_list, 'MBMol', dict_cutoffs)
    write_itp(mbmol)
```

---

## Vermouth Upgrade (>= 0.15.0)

CTGoMartini 1.0.0 requires vermouth >= 0.15.0:

### DSSP Handling

| Version | Command | Behavior |
|---------|---------|----------|
| 0.9.6 | `-dssp /path/to/dssp` | Requires external DSSP executable |
| 0.15.0 | `-dssp` (no argument) | Uses MDTraj (no external dependency) |

### Side-Chain Fix (`-scfix`)

| Version | Default | To Enable | To Disable |
|---------|---------|-----------|------------|
| 0.9.6 | Disabled | Add `-scfix` | (default) |
| 0.15.0 | **Enabled** | (default) | Add `-noscfix` |

### Go Contacts Workflow

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

### File Output Changes

| Old File (0.6.0) | New File (1.0.0) | Purpose |
|------------------|------------------|---------|
| `BB-part-def_VirtGoSites.itp` | `go_atomtypes.itp` | Go virtual site atom types |
| `go-table_VirtGoSites.itp` | `go_nbparams.itp` | Go interaction parameters |

---

## Version Information

- **This guide**: CTGoMartini >= 1.0.0
- **Latest stable**: 1.0.0
- **Python requirement**: >= 3.12

For questions or issues, please refer to the [GitHub repository](https://github.com/ComputBiophys/CTGoMartini).
