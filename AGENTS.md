# CTGoMartini Agent Guide

## Project Overview

CTGoMartini is a Python package for simulating protein conformational transitions using Go-Martini coarse-grained models. It provides tools for setting up and running molecular dynamics simulations with single-basin, switching, and multiple-basin Go-Martini potentials, with particular focus on membrane proteins.

**Version**: 1.0.0  
**Python Requirement**: >= 3.12  
**License**: MIT  
**Development Status**: Beta (4 - Beta)

### Key Features
- **New "Contacts" Interaction Type**: Replaces LJ-type contact interactions in classic GoMartini3 models
- **Automatic Topology Construction**: Generate GoMartini3 model topologies for different simulation methods (Single-basin, Multiple-basin, Switching)
- **Simplified MD Workflow**: Streamlined process for running coarse-grained molecular dynamics simulations with OpenMM
- **Modern Python**: Full type annotations using Python 3.12+ syntax

## Technology Stack

### Core Dependencies
- **OpenMM** (>=8.1): Molecular dynamics simulation engine
- **MDAnalysis** (>=2.4): Trajectory analysis
- **vermouth** (==0.14.0): Martinize2 integration for coarse-graining with native Go model support
- **NumPy**: Numerical computations
- **Pandas**: Data manipulation
- **Matplotlib** (>=3.3.3): Plotting and visualization

### Optional Dependencies
- **openmm-plumed** (>=2.0.0): PLUMED integration for enhanced sampling
- **dssp** (<4.0): Optional, for secondary structure assignment. Not required if using MDTraj (default)

### Development Dependencies
- **pytest** (>=7.0): Testing framework
- **pytest-cov**: Coverage reporting
- **black**: Code formatting (line-length: 100)
- **ruff**: Linting and import sorting
- **mypy**: Static type checking

## Project Structure

```
CTGoMartini/
├── ctgomartini/              # Main package directory
│   ├── __init__.py           # Package initialization
│   ├── _version.py           # Version information
│   ├── topology/             # Topology parsing and construction
│   │   ├── parser.py         # TopologyParser (GROMACS .top parser)
│   │   ├── models.py         # Molecule, ForceField classes
│   │   ├── builder.py        # MartiniTopFile (OpenMM system builder)
│   │   ├── generator/        # Topology generators
│   │   │   ├── single_basin.py   # create_sb_topology
│   │   │   ├── multi_basin.py    # create_mb_topology
│   │   │   └── combiner.py       # combine_* functions
│   │   └── interactions/     # Bonded and non-bonded interactions
│   │       ├── base.py
│   │       ├── bonds.py
│   │       ├── angles.py
│   │       ├── dihedrals.py
│   │       ├── contacts.py
│   │       ├── virtual_sites.py
│   │       ├── multi_state.py
│   │       ├── mixing.py
│   │       ├── registry.py
│   │       └── vsites.py
│   │
│   ├── simulation/           # Simulation execution
│   │   ├── config.py         # SimulationConfig, load_config
│   │   ├── base.py           # SimulationRunner base class
│   │   ├── md.py             # MDRunner (standard MD)
│   │   └── remd.py           # REMDRunner (replica exchange)
│   │
│   ├── cli/                  # Command-line interface
│   │   ├── run_ctgomartini.py   # Unified simulation runner
│   │   └── ctgomartinize.py     # Topology generation
│   │
│   ├── utils/                # Utility functions
│   │   ├── write_itp.py
│   │   ├── contacts.py
│   │   ├── convert_map_format.py
│   │   ├── convert_long_short_elastic_bonds.py
│   │   ├── constraints_to_bonds.py
│   │   └── pdb_validation.py
│   │
│   ├── analysis/             # Analysis tools
│   │   ├── drms_analysis.py
│   │   ├── remd_drms_analysis.py
│   │   ├── remd_exchange_ratio.py
│   │   ├── remd_free_energy.py
│   │   ├── remd_mbar_analysis.py
│   │   └── remd_replica_state.py
│   │
│   └── data/                 # Data files and templates
│       ├── ForceFields/Martini300/
│       ├── Membrane/
│       └── Soluble/
│
├── tutorials/                # Tutorial examples
│   ├── GlnBP/
│   └── TRAAK/
├── pyproject.toml
├── tests/                    # Test suite
├── MIGRATION_GUIDE.md
└── README.md
```

### Module Organization

| Module | Purpose | Key Classes/Functions |
|--------|---------|----------------------|
| `topology/` | Topology parsing, building, and interactions | `TopologyParser`, `MartiniTopFile`, `Molecule`, `ForceField`, `create_sb_topology`, `create_mb_topology` |
| `simulation/` | Simulation execution | `SimulationRunner`, `MDRunner`, `REMDRunner`, `SimulationConfig`, `load_config` |
| `cli/` | Command-line interface | `run_ctgomartini`, `ctgomartinize` |
| `utils/` | I/O utilities and helpers | `write_itp`, `convert_map_format`, `convert_long_short_elastic_bonds` |
| `analysis/` | Post-simulation analysis | `DRMSAnalyzer`, `remd_*` |
| `data/` | Data files and templates | Force fields, membrane/soluble templates |

## Build and Installation

### Installation Commands

```bash
# 1. Create conda environment
conda create -n ctgomartini python=3.12 -y
conda activate ctgomartini

# 2. Install core dependencies (note: NumPy < 2.4.0 is required)
pip install MDAnalysis "numpy>=2.0,<2.4.0" pandas matplotlib "vermouth==0.14.0"

# 3. Install OpenMM with CUDA support (adjust cuda-version as needed)
conda install -c conda-forge openmm cuda-version=12 -y

# 4. Install Jupyter for analysis notebooks
conda install jupyter -y

# 5. Install CTGoMartini
pip install -e .

# 6. Install REMD dependencies (optional, only for replica exchange simulations)
conda config --add channels omnia --add channels conda-forge
conda install openmmtools "numpy<2.4.0" "jax<=0.8.0" -y

# 7. Install testing framework
pip install pytest pytest-cov

# 8. For multi-GPU REMD, also install MPI support (optional)
# conda install -c conda-forge mpi4py mpich=3 -y
```

**Important Version Constraints**:
- **NumPy**: `>=2.0,<2.4.0` (REMD fails with 2.4.0+ due to openmmtools compatibility)
- **vermouth**: `==0.14.0` (pinned for force field consistency - includes proline helix fix)
- **CUDA**: Match your NVIDIA driver version (10.2, 11, or 12)

### Package Configuration (pyproject.toml)

The project uses modern Python packaging with `pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"
```

Key configurations:
- **Package discovery**: `ctgomartini*` packages included
- **Package data**: All files under `ctgomartini/data/` are included
- **Python version**: >=3.12 required

## Testing

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=ctgomartini
```

### Test Structure

| Category | Description |
|----------|-------------|
| `tests/unit/` | Unit tests for individual functions |
| `tests/integration/` | Integration tests |
| `tests/e2e/` | End-to-end workflow tests |

### NumPy Version Compatibility

**Required**: `numpy>=2.0, <2.4.0`

**Note**: NumPy 2.4.0+ causes REMD simulations to fail with:
```
TypeError: only 0-dimensional arrays can be converted to Python scalars
```
This is due to changes in NumPy's array scalar handling that affect `openmmtools`.

### Vermouth Version Compatibility

**Required**: `vermouth==0.14.0`

**Note**: vermouth 0.15.0+ uses MDTraj by default for secondary structure assignment.

#### Important: Version 0.13.0 vs 0.14.0 Force Field Differences

Starting from vermouth 0.14.0, the Martini 3.0 force field includes an important fix for proline (PRO) residues in helical conformations:

| Feature | vermouth 0.13.0 | vermouth >=0.14.0 |
|---------|-----------------|-------------------|
| **Angle type for helical PRO** | Harmonic (type 2) | Restricted bending (type 10) |
| **Physical behavior** | Allows angles approaching 180° | Prevents 180° angles (barrier at 180°) |
| **Impact on helices** | PRO may adopt unrealistic flat angles | PRO maintains proper helical geometry |

**Technical Details:**
The change affects backbone angle potentials (`-BB BB +BB`) when proline is present in helical secondary structures (`cgsecstruct: "H|1|2|3"`). The new restricted bending potential (function type 10) ensures that:
1. Proline's cyclic side chain does not conflict with helical dihedral constraints
2. Unphysical 180° backbone angles are prevented
3. Helical prolines maintain appropriate geometry

**For Publication Reproducibility:**
If your work was originally performed with vermouth 0.13.0 and you need to strictly reproduce those results:
- Create a dedicated conda environment with `vermouth==0.13.0`
- Note that 0.13.0 requires external DSSP (MDTraj integration was added in 0.15.0)

For new studies, we strongly recommend using `vermouth>=0.14.0` (or ideally >=0.15.0) as the proline fix provides more physically accurate behavior, particularly for proteins containing proline-rich helical regions.

## Code Style Guidelines

### Python Style

The project uses modern Python formatting and linting tools configured in `pyproject.toml`:

**Black Configuration**:
- Line length: 100 characters
- Target Python version: 3.12
- Tutorial directory excluded

**Ruff Configuration**:
- Target Python version: 3.12
- Line length: 100
- Enabled rules: pycodestyle (E, W), Pyflakes (F), isort (I), pep8-naming (N), pydocstyle (D), pyupgrade (UP), flake8-bugbear (B), flake8-comprehensions (C4), flake8-simplify (SIM)
- Ignored rules: D100, D104, D105, D107 (docstring requirements)
- Google-style docstrings enforced

### Naming Conventions

**Current Convention**:
- Module names: `snake_case.py`
- Class names: `PascalCase`
- Function names: `snake_case`
- Variable names: `snake_case`

**Key Classes**:
- `TopologyParser` (was `Topology`)
- `MartiniTopFile`
- `Molecule`
- `ForceField`
- `SimulationRunner`

**Key Functions**:
- `create_sb_topology` (was `GenSBPTop`)
- `create_mb_topology` (was `GenMBPTop`)
- `load_config` (was `read_inputs`)
- `write_itp` (was `WriteItp`)

### Type Annotations

All public functions should have full type annotations using Python 3.10+ syntax:

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

### Docstrings

Use Google-style docstrings:

```python
def create_system(
    self,
    nonbonded_cutoff: unit.Quantity = 1.1 * unit.nanometer,
    epsilon_r: float = 15.0,
    remove_com_motion: bool = True,
) -> mm.System:
    """
    Create an OpenMM System from the topology.

    Args:
        nonbonded_cutoff: Cutoff distance for nonbonded interactions.
        epsilon_r: Relative dielectric constant.
        remove_com_motion: Whether to remove center of mass motion.

    Returns:
        Configured OpenMM System object.
    """
```

## Usage Workflow

### 1. Prepare All-Atom Structures
- Clean PDB files (remove waters, trim residues)
- Ensure consistent residue numbering across states

### 2. Generate Contact Maps
- Upload PDB files to [rCSU web-server](http://info.ifpan.edu.pl/~rcsu/rcsu/index.html)
- Download and extract contact map files (.map)

### 3. Generate Multiple-Basin Topology

```bash
ctgomartinize \
    -s state1.pdb state2.pdb \
    -m state1.map state2.map \
    -mol State1Name State2Name \
    -mbmol CombinedName \
    -ff martini3001 \
    -method exp
```

### 4. Set Up Simulation System
- Use `insane.py` to solvate and add ions
- Edit `system.top` to include force field and molecule topology
- Prepare input files (`npt.inp` for equilibration, `md.inp` for production)

### Input File Format (.inp)

```ini
mini_nstep  = 10000         ; Minimization steps
nstep       = 50000000      ; Production steps
dt          = 0.020         ; Time step (ps)
temp        = 310           ; Temperature (K)
platform    = CUDA          ; CPU/CUDA/OpenCL/Reference
pcouple     = yes           ; Pressure coupling
p_type      = isotropic     ; isotropic/membrane/anisotropic
rest        = yes           ; Position restraints
plumed      = no            ; PLUMED support
```

### 5. Run Simulation

```bash
# Minimization and equilibration
run_ctgomartini -i npt.inp

# Production run
run_ctgomartini -i md.inp
```

### 6. Analysis

```bash
# REMD analysis
python -m ctgomartini.analysis.remd_replica_state -f output.nc
```

## Development Conventions

### Adding New Features

1. **Topology functions**: Add to `topology/generator/` or `topology/interactions/`
2. **Simulation features**: Add to `simulation/` directory
3. **CLI commands**: Update `cli/` directory
4. **Utilities**: Add to `utils/` directory
5. **Analysis tools**: Add to `analysis/` directory

### Import Patterns

```python
# Preferred: Import from package level
from ctgomartini.topology import TopologyParser, MartiniTopFile, create_mb_topology
from ctgomartini.simulation import load_config, MDRunner
from ctgomartini.utils import write_itp, convert_map_format

# For internal use
from ctgomartini.topology.interactions import HarmonicBonds, ContactsLJ
```

### Error Handling

- Use assertions for internal consistency checks
- Raise `ValueError` or `Exception` with descriptive messages
- Check file existence before operations

## Security Considerations

- The package uses `subprocess.run` with `shell=True` in some places (`ctgomartinize.py`)
- Input file paths are passed directly to shell commands - sanitize user inputs
- No network operations in the core package

## Known Limitations

1. **Global Parameters**: The `beta`, `C1`, `C2` parameters are global, supporting only one EXP system at a time
2. **Constraint Tolerance**: Default constraint tolerance may need adjustment for some systems
3. **CUDA Compatibility**: OpenMM 8.4 may have compatibility issues with CUDA 12.6

## Troubleshooting

### Common Issues

1. **DSSP version**: Must use dssp < 4.0 or mkdssp (different command-line interface in v4+)
2. **GPU Platform**: CUDA/OpenCL platform availability depends on hardware/drivers
3. **Ion naming**: insane.py may generate incorrect ion names (use sed to fix)

### Debug Tips

- Check input file syntax (semicolon comments, key=value format)
- Verify topology file paths in system.top
- Use `Reference` platform for debugging (slower but more precise)
- Enable PLUMED logging for enhanced sampling diagnostics

## Migration Notes

If migrating from older versions (<1.0.0):

```python
# Old imports (NO LONGER WORK)
from ctgomartini.api import MartiniTopFile, GenMBPTop
from ctgomartini.core import Topology, Molecule
from ctgomartini.util import read_inputs

# New imports
from ctgomartini.topology import TopologyParser, MartiniTopFile, create_mb_topology
from ctgomartini.simulation import load_config
```

See `MIGRATION_GUIDE.md` for complete migration instructions.

## Contact

**Author**: Song Yang  
**Email**: yangsong2015@pku.edu.cn  
**Affiliation**: Computational Biophysics Research Group  
**Repository**: https://github.com/ComputBiophys/CTGoMartini
