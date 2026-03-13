# CTGoMartini Agent Guide

## Project Overview

CTGoMartini is a Python package for simulating protein conformational transitions using Go-Martini coarse-grained models. It provides tools for setting up and running molecular dynamics simulations with single-basin, switching, and multiple-basin Go-Martini potentials, with particular focus on membrane proteins.

**Version**: 0.6.0  
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
- **vermouth** (==0.9.6): Martinize2 integration for coarse-graining (version locked due to API instability)
- **NumPy**: Numerical computations
- **Pandas**: Data manipulation
- **Matplotlib** (>=3.3.3): Plotting and visualization

### Optional Dependencies
- **openmm-plumed** (>=2.0.0): PLUMED integration for enhanced sampling
- **dssp** (<4.0): Secondary structure assignment (vermouth 0.9.6 only compatible with DSSP <4.0)

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
│   ├── _version.py           # Version information (0.5.0)
│   ├── api/                  # High-level API classes
│   │   ├── MartiniTopology.py    # MartiniTopFile class for topology handling
│   │   ├── MBMoleculeTop.py      # GenMBPTop for multiple-basin topologies
│   │   └── SBMoleculeTop.py      # GenSBPTop for single-basin topologies
│   ├── core/                 # Core functionality modules (formerly func/)
│   │   ├── Topology.py           # GROMACS topology parser
│   │   ├── Molecule.py           # Molecule representation
│   │   ├── ForceField.py         # Force field parameter handling
│   │   ├── Bonded_interaction.py # Bonded force terms
│   │   ├── Nonbonded_interaction.py  # Non-bonded force terms
│   │   ├── CombineMols.py        # Molecule combination utilities
│   │   └── vsites.py             # Virtual site handling
│   ├── utils/                # Utility functions (formerly util/)
│   │   ├── WriteItp.py           # ITP file generation
│   │   ├── ReadInp.py            # Input file parser
│   │   ├── ConvertLongShortElasticBonds.py
│   │   └── Create_goVirt_for_multimer.py
│   ├── analysis/             # Analysis tools
│   │   ├── QValueAnalysis.py     # Native contact (Q-value) analysis
│   │   └── MBAnalysis_v*.py      # Multiple-basin analysis scripts
│   └── data/                 # Data files and executable scripts
│       ├── run_ctgomartini.py    # Main simulation script
│       ├── ctgomartinize.py      # Topology generation script
│       ├── run_REMD.py           # Replica exchange MD script
│       ├── run_REMD_restart.py   # REMD restart script
│       ├── OVrCSU.py             # Contact map generation (OV+rCSU method)
│       ├── ForceFields/Martini300/   # Martini 3.0 force field
│       ├── Membrane/             # Membrane system input templates
│       └── Soluble/              # Soluble protein input templates
├── Tutorial/                 # Tutorial examples
│   ├── GlnBP/                # GlnBP conformational transition example
│   └── TRAAK/                # TRAAK channel example
├── pyproject.toml            # Package configuration (modern Python packaging)
├── setup.py                  # Minimal setup.py (backward compatibility only)
├── tests/                    # Test suite (unit/integration/e2e)
├── MIGRATION_GUIDE.md        # API migration guide (func/→core/, util/→utils/)
├── REFACTORING_STATUS.md     # Refactoring progress and known issues
└── README.md                 # Project documentation
```

### Module Organization

| Module | Purpose | Key Classes/Functions |
|--------|---------|----------------------|
| `api/` | High-level user API | `MartiniTopFile`, `GenMBPTop`, `GenSBPTop` |
| `core/` | Core algorithms and data structures | `Topology`, `Molecule`, `ForceField`, interaction classes |
| `utils/` | I/O utilities and helpers | `read_inputs`, `write_itp`, `convert_long_short_elastic_bonds` |
| `data/` | Executable scripts and templates | `run_ctgomartini.py`, `ctgomartinize.py` |
| `analysis/` | Post-simulation analysis | `NativeContanceAnalysis`, `MBAnalysis_v3_argparse.py` |

## Build and Installation

### Installation Commands

```bash
# Create conda environment
conda create -n ctgomartini python=3.12
conda activate ctgomartini

# Install core dependencies
conda install -c conda-forge openmm

# Install CTGoMartini (editable mode for development)
pip install -e .

# Install with optional dependencies
pip install -e ".[plumed,dssp]"

# Or install development dependencies
pip install -e ".[dev]"
```

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
# Run unit tests for ITP writing
pytest tests/

# Run pytest (if tests/ directory exists)
pytest

# Run with coverage
pytest --cov=ctgomartini
```

### Test Structure

1. **In-repo tests**: `tests/unit/` - Unit tests for ITP writing and other utilities
2. **Test suite**: `tests/` directory with unit, integration, and e2e tests

### Known Test Issues

From `REFACTORING_STATUS.md`:
- 3 out of 40 external tests fail due to incomplete test data (not code issues)
- GlnBP and Beta2AR test topologies missing `multi_angles`/`multi_contacts` sections
- MBP core functionality verified through `test_Beta2AR_HAM`

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

**MyPy Configuration**:
- Python version: 3.12
- `warn_return_any = true`
- `ignore_missing_imports = true` (for external libraries)

### Naming Conventions

**Current Convention** (post-refactoring):
- Module names: `snake_case.py`
- Class names: `PascalCase`
- Function names: `snake_case`
- Variable names: `snake_case`

**Note**: Some API functions retain PascalCase for backward compatibility:
- `GenMBPTop` (function)
- `GenSBPTop` (function)
- `MartiniTopFile` (class)

**Migration from old naming**:
| Old (camelCase) | New (snake_case) |
|-----------------|------------------|
| `WriteItp` | `write_itp` |
| `ConvertLongShortElasticBonds` | `convert_long_short_elastic_bonds` |
| `LoadStructure` | `load_structure` |
| `gen_restraints` | `generate_restraints` |

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
python ctgomartini/data/ctgomartinize.py \
    -s state1.pdb state2.pdb \
    -m state1.map state2.map \
    -mol State1Name State2Name \
    -mbmol CombinedName \
    -dssp path/to/dssp \
    -ff martini3001 \
    -method exp
```

### 4. Set Up Simulation System
- Use `insane.py` to solvate and add ions
- Edit `system.top` to include force field and molecule topology
- Prepare input files (`npt.inp` for equilibration, `md.inp` for production)

### Input File Format (.inp)

GROMACS-style parameter files:

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
python ctgomartini/data/run_ctgomartini.py -i npt.inp

# Production run
python ctgomartini/data/run_ctgomartini.py -i md.inp
```

### 6. Analysis

```bash
# Native contact analysis
python MBAnalysis_v3_argparse.py \
    -s npt.pdb -f md.xtc \
    -r state1.pdb state2.pdb \
    -sel "name BB" \
    -prefix output_prefix
```

## Development Conventions

### Adding New Features

1. API functions go in `api/` directory
2. Core algorithms go in `core/` directory
3. Utilities go in `utils/` directory
4. Analysis tools go in `analysis/` directory
5. Executable scripts go in `data/` directory

### Import Patterns

```python
# Preferred: Import from package level
from ctgomartini.core import Topology, Molecule
from ctgomartini.utils import read_inputs, write_itp
from ctgomartini.api import MartiniTopFile

# For data scripts
from ctgomartini.api import MartiniTopFile
from ctgomartini.utils import read_inputs
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
2. **Vermouth Compatibility**: Locked to vermouth == 0.9.6 due to API instability
3. **DSSP Compatibility**: vermouth 0.9.6 only supports DSSP <4.0 (use `mkdssp` instead of `dssp` v4+)
4. **Constraint Tolerance**: Default constraint tolerance may need adjustment for some systems
5. **CUDA Compatibility**: OpenMM 8.4 may have compatibility issues with CUDA 12.6

## Troubleshooting

### Common Issues

1. **DSSP version**: Must use dssp < 4.0 or mkdssp (different command-line interface in v4+)
2. **GPU Platform**: CUDA/OpenCL platform availability depends on hardware/drivers
3. **Ion naming**: insane.py may generate incorrect ion names (use sed to fix)
4. **Virtual sites**: Multimer systems may have CA/BB extraction issues

### Debug Tips

- Check input file syntax (semicolon comments, key=value format)
- Verify topology file paths in system.top
- Use `Reference` platform for debugging (slower but more precise)
- Enable PLUMED logging for enhanced sampling diagnostics

## Migration Notes

If migrating from older versions (<0.6.0):

```python
# Old imports (NO LONGER WORK)
from ctgomartini.func import Topology, Molecule
from ctgomartini.util import read_inputs

# New imports
from ctgomartini.core import Topology, Molecule
from ctgomartini.utils import read_inputs
```

See `MIGRATION_GUIDE.md` for complete migration instructions.

## Contact

**Author**: Song Yang  
**Email**: yangsong2015@pku.edu.cn  
**Affiliation**: Computational Biophysics Research Group  
**Repository**: https://github.com/ComputBiophys/CTGoMartini
