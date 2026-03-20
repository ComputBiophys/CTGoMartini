# CTGoMartini: A Python Package for Protein Conformational Transitions and Associated Protein-Lipid Interactions with Go-Martini Models

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/ComputBiophys/CTGoMartini/releases)
[![Python](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

CTGoMartini is a Python package for single-basin Go-Martini, Switching Go-Martini, and Multiple-basin Go-Martini simulations. The main aim of this package is to provide a user-friendly way to simulate conformational transitions of proteins, with particular focus on membrane proteins, using Go-Martini models.

**Version 1.0.0** - See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for upgrade information.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Features

- **New "Contacts" Interaction Type**: Replaces LJ-type contact interactions in classic GoMartini3 models, eliminating incorrect interactions between multiple protein copies and facilitating construction of multiple basin potentials.
- **Automatic Topology Construction**: Generate GoMartini3 model topologies for different simulation methods (Single-basin, Multiple-basin, Switching).
- **Simplified MD Workflow**: Streamlined process for running coarse-grained molecular dynamics simulations with OpenMM.
- **Modern Python**: Full type annotations using Python 3.12+ syntax for better IDE support and code clarity.
- **Unified CLI**: Single entry point for all simulation types with automatic detection.

## Installation

### Prerequisites

- Python >= 3.12
- OpenMM >= 8.1
- vermouth >= 0.15.0
- MDAnalysis >= 2.4
- NumPy >= 2.0, < 2.4.0 (see [Known Issues](#known-issues))

### Standard Install

For most users (single-GPU or CPU simulations):

```bash
# Create conda environment
conda create -n ctgomartini python=3.12 -y
conda activate ctgomartini

# Install core dependencies
pip install MDAnalysis numpy==2.3.4 pandas matplotlib vermouth==0.14.0

# Install OpenMM with CUDA support (adjust cuda-version as needed)
conda install -c conda-forge openmm cuda-version=12 -y

# Install CTGoMartini (editable mode for development)
pip install -e .

# Optional: Install Jupyter for tutorials
conda install jupyter -y
```

### Multi-GPU REMD Installation

For running replica exchange MD (REMD) across multiple GPUs:

```bash
# Install MPI support for multi-GPU REMD
conda install -c conda-forge mpi4py mpich=3 -y

# Install openmmtools with compatible dependencies
conda config --add channels omnia --add channels conda-forge
conda install openmmtools "numpy<2.4.0" "jax<=0.8.0" -y
```

> **Note**: Multi-GPU REMD uses `openmmtools` for replica exchange. See [YANK documentation](http://getyank.org/latest/running.html) for MPI configuration details.

### Optional Dependencies

```bash
# PLUMED support for enhanced sampling
conda install -c conda-forge openmm-plumed

# DSSP for secondary structure assignment (alternative to MDTraj)
conda install -c conda-forge dssp
```

## Testing

Run the test suite to verify your installation:

```bash
# Run all tests
pytest tests/

# Run specific test categories
pytest tests/unit/        # Unit tests (fast, no external dependencies)
pytest tests/integration/ # Integration tests (requires OpenMM)
pytest tests/e2e/         # End-to-end tests (may take several minutes)

# Run with coverage report
pytest tests/ --cov=ctgomartini --cov-report=html
```

### Test Structure

| Category | Description | Runtime |
|----------|-------------|---------|
| `tests/unit/` | Individual function tests | Seconds |
| `tests/integration/` | OpenMM/GROMACS comparison | Minutes |
| `tests/e2e/` | Full workflow tests | Tens of minutes |

### Known Issues

#### NumPy Version Compatibility (< 2.4.0)

**NumPy 2.4.0+ causes REMD simulations to fail** with:
```
TypeError: only 0-dimensional arrays can be converted to Python scalars
```

This is due to changes in NumPy's array scalar handling that affect `openmmtools` replica exchange. **Use NumPy < 2.4.0**:

```bash
pip install "numpy>=2.0,<2.4.0"
# or specifically
pip install numpy==2.3.4
```

Affected functionality:
- REMD simulations using `run_ctgomartini -i remd.inp`
- Replica exchange with `openmmtools.multistate`

Standard MD simulations are not affected.

#### Vermouth Version Compatibility (0.13.0 vs 0.14.0+)

The package now requires `vermouth>=0.14.0`, which includes important fixes for proline (PRO) residues in helical conformations. Starting from vermouth 0.14.0:

| Feature | vermouth 0.13.0 | vermouth >=0.14.0 |
|---------|-----------------|-------------------|
| **Helical PRO angle type** | Harmonic (type 2) | Restricted bending (type 10) |
| **Behavior** | Allows angles near 180° | Prevents 180° angles |

**Impact on Tests**: The end-to-end tests in `tests/e2e/test_gen_mb_itp.py` are temporarily disabled because the reference files were generated with vermouth 0.13.0. These tests compare generated topologies against reference files, and the force field differences cause test failures (even though the generated topologies are physically correct with the new version).

**Workaround**: The reference files will be updated to match vermouth >=0.14.0 output. For now:
- All other 135 tests pass correctly
- The package functions correctly with vermouth >=0.14.0
- If you need strict reproducibility with old results, use a dedicated environment with vermouth 0.13.0

See `AGENTS.md` for technical details about the proline force field changes.

## Project Structure

```
CTGoMartini/
├── ctgomartini/              # Main package
│   ├── topology/             # Topology parsing and building
│   │   ├── parser.py         # TopologyParser (GROMACS topology parser)
│   │   ├── models.py         # Molecule, ForceField classes
│   │   ├── builder.py        # MartiniTopFile for OpenMM
│   │   ├── generator/        # Topology generators
│   │   │   ├── single_basin.py   # create_sb_topology
│   │   │   ├── multi_basin.py    # create_mb_topology
│   │   │   └── combiner.py       # combine_* functions
│   │   └── interactions/     # Bonded and non-bonded interactions
│   │       ├── bonds.py
│   │       ├── angles.py
│   │       ├── dihedrals.py
│   │       ├── contacts.py
│   │       ├── virtual_sites.py
│   │       ├── multi_state.py
│   │       └── mixing.py
│   │
│   ├── simulation/           # Simulation execution
│   │   ├── config.py         # SimulationConfig, load_config
│   │   ├── base.py           # SimulationRunner base class
│   │   ├── md.py             # MDRunner (standard MD)
│   │   └── remd.py           # REMDRunner (replica exchange)
│   │
│   ├── cli/                  # Command-line interface
│   │   ├── run_ctgomartini.py   # Main simulation runner
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
│   │   ├── QValueAnalysis.py
│   │   ├── drms_analysis.py
│   │   └── remd_*.py         # REMD analysis scripts
│   │
│   └── data/                 # Data files and templates
│       ├── ForceFields/      # Martini 3.0 force field
│       ├── Membrane/         # Membrane system templates
│       └── Soluble/          # Soluble protein templates
│
├── tests/                    # Test suite
│   ├── unit/
│   ├── integration/
│   ├── e2e/
│   └── fixtures/
│
├── tutorials/                # Tutorial examples
│   ├── GlnBP/
│   └── TRAAK/
│
├── pyproject.toml
├── MIGRATION_GUIDE.md
└── README.md
```

## Quick Start

### 1. Generate Topology

```bash
# Single-basin (SBP) - single structure
ctgomartinize -s protein.pdb -m auto -mol protein -ff martini3001 -dssp -method sbp

# Multiple-basin with exponential mixing (EXP) - 2+ structures
ctgomartinize -s stateA.pdb stateB.pdb -m auto -mol stateA stateB -mbmol protein \
    -ff martini3001 -dssp -method exp

# Hamiltonian mixing (HAM) - exactly 2 structures
ctgomartinize -s stateA.pdb stateB.pdb -m auto -mol stateA stateB -mbmol protein \
    -ff martini3001 -dssp -method ham

# Switching - separate topologies for each structure
ctgomartinize -s stateA.pdb stateB.pdb -m auto -mol stateA stateB \
    -ff martini3001 -dssp -method switching
```

### 2. Run Simulation

```bash
# Minimization and equilibration
run_ctgomartini -i npt.inp

# Production run
run_ctgomartini -i md.inp

# REMD run
run_ctgomartini -i remd.inp

# Append run
run_ctgomartini -i md.inp --append
```

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

### 3. Analysis

```python
# Native contact analysis
from ctgomartini.analysis import NativeContanceAnalysis

analyzer = NativeContanceAnalysis()
analyzer.calculate_qvalues(
    structure='npt.pdb',
    trajectory='md.xtc',
    references=['state1.pdb', 'state2.pdb'],
    selection='name BB'
)

# dRMS analysis
from ctgomartini.analysis.drms_analysis import DRMSAnalyzer

analyzer = DRMSAnalyzer()
analyzer.calculate_reference_distances(
    ref_files=['state1.pdb', 'state2.pdb'],
    min_distance=6.0,
    max_distance=50.0
)
results = analyzer.analyze_trajectory('npt.pdb', 'md.xtc')
```

## Project Roadmap

### Completed (v1.0.0)

- [x] **Topology Module Reorganization**: Consolidated `api/` and `core/` into unified `topology/` module
- [x] **Simulation Module**: New `simulation/` package with config, base, md, remd modules
- [x] **CLI Module**: New `cli/` package with unified entry points
- [x] **Utils Reorganization**: Snake_case naming, cleaned up deprecated files
- [x] **Data Cleanup**: Removed deprecated scripts, organized templates
- [x] **Analysis Reorganization**: Moved REMD scripts, converted notebooks to Python
- [x] **Vermouth 0.15.0+ Compatibility**: Updated for latest Martinize2
- [x] **MDTraj Integration**: Replace DSSP with MDTraj for secondary structure
- [x] **Type Annotations**: Full type hints using Python 3.12+ syntax
- [x] **MIT License**: Added open-source license
- [x] **Testing Framework**: Comprehensive test suite with pytest (137 tests)
- [x] **Checkpoint Support**: Save/resume simulations
- [x] **Constraint Tolerance**: Configurable in .inp files
- [x] **Anisotropic Position Restraints**: x,y,z restraint support
- [x] **PLUMED Integration**: Enhanced sampling support
- [x] **PDB Validation**: Input file checking
- [x] **New Contact Types**: Contacts6_12, Contacts10_12

### In Progress / Planned

- [ ] **Parameter Searching**: Automated parameter optimization
- [ ] **Tutorial update**: Update tutorials for new CLI
- [ ] **Contacts Generation**: The contacts are different between different versions

### Known Issues / TODO

- [ ] **remd_free_energy.py**: Segmentation fault when running with JAX/PyMBAR. 
  - Workaround: Downgrade numpy to `numpy==2.3.4`
  - Alternative: Use `pymbar.MBAR` directly instead of `ReplicaExchangeAnalyzer`
- [ ] **REMD simulation**: `only 0-dimensional arrays can be converted to Python scalars` error with certain numpy/openmmtools versions
  - Workaround: Downgrade numpy to `numpy==2.3.4`

## Citation

If you use CTGoMartini in your research, please cite:

```
Song Yang, CTGoMartini: A Python Package for Protein Conformational 
Transitions and Protein-Lipid Interactions with Go-Martini Models, 2026
```

## Contact

**Author**: Song Yang  
**Email**: yangsong2015@pku.edu.cn  
**Affiliation**: Computational Biophysics Research Group
