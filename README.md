# CTGoMartini: A Python Package for Protein Conformational Transitions with Go-Martini Models

[![Version](https://img.shields.io/badge/version-0.6.0-blue.svg)](https://github.com/ComputBiophys/CTGoMartini/releases)
[![Python](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

CTGoMartini is a Python package for single-basin Go-Martini, Switching Go-Martini, and Multiple-basin Go-Martini simulations. The main aim of this package is to provide a user-friendly way to simulate conformational transitions of proteins, with particular focus on membrane proteins, using Go-Martini models.

**Version 0.6.0** - See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for upgrade information from 0.5.x.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Features

- **New "Contacts" Interaction Type**: Replaces LJ-type contact interactions in classic GoMartini3 models, eliminating incorrect interactions between multiple protein copies and facilitating construction of multiple basin potentials.
- **Automatic Topology Construction**: Generate GoMartini3 model topologies for different simulation methods (Single-basin, Multiple-basin, Switching).
- **Simplified MD Workflow**: Streamlined process for running coarse-grained molecular dynamics simulations with OpenMM.
- **Modern Python**: Full type annotations using Python 3.12+ syntax for better IDE support and code clarity.

## Installation

### Prerequisites

- Python >= 3.12
- OpenMM >= 8.1
- Vermouth == 0.9.6
- MDAnalysis >= 2.4

### Quick Install

```bash
# Create conda environment
conda create -n ctgomartini python=3.12
conda activate ctgomartini

# Install dependencies
conda install -c conda-forge openmm

# Install CTGoMartini
pip install -e .

# Optional: Install PLUMED support for enhanced sampling
conda install -c conda-forge openmm-plumed

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

## Project Structure

```
CTGoMartini/
├── ctgomartini/       # Main package
│   ├── api/           # High-level API classes
│   │   ├── MartiniTopology.py    # MartiniTopFile for topology handling
│   │   ├── MBMoleculeTop.py      # GenMBPTop for multiple-basin topologies
│   │   └── SBMoleculeTop.py      # GenSBPTop for single-basin topologies
│   ├── core/          # Core functionality modules (formerly func/)
│   │   ├── Topology.py           # GROMACS topology parser
│   │   ├── Molecule.py           # Molecule representation
│   │   ├── ForceField.py         # Force field parameter handling
│   │   ├── Bonded_interaction.py # Bonded force terms
│   │   ├── Nonbonded_interaction.py  # Non-bonded force terms
│   │   ├── CombineMols.py        # Molecule combination utilities
│   │   └── vsites.py             # Virtual site handling
│   ├── utils/         # Utility functions (formerly util/)
│   │   ├── ReadInp.py            # Input file parser
│   │   ├── WriteItp.py           # ITP file generation
│   │   ├── ConvertLongShortElasticBonds.py
│   │   └── Create_goVirt_for_multimer.py
│   ├── data/          # Data files and executable scripts
│   │   ├── run_ctgomartini.py    # Main simulation script
│   │   ├── ctgomartinize.py      # Topology generation script
│   │   ├── run_REMD.py           # Replica exchange MD script
│   │   └── ForceFields/          # Martini 3.0 force field files
│   └── analysis/      # Analysis tools
│       └── QValueAnalysis.py     # Native contact analysis
│
├── tests/             # Test suite
│   ├── unit/          # Unit tests for individual functions
│   ├── integration/   # Integration tests (OpenMM/GROMACS comparison)
│   ├── e2e/           # End-to-end workflow tests
│   └── fixtures/      # Test data and reference files
│
└── Tutorial/          # Tutorial examples
    ├── GlnBP/         # Glutamine Binding Protein example
    └── TRAAK/         # TRAAK channel example
```


### 2. Run Simulation

```bash
python ctgomartini/data/run_ctgomartini.py -i npt.inp
```

See the `Tutorial/` directory for complete examples.

## Project Roadmap

### Completed (v0.6.0)

- [x] **Force Class Refinement**: Transitioned to objective classification criteria
- [x] **Multi-basin Force Groups Overhaul**: Efficient CV-based approach
- [x] **Code Reconstruction**: Refactored with AI tools for improved readability
- [x] **Type Annotations**: Full type hints using Python 3.12+ syntax
- [x] **Modern Package Structure**: Migrated to `pyproject.toml`
- [x] **MIT License**: Added open-source license
- [x] **Module Renaming**: `func/` → `core/`, `util/` → `utils/` (see [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md))
- [x] **Testing Framework**: Comprehensive test suite with pytest (50 tests)
- [x] **Test Organization**: Tests organized into unit/integration/e2e categories
- [x] **Charge and Atom Number Checking**: Validation added
- [x] **XTC Reporter**: Compressed trajectory output
- [x] **Checkpoint Support**: Save/resume simulations
- [x] **Constraint Tolerance**: Configurable in .inp files
- [x] **Anisotropic Position Restraints**: x,y,z restraint support
- [x] **PLUMED Integration**: Enhanced sampling support
- [x] **Single-basin Contact Topology**: SBP generation
- [x] **Tutorials**: Examples for membrane and soluble proteins
- [x] **Bug Fixes**: Local_BondedInteraction_dict import, MBP energy expression, race conditions

### In Progress / Planned

#### v1.0.0 (In Development)
- [x] **Vermouth 0.15.0+ Compatibility**: Update for latest Vermouth/Martinize2
- [x] **MDTraj Integration**: Replace DSSP with MDTraj for secondary structure
- [x] **Contact Map Options**: Dual mode (OV+rCSU / martinize2 internal)

#### Future
- [ ] **Analysis Module**: Enhanced analysis tools
- [x] **Create_goVirt_for_multimer**: Fix CA/BB extraction bugs： Deprecated
- [x] **Global Parameters**: Support multiple EXP systems (beta, C1, C2)
- [ ] **Parameter Searching**: Automated parameter optimization
- [x] **ctgomartinize.py**: Fix `-other_params "-nt"` bug
- [x] **Contact Options**: Parameter to toggle contacts usage
- [x] **PDB Validation**: Input file checking
- [x] **Constraints to bonds**: Add flag to control
- [x] **Add new contact types**
- [ ] **Convert contact interactions into new types**
- [ ] **Tutorial update**

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
