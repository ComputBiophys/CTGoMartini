# CTGoMartini: A Python Package for Protein Conformational Transitions and Associated Protein-Lipid Interactions with Gō-Martini Models

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/ComputBiophys/CTGoMartini/releases)
[![Python](https://img.shields.io/badge/python-3.12+-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

CTGoMartini is a Python package for single-basin Gō-Martini, Switching Gō-Martini, and Multiple-basin Gō-Martini simulations. The main aim of this package is to provide a user-friendly way to simulate conformational transitions of proteins, with particular focus on membrane proteins, using Gō-Martini models.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Features

- **New "Contacts" Interaction Type**: Replaces LJ-type contact interactions in classic Gō-Martini3 models, eliminating incorrect interactions between multiple protein copies and facilitating construction of multiple basin potentials.
- **Automatic Topology Construction**: Generate Gō-Martini3 model topologies for different simulation methods (Single-basin, Multiple-basin, Switching).
- **Simplified MD Workflow**: Streamlined process for running coarse-grained molecular dynamics simulations with OpenMM.
- **Modern Python**: Full type annotations using Python 3.12+ syntax for better IDE support and code clarity.
- **Unified CLI**: Single entry point for all simulation types with automatic detection.

## Installation

### Prerequisites

- Python >= 3.12
- OpenMM >= 8.1
- vermouth == 0.14.0
- MDAnalysis >= 2.4
- NumPy >= 2.0, < 2.4.0

### Install

```bash
# 1. Create conda environment
conda create -n ctgomartini python=3.12 -y
conda activate ctgomartini

# 2. Install core dependencies (note: NumPy < 2.4.0 is required)
pip install MDAnalysis "numpy>=2.0,<2.4.0" pandas matplotlib "vermouth==0.14.0"

# 3. Install OpenMM with CUDA support (adjust cuda-version as needed)
conda install -c conda-forge openmm cuda-version=12 -y

# 4. Install CTGoMartini
pip install -e .

# 5. Install testing framework
pip install pytest pytest-cov

# 6. Install REMD dependencies (optional, only for replica exchange simulations)
conda config --add channels omnia --add channels conda-forge
conda install openmmtools "numpy<2.4.0" "jax<=0.8.0" -y

# 7. Install Jupyter for analysis notebooks (optional)
conda install jupyter -y

# 8. For multi-GPU REMD, also install MPI support (optional)
# conda install -c conda-forge mpi4py mpich=3 -y
```

> **Important Notes**:
> - `numpy<2.4.0` is **mandatory** for REMD - NumPy 2.4.0+ causes `TypeError` with `openmmtools`
> - `cuda-version=12` should match your NVIDIA driver; use `10.2`, `11`, or `12` as needed
> - Multi-GPU REMD uses `openmmtools` for replica exchange. See [YANK documentation](http://getyank.org/latest/running.html) for MPI configuration details.

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

## Quick Start

### 1. Generate Topology

```bash
# Single-basin (SBP) - single structure
ctgomartinize -s protein.pdb -m auto -mol protein -ff martini3001 -dssp -method sbp

# Multiple-basin with exponential mixing (EXP) - 2 structures
ctgomartinize -s stateA.pdb stateB.pdb -m auto -mol stateA stateB -mbmol protein -ff martini3001 -dssp -method exp

# Hamiltonian mixing (HAM) - 2 structures
ctgomartinize -s stateA.pdb stateB.pdb -m auto -mol stateA stateB -mbmol protein -ff martini3001 -dssp -method ham

# Switching - separate topologies for each structure
ctgomartinize -s stateA.pdb stateB.pdb -m auto -mol stateA stateB -ff martini3001 -dssp -method switching
```

### 2. Run Simulation

```bash
# Minimization and equilibration
run_ctgomartini -i npt.inp

# Production run
run_ctgomartini -i md.inp

# Append run
run_ctgomartini -i md.inp --append

# REMD run
run_ctgomartini -i remd.inp
```

## Roadmap

See [ROADMAP.md](ROADMAP.md) for project roadmap and future plans.

## Citation

If you use CTGoMartini in your research, please cite:

```
Song Yang, CTGoMartini: A Python Package for Protein Conformational 
Transitions and Protein-Lipid Interactions with Gō-Martini Models, 2026
```