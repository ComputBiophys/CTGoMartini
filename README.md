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

### Quick Install (One Command)

```bash
# Clone the repository
git clone https://github.com/ComputBiophys/CTGoMartini.git
cd CTGoMartini

# Create environment with all dependencies
conda env create -f environment.yml

# Activate and install CTGoMartini
conda activate ctgomartini
pip install -e .

```

### CUDA Version

The default `environment.yml` uses CUDA 12. If you need a different version, edit the file:

```yaml
cuda-version=11  # or 10.2
```

> **Important Notes**:
> - `numpy<2.4.0` is required for REMD compatibility
> - `cuda-version` should match your NVIDIA driver
> - Multi-GPU REMD uses `openmmtools` for replica exchange. See [YANK documentation](http://getyank.org/latest/running.html) for MPI configuration details.

## Testing

Run the test suite to verify your installation:

```bash
# Run all tests
pytest tests/
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
Yang, Song; Song, Chen*. CTGoMartini: A Python Package for Protein Conformational 
Transitions and Protein-Lipid Interactions with Gō-Martini Models, 2026.
```
