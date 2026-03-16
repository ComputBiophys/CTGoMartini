# CTGoMartini Example

This directory contains examples demonstrating CTGoMartini usage.

## Examples

### 1. Single-Basin (SB)
- **Directory**: `SB/`
- **Description**: Single-basin Go-Martini model for one conformation
- **Files**: `HAMP_Bicyclomycin.pdb`, `HAMP_Bicyclomycin.map`
- **Command**: `python generate_topology.py`

### 2. Switching
- **Directory**: `Switching/`
- **Description**: Switching potential for conformational transition
- **Files**: `State1.pdb`, `State1.map`, `State2.pdb`, `State2.map`
- **Command**: `python generate_topology.py`

### 3. MB-EXP (Multiple-Basin Exponential)
- **Directory**: `MB-EXP/`
- **Description**: Multiple-basin with exponential potential
- **Files**: `State1.pdb`, `State1.map`, `State2.pdb`, `State2.map`
- **Command**: `python generate_topology.py`

### 4. MB-HAM (Multiple-Basin Harmonic)
- **Directory**: `MB-HAM/`
- **Description**: Multiple-basin with harmonic approximation
- **Files**: `State1.pdb`, `State1.map`, `State2.pdb`, `State2.map`
- **Command**: `python generate_topology.py`

## Quick Test

Run all tests:
```bash
python run_all_tests.py
```

## Shared Files

- `State1.pdb`, `State2.pdb` - Sample protein structures
- `State1.map`, `State2.map` - Contact map files
- `martini_v3.0.0.itp` - Martini 3.0 force field

## Output Files

- `Protein.itp` / `Protein_SB.itp` / `Protein_MB.itp` - Generated topology
- `*_params.itp` - Force field parameters
