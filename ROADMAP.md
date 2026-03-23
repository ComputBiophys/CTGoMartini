# Project Roadmap

## Completed (v1.0.0)

- [x] **Topology Module Reorganization**: Consolidated `api/` and `core/` into unified `topology/` module
- [x] **Simulation Module**: New `simulation/` package with config, base, md, remd modules
- [x] **CLI Module**: New `cli/` package with unified entry points
- [x] **Utils Reorganization**: Snake_case naming, cleaned up deprecated files
- [x] **Data Cleanup**: Removed deprecated scripts, organized templates
- [x] **Analysis Reorganization**: Moved REMD scripts, converted notebooks to Python
- [x] **Vermouth 0.14.0+ Compatibility**: Updated for latest Martinize2
- [x] **MDTraj Integration**: Replace DSSP with MDTraj for secondary structure
- [x] **Type Annotations**: Full type hints using Python 3.12+ syntax
- [x] **MIT License**: Added open-source license
- [x] **Testing Framework**: Comprehensive test suite with pytest (137 tests)
- [x] **Checkpoint Support**: Save/resume simulations
- [x] **Constraint Tolerance**: Configurable in .inp files
- [x] **Anisotropic Position Restraints**: x,y,z restraint support
- [x] **PDB Validation**: Input structure compatibility checking with detailed error reporting
- [x] **New Contact Types**: Contacts6_12, Contacts10_12
- [x] **Parameter Searching**: Automated parameter optimization via FESAnalyzer.parameter_sweep()
- [x] **Tutorial update**: Update tutorials for new CLI
- [x] **Contacts Generation**: The contacts are different between different versions
- [x] **Analysis Module**: Comprehensive analysis tools including dRMS analysis, MBAR free energy calculation, REMD exchange ratio, and trajectory extraction
- [x] **Multi-EXP System Support**: Fixed global parameters bug (beta, C1, C2) by hardcoding values into energy expressions
- [x] **Command-line Fix**: Correct handling of dash-prefixed parameters in -other_params (e.g., "-nt")

## In Progress / Planned
