# Changelog

All notable changes to CTGoMartini will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.0] - 2026-03-13

### Breaking Changes

- **Module Restructure**: Renamed directories for consistency
  - `ctgomartini/func/` → `ctgomartini/core/`
  - `ctgomartini/util/` → `ctgomartini/utils/`
  - See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for upgrade instructions

- **Function Renaming**: Snake_case for Python naming conventions
  - `WriteItp` → `write_itp`
  - `ConvertLongShortElasticBonds` → `convert_long_short_elastic_bonds`
  - `Create_goVirt_for_multimer` → `create_go_virt_for_multimer`
  - `gen_restraints` → `generate_restraints`
  - `restraints` → `add_restraints`

### Added

- **Modern Package Structure**: Migrated to `pyproject.toml`
- **MIT License**: Added open-source license
- **Comprehensive Testing**: 50 tests organized into unit/integration/e2e categories
- **Type Annotations**: Full type hints using Python 3.12+ syntax
- **Test Infrastructure**: 
  - `WorkingDirectoryContext` for safe directory changes
  - Proper test isolation to prevent side effects
- **Documentation**:
  - `MIGRATION_GUIDE.md` for upgrade instructions
  - `AGENTS.md` for developer guidance
  - `MEMORY.md` for project context

### Fixed

- **Import Error**: `Local_BondedInteraction_dict` now correctly imported from `Bonded_interaction`
- **MBP Energy Expression**: Fixed empty energy expressions in multiple-basin potential
- **Race Condition**: Fixed parallel execution issues in `test_MBGoMartini.py`
- **Array Type Error**: Fixed `Compare_energy` to handle various array types properly

### Changed

- **Test Organization**: Moved from `CTGoMartini-tests/` to `tests/` with proper structure
- **Code Quality**: Improved error handling and type safety throughout

## [0.5.0] - Previous Version

- Initial stable release
- Single-basin and Multiple-basin Go-Martini support
- OpenMM integration
- PLUMED support for enhanced sampling

---

[0.6.0]: https://github.com/ComputBiophys/CTGoMartini/releases/tag/v0.6.0
