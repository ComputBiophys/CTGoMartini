# Changelog

All notable changes to CTGoMartini will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-03-15

### Breaking Changes

- **Vermouth Version**: Updated dependency to `vermouth>=0.10.0,<0.14.0` (was ==0.9.6)
  - Note: vermouth >= 0.14.0 modified Proline protein parameters; code supports it but tests use older reference data
  - martinize2 CLI changes: `-dssp <cmd>` → `-dssp` (uses MDTraj by default)
  - Side-chain fix (`-scfix`) is now **enabled by default** (use `-noscfix` to disable)
  - Output file naming changed: `{name}.itp` → `{name}_0.itp` (handled internally)
  - Virtual sites format changed in ITP output

- **Go Contacts Generation**: Replaced OV+rCSU workflow with martinize2 native `-go` support
  - Old: `create_go_virt_for_multimer` → `go-table_VirtGoSites.itp`
  - New: `martinize2 -go contact_map.out -go-low 0.3 -go-up 1.1 -go-eps 12.0`
  - Generates `go_atomtypes.itp` and `go_nbparams.itp` automatically

### Added

- **Command Line Tool**: `ctgomartinize` command available after installation
  - Entry point added to pyproject.toml
  - Can be invoked directly: `ctgomartinize -s protein.pdb ...`

- **contacts Module**: New `ctgomartini/utils/contacts.py` utility
  - `gen_go_contacts()`: Unified interface for contact generation (auto, .map, .out)
  - `run_ovrcsu()`: Automatic contact generation via OVrCSU
  - `detect_contact_file_format()`: File format detection by extension

- **convert_map_format Module**: New utility to convert rCSU web-server `.map` format to martinize2-compatible `contact_map.out`
  - Handles format differences (header, columns, Count/Model fields)
  - CLI and Python API support

- **E2E Tests for Vermouth 0.15.0**: Comprehensive tests for new Go contacts workflow
  - Validates contacts conversion from `go_nbparams.itp` to `[contacts]` in final topology
  - Verifies virtual site ID mapping (512 offset)
  - Checks sigma value preservation

### Changed

- **ctgomartinize.py**: Major refactoring for improved maintainability
  - Reduced code duplication via `_process_single_state()` helper function
  - Improved directory handling using `pathlib.Path` and proper context management
  - Enhanced CLI output with structured progress indicators and step markers
  - Removed shell=True from subprocess calls for better security
  - Added command-line entry point: `ctgomartinize`
  - Updated help text with `RawDescriptionHelpFormatter` for better formatting
  - New parameters: `-go-eps`, `-go-low`, `-go-up` for Go contacts control
  - `-m` parameter now supports 'auto' mode (default) for automatic contact generation

### Migration Notes

See [MIGRATION_GUIDE.md](MIGRATION_GUIDE.md) for detailed upgrade instructions from 0.6.0 to 1.0.0.

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
