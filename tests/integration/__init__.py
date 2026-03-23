"""
Integration tests for CTGoMartini.

These tests verify component interactions and require external dependencies:
    - OpenMM >= 8.1 for energy calculations
    - Pre-calculated reference data in fixtures/

Tested workflows:
    - Classic Martini topology parsing
    - Contact-based interactions
    - Multiple-basin potential (MBP) calculations
    - Energy and force comparisons (OpenMM vs GROMACS)
"""
