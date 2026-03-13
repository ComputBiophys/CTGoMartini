"""
End-to-end tests for CTGoMartini.

Full workflow tests that execute complete simulation pipelines.
These tests may take several minutes to complete.

Tested workflows:
    - Multi-basin ITP generation (ctgomartinize.py)
    - Molecular dynamics simulation runs (run_ctgomartini.py)

Requirements:
    - OpenMM >= 8.1
    - DSSP or mkdssp for secondary structure assignment
    - Pre-configured simulation templates in fixtures/
"""
