"""
Pytest configuration and shared fixtures for CTGoMartini tests.

This module provides common utilities and fixtures used across
the test suite, including safe directory change context managers.
"""

import os
from pathlib import Path

import pytest


class WorkingDirectoryContext:
    """
    Context manager to safely change and restore working directory.

    Ensures that the original working directory is restored even if
    an exception occurs within the context.

    Args:
        new_dir (str): Directory to change into.

    Example:
        >>> with WorkingDirectoryContext('/tmp'):
        ...     # Work in /tmp
        ...     pass
        >>> # Original directory is restored
    """

    def __init__(self, new_dir: str):
        """Initialize with target directory."""
        self.new_dir = new_dir
        self.original_dir = None

    def __enter__(self):
        """Save current directory and change to new directory."""
        try:
            self.original_dir = os.getcwd()
        except FileNotFoundError:
            # Current working directory was deleted (e.g., by previous test's
            # tempfile.TemporaryDirectory cleanup). Fall back to a safe directory.
            self.original_dir = os.path.expanduser("~")

        # Ensure target directory exists
        os.makedirs(self.new_dir, exist_ok=True)
        os.chdir(self.new_dir)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Restore original directory regardless of exceptions."""
        # Only restore if original_dir was set and still exists
        if self.original_dir is not None:
            try:
                os.chdir(self.original_dir)
            except FileNotFoundError:
                # Original directory no longer exists, fall back to home
                os.chdir(os.path.expanduser("~"))
        return False  # Don't suppress exceptions


# =============================================================================
# REMD Test Data Fixtures
# =============================================================================

@pytest.fixture(scope="session")
def remd_test_data_dir() -> Path:
    """
    Return the path to REMD test data directory.
    
    Returns:
        Path to tests/fixtures/remd/ directory containing:
        - output.nc: Analysis NetCDF file
        - output_checkpoint.nc: Checkpoint NetCDF with positions
        - npt.pdb: System topology
        - closed_cg.pdb: Closed state reference structure
        - open_cg.pdb: Open state reference structure
    """
    return Path(__file__).parent / "fixtures" / "remd"


@pytest.fixture(scope="session")
def remd_netcdf_file(remd_test_data_dir: Path) -> Path:
    """Path to REMD analysis NetCDF file."""
    return remd_test_data_dir / "output.nc"


@pytest.fixture(scope="session")
def remd_checkpoint_file(remd_test_data_dir: Path) -> Path:
    """Path to REMD checkpoint NetCDF file."""
    return remd_test_data_dir / "output_checkpoint.nc"


@pytest.fixture(scope="session")
def remd_topology_file(remd_test_data_dir: Path) -> Path:
    """Path to REMD system topology (PDB)."""
    return remd_test_data_dir / "npt.pdb"


@pytest.fixture(scope="session")
def remd_closed_ref(remd_test_data_dir: Path) -> Path:
    """Path to closed state reference structure."""
    return remd_test_data_dir / "closed_cg.pdb"


@pytest.fixture(scope="session")
def remd_open_ref(remd_test_data_dir: Path) -> Path:
    """Path to open state reference structure."""
    return remd_test_data_dir / "open_cg.pdb"
