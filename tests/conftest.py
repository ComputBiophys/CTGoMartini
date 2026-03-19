"""
Pytest configuration and shared fixtures for CTGoMartini tests.

This module provides common utilities and fixtures used across
the test suite, including safe directory change context managers.
"""

import os


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
