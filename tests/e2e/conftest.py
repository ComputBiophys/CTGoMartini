"""Pytest configuration for end-to-end tests.

Marks every test collected from this directory with the `e2e` marker so the
default run (`addopts` sets `-m "not e2e"`) skips them. Run them explicitly
with `pytest -m e2e`.
"""

from pathlib import Path

import pytest

_E2E_DIR = Path(__file__).parent.resolve()


def pytest_collection_modifyitems(items):
    """Apply the `e2e` marker to every test collected from this directory."""
    for item in items:
        if _E2E_DIR in Path(item.path).resolve().parents:
            item.add_marker(pytest.mark.e2e)
