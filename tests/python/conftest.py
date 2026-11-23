"""Shared set-up for the unit tests: import the scripts in scripts/ as modules, and locate the test data."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture
def test_bundle():
    """The synthetic bundle in release format (tests/data/make_test_bundle.py): genes G01-G06 changed."""
    return ROOT / "tests" / "data" / "bundle"
