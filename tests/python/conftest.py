"""Shared set-up for the unit tests: import the scripts as modules, and locate the test data.

Scripts run by the pipeline live with their module (workflow/nextflow/modules/<module>/scripts/); project scripts
(downloads, metadata) live in scripts/. Both are put on the import path.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
for folder in [ROOT / "scripts", *sorted((ROOT / "workflow/nextflow/modules").glob("*/scripts"))]:
    sys.path.insert(0, str(folder))


@pytest.fixture
def test_bundle():
    """The synthetic bundle in release format (tests/data/make_test_bundle.py): genes G01-G06 changed."""
    return ROOT / "tests" / "data" / "bundle"
