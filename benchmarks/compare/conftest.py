"""Fixtures for the scikit-fem comparison benchmarks.

The whole ``benchmarks/compare/`` tree is skipped unless scikit-fem is
installed (``uv run --group bench-compare pytest benchmarks/``).
"""

from __future__ import annotations

import pytest


@pytest.fixture(scope="session")
def skfem():
    return pytest.importorskip("skfem")
