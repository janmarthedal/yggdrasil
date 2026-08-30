"""Shared configuration for the benchmark suite.

Adds a ``--bench-size`` option that selects the mesh-resolution sweep.  The
chosen sizes are written into ``benchmarks.sizes`` before collection so the
parametrized test modules pick them up.
"""

from __future__ import annotations

import pytest

from benchmarks import sizes

_SIZES: dict[str, dict[str, list[int]]] = {
    "small": {"2d": [16, 32, 64], "3d": [4, 8]},
    "medium": {"2d": [32, 64, 128], "3d": [8, 16]},
    "large": {"2d": [64, 128, 256], "3d": [8, 16, 24]},
}


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--bench-size",
        action="store",
        default="small",
        choices=sorted(_SIZES),
        help="Mesh-resolution sweep for the benchmark suite (default: small).",
    )


def pytest_configure(config: pytest.Config) -> None:
    chosen = _SIZES[str(config.getoption("--bench-size"))]
    sizes.SIZES_2D[:] = chosen["2d"]
    sizes.SIZES_3D[:] = chosen["3d"]
