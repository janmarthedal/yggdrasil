"""Mesh-resolution sweep for the benchmark suite.

The lists are mutated in place by ``conftest.pytest_configure`` according to the
``--bench-size`` option, so test modules must import the list objects (not copy
them) and reference them at collection time.
"""

SIZES_2D: list[int] = [16, 32, 64]
SIZES_3D: list[int] = [4, 8]
