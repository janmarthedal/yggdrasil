"""Benchmark: structured mesh generation."""

from __future__ import annotations

import pytest

from benchmarks.sizes import SIZES_2D, SIZES_3D
from yggdrasil import unit_cube_tet_mesh, unit_square_tri_mesh


@pytest.mark.parametrize("n", SIZES_2D)
def test_unit_square_tri_mesh(benchmark, n):
    mesh = benchmark(unit_square_tri_mesh, n)
    assert mesh.num_elements == 2 * n * n
    benchmark.extra_info.update(num_nodes=mesh.num_nodes, num_elements=mesh.num_elements)


@pytest.mark.parametrize("n", SIZES_3D)
def test_unit_cube_tet_mesh(benchmark, n):
    mesh = benchmark(unit_cube_tet_mesh, n)
    assert mesh.num_elements == 6 * n**3
    benchmark.extra_info.update(num_nodes=mesh.num_nodes, num_elements=mesh.num_elements)
