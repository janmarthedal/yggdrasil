"""Benchmark: global matrix / vector assembly.

Assembly is a pure-Python per-element loop today (``yggdrasil/assemble.py``) and
is expected to dominate runtime at scale, so this is the primary optimization
target.
"""

from __future__ import annotations

import pytest

from benchmarks.problems import (
    make_laplace_eigenvalue_2d,
    make_poisson_2d_dirichlet,
    make_poisson_3d_dirichlet,
)
from benchmarks.sizes import SIZES_2D, SIZES_3D


def _record(benchmark, problem):
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("quadrature_order", [1, 2, 3])
@pytest.mark.parametrize("n", SIZES_2D)
def test_stiffness_assembly_2d(benchmark, n, quadrature_order):
    problem = make_poisson_2d_dirichlet(n, quadrature_order=quadrature_order)
    K = benchmark(problem.assemble_lhs)
    assert K.shape == ((n + 1) ** 2, (n + 1) ** 2)
    _record(benchmark, problem)


@pytest.mark.parametrize("n", SIZES_2D)
def test_load_vector_assembly_2d(benchmark, n):
    problem = make_poisson_2d_dirichlet(n)
    b = benchmark(problem.assemble_rhs)
    assert b.shape == ((n + 1) ** 2,)
    _record(benchmark, problem)


@pytest.mark.parametrize("n", SIZES_2D)
def test_mass_assembly_2d(benchmark, n):
    problem = make_laplace_eigenvalue_2d(n)
    K, M = benchmark(problem.assemble_lhs)
    assert M.shape == K.shape
    _record(benchmark, problem)


@pytest.mark.parametrize("n", SIZES_3D)
def test_stiffness_assembly_3d(benchmark, n):
    problem = make_poisson_3d_dirichlet(n)
    K = benchmark(problem.assemble_lhs)
    assert K.shape == ((n + 1) ** 3, (n + 1) ** 3)
    _record(benchmark, problem)
