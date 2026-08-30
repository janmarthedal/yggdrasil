"""Benchmark: linear solve and full end-to-end PDE solves."""

from __future__ import annotations

import numpy as np
import pytest

from benchmarks.problems import (
    make_laplace_eigenvalue_2d,
    make_poisson_2d_dirichlet,
    make_poisson_2d_mixed_bc,
    make_poisson_3d_dirichlet,
)
from benchmarks.sizes import SIZES_2D, SIZES_3D


@pytest.mark.parametrize("n", SIZES_2D)
def test_linear_solve_2d(benchmark, n):
    problem = make_poisson_2d_dirichlet(n)
    u = benchmark(problem.linear_solve)
    assert problem.compute_error(u) < 1e-2
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_full_poisson_2d_dirichlet(benchmark, n):
    problem = make_poisson_2d_dirichlet(n)
    u = benchmark(problem.run_full)
    assert problem.compute_error(u) < 1e-2
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_full_poisson_2d_mixed_bc(benchmark, n):
    problem = make_poisson_2d_mixed_bc(n)
    u = benchmark(problem.run_full)
    assert problem.compute_error(u) < 1e-10
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_3D)
def test_full_poisson_3d_dirichlet(benchmark, n):
    problem = make_poisson_3d_dirichlet(n)
    u = benchmark(problem.run_full)
    assert 0.04 < u.max() < 0.07
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_full_laplace_eigenvalue_2d(benchmark, n):
    problem = make_laplace_eigenvalue_2d(n)
    eigenvalues = benchmark(problem.run_full)
    np.testing.assert_allclose(eigenvalues[0], 2 * np.pi**2, rtol=5e-2)
    benchmark.extra_info.update(problem.meta)
