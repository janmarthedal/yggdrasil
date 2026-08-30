"""Benchmark: L2 error computation (per-element Python loop in yggdrasil/error.py)."""

from __future__ import annotations

import pytest

from benchmarks.problems import _sin_sin, make_poisson_2d_dirichlet
from benchmarks.sizes import SIZES_2D
from yggdrasil import l2_error


@pytest.mark.parametrize("quadrature_order", [2, 5])
@pytest.mark.parametrize("n", SIZES_2D)
def test_l2_error_2d(benchmark, n, quadrature_order):
    problem = make_poisson_2d_dirichlet(n)
    u = problem.linear_solve()
    err = benchmark(l2_error, problem.mesh, u, _sin_sin, quadrature_order)
    assert err < 1e-2
    benchmark.extra_info.update(problem.meta)
