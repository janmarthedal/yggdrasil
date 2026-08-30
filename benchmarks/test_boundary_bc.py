"""Benchmark: boundary extraction and boundary-condition application."""

from __future__ import annotations

import pytest

from benchmarks.problems import (
    make_mass_projection_2d,
    make_poisson_2d_dirichlet,
    make_poisson_2d_mixed_bc,
    make_poisson_3d_dirichlet,
)
from benchmarks.sizes import SIZES_2D, SIZES_3D


@pytest.mark.parametrize("n", SIZES_2D)
def test_extract_boundary_2d(benchmark, n):
    problem = make_poisson_2d_dirichlet(n)
    bnd = benchmark(problem.extract_bnd)
    assert bnd.num_elements == 4 * n
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_3D)
def test_extract_boundary_3d(benchmark, n):
    problem = make_poisson_3d_dirichlet(n)
    bnd = benchmark(problem.extract_bnd)
    assert bnd.num_elements == 12 * n**2
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_tag_and_select_boundary_faces(benchmark, n):
    problem = make_poisson_2d_mixed_bc(n)
    bnd = benchmark(problem.extract_bnd)
    assert "tag" in bnd.element_groups[0].cell_data
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_condense_dirichlet_bc_2d(benchmark, n):
    problem = make_poisson_2d_dirichlet(n)
    system = benchmark(problem.condense)
    assert system.K.shape[0] == (n + 1) ** 2 - 4 * n
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_assemble_neumann_bc(benchmark, n):
    problem = make_poisson_2d_mixed_bc(n)
    b = benchmark(problem.assemble_rhs)
    assert b.shape == ((n + 1) ** 2,)
    benchmark.extra_info.update(problem.meta)


@pytest.mark.parametrize("n", SIZES_2D)
def test_l2_project_2d(benchmark, n):
    problem = make_mass_projection_2d(n)
    u = benchmark(problem.run_full)
    assert u.shape == ((n + 1) ** 2,)
    benchmark.extra_info.update(problem.meta)
