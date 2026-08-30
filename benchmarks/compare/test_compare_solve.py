"""Compare yggdrasil vs scikit-fem end-to-end solves (timing + agreement)."""

from __future__ import annotations

import numpy as np
import pytest

from benchmarks.compare import problems_skfem
from benchmarks.problems import make_poisson_2d_dirichlet, make_poisson_3d_dirichlet
from benchmarks.sizes import SIZES_2D, SIZES_3D


def _agree(u_ygg, u_skf_on_ygg, rtol=1e-6):
    scale = max(np.abs(u_ygg).max(), 1e-30)
    np.testing.assert_allclose(u_ygg, u_skf_on_ygg, atol=rtol * scale)


@pytest.mark.parametrize("n", SIZES_2D)
@pytest.mark.parametrize("library", ["yggdrasil", "skfem"])
def test_full_poisson_2d_dirichlet(benchmark, skfem, n, library):  # noqa: ARG001
    ygg = make_poisson_2d_dirichlet(n)
    skf = problems_skfem.poisson_2d_dirichlet(ygg.mesh)

    u_ygg = ygg.run_full()
    m_skf, _ = problems_skfem._tri_basis(ygg.mesh)
    u_skf = problems_skfem.reorder_to_yggdrasil(m_skf, skf.solve_full(), ygg.mesh)
    _agree(u_ygg, u_skf)

    target = ygg.run_full if library == "yggdrasil" else skf.solve_full
    benchmark(target)
    benchmark.extra_info.update(ygg.meta, library=library)


@pytest.mark.parametrize("n", SIZES_3D)
@pytest.mark.parametrize("library", ["yggdrasil", "skfem"])
def test_full_poisson_3d_dirichlet(benchmark, skfem, n, library):  # noqa: ARG001
    ygg = make_poisson_3d_dirichlet(n)
    skf = problems_skfem.poisson_3d_dirichlet(ygg.mesh)

    u_ygg = ygg.run_full()
    m_skf, _ = problems_skfem._tet_basis(ygg.mesh)
    u_skf = problems_skfem.reorder_to_yggdrasil(m_skf, skf.solve_full(), ygg.mesh)
    _agree(u_ygg, u_skf)

    target = ygg.run_full if library == "yggdrasil" else skf.solve_full
    benchmark(target)
    benchmark.extra_info.update(ygg.meta, library=library)
