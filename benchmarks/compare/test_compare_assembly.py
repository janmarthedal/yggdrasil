"""Compare yggdrasil vs scikit-fem assembly (timing + agreement)."""

from __future__ import annotations

import numpy as np
import pytest

from benchmarks.compare import problems_skfem
from benchmarks.problems import make_laplace_eigenvalue_2d, make_poisson_2d_dirichlet
from benchmarks.sizes import SIZES_2D


@pytest.mark.parametrize("n", SIZES_2D)
@pytest.mark.parametrize("library", ["yggdrasil", "skfem"])
def test_stiffness_assembly_2d(benchmark, skfem, n, library):  # noqa: ARG001
    ygg = make_poisson_2d_dirichlet(n)
    skf = problems_skfem.poisson_2d_dirichlet(ygg.mesh)

    K_ygg = ygg.assemble_lhs()
    K_skf = skf.assemble_lhs()

    # Same discretisation -> identical operator (up to node permutation absorbed
    # by comparing eigenvalue-invariant quantities).
    assert K_ygg.shape == K_skf.shape
    np.testing.assert_allclose(K_ygg.diagonal().sum(), K_skf.diagonal().sum(), rtol=1e-10)
    np.testing.assert_allclose(abs(K_ygg).sum(), abs(K_skf).sum(), rtol=1e-10)

    target = ygg.assemble_lhs if library == "yggdrasil" else skf.assemble_lhs
    benchmark(target)
    benchmark.extra_info.update(ygg.meta, library=library)


@pytest.mark.parametrize("n", SIZES_2D)
@pytest.mark.parametrize("library", ["yggdrasil", "skfem"])
def test_mass_assembly_2d(benchmark, skfem, n, library):  # noqa: ARG001
    ygg = make_laplace_eigenvalue_2d(n)
    skf = problems_skfem.mass_2d(ygg.mesh)

    _, M_ygg = ygg.assemble_lhs()
    M_skf = skf.assemble_lhs()
    np.testing.assert_allclose(M_ygg.sum(), M_skf.sum(), rtol=1e-10)  # both = area = 1

    target = (lambda: ygg.assemble_lhs()[1]) if library == "yggdrasil" else skf.assemble_lhs
    benchmark(target)
    benchmark.extra_info.update(ygg.meta, library=library)
