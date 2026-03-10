"""Unit tests for yggdrasil.error.l2_error."""

import numpy as np
import scipy.sparse.linalg

from yggdrasil import (
    assemble_bilinear_form,
    assemble_load_vector,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    l2_error,
)
from yggdrasil.mesh_generators import unit_square_tri_mesh


def _solve_poisson(n):
    """Solve -∇²u = f on unit square, u=0 on boundary; return (mesh, u)."""
    mesh = unit_square_tri_mesh(n)

    def source(x):
        return 2 * np.pi**2 * np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])

    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, source, quadrature_order=5)
    bnd = extract_boundary(mesh)
    bc_nodes = bnd.point_data["original_node_index"]
    system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)
    u = system.reconstruct(scipy.sparse.linalg.spsolve(system.K, system.b))
    return mesh, u


def test_l2_error_exact_solution_is_zero():
    """uh = FE interpolant of a linear function gives L² error ≈ 0.

    A linear function is exactly representable by P1 elements, so
    evaluating l2_error with uh set to nodal values of u_exact must
    return essentially zero (up to floating-point rounding).
    """
    mesh = unit_square_tri_mesh(4)

    def u_exact(x):
        return x[:, 0] + 2 * x[:, 1]  # linear: exactly in P1 space

    uh = u_exact(mesh.nodes)
    err = l2_error(mesh, uh, u_exact, quadrature_order=2)
    assert err < 1e-13, f"Expected ~0, got {err}"


def test_l2_error_constant_one_on_unit_square():
    """uh=1, u_exact=0 on unit square: L² error should equal sqrt(area)=1."""
    mesh = unit_square_tri_mesh(8)
    uh = np.ones(mesh.num_nodes)
    err = l2_error(mesh, uh, lambda x: np.zeros(len(x)), quadrature_order=2)
    np.testing.assert_allclose(err, 1.0, rtol=1e-12)


def test_l2_error_monotone_in_resolution():
    """Coarser mesh has larger L² error than finer mesh."""

    def u_exact(x):
        return np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])

    mesh_coarse, u_coarse = _solve_poisson(8)
    mesh_fine, u_fine = _solve_poisson(16)

    err_coarse = l2_error(mesh_coarse, u_coarse, u_exact, quadrature_order=5)
    err_fine = l2_error(mesh_fine, u_fine, u_exact, quadrature_order=5)

    assert err_coarse > err_fine, f"Expected coarse ({err_coarse}) > fine ({err_fine})"
