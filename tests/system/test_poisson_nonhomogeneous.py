"""System test: Laplace equation with non-homogeneous Dirichlet BCs on the unit square.

Problem:
    -∇²u = 0   on Ω = [0,1]²
     u   = 0   on left, right, and bottom edges
     u   = sin(πx)  on the top edge (y = 1)

Analytical solution:
    u(x,y) = sin(πx) · sinh(πy) / sinh(π)

The boundary data is prescribed via L² projection (project_dirichlet_bc)
with g = u_exact evaluated on the full boundary; this is zero on three edges
and sin(πx)·sinh(πy)/sinh(π) on the top edge (where it equals sin(πx)).

P1 elements converge at O(h²) in the maximum nodal error.
"""

import numpy as np
import scipy.sparse.linalg

from yggdrasil import apply_dirichlet_bc, assemble_bilinear_form, assemble_load_vector, extract_boundary, grad_grad_form, project_dirichlet_bc
from yggdrasil.mesh_generators import unit_square_tri_mesh


def u_exact(x):
    return np.sin(np.pi * x[:, 0]) * np.sinh(np.pi * x[:, 1]) / np.sinh(np.pi)


def solve_laplace(n):
    """Return the max nodal error for a uniform mesh with n subdivisions."""
    mesh = unit_square_tri_mesh(n)

    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

    bnd = extract_boundary(mesh)
    bc_nodes, bc_vals = project_dirichlet_bc(bnd, g=u_exact, quadrature_order=2)
    K, b = apply_dirichlet_bc(K, b, bc_nodes, bc_val=bc_vals)

    u = scipy.sparse.linalg.spsolve(K, b)
    return np.max(np.abs(u - u_exact(mesh.nodes)))


def test_laplace_nonhomogeneous_accuracy():
    """Max nodal error on a 32×32 mesh matches the known value to 0.1%."""
    np.testing.assert_allclose(solve_laplace(32), 8.035128524226387e-4, rtol=1e-3)


def test_laplace_nonhomogeneous_convergence():
    """Halving h should reduce the max nodal error by a factor of ~4."""
    err_coarse = solve_laplace(16)
    err_fine = solve_laplace(32)
    np.testing.assert_allclose(err_coarse / err_fine, 4.004808646390466, rtol=1e-3)
