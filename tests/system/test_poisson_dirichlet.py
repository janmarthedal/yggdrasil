"""System test: Poisson equation with homogeneous Dirichlet BCs on the unit square.

Problem:
    -∇²u = f   on Ω = [0,1]²
     u   = 0   on ∂Ω

Analytical solution:
    u(x,y) = sin(πx) sin(πy)

Source term:
    f(x,y) = 2π² sin(πx) sin(πy)

P1 elements on a uniform triangular mesh converge at O(h²) in the maximum
nodal error.  The test verifies this rate by solving at two resolutions and
checking that the error drops by a factor of roughly 4 when h is halved.
"""

import numpy as np
import scipy.sparse.linalg

from yggdrasil import assemble_bilinear_form, assemble_load_vector, condense_dirichlet_bc, extract_boundary, grad_grad_form
from yggdrasil.mesh_generators import unit_square_tri_mesh


def u_exact(x):
    return np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])


def source(x):
    return 2 * np.pi**2 * np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])


def solve_poisson(n):
    """Return the max nodal error for a uniform mesh with n subdivisions."""
    mesh = unit_square_tri_mesh(n)

    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, source, quadrature_order=5)

    bnd = extract_boundary(mesh)
    bc_nodes = bnd.point_data["original_node_index"]
    system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)
    u = system.reconstruct(scipy.sparse.linalg.spsolve(system.K, system.b))
    return np.max(np.abs(u - u_exact(mesh.nodes)))


def test_poisson_dirichlet_accuracy():
    """Max nodal error on a 32×32 mesh matches the known value to 0.1%."""
    np.testing.assert_allclose(solve_poisson(32), 8.028035006669709e-4, rtol=1e-3)


def test_poisson_dirichlet_convergence():
    """Convergence ratio when halving h matches the known value of ~3.994 to 0.1%."""
    err_coarse = solve_poisson(16)
    err_fine = solve_poisson(32)
    np.testing.assert_allclose(err_coarse / err_fine, 3.994222254543219, rtol=1e-3)
