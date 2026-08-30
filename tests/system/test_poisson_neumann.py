"""System test: Poisson equation with mixed Dirichlet/Neumann BCs on [0,1]².

Problem:
    -∇²u = 0   on Ω = [0,1]²
     u   = 0   on Γ_D (left edge, x = 0)
    ∂u/∂n = 1  on Γ_N (right edge, x = 1)
    ∂u/∂n = 0  on top and bottom (natural, no assembly needed)

Exact solution:
    u(x, y) = x

P1 elements on a uniform triangular mesh recover the exact solution at all nodes
(u = x is in the FE space), so the test checks that the max nodal error is below
a tight tolerance for a modest mesh.
"""

import numpy as np
import scipy.sparse.linalg as spla

from yggdrasil import (
    assemble_bilinear_form,
    assemble_load_vector,
    assemble_neumann_bc,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    select_boundary_faces,
    tag_boundary_faces,
    unit_square_tri_mesh,
)


def solve_poisson_neumann(n):
    """Return the max nodal error for a uniform mesh with n subdivisions."""
    mesh = unit_square_tri_mesh(n)

    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] < 1e-10, tag=1)   # left: Dirichlet
    bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] > 1 - 1e-10, tag=2)  # right: Neumann

    neumann_mesh = select_boundary_faces(bnd, tag=2)
    b += assemble_neumann_bc(neumann_mesh, g=1.0, quadrature_order=1, dofs=mesh.num_nodes)

    dirichlet_mesh = select_boundary_faces(bnd, tag=1)
    bc_nodes = dirichlet_mesh.point_data["original_node_index"]
    system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    u = system.reconstruct(spla.spsolve(system.K, system.b))
    u_exact = mesh.nodes[:, 0]
    return np.max(np.abs(u - u_exact))


def test_poisson_neumann_accuracy():
    """Max nodal error on an 8×8 mesh is below 1e-12 (exact solution is in FE space)."""
    err = solve_poisson_neumann(8)
    assert err < 1e-12


def test_poisson_neumann_finer_mesh():
    """Max nodal error on a 16×16 mesh is also below 1e-12."""
    err = solve_poisson_neumann(16)
    assert err < 1e-12
