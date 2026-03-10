"""System test: Laplacian eigenvalue problem on the unit square.

Problem:
    -∇²u = λ u   on Ω = [0,1]²
         u = 0   on ∂Ω

Analytical eigenvalues:
    λ_{mn} = π²(m² + n²),   m, n = 1, 2, 3, ...

The smallest eigenvalue is λ_{11} = 2π² ≈ 19.739.

This test assembles K (stiffness) and M (mass) matrices and solves the
generalized eigenvalue problem  K u = λ M u  with homogeneous Dirichlet BCs.
"""

import numpy as np
import scipy.sparse.linalg

from yggdrasil import assemble_bilinear_form, condense_dirichlet_bc, extract_boundary, grad_grad_form, mass_form
from yggdrasil.mesh_generators import unit_square_tri_mesh


def test_smallest_eigenvalue():
    """Smallest eigenvalue of −Δ on the unit square is close to 2π²."""
    n = 16
    mesh = unit_square_tri_mesh(n)

    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)

    bnd = extract_boundary(mesh)
    bc_nodes = bnd.point_data["original_node_index"]

    # Condense BCs from both K and M: remove BC rows/cols
    system_K = condense_dirichlet_bc(K, np.zeros(mesh.num_nodes), bc_nodes, bc_val=0.0)
    system_M = condense_dirichlet_bc(M, np.zeros(mesh.num_nodes), bc_nodes, bc_val=0.0)

    K_free = system_K.K
    M_free = system_M.K

    # Solve for the 6 smallest eigenvalues
    eigenvalues, _ = scipy.sparse.linalg.eigsh(K_free, M=M_free, k=6, which="SM")
    eigenvalues = np.sort(eigenvalues)

    lambda_1_exact = 2 * np.pi**2  # ≈ 19.739

    # P1 on a 16×16 mesh: converges from above, so allow 2% relative tolerance
    np.testing.assert_allclose(eigenvalues[0], lambda_1_exact, rtol=2e-2)
