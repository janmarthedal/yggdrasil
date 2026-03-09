"""Solve the Poisson equation with mixed Dirichlet/Neumann BCs on [0,1]².

Problem:
    -∇²u = 0   on Ω = [0,1]²
     u   = 0   on left edge  (x = 0)   ← Dirichlet
    ∂u/∂n = 1  on right edge (x = 1)   ← Neumann
    ∂u/∂n = 0  on top/bottom           ← natural (no assembly)

Exact solution: u(x, y) = x
"""

import matplotlib.pyplot as plt
import matplotlib.tri
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


def main():
    n = 32
    mesh = unit_square_tri_mesh(n)

    # Assemble stiffness matrix and zero load vector (f = 0)
    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

    # Tag boundary faces
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] < 1e-10, tag=1)      # left: Dirichlet
    bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] > 1 - 1e-10, tag=2)  # right: Neumann

    # Add Neumann contribution (∂u/∂n = 1 on right edge)
    neumann_mesh = select_boundary_faces(bnd, tag=2)
    b += assemble_neumann_bc(neumann_mesh, g=1.0, quadrature_order=1, n_dofs=mesh.num_nodes)

    # Apply Dirichlet BC (u = 0 on left edge)
    dirichlet_mesh = select_boundary_faces(bnd, tag=1)
    bc_nodes = dirichlet_mesh.point_data["original_node_index"]
    system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    # Solve
    u = system.reconstruct(spla.spsolve(system.K, system.b))

    # Report error vs exact solution u = x
    u_exact = mesh.nodes[:, 0]
    print(f"Max nodal error vs u_exact = x: {np.max(np.abs(u - u_exact)):.2e}")

    # Plot numerical solution and error side by side
    triangulation = matplotlib.tri.Triangulation(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        mesh.element_groups[0].connectivity,
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    tpc0 = axes[0].tripcolor(triangulation, u, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc0, ax=axes[0], label="u")
    axes[0].set_title(f"Numerical solution (n={n})")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    tpc1 = axes[1].tripcolor(triangulation, np.abs(u - u_exact), shading="gouraud", cmap="magma")
    fig.colorbar(tpc1, ax=axes[1], label="|u - u_exact|")
    axes[1].set_title("Pointwise error |u - x|")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
