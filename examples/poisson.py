"""Solve the Poisson equation -Δu = 1 on [0,1]² with u=0 on the boundary.

Uses Tri3 elements on a structured triangular mesh.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla

from yggdrasil import apply_dirichlet_bc, assemble_bilinear_form, assemble_load_vector, extract_boundary, unit_square_tri_mesh


def main():
    n = 32
    mesh = unit_square_tri_mesh(n)

    # Assemble stiffness matrix
    def stiffness_form(N, grad_N):
        return np.einsum("qia,qja->qij", grad_N, grad_N)

    K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)

    # Assemble load vector (f = 1)
    b = assemble_load_vector(mesh, f_val=1.0, quadrature_order=1)

    # Apply Dirichlet BC (u = 0 on boundary)
    bc_nodes = extract_boundary(mesh).point_data["original_node_index"]
    K, b = apply_dirichlet_bc(K, b, bc_nodes)

    # Solve
    u = spla.spsolve(K, b)

    # Plot
    triangulation = plt.matplotlib.tri.Triangulation(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        mesh.element_groups[0].connectivity,
    )
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    tpc = ax.tripcolor(triangulation, u, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc, ax=ax, label="u")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Poisson: $-\\Delta u = 1$ on $[0,1]^2$ (n={n})")
    ax.set_aspect("equal")
    plt.tight_layout()
    plt.show()

    print(f"Max solution value: {u.max():.6f}")
    print("(Exact max ≈ 0.0737 at center)")


if __name__ == "__main__":
    main()
