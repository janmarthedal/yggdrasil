"""Solve the Poisson equation -Δu = 1 on [0,1]² with u=0 on the boundary.

Uses Tri3 elements on a structured triangular mesh.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from yggdrasil import ElementGroup, Mesh, assemble_bilinear_form, assemble_load_vector
from yggdrasil.elements import Tri3


# TODO: Make into general mesh creation function inside library and use that instead
def make_unit_square_tri_mesh(n: int) -> Mesh:
    """Create a structured triangular mesh on [0,1]² with (n+1)² nodes."""
    x = np.linspace(0, 1, n + 1)
    y = np.linspace(0, 1, n + 1)
    xx, yy = np.meshgrid(x, y)
    nodes = np.column_stack([xx.ravel(), yy.ravel()])

    triangles = []
    for j in range(n):
        for i in range(n):
            # Node indices for the quad cell
            n0 = j * (n + 1) + i
            n1 = n0 + 1
            n2 = n0 + (n + 1) + 1
            n3 = n0 + (n + 1)
            # Split into two triangles
            triangles.append([n0, n1, n2])
            triangles.append([n0, n2, n3])

    conn = np.array(triangles, dtype=np.intp)
    group = ElementGroup(element=Tri3(), connectivity=conn)
    return Mesh(nodes, [group])


# TODO: Make mesh library function that finds boundary of a mesh
def find_boundary_nodes(mesh: Mesh, tol: float = 1e-12) -> np.ndarray:
    """Find node indices on the boundary of [0,1]²."""
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    on_boundary = (x < tol) | (x > 1 - tol) | (y < tol) | (y > 1 - tol)
    return np.where(on_boundary)[0]


# TODO: Make library function that helps assemble the linear system
def apply_dirichlet_bc(K: sp.csr_matrix, b: np.ndarray, bc_nodes: np.ndarray, bc_val: float = 0.0):
    """Apply Dirichlet BC by zeroing rows/cols and setting diagonal to 1."""
    K = K.tolil()
    for node in bc_nodes:
        K[node, :] = 0
        K[:, node] = 0
        K[node, node] = 1.0
        b[node] = bc_val
    return K.tocsr(), b


def main():
    n = 32
    mesh = make_unit_square_tri_mesh(n)

    # Assemble stiffness matrix
    def stiffness_form(N, grad_N):
        return np.einsum("qia,qja->qij", grad_N, grad_N)

    K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)

    # Assemble load vector (f = 1)
    b = assemble_load_vector(mesh, f_val=1.0, quadrature_order=1)

    # Apply Dirichlet BC (u = 0 on boundary)
    bc_nodes = find_boundary_nodes(mesh)
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
