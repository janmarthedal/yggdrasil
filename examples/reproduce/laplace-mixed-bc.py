r"""Laplace with mixed boundary conditions on an annular sector.

This example corresponds to Example 13 of scikit-fem's Gallery of examples,
https://scikit-fem.readthedocs.io/en/latest/listofexamples.html#example-13-laplace-with-mixed-boundary-conditions,
with source code https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex13.py

Problem
-------
Solve the Laplace equation

    Δu = 0   in Ω

on the annular sector

    Ω = {(x, y) : 1 < x² + y² < 4,  0 < θ < π/2},

where θ = arctan(y/x), with mixed boundary conditions:

    u = 0         on Γ_ground   (y = 0, i.e. θ = 0)       ← Dirichlet (earthed)
    u = 1         on Γ_positive (x = 0, i.e. θ = π/2)     ← Dirichlet (charged)
    ∂u/∂n = 0     on inner/outer arcs (r = 1 and r = 2)   ← Neumann (insulated)

Exact solution
--------------
    u(x, y) = 2θ/π = (2/π) arctan(y/x)

The field strength is |∇u|² = 4 / (π² (x² + y²)), and integrating over Ω gives
the conductance (for unit potential difference and unit conductivity):

    C = ‖∇u‖²_{L²(Ω)} = 2 ln 2 / π ≈ 0.4413

Approach
--------
The mesh is first built as a rectangle in (r, θ)-space,

    r ∈ [1, 2],  θ ∈ [0, π/2],

where the two isopotential boundaries are aligned with the coordinate axes and
are easy to identify.  The mesh is then mapped to (x, y)-space via

    x = r cos θ,  y = r sin θ.

The remaining two boundaries (constant-r arcs) receive the homogeneous Neumann
condition automatically (no flux assembly is needed).
"""

import matplotlib.pyplot as plt
import matplotlib.tri
import numpy as np
import scipy.sparse.linalg as spla

from yggdrasil import (
    assemble_bilinear_form,
    assemble_load_vector,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    select_boundary_faces,
    tag_boundary_faces,
)
from yggdrasil.mesh_generators import map_mesh_points, rectangular_tri_mesh


def main():
    # Mesh: 21 nodes along r in [1,2], 31 nodes along θ in [0, π/2]
    nr, nt = 21, 31
    mesh_rtheta = rectangular_tri_mesh(1.0, 2.0, 0.0, np.pi / 2, nr, nt)

    # Map (r, θ) → (x, y)
    def polar_to_cartesian(pts):
        r, theta = pts[:, 0], pts[:, 1]
        return np.column_stack([r * np.cos(theta), r * np.sin(theta)])

    mesh = map_mesh_points(mesh_rtheta, polar_to_cartesian)

    # Assemble stiffness matrix; load vector is zero (Laplace, no source)
    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

    # Tag Dirichlet boundary faces (inner/outer arcs remain untagged → natural ∂u/∂n=0)
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda c: c[:, 1] < 1e-10, tag=1)  # y=0: ground
    bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] < 1e-10, tag=2)  # x=0: positive

    # Collect Dirichlet DOFs and values for both isopotential boundaries
    ground_mesh = select_boundary_faces(bnd, tag=1)
    positive_mesh = select_boundary_faces(bnd, tag=2)

    bc_dofs = np.concatenate([
        ground_mesh.point_data["original_node_index"],
        positive_mesh.point_data["original_node_index"],
    ])
    bc_vals = np.concatenate([
        np.zeros(ground_mesh.num_nodes),
        np.ones(positive_mesh.num_nodes),
    ])

    system = condense_dirichlet_bc(K, b, bc_dofs, bc_val=bc_vals)

    # Solve
    u = system.reconstruct(spla.spsolve(system.K, system.b))

    # Compare with exact solution u = 2θ/π
    x, y = mesh.nodes[:, 0], mesh.nodes[:, 1]
    u_exact = 2 * np.arctan2(y, x) / np.pi
    max_err = np.max(np.abs(u - u_exact))

    # Conductance = u^T K u (exact value: 2 ln 2 / π)
    conductance_fem = float(u @ K @ u)
    conductance_exact = 2 * np.log(2) / np.pi

    print(f"Max nodal error vs u_exact = 2θ/π: {max_err:.4e}")
    print(f"Conductance (FEM):   {conductance_fem:.6f}")
    print(f"Conductance (exact): {conductance_exact:.6f}  (2 ln 2 / π)")

    # Plot numerical solution and pointwise error side by side
    triangulation = matplotlib.tri.Triangulation(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        mesh.element_groups[0].connectivity,
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    tpc0 = axes[0].tripcolor(triangulation, u, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc0, ax=axes[0], label="u")
    axes[0].set_title("Numerical solution")
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("y")
    axes[0].set_aspect("equal")

    tpc1 = axes[1].tripcolor(triangulation, np.abs(u - u_exact), shading="gouraud", cmap="magma")
    fig.colorbar(tpc1, ax=axes[1], label="|u - u_exact|")
    axes[1].set_title("Pointwise error  |u − 2θ/π|")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_aspect("equal")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
