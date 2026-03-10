"""Solve the heat equation u_t = Δu on [0,1]² using Crank-Nicolson.

Problem:
    u_t = Δu   on Ω = [0,1]², t ∈ [0, T]
    u = 0      on ∂Ω  (homogeneous Dirichlet)
    u(x,y,0) = sin(πx) sin(πy)

Analytical solution:
    u(x,y,t) = exp(-2π²t) sin(πx) sin(πy)
"""

import matplotlib.pyplot as plt
import matplotlib.tri
import numpy as np
import scipy.sparse.linalg as spla

from yggdrasil import (
    assemble_bilinear_form,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    l2_error,
    mass_form,
    unit_square_tri_mesh,
)


def main():
    n = 32
    dt = 0.01
    T = 0.1
    n_steps = int(round(T / dt))

    mesh = unit_square_tri_mesh(n)

    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=2)
    M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)

    # All boundary nodes have homogeneous Dirichlet BC
    bnd = extract_boundary(mesh)
    bc_dofs = bnd.point_data["original_node_index"]

    # Condense both matrices to free DOFs
    dummy_b = np.zeros(mesh.num_nodes)
    sys_K = condense_dirichlet_bc(K, dummy_b, bc_dofs)
    sys_M = condense_dirichlet_bc(M, dummy_b, bc_dofs)
    K_free = sys_K.K
    M_free = sys_M.K

    # Pre-factor implicit system A = M_free + 0.5*dt*K_free (done once)
    A_lu = spla.splu((M_free + 0.5 * dt * K_free).tocsc())
    B = M_free - 0.5 * dt * K_free

    # Initial condition at free DOFs
    free_dofs = np.setdiff1d(np.arange(mesh.num_nodes), bc_dofs)
    x_free = mesh.nodes[free_dofs]
    u_free = np.sin(np.pi * x_free[:, 0]) * np.sin(np.pi * x_free[:, 1])

    u_initial = sys_K.reconstruct(u_free)

    # Crank-Nicolson time-stepping
    for _ in range(n_steps):
        u_free = A_lu.solve(B @ u_free)

    u_final = sys_K.reconstruct(u_free)

    t_final = n_steps * dt
    err = l2_error(
        mesh,
        u_final,
        lambda x: np.exp(-2 * np.pi**2 * t_final) * np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1]),
        quadrature_order=2,
    )
    print(f"T = {t_final:.3f},  dt = {dt},  n_steps = {n_steps}")
    print(f"L² error at T: {err:.3e}")

    # Plot initial and final state
    triang = matplotlib.tri.Triangulation(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        mesh.element_groups[0].connectivity,
    )
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    tpc0 = axes[0].tripcolor(triang, u_initial, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc0, ax=axes[0])
    axes[0].set_title("u(x,y,0) — initial condition")
    axes[0].set_aspect("equal")

    tpc1 = axes[1].tripcolor(triang, u_final, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc1, ax=axes[1])
    axes[1].set_title(f"u(x,y,{t_final}) — Crank-Nicolson")
    axes[1].set_aspect("equal")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
