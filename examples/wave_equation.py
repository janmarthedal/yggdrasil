"""Solve the scalar wave equation u_tt = Δu on [0,1]² using Newmark-β.

Problem:
    u_tt = Δu  on Ω = [0,1]², t ∈ [0, T]
    u = 0      on ∂Ω  (homogeneous Dirichlet)
    u(x,y,0)   = sin(πx) sin(πy)
    u_t(x,y,0) = 0

Analytical solution (free vibration at the lowest mode):
    ω₁ = π√2
    u(x,y,t) = cos(ω₁ t) sin(πx) sin(πy)

Scheme: Newmark-β with β=0.25, γ=0.5 (implicit trapezoidal, unconditionally stable).
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse.linalg as spla

from yggdrasil import (
    assemble_bilinear_form,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    mass_form,
    unit_square_tri_mesh,
)


def main():
    n = 32
    T = 1.0
    dt = 0.02
    n_steps = int(round(T / dt))
    beta = 0.25
    gamma = 0.5
    omega1 = np.pi * np.sqrt(2)

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

    # Pre-factor implicit system A = M_free + β·dt²·K_free (done once)
    A_lu = spla.splu((M_free + beta * dt**2 * K_free).tocsc())

    # Initial conditions at free DOFs
    free_dofs = np.setdiff1d(np.arange(mesh.num_nodes), bc_dofs)
    x_free = mesh.nodes[free_dofs]
    u_free = np.sin(np.pi * x_free[:, 0]) * np.sin(np.pi * x_free[:, 1])
    v_free = np.zeros_like(u_free)
    # a₀ = M_free⁻¹ (−K_free u₀) from M ü + K u = 0 at t=0
    a_free = spla.spsolve(M_free, -(K_free @ u_free))

    # Locate center node (0.5, 0.5) — interior node for n even
    center_node = int(np.argmin(np.sum((mesh.nodes - 0.5) ** 2, axis=1)))
    center_free_idx = int(np.searchsorted(free_dofs, center_node))
    assert free_dofs[center_free_idx] == center_node, "center node not in free DOFs"

    t_array = np.linspace(0, T, n_steps + 1)
    u_center = np.empty(n_steps + 1)
    u_center[0] = u_free[center_free_idx]

    # Newmark-β time-stepping
    for step in range(n_steps):
        # Predictor step
        u_pred = u_free + dt * v_free + dt**2 * (0.5 - beta) * a_free
        # Solve (M + β·dt²·K) a_{n+1} = −K·u_pred
        a_new = A_lu.solve(-(K_free @ u_pred))
        # Corrector step
        u_free = u_pred + beta * dt**2 * a_new
        v_free = v_free + dt * ((1 - gamma) * a_free + gamma * a_new)
        a_free = a_new
        u_center[step + 1] = u_free[center_free_idx]

    # Report error vs analytical at final time
    x_center = mesh.nodes[center_node]
    amplitude = np.sin(np.pi * x_center[0]) * np.sin(np.pi * x_center[1])
    u_center_exact = amplitude * np.cos(omega1 * t_array)
    max_err = float(np.max(np.abs(u_center - u_center_exact)))
    print(f"T = {T:.3f},  dt = {dt},  n_steps = {n_steps},  ω₁ = {omega1:.6f}")
    print(f"Max center-node error vs analytical: {max_err:.3e}")

    # Plot center-node displacement vs analytical
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t_array, u_center, label="Newmark-β (FEM)")
    ax.plot(t_array, u_center_exact, "--", label=r"Analytical $\cos(\omega_1 t)$")
    ax.set_xlabel("t")
    ax.set_ylabel("u(0.5, 0.5, t)")
    ax.set_title("Center-node displacement — wave equation")
    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
