# TODO

## Suggestions by Claude
Ordered feature suggestions (by applicability-to-effort ratio):

1. **Neumann boundary conditions**. The RHS surface integral ∫_Γ g·v dS is nearly free — assemble_load_vector already runs on boundary meshes (it's used inside project_dirichlet_bc). Just expose it. Unlocks: heat flux, traction BCs, mixed-BC Poisson.

2. **Boundary region splitting/tagging**. Without this, Neumann BCs are only usable on the entire boundary. Add tag-by-geometry to select_boundary_faces (e.g., "select faces where x=0"). Moderate effort, but it's a prerequisite for any realistic multi-BC problem (e.g., insulated vs. heated sides).

3. **Export mass form + time-stepping (heat equation)**. mass_form = lambda N, grad_N: np.einsum("qi,qj->qij", N, N) is one line — it's already used inside project_dirichlet_bc. Expose it, add a simple θ-method stepper (backward Euler first), and you can solve the heat equation ∂u/∂t - ∇²u = f. Big unlock for parabolic PDEs.

4. **L² and H¹ error norm utilities**. Low effort: a function compute_error(mesh, u_h, u_exact) returning the L² and H¹ seminorm errors over the mesh, using existing quadrature infrastructure. Enables proper convergence studies
  beyond the current "check a known value" approach, and is prerequisite for adaptive refinement later.

5. **Robin (mixed) boundary conditions**. Extends Neumann: ∂u/∂n + αu = g adds a boundary mass matrix term α∫_Γ u v dS to K plus a load term. Moderate effort. Unlocks: convection/radiation heat transfer BCs, which are ubiquitous in engineering.

6. **Eigenvalue problems**. Once the mass matrix is available, scipy.sparse.linalg.eigsh(K, M) is one call. Unlocks: natural frequencies of vibration, buckling loads, Laplacian eigenvalues. Very low additional implementation cost.

7. **Nonlinear problems via Newton–Raphson**. Add assemble_residual and assemble_tangent (the latter is just assemble_bilinear_form with a position-dependent integrand). A generic Newton loop with line search is ~30 lines. Unlocks: nonlinear diffusion, p-Laplacian, steady Navier–Stokes (later).

8. **Vector-valued DOFs (linear elasticity)**. This requires the biggest architectural change: multiple DOFs per node, block assembly. But it unlocks an enormous class of problems — linear elasticity, Stokes flow, beam/plate problems. Worth planning early so the architecture doesn't need to be reversed later.

9. **VTK/XDMF output for Paraview**. A writer that maps Mesh + solution vector to .vtu format. Moderate effort (use meshio as a dependency). Essential for visualizing 3D results; matplotlib doesn't scale. Particularly valuable once you solve elasticity or 3D problems.

10. **Adaptive mesh refinement (h-refinement)**. The hardest item: requires mesh refinement algorithms (bisection for triangles, hanging nodes or conforming refinement), error estimators, and marking strategies. But it's transformative for efficiency on problems with singularities (re-entrant corners, point loads). Worth considering as a longer-term goal that shapes architectural decisions made now.

My top recommendation: Do 1 + 2 + 3 together — they're tightly coupled, low effort individually, and together let you solve the full heat equation with realistic mixed BCs, which
dramatically expands the demonstrable capability of the library.
  
  
## Other
- Implement more examples that are comparable to the examples of scikit-fem,
  see https://scikit-fem.readthedocs.io/en/latest/listofexamples.html
- Add benchmarking of solving certain systems (and perhaps do equivalent
  benchmarking of other libraries for speed comparisons -- e.g., scikit-fem)
