# TODO

## Suggestions by Claude
Ordered feature suggestions (by applicability-to-effort ratio):

1. **Mass matrix form + built-in forms library** (DONE). forms.py currently has only grad_grad_form. Adding ∫ uv dx (mass form) is a one-liner. This alone unlocks eigenvalue problems (structural vibration, acoustic resonance, buckling) and time-dependent PDEs (heat equation, wave equation). Every subsequent feature below benefits from it.

2. **Robin boundary conditions**. A boundary bilinear form ∫_Γ α u v dS plus an optional load term. It requires ~20 lines reusing the existing boundary mesh infrastructure. This unlocks realistic heat transfer (convection BCs), absorbing boundary conditions for wave problems, and radiation conditions.

3. **L² error norm utility**. A function l2_error(mesh, uh, u_exact, order) that integrates (uh − u_exact)² over elements. Trivial to implement using the existing quadrature/mapping pipeline. Essential for convergence studies, which are the primary way to validate any new feature, so the compounding benefit is high.

4. **Transient solver (heat/wave equation)**. A thin time-stepping wrapper (Euler, Crank-Nicolson, Newmark) that takes K, M and advances the solution. No new assembly machinery is needed. Solves the heat equation, wave equation, and parabolic variants of Poisson. Medium effort for broad scope.

5. **Eigenvalue problem solver** (DONE). A function that takes K and M and calls scipy.sparse.linalg.eigsh. Essentially a one-page wrapper. Solves structural vibration, Helmholtz resonance, and stability/buckling problems. High payoff for near-zero effort once the mass matrix exists.

6. **Linear elasticity (vector-valued problems)**. The DOFMap abstraction already anticipates this. The remaining work is the Lamé bilinear form ∫ σ(u):ε(v) dx and vector load assembly. This is the most broadly useful mechanical PDE and a natural next step after the scalar infrastructure is solid. Moderate effort.

7. **Solution interpolation and gradient evaluation at arbitrary points**. A evaluate(mesh, uh, points) function: for each query point, find the containing element (spatial search via scipy.spatial.cKDTree on element centroids + reference-coordinate inversion), then evaluate N(ξ) · u_elem. Needed for post-processing, flux computation, and coupling. Moderate effort but high practical value.

8. **Unstructured 2D mesh generation**. Use scipy.spatial.Delaunay or wrap the triangle library to mesh arbitrary polygonal domains (defined by vertices + holes). This removes the restriction to structured grids and enables geometry-driven examples. Moderate effort; the rest of the pipeline already handles unstructured meshes.

9. **L² projection / nodal recovery of derived fields**. Project a non-nodal derived quantity (e.g. the stress tensor or heat flux −∇u) back onto nodal values by solving M α = b where b_i = ∫ q N_i dx. Improves post-processing accuracy and is needed before tackling adaptive refinement. Low-to-medium effort, very useful.

10. **Nonlinear problems via Newton iteration**. A general Newton loop that re-assembles a tangent stiffness matrix and residual vector at each iteration. This requires a new assembly path accepting user-supplied residual/tangent callables. High effort, but it is the gateway to nonlinear elasticity, incompressible flow (Stokes/Navier-Stokes), and phase-field models. 
  
## Other
- Implement more examples that are comparable to the examples of scikit-fem,
  see https://scikit-fem.readthedocs.io/en/latest/listofexamples.html
- Add benchmarking of solving certain systems (and perhaps do equivalent
  benchmarking of other libraries for speed comparisons -- e.g., scikit-fem)
- Is computing `grad_phys` of `compute_physical_gradients` sometimes
  unneccesary work (see `assemble_load_vector`)?
