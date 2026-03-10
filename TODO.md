# TODO

## Suggestions by Claude
Feature Suggestions (Highest Value-to-Effort First)

1. Extended built-in forms (reaction, advection, variable-coefficient diffusion)

Effort: Very low (~30–50 lines)
Add to forms.py:
- reaction_form(c): ∫ c u v dx — enables Helmholtz (-Δu + k²u = f), reaction-diffusion
- advection_form(b): ∫ (b·∇u) v dx — for convection-diffusion with given velocity
- diffusion_form(kappa): ∫ κ(x) ∇u·∇v dx — spatially varying material parameters

Unlocks: Helmholtz, convection-diffusion, reaction-diffusion, screened Poisson — a very large class of problems.

---
2. Robin boundary conditions

Effort: Low (~30 lines in assemble.py)
∫_Γ α u v dS — assembles a boundary mass matrix with coefficient α and adds it to the stiffness matrix. Currently you have Dirichlet and Neumann; Robin is the third classical type.

Unlocks: Natural thermal convection BCs, impedance BCs in acoustics/electromagnetics, radiation BCs, and enables well-posed problems that are ill-conditioned under pure Neumann.

---
3. Vector DOFMap (n_components > 1)

Effort: Moderate (~100–150 lines)
The stub already exists in dof_map.py with raise NotImplementedError. Implement interleaved DOF blocks for vector fields. Then add forms:
- linear_elasticity_form(lam, mu): the symmetric gradient bilinear form ∫ σ(u):ε(v) dx
- Body force load vector for vector RHS

Unlocks: Linear elasticity (structural mechanics), membrane/plate problems — a whole new application domain beyond scalars.

---
4. Crank-Nicolson and higher-order time integration

Effort: Low (~50–80 lines, likely a new timestepping.py)
Currently the heat/wave examples do raw time-stepping. A utility that takes (M, K, b_func, dt) and steps with θ-method (θ=0.5 → Crank-Nicolson, θ=1 → implicit Euler) or BDF2 would be
valuable. The assembled matrices only need to be factored once for constant-coefficient problems (cache the LU factorization).

Unlocks: Second-order accurate transient solvers without instability for stiff problems. Directly improves heat and wave equation accuracy.

---
5. Mesh import (Gmsh .msh format)

Effort: Moderate (~150–200 lines in io.py)
Support reading .msh v2/v4 ASCII files from Gmsh. Map Gmsh element type tags to your existing elements (Line2/Line3/Tri3/Tri6/Tet4/Hex8). Preserve physical group tags as cell data for
region/boundary selection.

Unlocks: Arbitrary complex geometries — curved boundaries, holes, real-world domains — without being limited to structured meshes. This is probably the single biggest unlock for practical
use.

---
6. Mixed finite elements / saddle-point systems

Effort: High (~300–400 lines, requires design work)
Support two different function spaces (e.g., Taylor-Hood P2/P1: velocity in Tri6, pressure in Tri3) sharing the same mesh. Assemble off-diagonal blocks ∫ div(u) q dx.

Unlocks: Stokes flow (incompressible viscous fluid), Darcy flow, mixed formulation of Poisson — opens up fluid mechanics.

---
7. Nonlinear solver support (Newton iteration)

Effort: Moderate (~150 lines in a new nonlinear.py)
A Newton loop that accepts residual and Jacobian assembly callbacks, with convergence monitoring. The user provides F(u) and dF/du as bilinear/linear form builders.

Unlocks: Nonlinear Poisson (e.g., p-Laplacian), nonlinear elasticity, steady Navier-Stokes (as Oseen iterations).

---
8. Error estimators and adaptive mesh refinement (AMR)

Effort: Very high (~500+ lines, requires refinement algorithms per element type)
Residual-based a posteriori error estimators (element residual + edge jump terms), followed by red-green or newest-vertex bisection refinement. Most of the complexity is in the mesh
refinement bookkeeping.

Unlocks: Efficient high-accuracy solutions near singularities (re-entrant corners, point loads, crack tips). High payoff for problems where uniform refinement is prohibitively expensive.

---
Summary Table

| # | Feature | Effort | New physics unlocked |
|---|---------|--------|----------------------|
| 1 | Extra forms (reaction, advection) | Very low | Helmholtz, convection-diffusion, reaction-diffusion |
| 2 | Robin BCs | Low | Thermal convection, impedance, radiation |
| 3 | Vector DOFMap + elasticity form | Moderate | Linear elasticity, structural mechanics |
| 4 | Crank-Nicolson / θ-method timestepping | Low | Stable 2nd-order transient solvers |
| 5 | Gmsh mesh import | Moderate | Complex real-world geometries |
| 6 | Mixed FEM (saddle-point) | High | Stokes/Darcy flow, incompressible fluids |
| 7 | Newton iteration | Moderate | Nonlinear PDEs |
| 8 | AMR + error estimators | Very high | Singular problems, efficiency |

## Other
- Add LICENCE
- Implement more examples that are comparable to the examples of scikit-fem,
  see https://scikit-fem.readthedocs.io/en/latest/listofexamples.html
- Add benchmarking of solving certain systems (and perhaps do equivalent
  benchmarking of other libraries for speed comparisons -- e.g., scikit-fem)
- Is computing `grad_phys` of `compute_physical_gradients` sometimes
  unneccesary work (see `assemble_load_vector`)?
