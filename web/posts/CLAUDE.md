# Posts — Finite Element Method

## Purpose
This folder contains a series of posts introducing the finite element method,
accompanying the `yggdrasil` library in this repository.
The posts explain the mathematical theory and code constructs and reference the
accompanying library code in ./yggdrasil/ and possibly also example code in
./examples/.
*Every* function and class in the library are referenced from at least one post.
Some posts may be purely theoretical or regarding specific library such as
NumPy or SciPy.
The posts are prefixed with 01, 02, and so on, such that they form a linear
narrative.
Each post is of small or medium size and as self-contained as possible.
Plenty of examples, illustrations and animations are included.
The file `index.md` provides an introduction to and overview of the post series.

## Writing conventions
- Use markdown.
- Reference library source files by path relative to the repo root, e.g.
  `yggdrasil/assemble.py`.
- Math is written in LaTeX delimited by `$...$` (inline) and `$$...$$` (display).
- Mark incomplete sections with `**TODO**`.
- Prefer SVG over PNG for inline images when the illustration is fit for vector
  graphics.
- Consider using [three.js](https://threejs.org/) for 3D illustrations.

## Overall post structure
- The continuous problem.
  - Present a general formulation, but start with the Poisson problem as an example.
  - Introduce Dirichlet and Neumann boundary conditions.
  - Mention basic results on well-posedness of such problems.
- The weak formulation.
  - Show how to derive the weak formulation in the Poisson and general case.
- Discrete formulation
  - Formulation using basis function.
  - Optimality of solutions.
- Elements
  - Reference domains
    - Parameterization
    - Quadrature
  - Shape functions
    - Mapping from reference domain to physical domain
    - Jacobians
- Meshes
  - Representation
  - Boundaries
- Assembly
  - Global stiffness matrix
  - Load vector
  - Dirichlet boundary conditions
    - Condensation
    - L² projection
  - Neumann contributions
- Solving time-dependent systems
  - Heat equation
  - Wave equation
- External tools
  - Gmsh
  - Matplotlib
  - MeshIO
  - Paraview

## Table of Contents

- `index.md` — **Introduction**
  - Overview of the post series, the finite element method, and the yggdrasil library.

- `01-poisson-problem.md` — **The Poisson Problem**
  - Introduces the Poisson equation as a prototypical elliptic PDE. Defines Dirichlet
    and Neumann boundary conditions and states basic results on existence, uniqueness,
    and well-posedness for this specific case.

- `02-elliptic-pdes.md` — **General Elliptic PDEs**
  - Generalises from Poisson to a broader class of second-order elliptic PDEs.
    Introduces the abstract operator and bilinear-form setting, and discusses
    well-posedness via the Lax–Milgram theorem.

- `03-poisson-weak-form.md` — **Weak Formulation of the Poisson Problem**
  - Derives the weak (variational) formulation of the Poisson equation by multiplying
    by a test function and integrating by parts. Shows how Neumann conditions enter
    naturally as boundary integrals.

- `04-elliptic-weak-form.md` — **Weak Formulation of Elliptic PDEs**
  - Extends the weak formulation to the general elliptic setting. Introduces the
    abstract bilinear form $a(u, v)$ and linear functional $\ell(v)$, and states
    coercivity and continuity conditions needed for well-posedness.

- `05-discrete-formulation.md` — **The Discrete Formulation**
  - Introduces the Galerkin method: restricting the weak problem to a finite-dimensional
    space spanned by basis functions. Derives the resulting linear system and discusses
    the optimality of the discrete solution (Céa's lemma).
