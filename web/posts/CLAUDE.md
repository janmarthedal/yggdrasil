# Posts — Finite Element Method

## Purpose

These markdown files form a series of posts introducing the finite element method,
developed alongside the yggdrasil library in this repository. Posts explain the
mathematical theory and reference the accompanying library code directly.

## Publication workflow

When a post is ready to publish, copy it to the companion repository:

```
/Users/jamr/repos/janmarthedal/janmr.com/content/posts/finite-element-method/
```

Update the index there if needed. Do not edit posts in the `janmr.com` repository
directly — always work here and copy over.

## Writing conventions

- Reference library source files by path relative to the repo root, e.g. `yggdrasil/assemble.py`.
- Math is written in LaTeX delimited by `$...$` (inline) and `$$...$$` (display).
- Mark incomplete sections with `**TODO**`.

## Table of Contents

- [Introduction](index.md)
- [The Poisson Problem](01-the-poisson-problem.md)
- [A More General Formulation](02-a-more-general-formulation.md)
- [Weak Formulation of the Poisson Problem](03-weak-formulation-of-poisson-problem.md)
- Weak Formulation of General Formulation
- A finite solution space
- Elements
  - 1D
    - Line2
    - Line3
  - 2D
    - Triangle (Tri3, Tri6)
    - Quad (Quad4, Quad9)
  - 3D
    - Tetrahedron (Tet4)
    - Hexahedron (Hex8)
- Meshes
  - Representation
  - Well-posedness
  - Boundaries
    - Normals
  - Third-party tools
    - gmsh
    - meshio
- Quadrature
  - Line
  - Triangle
- Assembly
  - Stiffness matrix
  - Right-hand side
  - Boundary conditions
    - Dirichlet
    - Neumann
- Solving
  - Direct methods
  - Iterative methods
- Time-dependent systems
  - Solvers
