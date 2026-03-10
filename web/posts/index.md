# Introduction

The finite element method (FEM) is one of the most widely used numerical techniques
for solving partial differential equations (PDEs). It underpins simulation software
across engineering, physics, and applied mathematics — from structural mechanics and
fluid dynamics to electrostatics and heat transfer.

This series of posts develops the mathematical foundation of FEM step by step, and
ties each concept directly to the [yggdrasil](https://github.com/janmarthedal/yggdrasil)
library implemented alongside these posts.

## What is the finite element method?

Given a PDE posed on some domain $\Omega \subset \mathbb{R}^d$, the finite element
method finds an approximate solution by:

1. **Reformulating** the PDE as an equivalent *weak* (variational) problem.
2. **Discretising** the domain into a *mesh* of small, simple elements (triangles,
   tetrahedra, hexahedra, etc.).
3. **Restricting** the infinite-dimensional function space to a
   finite-dimensional space of *piecewise polynomial* functions defined on the mesh.
4. **Assembling** a linear system whose solution gives the coefficients of the
   approximate solution in the chosen basis.

The method is extremely flexible: it handles complex geometries, higher-order
approximations, and a wide range of PDEs within a unified framework.

## The yggdrasil library

`yggdrasil` is a Python library for finite element analysis built on NumPy and SciPy.
It is intentionally kept small and transparent so that every piece of code corresponds
directly to a mathematical concept explained in these posts.

Key modules:

| Module | What it does |
|---|---|
| `yggdrasil/mesh.py` | Stores nodes, connectivity, and auxiliary mesh data |
| `yggdrasil/mesh_generators.py` | Creates structured meshes (triangles, tetrahedra) |
| `yggdrasil/elements/` | Reference element shape functions and gradients |
| `yggdrasil/domains/` | Reference domain definitions and quadrature rules |
| `yggdrasil/mapping.py` | Maps between reference and physical coordinates |
| `yggdrasil/boundary.py` | Extracts and tags boundary faces |
| `yggdrasil/assemble.py` | Assembles stiffness matrices, load vectors, and BCs |

## Post series overview

The posts follow the natural progression of the method:

1. **[The Poisson Problem](/01-poisson-problem/)** — a concrete PDE to build
   intuition: the strong form, boundary conditions, and well-posedness.
2. **[General Elliptic PDEs](/02-elliptic-pdes/)** — the abstract operator and
   bilinear-form setting that generalises Poisson.
3. **[Weak Formulation of the Poisson Problem](/03-poisson-weak-form/)** — deriving
   the variational form by multiplying by a test function and integrating by parts.
4. **[Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/)** — the abstract
   weak problem, coercivity, and the Lax–Milgram theorem.
5. **[The Discrete Formulation](/05-discrete-formulation/)** — the Galerkin method,
   basis functions, the linear system, and Céa's lemma.
6. Further posts cover elements, meshes, and assembly in detail.

## Prerequisites

The posts assume familiarity with:

- Multivariable calculus (gradient, divergence, integration).
- Basic linear algebra (matrices, linear systems).
- Elementary functional analysis is helpful but not required; key concepts are
  introduced as needed.

No prior knowledge of FEM is assumed.
