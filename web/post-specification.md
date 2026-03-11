# Post and Media Specification

## Media

- `poisson-2d-domain.svg`
  - An svg illustrating a 2d domain.
  - See `./media-generation/poisson_hole.py` for domain specification.
  - Use proper notation for the relevant parts of the illustration.
- `poisson-2d-solution.png`
  - A solution to the problem shown by `poisson-2d-domain.svg`.
- `poisson-1d-solution.svg`
  - A solution to $-u''(x) = 1, \quad x \in (0, 1), \qquad u(0) = u(1) = 0.$

## Posts

- `index.md` — **Introduction**
  - Overview of the post series, the finite element method, and the yggdrasil library.

- `notation.md` — **Notation**
  - Table of commonly used notation throughout the posts

- `01-poisson-problem.md` — **The Poisson Problem**
  - Introduce the Poisson equation as a prototypical elliptic PDE.
  - Refer to [applications in physics and engineering](https://en.wikipedia.org/wiki/Poisson%27s_equation#Applications_in_physics_and_engineering)
  - Define Dirichlet and Neumann boundary conditions.
  - Provide a 2D example of domain and boundary types and refer to the
    illustration `poisson-2d-domain.svg`. Avoid detailed description of the
    domain and boundaries.
  - Include a solution to the 2D example and show the solution
    `poisson-2d-solution.png`
  - Describe a 1D Poisson problem with an analytical solution.
    Illustrate using `poisson-1d-solution.svg`
  - Avoid headings apart from the main title

- `02-elliptic-pdes.md` — **General Elliptic PDEs**

- `03-poisson-weak-form.md` — **Weak Formulation of the Poisson Problem**

- `04-elliptic-weak-form.md` — **Weak Formulation of Elliptic PDEs**

- `05-discrete-formulation.md` — **The Discrete Formulation**
