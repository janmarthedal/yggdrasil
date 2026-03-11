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
  - Generalise from Poisson to a broader class of second-order elliptic PDEs.
  - Introduce the abstract operator $L$ and the formulation $L u = f$.
  - Outline results related to well-posedness without going into details,
    but refer to authoratative resources -- preferably accessible online.

- `03-poisson-weak-form.md` — **Weak Formulation of the Poisson Problem**
  - Derive the weak (variational) formulation of the Poisson equation by
    multiplying by a test function and integrating by parts.
  - Show how Neumann conditions enter naturally as boundary integrals.
  - Link to Wikipedia for the concepts Green's first identity, Sobolev spaces,
    the Lax-Milgram theorem, and the Poincaré inequality.

- `04-elliptic-weak-form.md` — **Weak Formulation of Elliptic PDEs**
  - Extend the weak formulation to the general elliptic setting.
  - Introduce the abstract bilinear form a(u,v) and linear functional ℓ(v).

- `05-discrete-formulation.md` — **The Discrete Formulation**
  - Introduce the Galerkin method: restricting the weak problem to a
    finite-dimensional space spanned by basis functions.
  - Derive the resulting linear system.
  - Discuss the optimality of the discrete solution (Céa's lemma).
  