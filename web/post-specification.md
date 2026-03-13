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
  - Short introduction to the finite element method
  - Short introduction to the yggdrasil library
  - Overview of the post series with a complete Table of Contents and links
    to all posts

- `notation.md` — **Notation**
  - Table of commonly used notation throughout the posts

- `continuous/poisson-problem.md` — **The Poisson Problem**
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

- `continuous/elliptic-pdes.md` — **General Elliptic PDEs**
  - Generalise from Poisson to a broader class of second-order elliptic PDEs.
  - Introduce the abstract operator $L$ and the formulation $L u = f$.
  - Outline results related to well-posedness without going into details,
    but refer to authoratative resources -- preferably accessible online.

- `continuous/poisson-weak-form.md` — **Weak Formulation of the Poisson Problem**
  - Derive the weak (variational) formulation of the Poisson equation by
    multiplying by a test function and integrating by parts.
  - Show how Neumann conditions enter naturally as boundary integrals.
  - Link to Wikipedia for the concepts Green's first identity, Sobolev spaces,
    the Lax-Milgram theorem, and the Poincaré inequality.

- `continuous/elliptic-weak-form.md` — **Weak Formulation of Elliptic PDEs**
  - Extend the weak formulation to the general elliptic setting.
  - Introduce the abstract bilinear form a(u,v) and linear functional ℓ(v).

- `discrete/discrete-formulation.md` — **The Discrete Formulation**
  - Introduce the Galerkin method: restricting the weak problem to a
    finite-dimensional space spanned by basis functions.
  - Derive the resulting linear system.
  - Discuss the optimality of the discrete solution (Céa's lemma).
  
- `discrete/basis-to-shape-functions.md` — **From Basis to Shape Functions**
  - Partition the domain into elements and define shape functions as locally
    supported functions on each element.
  - Explain how the global space $V_h$ is assembled from local shape functions,
    noting that both continuous and discontinuous Galerkin methods fit this
    framework.
  - Use bilinearity of $a(\cdot,\cdot)$ to decompose the stiffness matrix into
    element contributions, and note that sparsity follows from the local support
    of the basis functions.
  - Introduce element stiffness matrices and element load vectors as the objects
    computed during assembly.
  - Avoid details and introduce as little new notation as possible. The details
    will be described later in dedicated assembly posts.

- `domains/introduction.md` — **Reference Domains and Quadrature**
  - Introduce the need for local parametrization of standard elements (domains).
  - Introduce the concept of integration using sample points and weights.
  - Refer to the class ReferenceDomain.quadrature.
  
- `domains/line-domain.md` — **The Line Domain**
  - Link to [Gauss-Legendre quadrature](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html)
  - Refer to the class LineDomain
  
- `domains/triangle-domain.md` — **The Triangle Domain**
  - Refer to the class TriangleDomain

- `domains/quadrilateral-domain.md` — **The Quadrilateral Domain**
  - Refer to the class QuadrilateralDomain

- `domains/tetrahedron-domain.md` — **The Tetrahedron Domain**
  - Refer to the class TetrahedronDomain

- `domains/hexahedron-domain.md` — **The Hexahedron Domain**
  - Refer to the class HexahedronDomain
