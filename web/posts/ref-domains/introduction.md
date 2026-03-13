# Reference Domains and Quadrature

Computing the element stiffness matrices and load vectors from the
[previous post](POSTROOT/discrete/basis-to-shape-functions/) requires evaluating
integrals over each element $T_e$. In a general mesh the elements can differ
in position, size, and orientation, so it is impractical to design shape
functions and integration rules separately for each one. Instead, every element
of a given type is treated as the image of a single fixed **reference domain**
$\hat{T}$ under a mapping $F_e : \hat{T} \to T_e$. Shape functions and
quadrature rules are defined once on $\hat{T}$ and then pulled back to the
physical element through $F_e$.

The reference domain for each element type is a standard, convenient geometric
shape. A line element maps to the interval $[0, 1]$, a triangle to the right
triangle with vertices $(0,0)$, $(1,0)$, $(0,1)$, and so on. The class
`ReferenceDomain` in [`refdomains/refdomain.py`](LIBROOT/refdomains/refdomain.py) represents such a
domain and exposes a `quadrature` method; concrete subclasses implement it for
each element type.

Integrals over $\hat{T}$ are evaluated numerically by **quadrature**: a rule of
order $p$ consists of $Q$ sample points $\hat{x}_1, \ldots, \hat{x}_Q \in \hat{T}$
and associated weights $w_1, \ldots, w_Q$ such that

$$\int_{\hat{T}} f(\hat{x})\, \mathrm{d}\hat{x} \approx \sum_{q=1}^{Q} w_q\, f(\hat{x}_q),$$

with the approximation being exact whenever $f$ is a polynomial of degree at
most $p$. The `quadrature(order)` method returns the pair `(points, weights)`,
where `points` has shape `(Q, topological_dimension)` and `weights` has shape
`(Q,)`. The caller specifies the required polynomial degree, and the
implementation returns a rule with enough points to integrate it exactly.

Each element type has a dedicated post describing its reference domain and quadrature rules:

- 1D
  - [The Line Domain](POSTROOT/ref-domains/line/)
- 2D
  - [The Triangle Domain](POSTROOT/ref-domains/triangle/)
  - [The Quadrilateral Domain](POSTROOT/ref-domains/quadrilateral/)
- 3D
  - [The Tetrahedron Domain](POSTROOT/ref-domains/tetrahedron/)
  - [The Hexahedron Domain](POSTROOT/ref-domains/hexahedron/)
