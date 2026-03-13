# Shape Functions on Reference Elements

The [previous posts](POSTROOT/ref-domains/overview/) showed how each element
type has a fixed reference domain $\hat{T}$ and how integrals over $\hat{T}$ are
evaluated by quadrature. The next ingredient is the set of **shape functions**
defined on $\hat{T}$: polynomial functions $N_1, \ldots, N_n : \hat{T} \to \mathbb{R}$
that form a basis for the local approximation space on the element. Each shape
function is associated with one node, and is equal to one at that node and zero
at all others.

The abstract class `ReferenceElement` in
[`elements/element.py`](LIBROOT/elements/element.py) captures this structure and
exposes three core methods. `node_coords` is a property returning an array of
shape `(num_nodes, topo_dim)` with the reference-domain coordinates of each
node. `shape_functions(xi)` takes an array `xi` of shape
`(num_points, topo_dim)` and returns an array of shape `(num_points, num_nodes)`;
the $i$-th column of the result is the $i$-th shape function evaluated at every
query point. `shape_function_gradients(xi)` returns an array of shape
`(num_points, num_nodes, topo_dim)`, where entry `[q, i, j]` is
$\partial N_i / \partial \hat{x}_j$ evaluated at query point $q$.

A key property of shape functions is the **partition of unity**: at every point
$\hat{x} \in \hat{T}$ the shape functions sum to one,

$$\sum_{i=1}^{n} N_i(\hat{x}) = 1.$$

This ensures that a constant function is represented exactly, which is a
necessary condition for consistency of the approximation. Gradients are not
evaluated in physical space here — transforming them from the reference domain to
a physical element requires the Jacobian of the element mapping $F_e$, which is
covered in a later assembly post.

Each element type has a dedicated post:

- 1D
  - [The Line2 Element](POSTROOT/elements/line2/)
  - [The Line3 Element](POSTROOT/elements/line3/)
- 2D
  - [The Tri3 Element](POSTROOT/elements/tri3/)
  - [The Tri6 Element](POSTROOT/elements/tri6/)
  - [The Quad4 Element](POSTROOT/elements/quad4/)
  - [The Quad9 Element](POSTROOT/elements/quad9/)
- 3D
  - [The Tet4 Element](POSTROOT/elements/tet4/)
  - [The Hex8 Element](POSTROOT/elements/hex8/)
