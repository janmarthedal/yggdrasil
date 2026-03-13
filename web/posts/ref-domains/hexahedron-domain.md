# The Hexahedron Domain

The reference domain for a hexahedral element is the unit cube $\hat{T} = [0,1]^3$. Every hexahedral element in a mesh is the image of $\hat{T}$ under a trilinear map $F_e$.

Just as the [quadrilateral rule](POSTROOT/ref-domains/quadrilateral-domain/) is a tensor product of two 1D Gauss–Legendre rules, the hexahedron rule is a tensor product of three. The derivation follows exactly the same reasoning, extended to a third coordinate direction. The $n^3$ quadrature points are all triples $(\hat{x}_i, \hat{x}_j, \hat{x}_k)$ of 1D nodes, with weights $w_i w_j w_k$. All weights are positive, sum to $1$, and the rule integrates all monomials $x^a y^b z^c$ with $a, b, c \leq p$ exactly.

The figure below shows the order-5 rule, which places $3 \times 3 \times 3 = 27$ points in a regular grid inside the reference cube.

![Tensor-product quadrature points on the reference hexahedron, order 5](MEDIAROOT/hexahedron-domain-quadrature.svg)

This construction is implemented in [`refdomains/hexahedron.py`](LIBROOT/refdomains/hexahedron.py) by `HexahedronDomain`, which delegates to `LineDomain` and forms the tensor product using `np.meshgrid`. The `quadrature(order)` method supports any polynomial order for which `LineDomain` has a rule.
