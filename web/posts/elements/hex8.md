---
title: The Hex8 Element
---

The **Hex8** element is the simplest hexahedral element: 8 nodes, one at each
corner of the reference cube, with trilinear shape functions. The
reference domain is $\hat{H} = [0,1]^3$, matching the
[hexahedron reference domain](POSTROOT/ref-domains/hexahedron/). The 8 nodes
are ordered bottom face first, then top face, each traversed in the same
cyclic order:

$$\hat{x}_0 = (0,0,0), \quad \hat{x}_1 = (1,0,0), \quad
  \hat{x}_2 = (1,1,0), \quad \hat{x}_3 = (0,1,0),$$
$$\hat{x}_4 = (0,0,1), \quad \hat{x}_5 = (1,0,1), \quad
  \hat{x}_6 = (1,1,1), \quad \hat{x}_7 = (0,1,1).$$

In the library, it is implemented in [`elements/hex8.py`](LIBROOT/elements/hex8.py)
by the class `Hex8`.

![Nodes of the Hex8 element](MEDIAROOT/hex8-nodes.svg)

The shape functions extend the pattern of [Quad4](POSTROOT/elements/quad4/)
from two to three dimensions. Writing $\ell_0(t) = 1 - t$ and $\ell_1(t) = t$
for the linear Lagrange basis on $[0,1]$, each Hex8 shape function is a triple
tensor product

$$N_i(\hat{x}, \hat{y}, \hat{z}) =
  \ell_a(\hat{x})\,\ell_b(\hat{y})\,\ell_c(\hat{z}),$$

where $(a, b, c) \in \{0,1\}^3$ is the index triple that places node $i$ at
$(\hat{x}_a, \hat{y}_b, \hat{z}_c)$. Written out explicitly:

$$N_0 = (1-\hat{x})(1-\hat{y})(1-\hat{z}), \quad
  N_1 = \hat{x}(1-\hat{y})(1-\hat{z}), \quad
  N_2 = \hat{x}\,\hat{y}(1-\hat{z}), \quad
  N_3 = (1-\hat{x})\hat{y}(1-\hat{z}),$$
$$N_4 = (1-\hat{x})(1-\hat{y})\hat{z}, \quad
  N_5 = \hat{x}(1-\hat{y})\hat{z}, \quad
  N_6 = \hat{x}\,\hat{y}\,\hat{z}, \quad
  N_7 = (1-\hat{x})\hat{y}\,\hat{z}.$$

Each $N_i$ equals one at its own corner and zero at all other seven. The eight
functions form a partition of unity. The gradients follow from the product rule;
for node 0 for example,

$$\nabla N_0 = \begin{pmatrix}
  -(1-\hat{y})(1-\hat{z})\\
  -(1-\hat{x})(1-\hat{z})\\
  -(1-\hat{x})(1-\hat{y})
\end{pmatrix},$$

and similarly for the remaining nodes with appropriate sign changes on each
factor. As with Quad4, the gradients depend on position because of the cross
terms in the trilinear functions.

Hex8 spans the trilinear space $\mathbb{Q}_1$ — all polynomials of degree at
most one in each variable separately. This includes all linear functions, so
the element achieves first-order accuracy. The face element is
[Quad4](POSTROOT/elements/quad4/), since each of the six square faces carries
four corner nodes forming a bilinear quadrilateral element.

