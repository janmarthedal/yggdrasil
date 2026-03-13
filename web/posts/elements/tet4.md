# The Tet4 Element

The **Tet4** element is the simplest tetrahedral element: 4 nodes, one at each
vertex of the reference tetrahedron, with linear shape functions. It is
implemented in [`elements/tet4.py`](LIBROOT/elements/tet4.py) by the class
`Tet4`. The reference domain is the right tetrahedron $\hat{T}$ with vertices

$$\hat{x}_0 = (0,0,0), \quad \hat{x}_1 = (1,0,0), \quad
  \hat{x}_2 = (0,1,0), \quad \hat{x}_3 = (0,0,1),$$

matching the
[tetrahedron reference domain](POSTROOT/ref-domains/tetrahedron/).

![Nodes of the Tet4 element](MEDIAROOT/tet4-nodes.svg)

Just as the [Tri3](POSTROOT/elements/tri3/) element uses barycentric
coordinates on the triangle, the Tet4 element uses their three-dimensional
analogue. The four barycentric coordinates on $\hat{T}$ are

$$L_0 = 1 - \hat{x} - \hat{y} - \hat{z}, \qquad
  L_1 = \hat{x}, \qquad L_2 = \hat{y}, \qquad L_3 = \hat{z},$$

and they satisfy $L_0 + L_1 + L_2 + L_3 = 1$ everywhere on $\hat{T}$.
Each $L_i$ equals one at vertex $i$ and zero at the other three, so the Tet4
shape functions are simply

$$N_0 = 1 - \hat{x} - \hat{y} - \hat{z}, \qquad
  N_1 = \hat{x}, \qquad N_2 = \hat{y}, \qquad N_3 = \hat{z}.$$

Since all four functions are linear their reference-domain gradients are
constant:

$$\nabla N_0 = \begin{pmatrix}-1\\-1\\-1\end{pmatrix}, \quad
  \nabla N_1 = \begin{pmatrix}1\\0\\0\end{pmatrix}, \quad
  \nabla N_2 = \begin{pmatrix}0\\1\\0\end{pmatrix}, \quad
  \nabla N_3 = \begin{pmatrix}0\\0\\1\end{pmatrix}.$$

Any linear function $f$ on $\hat{T}$ is recovered exactly as $f =
\sum_{i=0}^3 f(\hat{x}_i)\,N_i$. When the physical element has vertices
$x_0, x_1, x_2, x_3 \in \mathbb{R}^3$, the isoparametric mapping
$F_e(\hat{x}) = \sum_{i=0}^3 x_i\,N_i$ is an affine map from $\hat{T}$ to
the physical tetrahedron, and the Jacobian $J = \partial F_e / \partial
\hat{x}$ is constant within each element — a useful property it shares with
Tri3 and Line2. The face element is [Tri3](POSTROOT/elements/tri3/), since
each of the four triangular faces carries three nodes that form a linear
triangular element.
