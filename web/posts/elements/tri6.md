# The Tri6 Element

The **Tri6** element is a 6-node quadratic triangle, implemented in
[`elements/tri6.py`](LIBROOT/elements/tri6.py) by the class `Tri6`. It lives on
the same reference triangle as [Tri3](POSTROOT/elements/tri3/) but adds a node
at the midpoint of each edge, giving three vertex nodes and three edge-midpoint
nodes:

$$\hat{x}_0 = (0,0), \quad \hat{x}_1 = (1,0), \quad \hat{x}_2 = (0,1),$$
$$\hat{x}_3 = \bigl(\tfrac{1}{2},0\bigr), \quad
  \hat{x}_4 = \bigl(\tfrac{1}{2},\tfrac{1}{2}\bigr), \quad
  \hat{x}_5 = \bigl(0,\tfrac{1}{2}\bigr).$$

The shape functions are built from the barycentric coordinates
$L_0 = 1 - \hat{x} - \hat{y}$, $L_1 = \hat{x}$, $L_2 = \hat{y}$ introduced in
the Tri3 post. The pattern mirrors the quadratic
[Lagrange basis polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial)
on $[0, 1]$: vertex functions use $L_i(2L_i - 1)$ and edge-midpoint functions
use $4L_i L_j$:

$$N_0 = L_0(2L_0 - 1), \quad N_1 = L_1(2L_1 - 1), \quad N_2 = L_2(2L_2 - 1),$$
$$N_3 = 4L_0 L_1, \quad N_4 = 4L_1 L_2, \quad N_5 = 4L_0 L_2.$$

Node $i$ of each edge is shared between the two adjacent vertex functions; the
edge-midpoint function reaches its maximum of one at the midpoint and vanishes at
both endpoints. The figure below shows all six nodes on the reference triangle.

![Nodes of the Tri6 element](MEDIAROOT/tri6-nodes.svg)

The gradients follow from the chain rule applied to the barycentric coordinates.
For the vertex functions,

$$\nabla N_0 = -(4L_0 - 1)\begin{pmatrix}1\\1\end{pmatrix}, \qquad
\nabla N_1 = (4L_1 - 1)\begin{pmatrix}1\\0\end{pmatrix}, \qquad
\nabla N_2 = (4L_2 - 1)\begin{pmatrix}0\\1\end{pmatrix},$$

and for the edge-midpoint functions,

$$\nabla N_3 = 4\begin{pmatrix}1 - 2\hat{x} - \hat{y}\\-\hat{x}\end{pmatrix}, \qquad
\nabla N_4 = 4\begin{pmatrix}\hat{y}\\\hat{x}\end{pmatrix}, \qquad
\nabla N_5 = 4\begin{pmatrix}-\hat{y}\\1 - \hat{x} - 2\hat{y}\end{pmatrix}.$$

The six shape functions form a partition of unity and span the full space of
quadratic polynomials on the triangle, so Tri6 can represent any quadratic
function exactly.
