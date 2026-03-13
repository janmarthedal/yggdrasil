# The Quad9 Element

The **Quad9** element is a 9-node biquadratic quadrilateral, implemented in
[`elements/quad9.py`](LIBROOT/elements/quad9.py) by the class `Quad9`. It
lives on the same reference square as [Quad4](POSTROOT/elements/quad4/) but
adds a node at the midpoint of each edge and one at the center, giving four
corner nodes, four edge-midpoint nodes, and one interior node:

$$\hat{x}_0 = (0,0), \quad \hat{x}_1 = (1,0), \quad
  \hat{x}_2 = (1,1), \quad \hat{x}_3 = (0,1),$$
$$\hat{x}_4 = \bigl(\tfrac{1}{2},0\bigr), \quad
  \hat{x}_5 = \bigl(1,\tfrac{1}{2}\bigr), \quad
  \hat{x}_6 = \bigl(\tfrac{1}{2},1\bigr), \quad
  \hat{x}_7 = \bigl(0,\tfrac{1}{2}\bigr), \quad
  \hat{x}_8 = \bigl(\tfrac{1}{2},\tfrac{1}{2}\bigr).$$

![Nodes of the Quad9 element](MEDIAROOT/quad9-nodes.svg)

The shape functions are built as tensor products of the 1D quadratic [Lagrange
basis polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial) on
$[0,1]$ introduced in the [Line3](POSTROOT/elements/line3/) post:

$$\ell_0(t) = (1-t)(1-2t), \qquad \ell_1(t) = t(2t-1), \qquad
  \ell_2(t) = 4t(1-t).$$

Each Quad9 shape function is $N_i(\hat{x},\hat{y}) =
\ell_a(\hat{x})\,\ell_b(\hat{y})$ for the index pair $(a,b)$ that places
node $i$ at $(\hat{x}_a, \hat{y}_b)$:

| Node | $(\hat{x},\hat{y})$ | $N_i$ |
|------|---------------------|-------|
| 0 | $(0,0)$ | $\ell_0(\hat{x})\,\ell_0(\hat{y})$ |
| 1 | $(1,0)$ | $\ell_1(\hat{x})\,\ell_0(\hat{y})$ |
| 2 | $(1,1)$ | $\ell_1(\hat{x})\,\ell_1(\hat{y})$ |
| 3 | $(0,1)$ | $\ell_0(\hat{x})\,\ell_1(\hat{y})$ |
| 4 | $(\tfrac{1}{2},0)$ | $\ell_2(\hat{x})\,\ell_0(\hat{y})$ |
| 5 | $(1,\tfrac{1}{2})$ | $\ell_1(\hat{x})\,\ell_2(\hat{y})$ |
| 6 | $(\tfrac{1}{2},1)$ | $\ell_2(\hat{x})\,\ell_1(\hat{y})$ |
| 7 | $(0,\tfrac{1}{2})$ | $\ell_0(\hat{x})\,\ell_2(\hat{y})$ |
| 8 | $(\tfrac{1}{2},\tfrac{1}{2})$ | $\ell_2(\hat{x})\,\ell_2(\hat{y})$ |

The gradients follow directly from the product rule. For any node with
$N_i = \ell_a(\hat{x})\,\ell_b(\hat{y})$:

$$\frac{\partial N_i}{\partial \hat{x}} = \ell_a'(\hat{x})\,\ell_b(\hat{y}),
\qquad
\frac{\partial N_i}{\partial \hat{y}} = \ell_a(\hat{x})\,\ell_b'(\hat{y}),$$

where $\ell_0'(t) = 4t - 3$, $\ell_1'(t) = 4t - 1$, $\ell_2'(t) = 4 - 8t$
are the derivatives from Line3.

The nine functions form a partition of unity and span the full biquadratic
space $\mathbb{Q}_2 = \mathrm{span}\{1, \hat{x}, \hat{y}, \hat{x}^2,
\hat{x}\hat{y}, \hat{y}^2, \hat{x}^2\hat{y}, \hat{x}\hat{y}^2,
\hat{x}^2\hat{y}^2\}$, so Quad9 can represent any biquadratic function
exactly and achieves second-order accuracy. The face element is
[Line3](POSTROOT/elements/line3/), since each edge carries three nodes (two
corners plus the edge midpoint) that together form a quadratic 1D element.
