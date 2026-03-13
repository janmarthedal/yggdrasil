# The Triangle Domain

The reference domain for a triangular element is the right triangle $\hat{T}$ with vertices $(0,0)$, $(1,0)$, and $(0,1)$. It has area $\tfrac{1}{2}$, and every triangular element in a mesh is the image of $\hat{T}$ under an affine map $F_e$.

Unlike the interval, where Gauss–Legendre quadrature provides an optimal family of rules for any order, the triangle has no single canonical quadrature construction. Instead, [`TriangleDomain`](LIBROOT/refdomains/triangle.py) uses tabulated **symmetric** rules whose points and weights respect the three-fold symmetry of the equilateral triangle. Rules for orders 1 through 5 are included, taken from [Dunavant (1985)](https://onlinelibrary.wiley.com/doi/10.1002/nme.1620210612) and other standard sources. The number of quadrature points for each order is 1, 3, 4, 6, and 7.

The figure below shows the point locations for each order. Marker size is proportional to the absolute weight, and colour indicates sign — blue for positive, red for negative.

![Quadrature points on the reference triangle for orders 1–5](MEDIAROOT/triangle-domain-quadrature.svg)

The order-3 rule is notable: it uses 4 points but one weight is negative. Quadrature rules with negative weights are mathematically valid — the signed sum still converges to the integral — but they can amplify rounding errors, so they are used only where no equally efficient all-positive rule exists. The weights of all rules sum to $\tfrac{1}{2}$, the area of the reference triangle.

The `quadrature(order)` method of `TriangleDomain` returns `(points, weights)` for the lowest-order rule that integrates polynomials of the requested degree exactly, up to a maximum of order 5.
