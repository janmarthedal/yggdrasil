# The Tetrahedron Domain

The reference domain for a tetrahedral element is the right tetrahedron $\hat{T}$ with vertices $(0,0,0)$, $(1,0,0)$, $(0,1,0)$, and $(0,0,1)$. It has volume $\tfrac{1}{6}$, and every tetrahedral element in a mesh is the image of $\hat{T}$ under an affine map $F_e$.

As with the triangle, there is no canonical quadrature family for the tetrahedron. [`TetrahedronDomain`](LIBROOT/refdomains/tetrahedron.py) uses tabulated symmetric rules from [Keast (1986)](https://doi.org/10.1016/0045-7825(86)90059-9) and other standard sources, supporting orders 1 through 3 with 1, 4, and 5 quadrature points respectively.

The figure below shows the point locations for each order. Marker size is proportional to the absolute weight, and colour indicates sign — blue for positive, red for negative.

![Quadrature points on the reference tetrahedron for orders 1–3](MEDIAROOT/tetrahedron-domain-quadrature.svg)

The order-1 rule places a single point at the centroid $(\tfrac{1}{4}, \tfrac{1}{4}, \tfrac{1}{4})$ with weight $\tfrac{1}{6}$. The order-2 rule uses 4 points arranged symmetrically near the four vertices, all with equal weight $\tfrac{1}{24}$. The order-3 rule adds a fifth point at the centroid with a negative weight, in the same spirit as the order-3 triangle rule. All weights sum to $\tfrac{1}{6}$, the volume of the reference tetrahedron.

The `quadrature(order)` method of `TetrahedronDomain` returns `(points, weights)` for the lowest-order rule that integrates polynomials of the requested degree exactly, up to a maximum of order 3.
