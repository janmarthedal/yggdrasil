---
title: The Quadrilateral Domain
---

The reference domain for a quadrilateral element is the unit square $\hat{T} = [0,1]^2$. Every quadrilateral element in a mesh is the image of $\hat{T}$ under a bilinear map $F_e$, so quadrature rules are again defined once on the reference domain.

Because $[0,1]^2$ is a product domain, a 2D quadrature rule can be built directly from the 1D Gauss–Legendre rule on $[0,1]$ by a **tensor-product** construction. Given a 1D rule with points $\hat{x}_1, \ldots, \hat{x}_n$ and weights $w_1, \ldots, w_n$, the 2D rule uses all $n^2$ pairs $(\hat{x}_i, \hat{x}_j)$ as quadrature points with weights $w_i w_j$:

$$\int_0^1\!\int_0^1 f(x, y)\, \mathrm{d}x\, \mathrm{d}y \approx \sum_{i=1}^n \sum_{j=1}^n w_i w_j\, f(\hat{x}_i, \hat{x}_j).$$

If the 1D rule integrates polynomials of degree $p$ exactly, the tensor-product rule integrates all polynomials of the form $x^a y^b$ with $a \leq p$ and $b \leq p$ exactly. This covers all bivariate polynomials of total degree at most $p$, and more. The weights $w_i w_j$ are all positive and sum to $1$, the area of the reference square.

The figure below shows the order-7 rule, which places a $4 \times 4$ grid of 16 points on the reference square. Marker size is proportional to the weight.

![Tensor-product quadrature points on the reference quadrilateral, order 7](MEDIAROOT/quadrilateral-domain-quadrature.svg)

In the library, this construction is implemented in [`refdomains/quadrilateral.py`](LIBROOT/refdomains/quadrilateral.py) by `QuadrilateralDomain`, which delegates to `LineDomain` internally and forms the tensor product using `np.meshgrid`. Its `quadrature(order)` method supports any polynomial order for which `LineDomain` has a rule.

