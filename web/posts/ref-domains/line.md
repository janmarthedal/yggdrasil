# The Line Domain

The reference domain for a 1D element is the unit interval $\hat{T} = [0, 1]$. Every line element in a mesh is the image of $\hat{T}$ under an affine map $F_e$, so quadrature rules need only be defined once on $[0, 1]$ and then applied to any element by a change of variables.

The standard choice for quadrature on $[0, 1]$ is **Gauss–Legendre quadrature**. An $n$-point Gauss–Legendre rule is the unique rule that integrates polynomials of degree up to $2n - 1$ exactly using only $n$ function evaluations — the maximum possible precision for $n$ points. The classical rule is defined on $[-1, 1]$ and computed by [`numpy.polynomial.legendre.leggauss`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html); the nodes and weights are then mapped to $[0, 1]$ by a linear change of variables.

The figure below shows the quadrature points on $[0, 1]$ for $n = 1, 2, 3, 4$ points, integrating polynomials of degree up to $1, 3, 5, 7$ exactly. The marker size is proportional to the weight.

![Gauss–Legendre quadrature points on [0,1] for n = 1, 2, 3, 4](MEDIAROOT/line-domain-quadrature.svg)

Given a required polynomial degree $p$, the number of points needed is $n = \lceil (p+1)/2 \rceil$. This is implemented in [`refdomains/line.py`](LIBROOT/refdomains/line.py) by the class `LineDomain`, whose `quadrature(order)` method returns `(points, weights)` for the rule that integrates polynomials of degree `order` exactly.
