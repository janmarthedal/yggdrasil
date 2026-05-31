---
title: The Line3 Element
---

The **Line3** element is a 3-node quadratic element on $[0, 1]$. The nodes
are placed at $\hat{x}_0 = 0$, $\hat{x}_1 = 1$, and $\hat{x}_2 = \tfrac{1}{2}$:
the first two sit at the endpoints, as in Line2; the third is the midpoint.

In the library, it is implemented in [`elements/line3.py`](LIBROOT/elements/line3.py)
by the class `Line3`.

The three shape functions are the quadratic [Lagrange basis polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial) associated
with these nodes:

$$N_0(\hat{x}) = (1 - \hat{x})(1 - 2\hat{x}), \qquad
N_1(\hat{x}) = \hat{x}(2\hat{x} - 1), \qquad
N_2(\hat{x}) = 4\hat{x}(1 - \hat{x}).$$

Each satisfies $N_i(\hat{x}_j) = \delta_{ij}$. The endpoint functions $N_0$ and
$N_1$ are negative in the interior of the interval — $N_0$ dips below zero near
$\hat{x} = 1$ and $N_1$ near $\hat{x} = 0$ — while the midpoint function $N_2$
has the characteristic bubble shape, reaching its maximum value of one at
$\hat{x} = \tfrac{1}{2}$. The three functions sum to one everywhere.

The reference-domain gradients are:

$$\frac{\mathrm{d}N_0}{\mathrm{d}\hat{x}} = 4\hat{x} - 3, \qquad
\frac{\mathrm{d}N_1}{\mathrm{d}\hat{x}} = 4\hat{x} - 1, \qquad
\frac{\mathrm{d}N_2}{\mathrm{d}\hat{x}} = 4 - 8\hat{x}.$$

![Shape functions of the Line3 element](MEDIAROOT/line3-shape-functions.svg)

Because Line3 uses the same quadratic Lagrange basis as the 1D building blocks
of the [Quad9](POSTROOT/elements/quad9/) and [Hex8](POSTROOT/elements/hex8/)
elements, its shape functions also appear as factors in the tensor-product
constructions for those higher-dimensional elements.

