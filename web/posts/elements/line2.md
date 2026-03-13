---
title: The Line2 Element
---

The simplest 1D element is the **Line2** element: a 2-node linear element on the
reference interval $[0, 1]$ with nodes at $\hat{x}_0 = 0$ and $\hat{x}_1 = 1$.
It is implemented in [`elements/line2.py`](LIBROOT/elements/line2.py) by the
class `Line2`.

The two shape functions are the linear [Lagrange basis polynomials](https://en.wikipedia.org/wiki/Lagrange_polynomial) on $[0, 1]$:

$$N_0(\hat{x}) = 1 - \hat{x}, \qquad N_1(\hat{x}) = \hat{x}.$$

Each is equal to one at its own node and zero at the other. Together they form a
partition of unity, $N_0 + N_1 = 1$, and any linear function on $[0, 1]$ can be
written as $f(\hat{x}) = f(0)\,N_0(\hat{x}) + f(1)\,N_1(\hat{x})$. The
reference-domain gradients are constant:

$$\frac{\mathrm{d}N_0}{\mathrm{d}\hat{x}} = -1, \qquad \frac{\mathrm{d}N_1}{\mathrm{d}\hat{x}} = 1.$$

![Shape functions of the Line2 element](MEDIAROOT/line2-shape-functions.svg)

When a physical line element has nodes at positions $x_0$ and $x_1$ in space,
the mapping $F_e(\hat{x}) = x_0 N_0(\hat{x}) + x_1 N_1(\hat{x}) = x_0 + (x_1 -
x_0)\hat{x}$, $0 \leq \hat{x} \leq 1$, carries the reference interval to the physical element. The same
shape functions therefore serve a dual role: they define both the approximation
space and the isoparametric geometry mapping.

