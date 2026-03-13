# The Tri3 Element

The **Tri3** element is the simplest triangular element: 3 nodes, one at each
vertex of the reference triangle, with linear shape functions. It is implemented
in [`elements/tri3.py`](LIBROOT/elements/tri3.py) by the class `Tri3`. The
reference domain is the right triangle $\hat{T}$ with vertices
$\hat{x}_0 = (0,0)$, $\hat{x}_1 = (1,0)$, $\hat{x}_2 = (0,1)$, matching the
[triangle reference domain](POSTROOT/ref-domains/triangle-domain/).

A natural coordinate system on $\hat{T}$ is provided by the **barycentric
coordinates** $L_0$, $L_1$, $L_2$, defined by

$$L_0(\hat{x}, \hat{y}) = 1 - \hat{x} - \hat{y}, \qquad
L_1(\hat{x}, \hat{y}) = \hat{x}, \qquad
L_2(\hat{x}, \hat{y}) = \hat{y}.$$

Each $L_i$ equals one at vertex $i$ and zero at the other two, and they satisfy
$L_0 + L_1 + L_2 = 1$ everywhere on $\hat{T}$. The Tri3 shape functions are
exactly the barycentric coordinates:

$$N_0 = L_0 = 1 - \hat{x} - \hat{y}, \qquad N_1 = L_1 = \hat{x}, \qquad N_2 = L_2 = \hat{y}.$$

The partition of unity $N_0 + N_1 + N_2 = 1$ follows immediately. Since all
three functions are linear, their reference-domain gradients are constant:

$$\nabla N_0 = \begin{pmatrix}-1\\-1\end{pmatrix}, \qquad
\nabla N_1 = \begin{pmatrix}1\\0\end{pmatrix}, \qquad
\nabla N_2 = \begin{pmatrix}0\\1\end{pmatrix}.$$

The figure below shows the three shape functions as filled contour plots over the
reference triangle.

![Shape functions of the Tri3 element](MEDIAROOT/tri3-shape-functions.svg)

Any linear function $f$ on $\hat{T}$ is recovered exactly as
$f = f(\hat{x}_0)\,N_0 + f(\hat{x}_1)\,N_1 + f(\hat{x}_2)\,N_2$. When the
physical element has vertices $x_0$, $x_1$, $x_2$ in $\mathbb{R}^d$, the same
combination $F_e(\hat{x}) = x_0 N_0 + x_1 N_1 + x_2 N_2$ gives the affine
isoparametric mapping from $\hat{T}$ to the physical triangle.
