---
title: From Basis to Shape Functions
---

The [discrete formulation](POSTROOT/discrete/discrete-formulation/) requires choosing
a basis $\{\phi_1, \ldots, \phi_N\}$ for the space $V_h$ and assembling the
stiffness matrix $K_{ij} = a(\phi_j, \phi_i)$ and load vector $f_i = \ell(\phi_i)$.
In practice these global basis functions are not constructed directly; they emerge
from simpler, locally defined functions called **shape functions**.

Partition $\bar{\Omega}$ into $M$ non-overlapping elements $T_1, \ldots, T_M$
whose interiors are disjoint and whose union covers $\bar{\Omega}$. On each
element $T_e$ a small set of **shape functions** is defined, each supported
entirely within $T_e$. The global space $V_h$ is then assembled from these local
pieces — either by identifying shape functions on adjacent elements along shared
boundaries (as in the standard continuous Galerkin method) or by keeping them
independent (as in the discontinuous Galerkin method). Either way, each global
basis function $\phi_i$ has support on at most a small number of elements.

Since $a(\cdot, \cdot)$ is bilinear, the stiffness matrix entry decomposes as

$$K_{ij} = a(\phi_j, \phi_i) = \sum_{e=1}^{M} a\big|_{T_e}(\phi_j, \phi_i),$$

where $a\big|_{T_e}$ denotes the restriction of $a$ to element $T_e$. The
contribution from $T_e$ is zero unless the supports of both $\phi_i$ and $\phi_j$
overlap with $T_e$, which can only happen for nearby basis functions. This is the
origin of the sparsity of the stiffness matrix. Assembly therefore proceeds
element by element: for each $T_e$ one computes a small **element stiffness
matrix** and scatters its entries into the global $K$. The load vector is
assembled the same way, each element contributing a small **element load vector**
to $\mathbf{f}$.

