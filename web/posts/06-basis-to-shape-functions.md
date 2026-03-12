# From Basis to Shape Functions

The [discrete formulation](POSTROOT/05-discrete-formulation/) requires choosing
a basis $\{\phi_1, \ldots, \phi_N\}$ for the space $V_h$ and assembling the
stiffness matrix $K_{ij} = a(\phi_j, \phi_i)$ and load vector
$f_i = \ell(\phi_i)$. The full basis is needed to represent solutions with
non-zero Dirichlet boundary values; the restriction to the subspace $V_{h,0}$
(functions vanishing on the Dirichlet boundary) is handled during assembly and
is the subject of a later post. In practice these global basis functions are not
constructed directly; they emerge from simpler, locally defined functions called
**shape functions**.

Let the domain $\overline{\Omega}$ be partitioned into $M$ non-overlapping
elements $T_1, \ldots, T_M$ whose union covers $\bar{\Omega}$.
Each shape function is now defined locally on an element and has support within that element.

Each global basis function can be considered a linear combination of shape functions.
This fact together with the fact that $a(\cdot,\cdot)$ is bilinear means that it is
sufficient to consider $a(\cdot,\cdot)$ only for shape function that *belong to
the same element*.

Similarly for the load vector $f_i = \ell(\phi_i)$, it is sufficient to consider
$\ell$ for each shape function.
