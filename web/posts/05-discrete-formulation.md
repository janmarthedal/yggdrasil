# The Discrete Formulation

The weak problem — find $u \in V$ such that $a(u, v) = \ell(v)$ for all $v \in V$ —
lives in an infinite-dimensional function space and cannot be solved directly on a
computer. The *Galerkin method* replaces $V$ with a finite-dimensional subspace $V_h$,
turning the abstract variational problem into a concrete linear system.

## The Galerkin method

Let $V_h \subset V$ be a finite-dimensional subspace of dimension $N$, with basis
$\{\phi_1, \phi_2, \ldots, \phi_N\}$.

The **discrete (Galerkin) problem** is:

> **Find $u_h \in V_h$ such that $a(u_h, v_h) = \ell(v_h)$ for all $v_h \in V_h$.**

This is the same problem as before, but restricted to $V_h$.

Since $u_h \in V_h$, it can be written as a linear combination of the basis functions:

$$
u_h = \sum_{j=1}^N u_j \phi_j.
$$

Similarly, it suffices to test against each basis function $\phi_i$ in turn (by
linearity, if the equation holds for each basis function it holds for all of $V_h$).

Substituting:

$$
a\!\left(\sum_{j=1}^N u_j \phi_j,\, \phi_i\right) = \ell(\phi_i)
\quad i = 1, \ldots, N.
$$

Using bilinearity of $a$:

$$
\sum_{j=1}^N u_j\, a(\phi_j, \phi_i) = \ell(\phi_i)
\quad i = 1, \ldots, N.
$$

## The linear system

Define the **stiffness matrix** $K \in \mathbb{R}^{N \times N}$ and the
**load vector** $\mathbf{f} \in \mathbb{R}^N$ by

$$
K_{ij} = a(\phi_j, \phi_i),
\qquad
f_i = \ell(\phi_i).
$$

The discrete problem becomes the linear system:

$$
K \mathbf{u} = \mathbf{f},
$$

where $\mathbf{u} = (u_1, \ldots, u_N)^\top$ is the vector of unknown coefficients.

In `yggdrasil`, $K$ is assembled by `assemble_bilinear_form` and $\mathbf{f}$ by
`assemble_load_vector` (both in `yggdrasil/assemble.py`). The matrix is stored as a
sparse CSR matrix (`scipy.sparse.csr_matrix`) since, for FEM basis functions,
$K_{ij} \neq 0$ only when the supports of $\phi_i$ and $\phi_j$ overlap — i.e., when
nodes $i$ and $j$ belong to the same element.

## Finite element basis functions

The choice of $V_h$ and its basis $\{\phi_i\}$ defines the *finite element method*
proper.

Given a mesh of $\Omega$ — a partition into non-overlapping elements $\{K_e\}$ — the
FEM basis functions are:

- **Supported on a small patch of elements** (typically those sharing node $i$).
- **Polynomial on each element** of a chosen degree $p$.
- **Continuous across element boundaries** (for conforming $H^1$ elements).

The canonical choice for the Poisson problem is *continuous piecewise linear* functions
($p = 1$):

$$
\phi_i(\mathbf{x}) = \begin{cases}
  \text{linear on each element}, \\
  1 \text{ at node } i, \\
  0 \text{ at all other nodes.}
\end{cases}
$$

These are the *hat functions* (in 1D) or *tent functions* (in 2D/3D).

In `yggdrasil` the supported elements are:

| Element | Type | Degree |
|---|---|---|
| `Line2` | 1D linear | 1 |
| `Line3` | 1D quadratic | 2 |
| `Tri3` | 2D linear triangle | 1 |
| `Tri6` | 2D quadratic triangle | 2 |
| `Quad4` | 2D bilinear quadrilateral | 1 |
| `Quad9` | 2D biquadratic quadrilateral | 2 |
| `Tet4` | 3D linear tetrahedron | 1 |
| `Hex8` | 3D trilinear hexahedron | 1 |

All element classes live in `yggdrasil/elements/` and implement a common
`ReferenceElement` interface.

## Boundary conditions in the discrete system

### Dirichlet conditions

Dirichlet conditions $u = g_D$ on $\Gamma_D$ constrain specific DOFs directly.
Two strategies are supported in `yggdrasil/assemble.py`:

**Condensation** (`condense_dirichlet_bc`): partition the DOFs into free (interior)
and constrained (Dirichlet) sets. The constrained values are moved to the right-hand
side and the free DOFs are solved for in a reduced system.

**L² projection** (`project_dirichlet_bc`): find $g_h$ in the FE trace space on
$\Gamma_D$ that best approximates $g_D$ in the $L^2$ sense. Useful when $g_D$ is not
nodally interpolable (e.g., curved boundaries or high-order elements).

### Neumann conditions

As established in [Weak Formulation of the Poisson Problem](/03-poisson-weak-form/),
Neumann conditions enter naturally as boundary integrals on $\Gamma_N$:

$$
\ell_N(v_h) = \int_{\Gamma_N} g_N\, v_h \, ds.
$$

This is assembled by `assemble_neumann_bc` in `yggdrasil/assemble.py` using a boundary
mesh obtained from `extract_boundary` (in `yggdrasil/boundary.py`).

**Key ordering rule:** add Neumann contributions to $\mathbf{f}$ *before* calling
`condense_dirichlet_bc`.

## A complete assembly workflow

```python
from yggdrasil.mesh_generators import unit_square_tri_mesh
from yggdrasil.boundary import extract_boundary, tag_boundary_faces, select_boundary_faces
from yggdrasil.assemble import (
    assemble_bilinear_form, assemble_load_vector,
    assemble_neumann_bc, condense_dirichlet_bc,
)
import scipy.sparse.linalg as spla
import numpy as np

mesh = unit_square_tri_mesh(20)

def poisson_bilinear_form(N, grad_N):
    return grad_N @ grad_N.transpose(0, 2, 1)

def f(x):
    return 2 * np.pi**2 * np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])

K = assemble_bilinear_form(mesh, poisson_bilinear_form, quadrature_order=2)
b = assemble_load_vector(mesh, f, quadrature_order=2)

# Dirichlet BC: u = 0 on entire boundary
boundary = extract_boundary(mesh)
bc_dofs = boundary.point_data["original_node_index"]
system = condense_dirichlet_bc(K, b, bc_dofs, bc_val=0.0)

u_free = spla.spsolve(system.K, system.b)
u = system.reconstruct(u_free)
```

## Properties of the discrete solution

### Non-singularity of $K$

Under the same conditions that guarantee well-posedness of the continuous problem
(coercivity and continuity of $a$, see
[Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/)), the stiffness matrix
$K$ is non-singular. For the symmetric Poisson problem, $K$ is symmetric positive
definite after applying Dirichlet conditions.

### Optimality: Céa's lemma

The Galerkin solution $u_h$ is the *best approximation* to $u$ in $V_h$ with respect
to the *energy norm* $\|v\|_a = \sqrt{a(v,v)}$.

**Lemma (Céa, 1964).** Let $u$ be the exact solution and $u_h$ the Galerkin solution.
Then

$$
\|u - u_h\|_a \leq \frac{M}{\alpha}\, \inf_{w_h \in V_h} \|u - w_h\|_a,
$$

where $M$ is the continuity constant and $\alpha$ is the coercivity constant of $a$.

For symmetric $a$ (i.e., $a(u,v) = a(v,u)$), $M = \alpha$ and the Galerkin solution
is the *exact* best approximation (orthogonal projection) in the energy norm:

$$
\|u - u_h\|_a = \inf_{w_h \in V_h} \|u - w_h\|_a.
$$

### Galerkin orthogonality

A direct consequence: the error $e_h = u - u_h$ is orthogonal to $V_h$ in the $a$-inner
product:

$$
a(u - u_h, v_h) = 0 \quad \forall\, v_h \in V_h.
$$

*Proof:* $a(u, v_h) = \ell(v_h) = a(u_h, v_h)$ for all $v_h \in V_h$. $\square$

Galerkin orthogonality is the central identity used to prove convergence rates.

## Convergence

Céa's lemma reduces the question of how well $u_h$ approximates $u$ to the question
of how well $V_h$ can *represent* $u$ (the approximation problem). For degree-$p$
elements on a mesh with maximum element diameter $h$:

$$
\inf_{w_h \in V_h} \|u - w_h\|_{H^1} \leq C h^{\min(p, s)}\, |u|_{H^{s+1}(\Omega)},
$$

where $s$ characterises the smoothness of $u$. Combined with Céa's lemma this gives:

$$
\|u - u_h\|_{H^1} = O(h^p) \quad \text{as } h \to 0,
$$

provided $u$ is sufficiently smooth ($u \in H^{p+1}(\Omega)$). Refining the mesh
(smaller $h$) or increasing the polynomial degree $p$ both improve accuracy.

## What comes next

With the linear system derived, the remaining posts cover how the integrals defining
$K$ and $\mathbf{f}$ are evaluated in practice:

- **Elements** — reference elements, shape functions, and quadrature.
- **Mapping** — transforming integrals from reference to physical elements using
  Jacobians (`yggdrasil/mapping.py`).
- **Meshes** — how meshes are stored and boundaries extracted
  (`yggdrasil/mesh.py`, `yggdrasil/boundary.py`).
- **Assembly** — the full loop over elements that builds $K$ and $\mathbf{f}$
  (`yggdrasil/assemble.py`).
