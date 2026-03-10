# Weak Formulation of the Poisson Problem

The *strong form* of the Poisson problem requires the solution $u$ to be twice
differentiable everywhere in $\Omega$. This is a restrictive requirement: many source
terms $f$ and domain geometries do not yield such smooth solutions. The *weak*
(or *variational*) formulation is an equivalent reformulation that requires only one
derivative from $u$ — exactly the right amount for finite element approximation.

## Multiplying by a test function

Start with the strong form

$$
-\Delta u = f \quad \text{in } \Omega.
$$

Choose an arbitrary *test function* $v \in H^1_{0,\Gamma_D}(\Omega)$ — any function
that is in $H^1(\Omega)$ and vanishes on the Dirichlet boundary $\Gamma_D$. Multiply
both sides of the PDE by $v$ and integrate over $\Omega$:

$$
-\int_\Omega (\Delta u)\, v \, dx = \int_\Omega f\, v \, dx.
$$

If this equation holds for *all* such $v$, it is equivalent to the original PDE (in an
appropriate function-space sense).

## Integration by parts (Green's first identity)

The left-hand side contains second derivatives of $u$, which we want to avoid.
Green's first identity states:

$$
-\int_\Omega (\Delta u)\, v \, dx
= \int_\Omega \nabla u \cdot \nabla v \, dx
  - \int_{\partial\Omega} (\nabla u \cdot \mathbf{n})\, v \, ds.
$$

This is exactly the $d$-dimensional integration-by-parts formula
$\int_\Omega (\nabla \cdot \mathbf{F}) v = -\int_\Omega \mathbf{F} \cdot \nabla v + \int_{\partial\Omega} \mathbf{F} \cdot \mathbf{n} \, v$
with $\mathbf{F} = \nabla u$.

Substituting:

$$
\int_\Omega \nabla u \cdot \nabla v \, dx
- \int_{\partial\Omega} (\nabla u \cdot \mathbf{n})\, v \, ds
= \int_\Omega f\, v \, dx.
$$

## Incorporating boundary conditions

Split the boundary integral over $\Gamma_D$ and $\Gamma_N$:

$$
\int_{\partial\Omega} (\nabla u \cdot \mathbf{n})\, v \, ds
= \int_{\Gamma_D} (\nabla u \cdot \mathbf{n})\, v \, ds
+ \int_{\Gamma_N} (\nabla u \cdot \mathbf{n})\, v \, ds.
$$

Since $v = 0$ on $\Gamma_D$ by the choice of test space, the first integral vanishes.
On $\Gamma_N$, the Neumann condition says $\nabla u \cdot \mathbf{n} = g_N$, so the
second integral becomes $\int_{\Gamma_N} g_N v \, ds$.

The equation becomes:

$$
\int_\Omega \nabla u \cdot \nabla v \, dx
= \int_\Omega f\, v \, dx + \int_{\Gamma_N} g_N\, v \, ds.
$$

This is the **weak formulation** of the Poisson problem.

## The weak problem

**Find** $u \in H^1(\Omega)$ with $u = g_D$ on $\Gamma_D$ such that

$$
\int_\Omega \nabla u \cdot \nabla v \, dx
= \int_\Omega f\, v \, dx + \int_{\Gamma_N} g_N\, v \, ds
\quad \forall\, v \in H^1_{0,\Gamma_D}(\Omega).
$$

### Key observations

1. **Reduced regularity.** Both $u$ and $v$ now need only one weak derivative — they
   only need to be in $H^1(\Omega)$, not $H^2(\Omega)$.

2. **Neumann conditions appear naturally.** They enter as a boundary integral on the
   right-hand side. This is why Neumann conditions are called *natural* boundary
   conditions in the FEM context.

3. **Dirichlet conditions are essential.** They constrain the function space itself
   and must be enforced explicitly. This is why they are called *essential* boundary
   conditions.

4. **Equivalence.** If $u$ is a strong solution, it is also a weak solution. If a
   weak solution happens to be twice differentiable, it is a strong solution. For
   smooth data and domains, the two formulations are equivalent.

## Homogeneous Dirichlet conditions

When $g_D = 0$ on all of $\partial\Omega$ (i.e., $\Gamma_N = \emptyset$ and $g_D = 0$),
the trial and test spaces coincide: both equal $H^1_0(\Omega)$. The weak problem
simplifies to:

**Find** $u \in H^1_0(\Omega)$ such that

$$
\int_\Omega \nabla u \cdot \nabla v \, dx = \int_\Omega f\, v \, dx
\quad \forall\, v \in H^1_0(\Omega).
$$

This is the form most commonly seen in textbooks and the starting point for
well-posedness analysis.

## Bilinear and linear forms

Introduce the notation:

$$
a(u, v) = \int_\Omega \nabla u \cdot \nabla v \, dx,
\qquad
\ell(v) = \int_\Omega f\, v \, dx + \int_{\Gamma_N} g_N\, v \, ds.
$$

Then the weak problem reads:

$$
\text{find } u \in V_g \text{ such that } a(u, v) = \ell(v) \quad \forall\, v \in V_0,
$$

where $V_g = \{ w \in H^1(\Omega) : w = g_D \text{ on } \Gamma_D \}$ and
$V_0 = H^1_{0,\Gamma_D}(\Omega)$.

$a(\cdot,\cdot)$ is a *bilinear form* — linear in each argument separately.
$\ell(\cdot)$ is a *linear functional* (or *linear form*).

This abstract structure is exactly what is generalised in
[Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/) and implemented in
`yggdrasil/assemble.py`.

## Implementation in yggdrasil

In `yggdrasil/assemble.py`, `assemble_bilinear_form` evaluates the discrete version of
$a(u, v)$ and `assemble_load_vector` evaluates the discrete version of the volume part
of $\ell(v)$. The Neumann term $\int_{\Gamma_N} g_N v \, ds$ is assembled separately
by `assemble_neumann_bc` using a boundary mesh extracted with `extract_boundary` from
`yggdrasil/boundary.py`.

**Important:** Neumann contributions must be added to the right-hand side *before*
calling `condense_dirichlet_bc`, which modifies the system to enforce Dirichlet
conditions.

For the Poisson bilinear form, the integrand is $\nabla u \cdot \nabla v$, which in
`yggdrasil` is expressed as:

```python
def poisson_bilinear_form(N, grad_N):
    # grad_N: (num_quad, nodes_per_elem, spatial_dim)
    # returns: (num_quad, nodes_per_elem, nodes_per_elem)
    return grad_N @ grad_N.transpose(0, 2, 1)
```

## What comes next

The weak formulation of Poisson extends naturally to the general elliptic operator
from [General Elliptic PDEs](/02-elliptic-pdes/). This is the subject of
[Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/), where we also state
the conditions — coercivity and continuity — that guarantee well-posedness via the
Lax–Milgram theorem.
