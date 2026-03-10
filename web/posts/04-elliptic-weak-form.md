# Weak Formulation of Elliptic PDEs

In [Weak Formulation of the Poisson Problem](/03-poisson-weak-form/) we derived the
weak form of $-\Delta u = f$ by multiplying by a test function and integrating by
parts. Here we carry out the same procedure for the general elliptic operator
introduced in [General Elliptic PDEs](/02-elliptic-pdes/), and state the abstract
conditions — coercivity and continuity — that guarantee a unique solution via the
Lax–Milgram theorem.

## Deriving the weak form

Recall the strong problem:

$$
\begin{cases}
-\nabla \cdot (A \nabla u) + \mathbf{b} \cdot \nabla u + c\, u = f & \text{in } \Omega, \\
u = g_D & \text{on } \Gamma_D, \\
(A \nabla u) \cdot \mathbf{n} = g_N & \text{on } \Gamma_N.
\end{cases}
$$

Multiply by a test function $v \in H^1_{0,\Gamma_D}(\Omega)$ and integrate over
$\Omega$:

$$
-\int_\Omega \nabla \cdot (A \nabla u)\, v \, dx
+ \int_\Omega (\mathbf{b} \cdot \nabla u)\, v \, dx
+ \int_\Omega c\, u\, v \, dx
= \int_\Omega f\, v \, dx.
$$

Apply Green's first identity to the first term:

$$
-\int_\Omega \nabla \cdot (A \nabla u)\, v \, dx
= \int_\Omega (A \nabla u) \cdot \nabla v \, dx
  - \int_{\partial\Omega} (A \nabla u) \cdot \mathbf{n}\, v \, ds.
$$

Since $v = 0$ on $\Gamma_D$, the boundary integral over $\Gamma_D$ vanishes.
On $\Gamma_N$, the conormal flux condition $(A\nabla u) \cdot \mathbf{n} = g_N$ applies.

Assembling all terms gives the **weak formulation**:

**Find** $u \in H^1(\Omega)$ with $u = g_D$ on $\Gamma_D$ such that

$$
a(u, v) = \ell(v) \quad \forall\, v \in H^1_{0,\Gamma_D}(\Omega),
$$

where

$$
a(u, v)
= \int_\Omega (A \nabla u) \cdot \nabla v \, dx
+ \int_\Omega (\mathbf{b} \cdot \nabla u)\, v \, dx
+ \int_\Omega c\, u\, v \, dx,
$$

$$
\ell(v)
= \int_\Omega f\, v \, dx
+ \int_{\Gamma_N} g_N\, v \, ds.
$$

For Poisson ($A = I$, $\mathbf{b} = 0$, $c = 0$) this reduces exactly to the form
from [Weak Formulation of the Poisson Problem](/03-poisson-weak-form/).

## The abstract variational problem

Let $V = H^1_{0,\Gamma_D}(\Omega)$ (assuming homogeneous Dirichlet for simplicity;
non-homogeneous conditions are handled by a lifting). The abstract problem is:

> **Find $u \in V$ such that $a(u, v) = \ell(v)$ for all $v \in V$.**

This is a statement purely in terms of the bilinear form $a : V \times V \to \mathbb{R}$
and the linear functional $\ell : V \to \mathbb{R}$. The specific PDE and geometry
enter only through these two objects.

## Well-posedness: the Lax–Milgram theorem

To guarantee a unique solution, we need two properties of $a(\cdot,\cdot)$ on $V$.

### Continuity (boundedness)

There exists a constant $M > 0$ such that

$$
|a(u, v)| \leq M \|u\|_V \|v\|_V \quad \forall\, u, v \in V.
$$

Continuity says the bilinear form does not blow up as $u$ and $v$ grow.

For the general elliptic operator, continuity holds if $A$ is bounded (i.e.,
$\|A(\mathbf{x})\|_2 \leq \Lambda$ for some $\Lambda < \infty$), $\mathbf{b} \in
[L^\infty(\Omega)]^d$, and $c \in L^\infty(\Omega)$.

### Coercivity (V-ellipticity)

There exists a constant $\alpha > 0$ such that

$$
a(v, v) \geq \alpha \|v\|_V^2 \quad \forall\, v \in V.
$$

Coercivity says that $a(v,v)$ controls the full norm of $v$ from below. It is the
variational counterpart of the ellipticity condition on $A$.

For the Poisson problem on $V = H^1_0(\Omega)$:

$$
a(v, v) = \int_\Omega |\nabla v|^2 \, dx = \|\nabla v\|_{L^2}^2.
$$

By the *Poincaré inequality*, $\|\nabla v\|_{L^2} \geq C_P \|v\|_{H^1}$ for some
constant $C_P > 0$ depending only on $\Omega$. Hence coercivity holds with $\alpha = C_P^2$.

For the general operator, coercivity requires additional assumptions on $\mathbf{b}$
and $c$ (e.g., $c \geq 0$ and $\mathbf{b}$ not too large relative to $A$).

### The theorem

**Theorem (Lax–Milgram).** Let $V$ be a Hilbert space, $a : V \times V \to \mathbb{R}$
a continuous and coercive bilinear form, and $\ell : V \to \mathbb{R}$ a continuous
linear functional. Then there exists a unique $u \in V$ such that

$$
a(u, v) = \ell(v) \quad \forall\, v \in V,
$$

and the solution satisfies the stability estimate

$$
\|u\|_V \leq \frac{1}{\alpha} \|\ell\|_{V^*}.
$$

For symmetric $a$ (i.e., $a(u,v) = a(v,u)$, which holds when $\mathbf{b} = 0$ and
$A = A^\top$), the Lax–Milgram theorem reduces to the classical *Riesz representation
theorem*, and $u$ is the unique minimiser of the energy functional
$J(v) = \frac{1}{2} a(v,v) - \ell(v)$.

## The Poincaré inequality

The Poincaré inequality deserves special mention as it is the key estimate that allows
coercivity on $H^1_0(\Omega)$:

$$
\|v\|_{L^2(\Omega)} \leq C_P(\Omega)\, \|\nabla v\|_{L^2(\Omega)}
\quad \forall\, v \in H^1_0(\Omega).
$$

The constant $C_P$ depends on the domain $\Omega$ (roughly, on its diameter) but not
on $v$. This inequality holds precisely because functions in $H^1_0(\Omega)$ vanish on
the boundary — if there were no Dirichlet conditions, one could add a non-zero constant
to $v$ without changing $\|\nabla v\|$ but increasing $\|v\|$, violating the inequality.

## Symmetry and the energy interpretation

When $a$ is symmetric and coercive, the weak problem is equivalent to minimising the
*potential energy*:

$$
u = \operatorname*{arg\,min}_{w \in V_g} J(w),
\qquad
J(w) = \frac{1}{2} a(w, w) - \ell(w).
$$

For Poisson, this is

$$
J(w) = \frac{1}{2} \int_\Omega |\nabla w|^2 \, dx - \int_\Omega f w \, dx
        - \int_{\Gamma_N} g_N w \, ds,
$$

which has a clear physical interpretation (e.g., elastic potential energy minus work
done by external forces).

## Implementation in yggdrasil

The general bilinear form $a(u, v)$ for the elliptic operator is assembled by
`assemble_bilinear_form` in `yggdrasil/assemble.py`. The user provides a callable that
computes the element-level integrand.

For example, with $A = \kappa I$, $\mathbf{b} = 0$, $c = 0$:

```python
def bilinear_form(N, grad_N):
    # N:      (num_quad, nodes_per_elem)
    # grad_N: (num_quad, nodes_per_elem, spatial_dim)
    kappa = 2.0  # diffusion coefficient
    return kappa * (grad_N @ grad_N.transpose(0, 2, 1))
```

Adding a reaction term $c \cdot u \cdot v$:

```python
def bilinear_form(N, grad_N):
    c = 1.0
    diffusion = grad_N @ grad_N.transpose(0, 2, 1)
    reaction  = c * N[:, :, None] * N[:, None, :]
    return diffusion + reaction
```

The shape function values `N` have shape `(num_quad, nodes_per_elem)` and the outer
product `N[:, :, None] * N[:, None, :]` gives the `(num_quad, nodes_per_elem,
nodes_per_elem)` array representing $N_i N_j$ at each quadrature point.

## What comes next

With the abstract weak problem stated and well-posedness guaranteed by Lax–Milgram,
the next step is *discretisation*. The [Discrete Formulation](/05-discrete-formulation/)
restricts the infinite-dimensional problem to a finite-dimensional subspace spanned by
piecewise polynomial basis functions — this is the Galerkin method.
