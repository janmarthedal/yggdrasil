# General Elliptic PDEs

The Poisson equation $-\Delta u = f$ is the prototypical elliptic PDE, but many
physical models require a more general operator. This post introduces the class of
*second-order linear elliptic PDEs* that the finite element framework in `yggdrasil`
is designed to handle.

## A more general operator

Replace the scalar Laplacian with a general second-order operator acting on
$u : \Omega \to \mathbb{R}$:

$$
Lu = -\nabla \cdot (A \nabla u) + \mathbf{b} \cdot \nabla u + c\, u,
$$

where:

- $A : \Omega \to \mathbb{R}^{d \times d}$ is a *diffusion tensor* (matrix-valued),
- $\mathbf{b} : \Omega \to \mathbb{R}^d$ is an *advection velocity* (vector-valued),
- $c : \Omega \to \mathbb{R}$ is a *reaction coefficient* (scalar).

The corresponding PDE is

$$
-\nabla \cdot (A \nabla u) + \mathbf{b} \cdot \nabla u + c\, u = f \quad \text{in } \Omega,
$$

with the same types of boundary conditions as for Poisson:

$$
\begin{cases}
u = g_D & \text{on } \Gamma_D, \\
(A \nabla u) \cdot \mathbf{n} = g_N & \text{on } \Gamma_N.
\end{cases}
$$

Note that the Neumann condition now involves the *conormal flux* $(A\nabla u) \cdot \mathbf{n}$
rather than the plain normal derivative.

## Special cases

| $A$ | $\mathbf{b}$ | $c$ | Name |
|---|---|---|---|
| $I$ (identity) | $0$ | $0$ | Poisson / Laplace |
| $\kappa(x) I$ | $0$ | $0$ | Variable-coefficient diffusion |
| diagonal $\kappa_{ij}$ | $0$ | $0$ | Anisotropic diffusion |
| $I$ | $\mathbf{b}(x)$ | $0$ | Advection-diffusion |
| $I$ | $0$ | $c(x) \geq 0$ | Reaction-diffusion |
| $I$ | $0$ | $-\lambda$ | Helmholtz (eigenvalue problem) |

The yggdrasil assembly routines handle the general case; the user supplies the
bilinear form as a Python callable (see `yggdrasil/assemble.py`).

## Ellipticity condition

The key property that makes $L$ well-behaved is *ellipticity* of the diffusion tensor
$A$.

**Uniform ellipticity.** There exists a constant $\alpha > 0$ such that for all
$\mathbf{x} \in \Omega$ and all $\boldsymbol{\xi} \in \mathbb{R}^d$:

$$
\boldsymbol{\xi}^\top A(\mathbf{x})\, \boldsymbol{\xi} \geq \alpha \|\boldsymbol{\xi}\|^2.
$$

This says the smallest eigenvalue of $A(\mathbf{x})$ is bounded away from zero
uniformly over $\Omega$. For $A = I$ (Poisson), $\alpha = 1$.

Ellipticity is the critical assumption that ensures:
- The associated bilinear form is *coercive* (see
  [Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/)).
- The discrete linear system is non-singular.
- The FEM approximation converges to the true solution.

Without ellipticity — for example in the pure advection problem $\mathbf{b} \cdot \nabla u = f$
— the analysis breaks down and completely different numerical methods are required.

## The abstract operator form

It is useful to think of the PDE abstractly as

$$
Lu = f \quad \text{in } \Omega,
$$

where $L : V \to V^*$ is a bounded linear operator between a Hilbert space $V$
(typically $H^1_{0,\Gamma_D}(\Omega)$) and its dual $V^*$.

The question "does the PDE have a unique solution?" becomes: "is $L$ an isomorphism?"
The Lax–Milgram theorem (see [Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/))
provides sufficient conditions — coercivity and boundedness of the associated bilinear
form — that guarantee this.

## Symmetry

When $\mathbf{b} = 0$ and $A$ is symmetric ($A = A^\top$), the operator $L$ is
*self-adjoint*:

$$
\int_\Omega (Lu) v \, dx = \int_\Omega u (Lv) \, dx \quad \forall u, v.
$$

In this case the associated bilinear form is symmetric and the global stiffness matrix
assembled by `assemble_bilinear_form` in `yggdrasil/assemble.py` is symmetric positive
definite. This allows the use of efficient solvers such as the conjugate gradient method.

The Poisson problem is symmetric; adding advection ($\mathbf{b} \neq 0$) breaks
symmetry.

## A two-dimensional example

On $\Omega = (0,1)^2$ with anisotropic diffusion

$$
A = \begin{pmatrix} 2 & 0 \\ 0 & 1 \end{pmatrix}, \quad \mathbf{b} = 0, \quad c = 0,
$$

the PDE becomes

$$
-2\frac{\partial^2 u}{\partial x^2} - \frac{\partial^2 u}{\partial y^2} = f.
$$

Diffusion is twice as fast in the $x$-direction. If we choose
$f(x,y) = (2\pi^2 + \pi^2)\sin(\pi x)\sin(\pi y) = 3\pi^2 \sin(\pi x)\sin(\pi y)$,
the exact solution is

$$
u(x,y) = \sin(\pi x)\sin(\pi y).
$$

## What comes next

With the strong form of both Poisson and the general elliptic PDE in hand, the next
step is to derive the *weak formulation*. The Poisson case is worked out first in
[Weak Formulation of the Poisson Problem](/03-poisson-weak-form/), and the general
case follows in [Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/).
