# The Poisson Problem

The Poisson equation is one of the simplest and most studied PDEs. It is the natural
starting point for finite element analysis because it is rich enough to illustrate all
the key ideas — boundary conditions, well-posedness, weak formulation — while keeping
the algebra transparent.

## The strong form

Let $\Omega \subset \mathbb{R}^d$ ($d = 1, 2,$ or $3$) be a bounded domain with
boundary $\partial\Omega$. The *Poisson problem* is: find $u : \Omega \to \mathbb{R}$
such that

$$
-\Delta u = f \quad \text{in } \Omega,
$$

where $\Delta u = \nabla \cdot \nabla u = \sum_{i=1}^{d} \frac{\partial^2 u}{\partial x_i^2}$
is the Laplacian and $f : \Omega \to \mathbb{R}$ is a given source term.

The minus sign is conventional; it makes the associated operator positive (semi-)definite.

### Physical interpretations

The Poisson equation models many physical phenomena:

| $u$ | $f$ | Setting |
|---|---|---|
| Electrostatic potential | Charge density / $\varepsilon_0$ | Electrostatics |
| Temperature | Heat source | Steady-state heat conduction |
| Pressure head | Recharge rate | Groundwater flow (Darcy) |
| Displacement | Body force (1D rod) | Linear elasticity (1D) |

## Boundary conditions

The PDE alone does not uniquely determine $u$; boundary conditions are needed on
$\partial\Omega$.

The boundary $\partial\Omega$ is split into two disjoint parts,
$\partial\Omega = \Gamma_D \cup \Gamma_N$, on which different conditions are imposed.

### Dirichlet boundary condition

On $\Gamma_D$ (also called the *essential* boundary), the solution value is prescribed:

$$
u = g_D \quad \text{on } \Gamma_D,
$$

for a given function $g_D : \Gamma_D \to \mathbb{R}$.

The simplest case is the *homogeneous Dirichlet condition* $g_D = 0$, which models,
for example, a grounded conductor or a surface held at zero temperature.

### Neumann boundary condition

On $\Gamma_N$ (also called the *natural* boundary), the outward normal flux is
prescribed:

$$
\frac{\partial u}{\partial n} \equiv \nabla u \cdot \mathbf{n} = g_N \quad \text{on } \Gamma_N,
$$

where $\mathbf{n}$ is the outward unit normal to $\partial\Omega$ and
$g_N : \Gamma_N \to \mathbb{R}$ is given.

A homogeneous Neumann condition $g_N = 0$ means no flux through that part of the
boundary — a natural symmetry or insulation condition.

### The complete problem

Combining these gives the Poisson problem in strong form:

$$
\begin{cases}
-\Delta u = f & \text{in } \Omega, \\
u = g_D & \text{on } \Gamma_D, \\
\nabla u \cdot \mathbf{n} = g_N & \text{on } \Gamma_N.
\end{cases}
$$

We assume $\Gamma_D \neq \emptyset$ throughout (at least part of the boundary carries a
Dirichlet condition), which is the standard case for well-posedness.

## Function spaces

To state well-posedness precisely, we need to specify in which space solutions are
sought.

The natural space is the *Sobolev space*

$$
H^1(\Omega) = \{ v \in L^2(\Omega) : \nabla v \in [L^2(\Omega)]^d \},
$$

equipped with the inner product

$$
(u, v)_{H^1} = \int_\Omega u v \, dx + \int_\Omega \nabla u \cdot \nabla v \, dx.
$$

Functions in $H^1(\Omega)$ are square-integrable and have square-integrable first
derivatives — the minimal regularity needed to make the problem meaningful.

The subspace with homogeneous Dirichlet conditions on $\Gamma_D$ is

$$
H^1_{0,\Gamma_D}(\Omega) = \{ v \in H^1(\Omega) : v = 0 \text{ on } \Gamma_D \}.
$$

## Well-posedness

**Existence and uniqueness.** If $\Gamma_D \neq \emptyset$, $f \in L^2(\Omega)$,
$g_D \in H^{1/2}(\Gamma_D)$, and $g_N \in L^2(\Gamma_N)$, then the Poisson problem
has a unique solution $u \in H^1(\Omega)$.

This follows from the Lax–Milgram theorem applied to the weak formulation (see
[Weak Formulation of the Poisson Problem](/03-poisson-weak-form/) and
[Weak Formulation of Elliptic PDEs](/04-elliptic-weak-form/)).

**Regularity.** If the domain $\Omega$ and data are sufficiently smooth, the solution
enjoys higher regularity. On a convex domain with $f \in L^2(\Omega)$ and homogeneous
Dirichlet conditions, one has $u \in H^2(\Omega)$. Corners and re-entrant edges in
the domain can reduce regularity, which has important consequences for convergence of
the finite element method.

**Pure Neumann problem.** When $\Gamma_D = \emptyset$, a solution exists only if
$\int_\Omega f \, dx + \int_{\partial\Omega} g_N \, ds = 0$ (compatibility condition),
and it is then unique only up to an additive constant. This case requires special
treatment and is not the focus of the initial posts.

## A one-dimensional example

On $\Omega = (0, 1)$ with $f(x) = \pi^2 \sin(\pi x)$, $u(0) = 0$, $u(1) = 0$
(homogeneous Dirichlet on both ends), the exact solution is

$$
u(x) = \sin(\pi x).
$$

**Verification:**
$$
-u''(x) = -(-\pi^2 \sin(\pi x)) = \pi^2 \sin(\pi x) = f(x). \checkmark
$$

This example is used in the system tests in `tests/system/` to validate the full FEM
pipeline.

## What comes next

The strong form requires $u$ to be twice differentiable — a restrictive condition
that excludes many practically important situations (non-smooth data, domains with
corners). The *weak formulation*, derived in
[Weak Formulation of the Poisson Problem](/03-poisson-weak-form/), relaxes this
requirement and is the actual foundation on which finite element methods are built.

First, though, we generalise the problem: [General Elliptic PDEs](/02-elliptic-pdes/)
extends Poisson to a broader class of second-order operators.
