# Weak Formulation of the Poisson Problem

Let us consider the Poisson problem again, but in a slightly more general version:

$$
-\nabla \cdot \left( a(x) \nabla u(x) \right) = f(x) \quad \text{for $x \in \Omega$}
$$

with boundary conditions
- **Dirichlet**: $u=g$ on $\Gamma_D$
- **Neumann**: $a(x) \nabla u \cdot \mathbf{n} = h$ on $\Gamma_N$

and where:
- $\Omega \subset \mathbb{R}^d$ is the domain
- $a(x) > 0$ is a spatially varying scalar coefficient
- $f(x)$ is a known source
- $\mathbf{n}$ is the outward normal on the boundary

We seek a weak (integral) form of this PDE.

Let $v$ be a **test function**, which is a function $\Omega \to \mathbb{R}$ which is
- Smooth (sufficiently differentiable on the domain)
- $v=0$ on $\Gamma_D$

By multiplying both sides of the PDE by $v(x)$ we get:

$$
-\nabla \cdot \left( a(x) \nabla u \right) v = f v \quad \text{in } \Omega
$$

and integrate over the domain

$$
\int_\Omega - \nabla \cdot \left( a(x) \nabla u \right) v \, dx = \int_\Omega f v \, dx
$$

Apply integration by parts (divergence theorem) to the left-hand side:

$$
\int_\Omega - \nabla \cdot \left( a(x) \nabla u \right) v \, dx = \int_\Omega a(x) \nabla u \cdot \nabla v \, dx - \int_{\partial\Omega} a(x) (\nabla u \cdot \mathbf{n}) v \, ds
$$

So we have:

$$
\int_\Omega a(x) \nabla u \cdot \nabla v \, dx = \int_\Omega f v \, dx + \int_{\partial\Omega} a(x) (\nabla u \cdot \mathbf{n}) v \, ds
$$

When applying the boundary conditions we split the boundary integral:
- On $\Gamma_D$ we have $v=0$, so it vanishes.
- On $\Gamma_N$ we have $a(x) \nabla u \cdot \mathbf{n} = h$, so:

$$
\int_{\partial\Omega} a(x) (\nabla u \cdot \mathbf{n}) v \, ds = \int_{\Gamma_N} h v \, ds
$$

The final Weak Form is: Find $u \in V$ such that

$$
\int_\Omega a(x) \nabla u \cdot \nabla v \, dx = \int_\Omega f v \, dx + \int_{\Gamma_N} h v \, ds \quad \forall v \in V_0
$$

where:
- $V$: functions satisfying $u=g$ on $\Gamma_D$
- $V_0$: test functions with $v=0$ on $\Gamma_D$
