---
title: Weak Formulation of the Poisson Problem
---

The strong form of the [Poisson problem](POSTROOT/continuous/poisson-problem/) — find $u$ satisfying $-\Delta u = f$ pointwise in $\Omega$ together with boundary conditions on $\partial\Omega$ — requires $u$ to be twice continuously differentiable. This regularity is often unavailable for domains or data that arise in practice, and it is also inconvenient for the finite element method. The **weak formulation** (also called the variational formulation) relaxes this requirement and provides the correct mathematical framework for approximation.

The idea is to multiply the equation $-\Delta u = f$ by a smooth **test function** $v$ and integrate over $\Omega$:

$$-\int_\Omega (\Delta u)\, v \, \mathrm{d}x = \int_\Omega f\, v \, \mathrm{d}x.$$

The left-hand side involves second derivatives of $u$, which is precisely what we want to avoid. Applying [Green's first identity](https://en.wikipedia.org/wiki/Green%27s_identities) (integration by parts in multiple dimensions),

$$-\int_\Omega (\Delta u)\, v \, \mathrm{d}x = \int_\Omega \nabla u \cdot \nabla v \, \mathrm{d}x - \int_{\partial\Omega} \frac{\partial u}{\partial n}\, v \, \mathrm{d}s,$$

where $\mathbf{n}$ is the unit outward normal and $\partial u / \partial n = \nabla u \cdot \mathbf{n}$ is the outward normal derivative. Combining with the right-hand side,

$$\int_\Omega \nabla u \cdot \nabla v \, \mathrm{d}x - \int_{\partial\Omega} \frac{\partial u}{\partial n}\, v \, \mathrm{d}s = \int_\Omega f\, v \, \mathrm{d}x.$$

Now the boundary conditions come into play. We choose the test function $v$ to vanish on the Dirichlet boundary: $v = 0$ on $\Gamma_D$. This makes the Dirichlet part of the boundary integral disappear. On the Neumann boundary $\Gamma_N$ the flux is prescribed as $\partial u / \partial n = g_N$, so that part of the boundary integral is known. Moving it to the right-hand side gives

$$\int_\Omega \nabla u \cdot \nabla v \, \mathrm{d}x = \int_\Omega f\, v \, \mathrm{d}x + \int_{\Gamma_N} g_N\, v \, \mathrm{d}s.$$

This is the weak form of the Poisson problem. Two things are worth noting. First, the equation now involves only first derivatives of both $u$ and $v$, which considerably widens the class of admissible functions. Second, the Neumann condition has entered **naturally**: it was not imposed explicitly but appeared as a boundary integral from the integration-by-parts step. For this reason Neumann conditions are often called *natural boundary conditions* in the finite element literature.

To state the weak problem precisely we need to specify the function spaces. The natural space is the **[Sobolev space](https://en.wikipedia.org/wiki/Sobolev_space)** $H^1(\Omega)$, consisting of all square-integrable functions whose first-order partial derivatives are also square-integrable. We define the trial space

$$V = \{ v \in H^1(\Omega) : v = g_D \text{ on } \Gamma_D \}$$

and the test space

$$V_0 = \{ v \in H^1(\Omega) : v = 0 \text{ on } \Gamma_D \}.$$

The weak formulation then reads: find $u \in V$ such that

$$\int_\Omega \nabla u \cdot \nabla v \, \mathrm{d}x = \int_\Omega f\, v \, \mathrm{d}x + \int_{\Gamma_N} g_N\, v \, \mathrm{d}s \quad \text{for all } v \in V_0.$$

The boundary values $g_D$ are said to be imposed **essentially** (they restrict the space $V$), while the Neumann data $g_N$ appear in the right-hand side and are imposed naturally. When $\Gamma_D = \partial\Omega$ (pure Dirichlet problem) and $g_D = 0$, the two spaces coincide: $V = V_0 = H^1_0(\Omega)$.

Well-posedness of this problem — existence and uniqueness of $u \in V$ and continuous dependence on $f$ and $g_N$ — follows from the **[Lax–Milgram theorem](https://en.wikipedia.org/wiki/Lax%E2%80%93Milgram_theorem)**, a fundamental result in functional analysis. The left-hand side is a bounded, coercive bilinear form on $H^1(\Omega)$ and the right-hand side is a bounded linear functional, so the theorem applies. The coercivity here ultimately rests on the **[Poincaré inequality](https://en.wikipedia.org/wiki/Poincar%C3%A9_inequality)**, which guarantees that $\|\nabla v\|_{L^2(\Omega)}$ controls $\|v\|_{H^1(\Omega)}$ for functions vanishing on $\Gamma_D$. A concise account of these results can be found in Evans, *Partial Differential Equations* (AMS, 2010), Chapter 6.

The weak formulation is the starting point for the finite element method. Instead of seeking $u$ in the infinite-dimensional space $V$, we will restrict attention to a finite-dimensional subspace — spanned by carefully chosen basis functions defined on a mesh — and solve the resulting linear system. This is the subject of the next posts.

