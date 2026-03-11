# Weak Formulation of Elliptic PDEs

The [weak formulation of the Poisson problem](POSTROOT/03-poisson-weak-form/) extends naturally to the [general elliptic operator](POSTROOT/02-elliptic-pdes/)

$$Lu = -\nabla \cdot (A \nabla u) + \mathbf{b} \cdot \nabla u + c\, u.$$

The derivation follows the same steps: multiply $Lu = f$ by a test function $v \in V_0$ and integrate over $\Omega$,

$$\int_\Omega \bigl(-\nabla \cdot (A \nabla u) + \mathbf{b} \cdot \nabla u + c\, u\bigr)\, v \, \mathrm{d}x = \int_\Omega f\, v \, \mathrm{d}x.$$

Applying [Green's first identity](https://en.wikipedia.org/wiki/Green%27s_identities) to the divergence term,

$$-\int_\Omega \nabla \cdot (A \nabla u)\, v \, \mathrm{d}x = \int_\Omega (A \nabla u) \cdot \nabla v \, \mathrm{d}x - \int_{\partial\Omega} (A \nabla u) \cdot \mathbf{n}\, v \, \mathrm{d}s.$$

The boundary integral vanishes on $\Gamma_D$ because $v = 0$ there, and on $\Gamma_N$ the generalised Neumann condition prescribes $(A \nabla u) \cdot \mathbf{n} = g_N$. Moving the known flux to the right-hand side and collecting all terms gives

$$\int_\Omega \bigl[(A \nabla u) \cdot \nabla v + (\mathbf{b} \cdot \nabla u)\, v + c\, u\, v\bigr] \mathrm{d}x = \int_\Omega f\, v \, \mathrm{d}x + \int_{\Gamma_N} g_N\, v \, \mathrm{d}s.$$

It is conventional to name the two sides separately. Define the **bilinear form** $a : V_0 \times V_0 \to \mathbb{R}$ by

$$a(u, v) = \int_\Omega \bigl[(A \nabla u) \cdot \nabla v + (\mathbf{b} \cdot \nabla u)\, v + c\, u\, v\bigr] \mathrm{d}x,$$

and the **linear functional** $\ell : V_0 \to \mathbb{R}$ by

$$\ell(v) = \int_\Omega f\, v \, \mathrm{d}x + \int_{\Gamma_N} g_N\, v \, \mathrm{d}s.$$

The weak problem then takes the compact abstract form: find $u \in V$ such that

$$a(u, v) = \ell(v) \quad \text{for all } v \in V_0.$$

The Poisson problem is recovered by setting $A = I$, $\mathbf{b} = 0$, $c = 0$, which gives $a(u, v) = \int_\Omega \nabla u \cdot \nabla v \, \mathrm{d}x$. The abstract notation $a(u,v) = \ell(v)$ is standard throughout the finite element literature and applies equally to far more general problems.

Well-posedness again follows from the [Lax–Milgram theorem](https://en.wikipedia.org/wiki/Lax%E2%80%93Milgram_theorem), provided $a$ is *bounded* and *coercive* on $V_0$ and $\ell$ is bounded. Boundedness is straightforward under mild regularity assumptions on $A$, $\mathbf{b}$, and $c$. Coercivity is more delicate: for the purely diffusive case ($\mathbf{b} = 0$, $c \geq 0$) it follows directly from the [Poincaré inequality](https://en.wikipedia.org/wiki/Poincar%C3%A9_inequality) and the ellipticity condition on $A$. When convection is present, additional assumptions on $\mathbf{b}$ and $c$ are required; a sufficient condition is that $c - \tfrac{1}{2}\nabla \cdot \mathbf{b} \geq 0$ almost everywhere in $\Omega$. A full treatment is given in Brenner and Scott, *The Mathematical Theory of Finite Element Methods* (Springer, 2008), Chapter 5.

The abstract formulation $a(u,v) = \ell(v)$ is the foundation on which the discrete approximation is built. In the next post we replace the infinite-dimensional space $V$ with a finite-dimensional subspace and derive the linear system that must be solved.
