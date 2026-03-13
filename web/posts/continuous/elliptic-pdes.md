---
title: General Elliptic PDEs
---

The [Poisson problem](POSTROOT/continuous/poisson-problem/) is a special case of a much broader family of PDEs. A general second-order linear elliptic PDE on a domain $\Omega \subset \mathbb{R}^n$ can be written as

$$Lu = f \quad \text{in } \Omega,$$

where $L$ is a second-order linear differential operator and $f : \Omega \to \mathbb{R}$ is a given source term. A canonical form of $L$ is

$$Lu = -\nabla \cdot (A \nabla u) + \mathbf{b} \cdot \nabla u + c\, u,$$

with $A : \Omega \to \mathbb{R}^{n \times n}$ a matrix-valued diffusion coefficient, $\mathbf{b} : \Omega \to \mathbb{R}^n$ a convection field, and $c : \Omega \to \mathbb{R}$ a reaction coefficient. The Poisson equation corresponds to $A = I$, $\mathbf{b} = 0$, $c = 0$.

The term *elliptic* refers to a condition on $A$ that rules out wave-like or hyperbolic degeneracies. The operator $L$ is **uniformly elliptic** if there exists a constant $\alpha > 0$ such that

$$\xi^T A(x)\, \xi \geq \alpha |\xi|^2 \quad \text{for all } \xi \in \mathbb{R}^n \text{ and almost all } x \in \Omega.$$

This ensures the problem has a diffusive character in every spatial direction. Boundary conditions are of the same Dirichlet and Neumann types introduced for the Poisson problem, with the Neumann condition now prescribing the flux $(A \nabla u) \cdot \mathbf{n} = g_N$ to account for anisotropic diffusion.

Well-posedness — existence, uniqueness, and continuous dependence of $u$ on $f$ — is the central question before attempting a numerical solution. Under sufficient conditions on $A$, $\mathbf{b}$, $c$, and the domain, classical results guarantee that the boundary value problem $Lu = f$ has a unique solution that depends stably on the data. The precise hypotheses and proofs are given in Evans, *Partial Differential Equations* (AMS, 2010), Chapter 6, which is the standard graduate reference. A more numerically oriented treatment can be found in Brenner and Scott, *The Mathematical Theory of Finite Element Methods*, Chapter 5.

The precise mechanism by which well-posedness is established — via a weak formulation and the Lax–Milgram theorem — is the subject of the next posts.

