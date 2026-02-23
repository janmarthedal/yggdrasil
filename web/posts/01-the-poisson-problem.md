---
title: The Poisson Problem
layout: post
---
Let us begin with one of the most fundamental PDEs: the Poisson problem.

Consider a bounded domain $\Omega \subset \mathbb{R}^d$ with boundary $\partial\Omega$,
where $d$ is the spatial dimension.
The Poisson problem seeks a function $u : \Omega \to \mathbb{R}$ satisfying

$$
-\Delta u = f \quad \text{in } \Omega,
$$

where $\Delta$ denotes the Laplace operator,

$$
\Delta u = \sum_{i=1}^d \frac{\partial^2 u}{\partial x_i^2},
$$

and $f : \Omega \to \mathbb{R}$ is a given source term.

To make the problem well-posed, we must specify boundary conditions.
Let us split the boundary $\partial\Omega$ into two, non-overlapping, parts:
$\Gamma_D$ and $\Gamma_N$, where $\Gamma_D \cup \Gamma_N = \partial\Omega$
and $\Gamma_D \cap \Gamma_N = \emptyset$.

We now impose Dirichlet boundary conditions on $\Gamma_D$:

$$
u = g \quad \text{on } \Gamma_D.
$$

where $g : \Gamma_D \to \mathbb{R}$ is a given boundary value.

We also impose Neumann boundary conditions on $\Gamma_N$:

$$
\frac{\partial u}{\partial n} = h \quad \text{on } \Gamma_N,
$$

where $h : \Gamma_N \to \mathbb{R}$ is a given boundary flux.

As an example, consider **TODO**
