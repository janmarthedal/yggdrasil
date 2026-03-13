---
title: Notation
---

| Symbol | Description |
|--------|-------------|
| $n$ | Spatial dimension |
| $\Omega \subset \mathbb{R}^n$ | Open bounded domain |
| $\partial\Omega$ | Boundary of $\Omega$ |
| $\Gamma_D$ | Dirichlet part of the boundary |
| $\Gamma_N$ | Neumann part of the boundary |
| $\mathbf{n}$ | Unit outward normal to $\partial\Omega$ |
| $u$ | Unknown scalar field (solution) |
| $f$ | Source term (right-hand side) |
| $g_D$ | Prescribed boundary value (Dirichlet data) |
| $g_N$ | Prescribed normal flux (Neumann data) |
| $\Delta u$ | Laplacian of $u$, $\sum_{i=1}^n \partial^2 u / \partial x_i^2$ |
| $\nabla u$ | Gradient of $u$ |
| $\nabla \cdot \mathbf{v}$ | Divergence of a vector field $\mathbf{v}$ |
| $\partial u / \partial n$ | Outward normal derivative of $u$, equal to $\nabla u \cdot \mathbf{n}$ |
| $L$ | Second-order linear elliptic differential operator |
| $A$ | Matrix-valued diffusion coefficient, $A : \Omega \to \mathbb{R}^{n \times n}$ |
| $\mathbf{b}$ | Convection field, $\mathbf{b} : \Omega \to \mathbb{R}^n$ |
| $c$ | Reaction coefficient, $c : \Omega \to \mathbb{R}$ |
| $\alpha$ | Ellipticity constant |
| $H^1(\Omega)$ | Sobolev space of square-integrable functions with square-integrable first derivatives |
| $H^1_0(\Omega)$ | Functions in $H^1(\Omega)$ vanishing on $\partial\Omega$ |
| $V$ | Trial function space |
| $V_0$ | Test function space (functions vanishing on $\Gamma_D$) |
| $v$ | Test function |
| $\mathrm{d}x$ | Volume integration element in $\Omega$ |
| $\mathrm{d}s$ | Surface integration element on $\partial\Omega$ |
| $a(u, v)$ | Bilinear form |
| $\ell(v)$ | Linear functional (right-hand side) |
| $V_h$, $V_{h,0}$ | Discrete trial and test spaces |
| $N$ | Dimension of the discrete space |
| $\phi_i$ | Basis functions |
| $u_h$ | Discrete solution |
| $K$ | Stiffness matrix, $K_{ij} = a(\phi_j, \phi_i)$ |
| $\mathbf{u}$ | Vector of unknown coefficients |
| $\mathbf{f}$ | Load vector, $f_i = \ell(\phi_i)$ |

