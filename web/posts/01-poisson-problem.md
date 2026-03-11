# The Poisson Problem

The **Poisson equation** is one of the most important partial differential equations (PDEs) in applied mathematics. Given a domain $\Omega \subset \mathbb{R}^n$ and a source term $f : \Omega \to \mathbb{R}$, the Poisson equation seeks a function $u : \Omega \to \mathbb{R}$ satisfying

$$-\Delta u = f \quad \text{in } \Omega,$$

where $\Delta u = \sum_{i=1}^n \frac{\partial^2 u}{\partial x_i^2}$ is the Laplacian of $u$. The negative sign is a convention that makes the operator positive definite, as we will see later. It arises in [numerous areas of physics and engineering](https://en.wikipedia.org/wiki/Poisson%27s_equation#Applications_in_physics_and_engineering), including electrostatics, gravitational potential, heat conduction, and fluid mechanics.

To obtain a unique solution, the PDE must be supplemented with **boundary conditions** on $\partial\Omega$. The boundary is typically partitioned as $\partial\Omega = \Gamma_D \cup \Gamma_N$ (with $\Gamma_D \cap \Gamma_N = \emptyset$), and two types of conditions are imposed.

**Dirichlet boundary conditions** prescribe the value of the solution directly:

$$u = g_D \quad \text{on } \Gamma_D.$$

**Neumann boundary conditions** prescribe the outward normal derivative (flux):

$$\frac{\partial u}{\partial n} = g_N \quad \text{on } \Gamma_N,$$

where $\mathbf{n}$ is the unit outward normal to $\partial\Omega$. When $g_N = 0$, the condition is called a **zero-flux** or **natural** condition; it arises automatically in the weak formulation and requires no special treatment.

Consider a square domain $\Omega = [-1, 1]^2$ with an off-centre circular hole of radius $r = 0.5$ centered at $(0.3, 0.2)$. The outer square boundary $\Gamma_D$ carries a homogeneous Dirichlet condition $u = 0$, while the circular boundary $\Gamma_N$ carries a zero-flux Neumann condition $\partial u / \partial n = 0$. The source term is $f = 1$.

![Domain with Dirichlet boundary (square sides) and Neumann boundary (circle)](MEDIAROOT/poisson-2d-domain.svg)

Solving this problem — using the finite element method developed throughout these posts — gives the solution shown below.

![Solution to the Poisson problem on the square with circular hole](MEDIAROOT/poisson-2d-solution.png)

In one dimension the Poisson equation on $\Omega = (0, 1)$ with $f = 1$ and homogeneous Dirichlet conditions $u(0) = u(1) = 0$ reads

$$-u''(x) = 1, \quad x \in (0, 1), \qquad u(0) = u(1) = 0.$$

Integrating twice and applying the boundary conditions gives the exact solution

$$u(x) = \frac{x(1 - x)}{2}.$$

![Analytical solution to the 1D Poisson problem](MEDIAROOT/poisson-1d-solution.svg)

This simple example will serve as a useful reference throughout the series: it is easy to verify computed approximations against the exact parabolic profile.
