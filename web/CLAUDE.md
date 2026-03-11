# Posts and media on the Finite Element Method and the Yggdrasil project

## Purpose
The folder ./posts/ contains a series of posts introducing the finite element
method, accompanying the `yggdrasil` library in this repository.
The posts explain the mathematical foundation of the finite element method
and the code constructs that make up the library in relation to the theory.
There will be little abstract mathematical theory, but references to the
underlying theory should be included. Resources accessible online are preferred,
but authoratative books and papars may also be mentioned.
It is an important goal that every function and class from the library are
referenced from at least one post at some point.

## Post guidelines
- Plenty of examples, illustrations and animations should be included.
- Each post is of small size and as self-contained as possible.
- Use a shared set of notation, which is listed in `notation.md`.
  Update this file as posts are written.

## Writing conventions
- Use markdown syntax.
- Reference library source files by path relative to `LIBROOT`, e.g.
  `LIBROOT/assemble.py`.
- Reference image files by path relative to `MEDIAROOT`.
- Reference other posts by path relative to `POSTROOT`. The post link for file
  `some-post.md` should be `POSTROOT/some-post/`.
- Math is written in LaTeX delimited by `$...$` (inline) and `$$...$$` (display).
- Mark incomplete sections with `**TODO**`.
- Prefer SVG over PNG for inline images when the illustration is fit for vector
  graphics.
- Consider using [three.js](https://threejs.org/) for 3D illustrations.
- The posts are prefixed with 01, 02, and so on, such that they form a linear
  narrative.
- Each post has a main heading, but headings should otherwise be avoided.
  Let the text flow as continuous prose.

## Overall post structure
- The continuous problem.
  - Present a general formulation, but start with the Poisson problem as an example.
  - Introduce Dirichlet and Neumann boundary conditions.
  - Mention basic results on well-posedness of such problems.
- The weak formulation.
  - Show how to derive the weak formulation in the Poisson and general case.
- Discrete formulation
  - Formulation using basis function.
  - Optimality of solutions.
- Elements
  - Reference domains
    - Parameterization
    - Quadrature
  - Shape functions
    - Mapping from reference domain to physical domain
    - Jacobians
- Meshes
  - Representation
  - Boundaries
- Assembly
  - Global stiffness matrix
  - Load vector
  - Dirichlet boundary conditions
    - Condensation
    - L² projection
  - Neumann contributions
- Solving time-dependent systems
  - Heat equation
  - Wave equation
- External tools
  - Gmsh
  - Matplotlib
  - MeshIO
  - Paraview

See [post specification](./post-specification.md) for a list of media and post
files and what they depict or contain.
