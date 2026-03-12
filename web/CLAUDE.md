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
- Reference library source files using links, e.g., reference
  yggdrasil/assemble.py as "[`assemble.py`](LIBROOT/assemble.py)".
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
- From Continuous to Discrete
  - The Poisson Problem
  - General Elliptic PDEs
  - Weak Formulation of the Poisson Problem
  - Weak Formulation of Elliptic PDEs
  - The Discrete Formulation
- Elements
  - From Basis to Shape Functions
  - Reference Domains and Quadrature. Code: ReferenceDomain
    - 1D
      - Line. [Gauss-Legendre quadrature](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html). Code: LineDomain
    - 2D
      - Triangle. Code: TriangleDomain
      - Quadrilateral. Code: QuadrilateralDomain
    - 3D
      - Tetrahedron. Code: TetrahedronDomain
      - Hexahedron. Code: HexahedronDomain
  - Shape Functions. Code: ReferenceElement
    - 1D
      - Line-2. Code: Line2
      - Line-3. Code: Line3
    - 2D
      - Triangle-3. Code: Tri3
      - Triangle-6. Code: Tri6
      - Quadrature-4. Code: Quad4
      - Quadrature-9. Code: Quad9
    - 3D
      - Tetrahedron-4. Code: Tet4
      - Hexahedron-8. Code: Hex8
- Meshes
  - Representation. Code: Mesh, ElementGroup
  - Boundaries. Code: extract_boundary
  - Loading and Saving Meshes Using MeshIO
- Assembly
  - Computing the Jacobian. Code: compute_physical_gradients
  - Global stiffness matrix. Code: assemble_bilinear_form
  - Load vector. Code: assemble_load_vector
  - Dirichlet boundary conditions
    - Condensation. Code: CondensedSystem, condense_dirichlet_bc
    - L² projection. Code: project_dirichlet_bc
  - Neumann contributions. Code: assemble_neumann_bc

See [post specification](./post-specification.md) for a list of media and post
files and what they depict or contain.
