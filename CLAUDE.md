# Yggdrasil

## Project outline
The project aims to be a library for finite element analysis.

It should offer functionality such as:
- Loading or creating meshes.
- A wide range of elements of different orders (initially restricted to the most used).
- Quadrature for the elements of different orders.
- Specifying PDEs from their weak (perhaps strong?) form.
- Applying Dirichlet and/or Neumann boundary conditions to the underlying geometry and/or mesh.
- Assembling the discrete systems needed for approximate solutions.
- Making it easy to visualize elements, meshes, solutions, etc.

## Environment
- Use uv for package management

## Running tests and examples
- Tests: `uv run pytest`
- Examples: `uv run --group examples python examples/<path>.py`

## File structure
- `yggdrasil/` — library source
  - `assemble.py` — assembling matrices and vectors from (bilinear) functions
  - `domains/` — reference domain definitions (line, triangle, quadrilateral, tetrahedron, hexahedron) with quadrature rules
  - `elements/` — reference element implementations (shape functions and gradients) for each supported element type
  - `mapping.py` — Jacobian, physical-space gradients, and mapped quadrature utilities
  - `mesh.py` — `Mesh` class for storing nodes, connectivity, and auxiliary data
  - `mesh_generators.py` — helper functions for creating structured meshes
- `tests/` — pytest test suite
- `examples/` — example scripts (require `--group examples` for extra dependencies like matplotlib)
- `web/`
  - `pages/` — interactive web-based visualizations (standalone HTML using Three.js)
  - `posts/` — markdown posts describing the theoretical foundation of the library code

## Packages
SciPy and NumPy are used as the foundation for all computations.
The functionality is otherwise self-contained.

## Technical specifications

### Supported elements
- 1D
	- Line2 (2 nodes, linear)
	- Line3 (3 nodes, quadratic)
- 2D
	- Tri3 (3 nodes, linear)
	- Tri6 (6 nodes, quadratic)
	- Quad4 (4 nodes, bilinear)
	- Quad9 (9 nodes, biquadratic)
- 3D
	- Tet4 (4 nodes, linear)
	- Hex8 (8 nodes, trilinear)
