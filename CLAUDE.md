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
- All tests: `uv run pytest`
- Unit tests: `uv run pytest tests/unit/`
- System tests: `uv run pytest tests/system/`
- Linting: `uv run ruff check`
- Type checking: `uv run --group dev ty check`
- Examples: `uv run --group examples python examples/<path>.py`
- Web build: `cd web && bun run build` (outputs to `web/_site/`)
- Web dev server: `cd web && bun run dev` (live reload via browser-sync)

## File structure
- `yggdrasil/` — library source
  - `assemble.py` — assembling matrices and vectors from (bilinear) functions
  - `forms.py` — bilinear and linear form definitions (mass, stiffness, load)
  - `dof_map.py` — DOF mapping for multi-field problems
  - `error.py` — error computation utilities (L² error, H¹ error)
  - `io.py` — mesh I/O (reading/writing mesh files)
  - `refdomains/` — reference domain definitions (line, triangle, quadrilateral, tetrahedron, hexahedron) with quadrature rules
  - `elements/` — reference element implementations (shape functions and gradients) for each supported element type
  - `mapping.py` — Jacobian, physical-space gradients, and mapped quadrature utilities
  - `boundary.py` — boundary mesh extraction (`extract_boundary`) and face selection (`select_boundary_faces`)
  - `mesh.py` — `Mesh` class for storing nodes, connectivity, and auxiliary data
  - `mesh_generators.py` — helper functions for creating structured meshes
- `tests/`
  - `unit/` — unit tests for individual modules (mesh, mapping, boundary, quadrature, elements, assembly)
  - `system/` — full PDE solutions compared against analytical solutions where possible
- `examples/` — example scripts (require `--group examples` for extra dependencies like matplotlib)
- `web/`
  - `posts/` — markdown posts describing the theoretical foundation of the library code (see `web/CLAUDE.md` for writing conventions)
  - `media/` — images and media files for posts
  - `media-generation/` — Python scripts that generate SVG/PNG assets for posts (shape functions, quadrature points, domain illustrations, etc.)

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

## Key gotchas
- **Dirichlet/Neumann ordering**: add Neumann contributions to RHS *before* applying `apply_dirichlet_bc`
- **Non-square Jacobians**: `mapping.py` handles `topo_dim < spatial_dim` via metric tensor `g = J^T J` — used for surface integrals in `boundary.py`
- **`original_node_index`**: boundary meshes carry this in `point_data` to map local DOF indices back to global
