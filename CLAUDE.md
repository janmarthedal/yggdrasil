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
- Benchmarks: `uv run --group bench pytest benchmarks/` (not collected by `uv run pytest`)
- Benchmarks + scikit-fem comparison (opt-in): `uv run --group bench-compare pytest benchmarks/compare/`

## File structure
- `yggdrasil/` — library source
  - `assemble.py` — assembling matrices and vectors from forms; Dirichlet condensation/projection
  - `forms.py` — bilinear form definitions (`mass_form`, `grad_grad_form`)
  - `dof_map.py` — node → global DOF mapping; scalar only today, seam for future vector/multi-field problems
  - `error.py` — error computation utilities (L² error)
  - `io.py` — mesh I/O (reading/writing mesh files)
  - `refdomains/` — reference domain definitions (line, triangle, quadrilateral, tetrahedron, hexahedron) with quadrature rules
  - `elements/` — reference element implementations (shape functions and gradients) for each supported element type
  - `mapping.py` — Jacobian, physical-space gradients, and mapped quadrature utilities
  - `boundary.py` — boundary mesh extraction (`extract_boundary`), tagging (`tag_boundary_faces`), and face selection (`select_boundary_faces`)
  - `mesh.py` — immutable `Mesh` and `ElementGroup` classes for storing nodes, connectivity, and auxiliary data
  - `mesh_generators.py` — structured mesh constructors and `map_mesh_points` (remap node coordinates)
- `tests/`
  - `unit/` — unit tests for individual modules (mesh, mapping, boundary, assembly, forms, error, dof_map, elements, refdomains)
  - `system/` — full PDE solutions compared against analytical solutions where possible
- `examples/` — example scripts (require `--group examples` for extra dependencies like matplotlib)
- `benchmarks/` — `pytest-benchmark` timing suite (`--group bench`); `problems.py` holds shared solve-stage factories mirroring `tests/system/`; `compare/` is the opt-in scikit-fem comparison. See `benchmarks/README.md`.

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
- **Dirichlet/Neumann ordering**: add Neumann contributions (`assemble_neumann_bc`) to the RHS *before* calling `condense_dirichlet_bc`
- **Non-square Jacobians**: `mapping.py` handles `topo_dim < spatial_dim` via metric tensor `g = J^T J` — used for surface integrals in `boundary.py`
- **`original_node_index`**: boundary meshes carry this in `point_data` to map local node indices back to global
- **Mesh immutability**: `Mesh`/`ElementGroup` are immutable — `nodes`/`connectivity` arrays are non-writeable, `element_groups` is a tuple, `point_data` is a `MappingProxyType`; operations like `tag_boundary_faces` and `map_mesh_points` return new meshes
