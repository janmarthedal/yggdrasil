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

## Neumann Boundary Condition Plan

### Context
Dirichlet BCs are applied manually in `examples/poisson.py` by zeroing rows/columns. Neumann BCs require integrating a flux function `g(x)` over boundary faces: `b_i += integral_Gamma g(x) N_i(x) dS`.

### What already works
1. **Boundary mesh extraction** (`yggdrasil/boundary.py`) — `extract_boundary(mesh)` returns a `Mesh` with lower-dimensional boundary elements and `point_data["original_node_index"]` mapping back to global DOF indices.
2. **Non-square Jacobian mapping** (`yggdrasil/mapping.py`) — handles `topo_dim < spatial_dim` via the metric tensor `g = J^T J`, giving correct surface `det_J`.
3. **Mesh/ElementGroup infrastructure** (`yggdrasil/mesh.py`) — boundary mesh is a valid `Mesh` iterable with `iter_element_groups()`.

### Remaining work

#### 1. `assemble_boundary_load_vector` (core piece) — `yggdrasil/assemble.py`
```python
def assemble_boundary_load_vector(
    boundary_mesh: Mesh,
    g: float | Callable,
    quadrature_order: int,
    num_global_nodes: int,
) -> NDArray[np.float64]:
```
Loop over boundary elements, compute N and det_J (non-square Jacobian path), evaluate g at physical quadrature points, compute `be_i = sum_q g(x_q) * N_i(x_q) * |det_J_q| * w_q`, scatter into global vector using `original_node_index`.

#### 2. Boundary sub-mesh selection — `yggdrasil/boundary.py`
```python
def select_boundary_faces(boundary_mesh: Mesh, predicate: Callable[[NDArray], bool]) -> Mesh:
```
Filter boundary_mesh to faces whose centroid satisfies the predicate, renumber nodes, preserve `original_node_index`.

#### 3. Generalize `assemble_load_vector` to accept callables — DONE

#### 4. Library-level `apply_dirichlet_bc`
Move `apply_dirichlet_bc` from `examples/poisson.py` into the library. Dirichlet elimination must happen *after* Neumann contributions are added to the RHS.

#### 5. Outward normal computation (optional)
For flux-type BCs where g depends on the normal direction. 2D: rotate tangent 90°. 3D: cross product of Jacobian columns.

### Minimum viable implementation
Items 1–2 are sufficient:
```python
bnd = extract_boundary(mesh)
neumann_faces = select_boundary_faces(bnd, lambda x: x[0] > 0.99)
b += assemble_boundary_load_vector(neumann_faces, g=1.0, order=1, num_global_nodes=mesh.num_nodes)
```
Items 4–5 can be added incrementally.
