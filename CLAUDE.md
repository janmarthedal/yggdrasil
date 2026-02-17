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

## Packages
Use SciPy and NumPy as the foundation for all computations.
The functionality should otherwise be self-contained.

Feel free to suggest packages for generating visualizations of elements, meshes, solutions, etc.

Feel free to suggest packages for testing and other development-related parts of the project.

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
