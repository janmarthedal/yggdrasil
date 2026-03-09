r""" Three-dimensional Poisson equation

This example corresponds to Example 9 of scikit-fem's Gallery of examples,
https://scikit-fem.readthedocs.io/en/latest/listofexamples.html#example-9-three-dimensional-poisson-equation
with source code https://github.com/kinnala/scikit-fem/blob/master/docs/examples/ex09.py

Problem
-------
Find u : Ω → ℝ such that

    -Δu = 1  in Ω = [0, 1]³
      u = 0  on ∂Ω

The source term is f = 1 throughout the unit cube and homogeneous Dirichlet
boundary conditions are imposed on all six faces.

There is no simple closed-form solution.  The Fourier sine series gives

    u(x, y, z) = (8/π⁵) Σ_{m, n, p odd}
                     sin(mπx) sin(nπy) sin(pπz) / (mnp(m² + n² + p²))

The maximum is attained at the centre (½, ½, ½) and equals approximately
0.0562 (the alternating signs of sin(mπ/2) cause significant cancellation
beyond the m=n=p=1 term).

Weak formulation
----------------
Find u ∈ H¹₀(Ω) such that for all v ∈ H¹₀(Ω)

    ∫_Ω ∇u · ∇v dΩ = ∫_Ω v dΩ

Discretisation
--------------
A structured tetrahedral mesh of [0, 1]³ is constructed by subdividing each
axis-aligned cube into 6 tetrahedra via the main-diagonal decomposition (all
tetrahedra share the vertex diagonal from (0,0,0) to (1,1,1) of each cube).
Linear Tet4 elements are used with quadrature order 1.
"""

import numpy as np
import scipy.sparse.linalg as spla

from yggdrasil import (
    assemble_bilinear_form,
    assemble_load_vector,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    unit_cube_tet_mesh,
)
from yggdrasil.io import to_meshio


def main() -> None:
    n = 8  # 8³ × 6 = 3 072 tetrahedra, (n+1)³ = 729 nodes
    mesh = unit_cube_tet_mesh(n)

    # Assemble stiffness matrix and unit load vector (f = 1)
    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, 1.0, quadrature_order=1)

    # Identify all boundary nodes and impose u = 0 on ∂Ω (full Dirichlet)
    bnd = extract_boundary(mesh)
    bc_nodes = np.unique(bnd.point_data["original_node_index"])
    system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    # Solve the condensed linear system
    u = system.reconstruct(spla.spsolve(system.K, system.b))

    # Report maximum value; Fourier series gives ≈ 0.0562 in the fine-mesh limit
    print(f"Max u = {u.max():.6f}  (expected ≈ 0.0562 as mesh is refined)")

    # Save solution to VTK for post-processing
    output_path = "poisson-3d.vtk"
    meshio_mesh = to_meshio(mesh)
    meshio_mesh.point_data["u"] = u
    meshio_mesh.write(output_path)
    print(f"Solution saved to {output_path}")


if __name__ == "__main__":
    main()
