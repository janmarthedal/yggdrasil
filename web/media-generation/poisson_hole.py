"""Poisson problem on a square with an off-centre circular hole.

Domain:  Ω = [-1, 1]² minus the open disk (x-cx)² + (y-cy)² < r²
         with hole centre (cx, cy) = (0.3, 0.2) and r = 0.5

Boundary conditions:
    u = 0         on the four sides of the square   ← Dirichlet (homogeneous)
    ∂u/∂n = 0    on the circular boundary            ← Neumann (zero flux, natural)

Source term: f = 1

The mesh is built from scratch using SciPy's Delaunay triangulation:
  1. Place boundary nodes on the four square sides and on the circle.
  2. Seed random interior nodes outside the circle.
  3. Triangulate with scipy.spatial.Delaunay.
  4. Remove any triangle whose centroid falls inside the circle.
"""

from pathlib import Path

import meshio
import matplotlib.pyplot as plt
import matplotlib.tri
import numpy as np
import scipy.sparse.linalg as spla
from scipy.spatial import Delaunay

from yggdrasil import (
    ElementGroup,
    Mesh,
    assemble_bilinear_form,
    assemble_load_vector,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    select_boundary_faces,
    tag_boundary_faces,
)
from yggdrasil.elements import Tri3


def make_mesh(
    half_side: float = 1.0,
    hole_center: tuple[float, float] = (0.3, 0.2),
    radius: float = 0.5,
    n_side: int = 24,
    n_circle: int = 36,
    n_interior: int = 1000,
    rng_seed: int = 0,
) -> Mesh:
    """Build a triangular mesh of a square with a circular hole.

    Parameters
    ----------
    half_side : float
        Half-side length of the square (domain is [-half_side, half_side]²).
    hole_center : tuple[float, float]
        (cx, cy) centre of the circular hole.
    radius : float
        Radius of the circular hole.
    n_side : int
        Number of nodes per side of the square (corners shared between sides).
    n_circle : int
        Number of nodes on the circular boundary.
    n_interior : int
        Number of random interior seed nodes.
    rng_seed : int
        Random seed for reproducibility.
    """
    cx, cy = hole_center

    # Square boundary: traverse all four sides, excluding repeated corners
    t = np.linspace(-half_side, half_side, n_side, endpoint=False)
    square_pts = np.vstack([
        np.column_stack([t,                          -half_side * np.ones(n_side)]),  # bottom
        np.column_stack([half_side * np.ones(n_side), t]),                            # right
        np.column_stack([-t,                          half_side * np.ones(n_side)]),  # top
        np.column_stack([-half_side * np.ones(n_side), -t]),                          # left
    ])

    # Circular boundary
    theta = np.linspace(0.0, 2.0 * np.pi, n_circle, endpoint=False)
    circle_pts = np.column_stack([cx + radius * np.cos(theta), cy + radius * np.sin(theta)])

    # Random interior nodes (outside the circle, with a small clearance)
    rng = np.random.default_rng(rng_seed)
    interior_pts: list[np.ndarray] = []
    while len(interior_pts) < n_interior:
        candidates = rng.uniform(-half_side, half_side, size=(n_interior * 4, 2))
        dx = candidates[:, 0] - cx
        dy = candidates[:, 1] - cy
        keep = dx ** 2 + dy ** 2 > (radius * 1.1) ** 2
        interior_pts.extend(candidates[keep])
    interior_pts_arr = np.array(interior_pts[:n_interior])

    # Combine and triangulate
    points = np.vstack([square_pts, circle_pts, interior_pts_arr])
    tri = Delaunay(points)

    # Drop triangles whose centroid lies inside the circle
    centroids = points[tri.simplices].mean(axis=1)
    dx = centroids[:, 0] - cx
    dy = centroids[:, 1] - cy
    outside = dx ** 2 + dy ** 2 >= radius ** 2
    connectivity = tri.simplices[outside].astype(np.intp)

    return Mesh(points, [ElementGroup(element=Tri3(), connectivity=connectivity)])


def main():
    mesh = make_mesh()

    # Assemble stiffness matrix and load vector (f = 1)
    K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
    b = assemble_load_vector(mesh, 1.0, quadrature_order=1)

    # Tag square-side faces as Dirichlet (tag=1).
    # Circle faces keep tag=0 — zero Neumann is the natural BC, no assembly needed.
    bnd = extract_boundary(mesh)
    half_side = 1.0
    bnd = tag_boundary_faces(
        bnd,
        lambda c: (np.abs(c[:, 0]) > half_side - 1e-8) | (np.abs(c[:, 1]) > half_side - 1e-8),
        tag=1,
    )

    # Apply homogeneous Dirichlet BC on the square boundary
    dirichlet_mesh = select_boundary_faces(bnd, tag=1)
    bc_nodes = dirichlet_mesh.point_data["original_node_index"]
    system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    # Solve
    u = system.reconstruct(spla.spsolve(system.K, system.b))

    print(f"Nodes: {mesh.num_nodes},  Elements: {mesh.num_elements}")
    print(f"Max u: {u.max():.6f}  (node {u.argmax()})")

    # Save mesh and solution to disk
    out_path = "poisson_hole.vtu"
    points_3d = np.column_stack([mesh.nodes, np.zeros(mesh.num_nodes)])
    meshio.write(
        out_path,
        meshio.Mesh(
            points=points_3d,
            cells=[("triangle", mesh.element_groups[0].connectivity)],
            point_data={"u": u},
        ),
    )
    print(f"Saved mesh and solution to {out_path}")

    # Plot
    triangulation = matplotlib.tri.Triangulation(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        mesh.element_groups[0].connectivity,
    )
    fig, ax = plt.subplots(figsize=(9, 6))
    tpc = ax.tripcolor(triangulation, u, shading="gouraud", cmap="viridis")
    fig.colorbar(tpc, ax=ax, label="u")
    ax.set_title("Poisson on square with off-centre circular hole  (f = 1, u = 0 on sides)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.tight_layout()
    out_img = Path(__file__).parent.parent / "media" / "poisson-2d-solution.png"
    out_img.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_img, dpi=150)
    print(f"Saved solution image to {out_img}")


if __name__ == "__main__":
    main()
