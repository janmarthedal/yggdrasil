import numpy as np

from .elements import Tri3
from .mesh import ElementGroup, Mesh


def map_mesh_points(mesh: Mesh, f) -> Mesh:
    """Return a new Mesh with nodes mapped through f, keeping connectivity.

    Parameters
    ----------
    mesh : Mesh
    f : callable
        Vectorized function ``(N, d) -> (N, d')`` applied to the node array.

    Returns
    -------
    Mesh
    """
    return Mesh(f(mesh.nodes), mesh.element_groups, mesh.point_data)


def rectangular_tri_mesh(x0: float, x1: float, y0: float, y1: float, nx: int, ny: int) -> Mesh:
    """Create a structured triangular mesh on [x0,x1]×[y0,y1].

    Parameters
    ----------
    x0, x1 : float
        x-axis bounds.
    y0, y1 : float
        y-axis bounds.
    nx, ny : int
        Number of nodes in each direction; produces (nx-1)×(ny-1) cells,
        each split into 2 triangles, giving 2*(nx-1)*(ny-1) elements total.

    Returns
    -------
    Mesh with a single Tri3 element group.
    """
    xx, yy = np.meshgrid(np.linspace(x0, x1, nx), np.linspace(y0, y1, ny))
    nodes = np.column_stack([xx.ravel(), yy.ravel()])

    triangles = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            n0 = j * nx + i
            n1 = n0 + 1
            n2 = n0 + nx + 1
            n3 = n0 + nx
            triangles.append([n0, n1, n2])
            triangles.append([n0, n2, n3])

    conn = np.array(triangles, dtype=np.intp)
    return Mesh(nodes, [ElementGroup(element=Tri3(), connectivity=conn)])


def unit_square_tri_mesh(n: int) -> Mesh:
    """Create a structured triangular mesh on [0,1]² with (n+1)² nodes.

    The unit square is divided into n×n quadrilateral cells, each split
    into 2 triangles (Tri3 elements), giving 2n² elements total.

    Parameters
    ----------
    n : int
        Number of subdivisions along each axis.

    Returns
    -------
    Mesh with a single Tri3 element group.
    """
    return rectangular_tri_mesh(0.0, 1.0, 0.0, 1.0, n + 1, n + 1)
