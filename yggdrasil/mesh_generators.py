import numpy as np

from .elements import Tri3
from .mesh import ElementGroup, Mesh


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
    x = np.linspace(0, 1, n + 1)
    y = np.linspace(0, 1, n + 1)
    xx, yy = np.meshgrid(x, y)
    nodes = np.column_stack([xx.ravel(), yy.ravel()])

    triangles = []
    for j in range(n):
        for i in range(n):
            n0 = j * (n + 1) + i
            n1 = n0 + 1
            n2 = n0 + (n + 1) + 1
            n3 = n0 + (n + 1)
            triangles.append([n0, n1, n2])
            triangles.append([n0, n2, n3])

    conn = np.array(triangles, dtype=np.intp)
    group = ElementGroup(element=Tri3(), connectivity=conn)
    return Mesh(nodes, [group])
