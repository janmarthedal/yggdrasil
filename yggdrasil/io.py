"""MeshIO integration: import and export yggdrasil meshes via meshio.

Requires the optional ``meshio`` package.  Install it with::

    pip install meshio
    # or, with uv:
    uv add --optional meshio meshio

Node orderings for all supported element types are identical between
yggdrasil and the VTK convention used by meshio, so no index remapping
is performed.

Supported cell types
--------------------
==============  ================  =====================
meshio type     yggdrasil class   description
==============  ================  =====================
line            Line2             2-node linear line
line3           Line3             3-node quadratic line
triangle        Tri3              3-node linear triangle
triangle6       Tri6              6-node quadratic triangle
quad            Quad4             4-node bilinear quad
quad9           Quad9             9-node biquadratic quad
tetra           Tet4              4-node linear tetrahedron
hexahedron      Hex8              8-node trilinear hexahedron
==============  ================  =====================
"""

import numpy as np

from .elements import Hex8, Line2, Line3, Quad4, Quad9, Tet4, Tri3, Tri6
from .elements.element import ReferenceElement
from .mesh import ElementGroup, Mesh

_MESHIO_TO_ELEMENT: dict[str, type[ReferenceElement]] = {
    "line": Line2,
    "line3": Line3,
    "triangle": Tri3,
    "triangle6": Tri6,
    "quad": Quad4,
    "quad9": Quad9,
    "tetra": Tet4,
    "hexahedron": Hex8,
}

_ELEMENT_TO_MESHIO: dict[type[ReferenceElement], str] = {
    v: k for k, v in _MESHIO_TO_ELEMENT.items()
}


def _require_meshio():
    try:
        import meshio  # noqa: F401
    except ImportError:
        raise ImportError(
            "meshio is required for this function. "
            "Install it with: pip install meshio"
        ) from None


def from_meshio(mesh) -> Mesh:
    """Convert a meshio Mesh to a yggdrasil Mesh.

    Cell blocks whose type is not in the supported list (e.g. ``vertex``,
    ``line`` when loading a 2-D surface mesh that includes boundary markers)
    are silently skipped.

    Parameters
    ----------
    mesh : meshio.Mesh
        The source mesh.

    Returns
    -------
    Mesh
        A yggdrasil Mesh with one ElementGroup per supported cell block.

    Raises
    ------
    ImportError
        If meshio is not installed.
    """
    _require_meshio()

    nodes = np.asarray(mesh.points, dtype=np.float64)

    element_groups = []
    for block in mesh.cells:
        elem_class = _MESHIO_TO_ELEMENT.get(block.type)
        if elem_class is None:
            continue
        conn = np.asarray(block.data, dtype=np.intp)
        element_groups.append(ElementGroup(element=elem_class(), connectivity=conn))

    return Mesh(nodes, element_groups)


def to_meshio(mesh: Mesh):
    """Convert a yggdrasil Mesh to a meshio Mesh.

    Parameters
    ----------
    mesh : Mesh
        The yggdrasil mesh to export.

    Returns
    -------
    meshio.Mesh
        A meshio Mesh that can be written to disk with ``mesh.write(...)``.

    Raises
    ------
    ImportError
        If meshio is not installed.
    TypeError
        If any element group uses an element type that has no meshio
        equivalent.
    """
    _require_meshio()
    import meshio

    cells = []
    for group in mesh.element_groups:
        cell_type = _ELEMENT_TO_MESHIO.get(type(group.element))
        if cell_type is None:
            raise TypeError(
                f"No meshio cell type for {type(group.element).__name__}"
            )
        cells.append(meshio.CellBlock(cell_type, group.connectivity))

    return meshio.Mesh(
        points=np.array(mesh.nodes),
        cells=cells,
        point_data={k: np.array(v) for k, v in mesh.point_data.items()},
    )
