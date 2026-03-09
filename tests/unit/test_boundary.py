import numpy as np
import numpy.testing as npt

from yggdrasil.boundary import extract_boundary, select_boundary_faces, tag_boundary_faces
from yggdrasil.elements import Hex8, Line2, Point1, Quad4, Tri3, Tet4
from yggdrasil.mesh import ElementGroup, Mesh
from yggdrasil.mesh_generators import unit_square_tri_mesh


def test_single_line2():
    """Single Line2 element: 2 boundary Point1 elements, 2 nodes."""
    nodes = np.array([[0.0], [1.0]])
    conn = np.array([[0, 1]], dtype=np.intp)
    mesh = Mesh(nodes, [ElementGroup(element=Line2(), connectivity=conn)])

    bnd = extract_boundary(mesh)

    assert bnd.num_nodes == 2
    assert len(bnd.element_groups) == 1
    group = bnd.element_groups[0]
    assert isinstance(group.element, Point1)
    assert group.num_elements == 2


def test_single_tri3():
    """Single Tri3 element: 3 boundary Line2 edges, 3 nodes."""
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    conn = np.array([[0, 1, 2]], dtype=np.intp)
    mesh = Mesh(nodes, [ElementGroup(element=Tri3(), connectivity=conn)])

    bnd = extract_boundary(mesh)

    assert bnd.num_nodes == 3
    assert len(bnd.element_groups) == 1
    group = bnd.element_groups[0]
    assert isinstance(group.element, Line2)
    assert group.num_elements == 3


def test_two_tri3_shared_edge():
    """Two Tri3 sharing an edge: shared edge is interior -> 4 boundary edges."""
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ])
    conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
    mesh = Mesh(nodes, [ElementGroup(element=Tri3(), connectivity=conn)])

    bnd = extract_boundary(mesh)

    assert bnd.num_nodes == 4
    assert len(bnd.element_groups) == 1
    group = bnd.element_groups[0]
    assert isinstance(group.element, Line2)
    assert group.num_elements == 4


def test_unit_square_mesh():
    """Unit square mesh with n=2: 8 boundary edges forming perimeter."""
    mesh = unit_square_tri_mesh(2)

    bnd = extract_boundary(mesh)

    assert len(bnd.element_groups) == 1
    group = bnd.element_groups[0]
    assert isinstance(group.element, Line2)
    assert group.num_elements == 8


def test_single_tet4():
    """Single Tet4 element: 4 boundary Tri3 faces, 4 nodes."""
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    conn = np.array([[0, 1, 2, 3]], dtype=np.intp)
    mesh = Mesh(nodes, [ElementGroup(element=Tet4(), connectivity=conn)])

    bnd = extract_boundary(mesh)

    assert bnd.num_nodes == 4
    assert len(bnd.element_groups) == 1
    group = bnd.element_groups[0]
    assert isinstance(group.element, Tri3)
    assert group.num_elements == 4


def test_single_hex8():
    """Single Hex8 element: 6 boundary Quad4 faces, 8 nodes."""
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ])
    conn = np.array([[0, 1, 2, 3, 4, 5, 6, 7]], dtype=np.intp)
    mesh = Mesh(nodes, [ElementGroup(element=Hex8(), connectivity=conn)])

    bnd = extract_boundary(mesh)

    assert bnd.num_nodes == 8
    assert len(bnd.element_groups) == 1
    group = bnd.element_groups[0]
    assert isinstance(group.element, Quad4)
    assert group.num_elements == 6


def test_original_node_index():
    """point_data['original_node_index'] correctly maps new -> original coords."""
    mesh = unit_square_tri_mesh(2)

    bnd = extract_boundary(mesh)

    orig_idx = bnd.point_data["original_node_index"]
    assert orig_idx.shape == (bnd.num_nodes,)
    npt.assert_array_equal(bnd.nodes, mesh.nodes[orig_idx])


def test_boundary_coords_match():
    """Boundary node coords match mesh.nodes[original_indices] for a tet."""
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    conn = np.array([[0, 1, 2, 3]], dtype=np.intp)
    mesh = Mesh(nodes, [ElementGroup(element=Tet4(), connectivity=conn)])

    bnd = extract_boundary(mesh)

    orig_idx = bnd.point_data["original_node_index"]
    npt.assert_array_equal(bnd.nodes, mesh.nodes[orig_idx])


# --- tag_boundary_faces / select_boundary_faces ---

def test_tag_boundary_faces_count():
    """Bottom edges (y≈0) of 2×2 unit-square mesh: exactly 2 tagged faces."""
    mesh = unit_square_tri_mesh(2)
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda x: np.isclose(x[:, 1], 0.0), tag=1)
    bottom = select_boundary_faces(bnd, tag=1)
    assert bottom.num_elements == 2


def test_tag_boundary_faces_coords():
    """All nodes in the selected sub-mesh lie on the tagged boundary."""
    mesh = unit_square_tri_mesh(2)
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda x: np.isclose(x[:, 1], 0.0), tag=1)
    bottom = select_boundary_faces(bnd, tag=1)
    assert np.allclose(bottom.nodes[:, 1], 0.0)


def test_select_boundary_faces_original_node_index():
    """original_node_index of the sub-mesh maps correctly to volume-mesh coords."""
    mesh = unit_square_tri_mesh(2)
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda x: np.isclose(x[:, 1], 0.0), tag=1)
    bottom = select_boundary_faces(bnd, tag=1)
    orig_idx = bottom.point_data["original_node_index"]
    npt.assert_array_equal(bottom.nodes, mesh.nodes[orig_idx])


def test_select_boundary_faces_unmatched_tag_is_empty():
    """Selecting a tag with no matches returns an empty mesh."""
    mesh = unit_square_tri_mesh(2)
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda x: np.isclose(x[:, 1], 0.0), tag=1)
    empty = select_boundary_faces(bnd, tag=99)
    assert empty.num_elements == 0
    assert empty.num_nodes == 0


def test_multiple_tags_disjoint():
    """Bottom (tag=1) and top (tag=2) together cover exactly the tagged faces; left/right remain untagged."""
    mesh = unit_square_tri_mesh(2)
    bnd = extract_boundary(mesh)
    bnd = tag_boundary_faces(bnd, lambda x: np.isclose(x[:, 1], 0.0), tag=1)
    bnd = tag_boundary_faces(bnd, lambda x: np.isclose(x[:, 1], 1.0), tag=2)

    bottom = select_boundary_faces(bnd, tag=1)
    top = select_boundary_faces(bnd, tag=2)
    untagged = select_boundary_faces(bnd, tag=0)

    assert bottom.num_elements == 2
    assert top.num_elements == 2
    assert untagged.num_elements == 4  # left + right sides

    assert np.allclose(bottom.nodes[:, 1], 0.0)
    assert np.allclose(top.nodes[:, 1], 1.0)
