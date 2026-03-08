import numpy as np
import pytest

from yggdrasil.elements import Tri3, Quad4
from yggdrasil.mesh import ElementGroup, Mesh


class TestElementGroup:
    def test_wrong_connectivity_width(self):
        with pytest.raises(ValueError, match="expected 3"):
            ElementGroup(
                element=Tri3(),
                connectivity=np.array([[0, 1, 2, 3]], dtype=np.intp),
            )

    def test_non_2d_connectivity(self):
        with pytest.raises(ValueError, match="must be 2D"):
            ElementGroup(
                element=Tri3(),
                connectivity=np.array([0, 1, 2], dtype=np.intp),
            )

    def test_num_elements(self):
        group = ElementGroup(
            element=Tri3(),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.intp),
        )
        assert group.num_elements == 2

    def test_physical_coords(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [1.5, 1.0],
        ])
        group = ElementGroup(
            element=Tri3(),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.intp),
        )
        coords = group.physical_coords(nodes)
        assert coords.shape == (2, 3, 2)
        np.testing.assert_array_equal(coords[0], nodes[[0, 1, 2]])
        np.testing.assert_array_equal(coords[1], nodes[[1, 3, 2]])

    def test_cell_data(self):
        group = ElementGroup(
            element=Tri3(),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.intp),
            cell_data={"material": np.array([1, 2])},
        )
        np.testing.assert_array_equal(group.cell_data["material"], [1, 2])


class TestMesh:
    @pytest.fixture
    def tri_mesh(self):
        """Two Tri3 elements sharing an edge."""
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, 1.0],
            [1.5, 1.0],
        ])
        group = ElementGroup(
            element=Tri3(),
            connectivity=np.array([[0, 1, 2], [1, 3, 2]], dtype=np.intp),
        )
        return Mesh(nodes, [group])

    @pytest.fixture
    def mixed_mesh(self):
        """Tri3 + Quad4 sharing nodes."""
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
        ])
        tri_group = ElementGroup(
            element=Tri3(),
            connectivity=np.array([[1, 4, 5]], dtype=np.intp),
        )
        quad_group = ElementGroup(
            element=Quad4(),
            connectivity=np.array([[0, 1, 2, 3]], dtype=np.intp),
        )
        return Mesh(nodes, [tri_group, quad_group])

    def test_num_nodes(self, tri_mesh):
        assert tri_mesh.num_nodes == 4

    def test_num_elements(self, tri_mesh):
        assert tri_mesh.num_elements == 2

    def test_spatial_dim(self, tri_mesh):
        assert tri_mesh.spatial_dim == 2

    def test_mixed_num_elements(self, mixed_mesh):
        assert mixed_mesh.num_elements == 2

    def test_mixed_num_nodes(self, mixed_mesh):
        assert mixed_mesh.num_nodes == 6

    def test_element_groups_by_type(self, mixed_mesh):
        tri_groups = mixed_mesh.element_groups_by_type(Tri3)
        quad_groups = mixed_mesh.element_groups_by_type(Quad4)
        assert len(tri_groups) == 1
        assert len(quad_groups) == 1
        assert tri_groups[0].num_elements == 1
        assert quad_groups[0].num_elements == 1

    def test_element_groups_by_type_no_match(self, tri_mesh):
        assert tri_mesh.element_groups_by_type(Quad4) == []

    def test_iter_element_groups(self, mixed_mesh):
        groups = list(mixed_mesh.iter_element_groups())
        assert len(groups) == 2

    def test_physical_coords_via_mesh(self, tri_mesh):
        group = tri_mesh.element_groups[0]
        coords = group.physical_coords(tri_mesh.nodes)
        assert coords.shape == (2, 3, 2)
        np.testing.assert_array_equal(coords[0, 0], [0.0, 0.0])
        np.testing.assert_array_equal(coords[0, 1], [1.0, 0.0])

    def test_point_data(self):
        nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, 1.0]])
        group = ElementGroup(
            element=Tri3(),
            connectivity=np.array([[0, 1, 2]], dtype=np.intp),
        )
        temps = np.array([100.0, 200.0, 150.0])
        mesh = Mesh(nodes, [group], point_data={"temperature": temps})
        np.testing.assert_array_equal(mesh.point_data["temperature"], temps)

    def test_non_2d_nodes_raises(self):
        with pytest.raises(ValueError, match="must be 2D"):
            Mesh(
                nodes=np.array([0.0, 1.0, 2.0]),
                element_groups=[],
            )

    def test_3d_nodes_raises_no_error(self):
        nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        group = ElementGroup(
            element=Tri3(),
            connectivity=np.zeros((0, 3), dtype=np.intp),
        )
        mesh = Mesh(nodes, [group])
        assert mesh.spatial_dim == 3
