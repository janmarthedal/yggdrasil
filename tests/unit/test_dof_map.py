import numpy as np
import pytest

from yggdrasil import DOFMap
from yggdrasil.mesh_generators import unit_square_tri_mesh


@pytest.fixture
def mesh():
    return unit_square_tri_mesh(2)


def test_n_dofs_equals_num_nodes(mesh):
    dof_map = DOFMap(mesh)
    assert dof_map.n_dofs == mesh.num_nodes


def test_n_components_default(mesh):
    dof_map = DOFMap(mesh)
    assert dof_map.n_components == 1


def test_element_dofs_identity(mesh):
    dof_map = DOFMap(mesh)
    arr = np.array([0, 2, 5], dtype=np.intp)
    result = dof_map.element_dofs(arr)
    np.testing.assert_array_equal(result, arr)


def test_boundary_dofs_identity(mesh):
    dof_map = DOFMap(mesh)
    arr = np.array([1, 3, 7], dtype=np.intp)
    result = dof_map.boundary_dofs(arr)
    np.testing.assert_array_equal(result, arr)


def test_n_components_2_raises(mesh):
    with pytest.raises(NotImplementedError):
        DOFMap(mesh, n_components=2)
