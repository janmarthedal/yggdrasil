"""Simple example: create a small 2D mesh with triangles and a quad."""

import numpy as np

from yggdrasil import ElementGroup, Mesh
from yggdrasil.elements import Quad4, Tri3

# Define nodes for a small domain:
#
#  3-------2-------5
#  |       | T0  / |
#  |  Q0   |   /   |
#  |       | /  T1 |
#  0-------1-------4
#
nodes = np.array([
    [0.0, 0.0],  # 0
    [1.0, 0.0],  # 1
    [1.0, 1.0],  # 2
    [0.0, 1.0],  # 3
    [2.0, 0.0],  # 4
    [2.0, 1.0],  # 5
])

quad_group = ElementGroup(
    element=Quad4(),
    connectivity=np.array([[0, 1, 2, 3]], dtype=np.intp),
    cell_data={"material_id": np.array([1])},
)

tri_group = ElementGroup(
    element=Tri3(),
    connectivity=np.array([[1, 5, 2], [4, 5, 1]], dtype=np.intp),
    cell_data={"material_id": np.array([2, 2])},
)

mesh = Mesh(nodes, [quad_group, tri_group])

print(f"Nodes: {mesh.num_nodes}")
print(f"Elements: {mesh.num_elements}")
print(f"Spatial dimension: {mesh.spatial_dim}")
print(f"Element groups: {len(mesh.element_groups)}")

for group in mesh.iter_element_groups():
    name = type(group.element).__name__
    coords = group.physical_coords(mesh.nodes)
    print(f"\n  {name}: {group.num_elements} element(s)")
    print(f"    Physical coords shape: {coords.shape}")
    print(f"    Material IDs: {group.cell_data['material_id']}")
