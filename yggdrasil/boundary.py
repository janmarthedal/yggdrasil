from collections import Counter

import numpy as np

from .elements.element import ReferenceElement
from .mesh import ElementGroup, Mesh


def extract_boundary(mesh: Mesh) -> Mesh:
    """Extract the boundary of a mesh as a new lower-dimensional mesh.

    Boundary faces are those that appear exactly once across all elements.
    The returned mesh contains only boundary nodes (renumbered), with
    ``point_data["original_node_index"]`` mapping back to the input mesh.
    """
    # Collect all faces: (canonical_key, full_global_nodes, face_element)
    face_records: list[tuple[tuple[int, ...], tuple[int, ...], ReferenceElement]] = []
    canonical_keys: list[tuple[int, ...]] = []

    for group in mesh.element_groups:
        elem = group.element
        face_elem = elem.face_element
        if face_elem is None:
            continue

        faces = elem.faces
        # Number of corner nodes used for the canonical key
        num_corner = face_elem.domain.topological_dimension + 1
        if num_corner < 1:
            num_corner = 1

        for cell_nodes in group.connectivity:
            for face in faces:
                global_nodes = tuple(int(cell_nodes[i]) for i in face)
                key = tuple(sorted(global_nodes[:num_corner]))
                canonical_keys.append(key)
                face_records.append((key, global_nodes, face_elem))

    # Count occurrences of each canonical key
    key_counts = Counter(canonical_keys)

    # Boundary faces are those appearing exactly once
    boundary_faces: list[tuple[tuple[int, ...], ReferenceElement]] = []
    for key, global_nodes, face_elem in face_records:
        if key_counts[key] == 1:
            boundary_faces.append((global_nodes, face_elem))

    if not boundary_faces:
        return Mesh(
            np.empty((0, mesh.spatial_dim), dtype=np.float64),
            [],
            point_data={"original_node_index": np.array([], dtype=np.intp)},
        )

    # Collect unique boundary node indices
    all_boundary_nodes: set[int] = set()
    for global_nodes, _ in boundary_faces:
        all_boundary_nodes.update(global_nodes)

    sorted_boundary_nodes = sorted(all_boundary_nodes)
    old_to_new = {old: new for new, old in enumerate(sorted_boundary_nodes)}

    new_nodes = mesh.nodes[sorted_boundary_nodes]

    # Group boundary faces by face element type
    groups_by_type: dict[type, list[tuple[tuple[int, ...], ReferenceElement]]] = {}
    for global_nodes, face_elem in boundary_faces:
        elem_type = type(face_elem)
        if elem_type not in groups_by_type:
            groups_by_type[elem_type] = []
        groups_by_type[elem_type].append((global_nodes, face_elem))

    element_groups = []
    for elem_type, faces_list in groups_by_type.items():
        face_elem = faces_list[0][1]
        conn = np.array(
            [[old_to_new[n] for n in global_nodes] for global_nodes, _ in faces_list],
            dtype=np.intp,
        )
        element_groups.append(ElementGroup(element=face_elem, connectivity=conn))

    return Mesh(
        new_nodes,
        element_groups,
        point_data={
            "original_node_index": np.array(sorted_boundary_nodes, dtype=np.intp)
        },
    )
