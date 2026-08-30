from collections import Counter
from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

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

        for cell_nodes in group.connectivity:
            for face in faces:
                global_nodes = tuple(int(cell_nodes[i]) for i in face)
                key = tuple(sorted(global_nodes))
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


def tag_boundary_faces(
    boundary_mesh: Mesh,
    predicate: Callable[[NDArray[np.float64]], NDArray[np.bool_]],
    tag: int,
) -> Mesh:
    """Assign an integer tag to boundary faces whose centroid satisfies a predicate.

    Parameters
    ----------
    boundary_mesh : Mesh
        A boundary mesh as returned by ``extract_boundary``.
    predicate : callable
        Signature: ``predicate(centroids) -> bool_array`` where ``centroids``
        has shape ``(num_faces, spatial_dim)`` and the return value has shape
        ``(num_faces,)``.
    tag : int
        The integer tag to assign to matching faces. Faces not matched retain
        their existing tag (default 0 if never tagged).

    Returns
    -------
    Mesh
        A new boundary mesh with ``cell_data["tag"]`` updated for each element
        group. Existing tags on non-matching faces are preserved.
    """
    new_groups = []
    for group in boundary_mesh.element_groups:
        coords = boundary_mesh.nodes[group.connectivity]  # (num_elems, npe, spatial_dim)
        centroids = coords.mean(axis=1)  # (num_elems, spatial_dim)
        mask = np.asarray(predicate(centroids), dtype=bool)

        existing = group.cell_data.get("tag", np.zeros(group.num_elements, dtype=np.intp))
        new_tags = existing.copy()
        new_tags[mask] = tag

        new_cell_data = dict(group.cell_data)
        new_cell_data["tag"] = new_tags
        new_groups.append(ElementGroup(element=group.element, connectivity=group.connectivity, cell_data=new_cell_data))

    return Mesh(boundary_mesh.nodes, new_groups, point_data=dict(boundary_mesh.point_data))


def select_boundary_faces(boundary_mesh: Mesh, tag: int) -> Mesh:
    """Extract a sub-mesh of boundary faces with the given tag.

    Parameters
    ----------
    boundary_mesh : Mesh
        A tagged boundary mesh, e.g. from ``tag_boundary_faces``. Faces with
        no ``cell_data["tag"]`` entry are treated as having tag 0.
    tag : int
        The tag value to select.

    Returns
    -------
    Mesh
        A sub-mesh containing only the faces with the given tag, with
        ``point_data["original_node_index"]`` mapping local nodes to global
        volume-mesh node indices.
    """
    bnd_original = boundary_mesh.point_data["original_node_index"]

    selected_local_nodes: set[int] = set()
    selected_groups_data: list[tuple[ReferenceElement, NDArray[np.intp]]] = []

    for group in boundary_mesh.element_groups:
        tags = group.cell_data.get("tag", np.zeros(group.num_elements, dtype=np.intp))
        mask = np.asarray(tags) == tag
        if not mask.any():
            continue
        selected_conn = group.connectivity[mask]
        selected_local_nodes.update(selected_conn.ravel().tolist())
        selected_groups_data.append((group.element, selected_conn))

    if not selected_groups_data:
        return Mesh(
            np.empty((0, boundary_mesh.spatial_dim), dtype=np.float64),
            [],
            point_data={"original_node_index": np.array([], dtype=np.intp)},
        )

    sorted_local = sorted(selected_local_nodes)
    local_to_sub = {old: new for new, old in enumerate(sorted_local)}

    sub_nodes = boundary_mesh.nodes[sorted_local]
    sub_original = bnd_original[sorted_local]

    new_groups = []
    for element, old_conn in selected_groups_data:
        new_conn = np.array(
            [[local_to_sub[n] for n in row] for row in old_conn],
            dtype=np.intp,
        )
        new_groups.append(ElementGroup(element=element, connectivity=new_conn))

    return Mesh(sub_nodes, new_groups, point_data={"original_node_index": sub_original})
