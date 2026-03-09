import numpy as np
from numpy.typing import NDArray

from .mesh import Mesh


class DOFMap:
    """Maps mesh nodes to global DOF indices.

    For n_components=1 (scalar), DOF index equals node index.
    For n_components=d (vector), node i maps to DOF block [d*i, ..., d*i+d-1].
    Only n_components=1 is currently supported.
    """

    def __init__(self, mesh: Mesh, n_components: int = 1) -> None:
        if n_components != 1:
            raise NotImplementedError("n_components > 1 not yet implemented")
        self.n_components = n_components
        self.n_dofs: int = mesh.num_nodes * n_components

    def element_dofs(self, elem_nodes: NDArray[np.intp]) -> NDArray[np.intp]:
        """Global DOF indices for the given element node indices.

        For n_components=1: identity. For n_components=d: interleaved blocks.
        """
        return elem_nodes  # n_components=1 fast path

    def boundary_dofs(self, original_node_index: NDArray[np.intp]) -> NDArray[np.intp]:
        """Convert original_node_index to global DOF indices.

        For n_components=1: identity.
        """
        return original_node_index
