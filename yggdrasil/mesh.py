from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from .elements.element import ReferenceElement


@dataclass
class ElementGroup:
    """A group of elements sharing the same type."""

    element: ReferenceElement
    connectivity: NDArray[np.intp]
    cell_data: dict[str, NDArray] = field(default_factory=dict)

    def __post_init__(self):
        if self.connectivity.ndim != 2:
            raise ValueError(
                f"connectivity must be 2D, got shape {self.connectivity.shape}"
            )
        if self.connectivity.shape[1] != self.element.num_nodes:
            raise ValueError(
                f"connectivity has {self.connectivity.shape[1]} columns, "
                f"expected {self.element.num_nodes} for {type(self.element).__name__}"
            )

    @property
    def num_elements(self) -> int:
        return self.connectivity.shape[0]

    def physical_coords(self, nodes: NDArray[np.float64]) -> NDArray[np.float64]:
        """Return physical coordinates for each element's nodes.

        Parameters
        ----------
        nodes : (num_nodes, spatial_dim)

        Returns
        -------
        coords : (num_elements, nodes_per_elem, spatial_dim)
        """
        return nodes[self.connectivity]


class Mesh:
    """A finite element mesh with nodes, element groups, and auxiliary data."""

    def __init__(
        self,
        nodes: NDArray[np.float64],
        element_groups: list[ElementGroup],
        point_data: dict[str, NDArray] | None = None,
    ):
        nodes = np.asarray(nodes, dtype=np.float64)
        if nodes.ndim != 2:
            raise ValueError(f"nodes must be 2D, got shape {nodes.shape}")
        self._nodes = nodes
        self._element_groups = list(element_groups)
        self._point_data = dict(point_data) if point_data else {}

    @property
    def nodes(self) -> NDArray[np.float64]:
        return self._nodes

    @property
    def num_nodes(self) -> int:
        return self._nodes.shape[0]

    @property
    def spatial_dim(self) -> int:
        return self._nodes.shape[1]

    @property
    def element_groups(self) -> list[ElementGroup]:
        return self._element_groups

    @property
    def num_elements(self) -> int:
        return sum(g.num_elements for g in self._element_groups)

    @property
    def point_data(self) -> dict[str, NDArray]:
        return self._point_data

    def iter_element_groups(self):
        """Iterate over element groups."""
        return iter(self._element_groups)

    def element_groups_by_type(self, element_type: type) -> list[ElementGroup]:
        """Return element groups matching the given element class."""
        return [g for g in self._element_groups if isinstance(g.element, element_type)]
