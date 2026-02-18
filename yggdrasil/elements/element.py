from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray

from ..domains.domain import ReferenceDomain


class ReferenceElement(ABC):
    @property
    @abstractmethod
    def domain(self) -> ReferenceDomain: ...

    @property
    @abstractmethod
    def num_nodes(self) -> int: ...

    @property
    @abstractmethod
    def node_coords(self) -> NDArray[np.float64]:
        """Node coordinates on the reference domain. Shape: (num_nodes, topo_dim)."""
        ...

    @abstractmethod
    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluate shape functions at given points.

        xi: (num_points, topo_dim) -> returns (num_points, num_nodes)
        """
        ...

    @abstractmethod
    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        """Evaluate shape function gradients at given points.

        xi: (num_points, topo_dim) -> returns (num_points, num_nodes, topo_dim)
        """
        ...

    @property
    @abstractmethod
    def faces(self) -> tuple[tuple[int, ...], ...]:
        """Local node indices for each face/edge of the element."""
        ...

    @property
    @abstractmethod
    def face_element(self) -> ReferenceElement | None:
        """The element type used for each face, or None for 0D elements."""
        ...
