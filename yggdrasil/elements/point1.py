from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..domains.point import PointDomain
from .element import ReferenceElement


class Point1(ReferenceElement):
    """1-node point element (0-dimensional)."""

    @property
    def domain(self) -> PointDomain:
        return PointDomain()

    @property
    def num_nodes(self) -> int:
        return 1

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.empty((1, 0), dtype=np.float64)

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        n = xi.shape[0]
        return np.ones((n, 1))

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        n = xi.shape[0]
        return np.empty((n, 1, 0), dtype=np.float64)

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return ()

    @property
    def face_element(self) -> ReferenceElement | None:
        return None
