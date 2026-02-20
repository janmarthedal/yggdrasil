from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..domains.line import LineDomain
from .element import ReferenceElement


class Line2(ReferenceElement):
    """2-node linear line element on [0, 1]."""

    @property
    def domain(self) -> LineDomain:
        return LineDomain()

    @property
    def num_nodes(self) -> int:
        return 2

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([[0.0], [1.0]])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x = xi[:, 0]
        N0 = 1.0 - x
        N1 = x
        return np.column_stack([N0, N1])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        n = xi.shape[0]
        grad = np.empty((n, 2, 1))
        grad[:, 0, 0] = -1.0
        grad[:, 1, 0] = 1.0
        return grad

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return ((0,), (1,))

    @property
    def face_element(self) -> ReferenceElement:
        from .point1 import Point1

        return Point1()

    @property
    def polynomial_degree(self) -> int:
        return 1
