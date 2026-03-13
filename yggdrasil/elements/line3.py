from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..refdomains.line import LineDomain
from .element import ReferenceElement


class Line3(ReferenceElement):
    """3-node quadratic line element on [0, 1]."""

    @property
    def domain(self) -> LineDomain:
        return LineDomain()

    @property
    def num_nodes(self) -> int:
        return 3

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([[0.0], [1.0], [0.5]])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x = xi[:, 0]
        # N0 = (1 - x)(1 - 2x) = 2x^2 - 3x + 1
        N0 = (1.0 - x) * (1.0 - 2.0 * x)
        # N1 = x(2x - 1)
        N1 = x * (2.0 * x - 1.0)
        # N2 = 4x(1 - x)
        N2 = 4.0 * x * (1.0 - x)
        return np.column_stack([N0, N1, N2])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x = xi[:, 0]
        n = xi.shape[0]
        grad = np.empty((n, 3, 1))
        grad[:, 0, 0] = 4.0 * x - 3.0
        grad[:, 1, 0] = 4.0 * x - 1.0
        grad[:, 2, 0] = 4.0 - 8.0 * x
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
        return 2
