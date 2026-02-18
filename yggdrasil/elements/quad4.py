from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..domains.quadrilateral import QuadrilateralDomain
from .element import ReferenceElement


class Quad4(ReferenceElement):
    """4-node bilinear quadrilateral on [0,1]^2.

    Node ordering:
      0: (0,0), 1: (1,0), 2: (1,1), 3: (0,1)
    """

    @property
    def domain(self) -> QuadrilateralDomain:
        return QuadrilateralDomain()

    @property
    def num_nodes(self) -> int:
        return 4

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        N0 = (1.0 - x) * (1.0 - y)
        N1 = x * (1.0 - y)
        N2 = x * y
        N3 = (1.0 - x) * y
        return np.column_stack([N0, N1, N2, N3])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        n = xi.shape[0]
        grad = np.empty((n, 4, 2))

        grad[:, 0, 0] = -(1.0 - y)
        grad[:, 0, 1] = -(1.0 - x)

        grad[:, 1, 0] = 1.0 - y
        grad[:, 1, 1] = -x

        grad[:, 2, 0] = y
        grad[:, 2, 1] = x

        grad[:, 3, 0] = -y
        grad[:, 3, 1] = 1.0 - x

        return grad

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return ((0, 1), (1, 2), (2, 3), (3, 0))

    @property
    def face_element(self) -> ReferenceElement:
        from .line2 import Line2

        return Line2()
