from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..domains.tetrahedron import TetrahedronDomain
from .element import ReferenceElement


class Tet4(ReferenceElement):
    """4-node linear tetrahedron.

    Vertices: (0,0,0), (1,0,0), (0,1,0), (0,0,1).
    """

    @property
    def domain(self) -> TetrahedronDomain:
        return TetrahedronDomain()

    @property
    def num_nodes(self) -> int:
        return 4

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y, z = xi[:, 0], xi[:, 1], xi[:, 2]
        N0 = 1.0 - x - y - z
        N1 = x
        N2 = y
        N3 = z
        return np.column_stack([N0, N1, N2, N3])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        n = xi.shape[0]
        grad = np.empty((n, 4, 3))
        grad[:, 0, :] = -1.0
        grad[:, 1, 0] = 1.0
        grad[:, 1, 1] = 0.0
        grad[:, 1, 2] = 0.0
        grad[:, 2, 0] = 0.0
        grad[:, 2, 1] = 1.0
        grad[:, 2, 2] = 0.0
        grad[:, 3, 0] = 0.0
        grad[:, 3, 1] = 0.0
        grad[:, 3, 2] = 1.0
        return grad

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))

    @property
    def face_element(self) -> ReferenceElement:
        from .tri3 import Tri3

        return Tri3()

    @property
    def polynomial_degree(self) -> int:
        return 1
