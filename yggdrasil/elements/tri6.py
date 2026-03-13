from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..refdomains.triangle import TriangleDomain
from .element import ReferenceElement


class Tri6(ReferenceElement):
    """6-node quadratic triangle on reference triangle (0,0)-(1,0)-(0,1).

    Node ordering: vertices first, then edge midpoints.
      0: (0,0), 1: (1,0), 2: (0,1), 3: (0.5,0), 4: (0.5,0.5), 5: (0,0.5)
    """

    @property
    def domain(self) -> TriangleDomain:
        return TriangleDomain()

    @property
    def num_nodes(self) -> int:
        return 6

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [0.5, 0.5],
            [0.0, 0.5],
        ])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        L0 = 1.0 - x - y  # barycentric coord for vertex 0
        L1 = x
        L2 = y

        N0 = L0 * (2.0 * L0 - 1.0)
        N1 = L1 * (2.0 * L1 - 1.0)
        N2 = L2 * (2.0 * L2 - 1.0)
        N3 = 4.0 * L0 * L1
        N4 = 4.0 * L1 * L2
        N5 = 4.0 * L0 * L2

        return np.column_stack([N0, N1, N2, N3, N4, N5])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        n = xi.shape[0]
        L0 = 1.0 - x - y

        grad = np.empty((n, 6, 2))

        # dN0/dx, dN0/dy: N0 = L0(2L0-1), dL0/dx = -1, dL0/dy = -1
        grad[:, 0, 0] = -(4.0 * L0 - 1.0)
        grad[:, 0, 1] = -(4.0 * L0 - 1.0)

        # dN1/dx, dN1/dy: N1 = L1(2L1-1), dL1/dx = 1, dL1/dy = 0
        grad[:, 1, 0] = 4.0 * x - 1.0
        grad[:, 1, 1] = 0.0

        # dN2/dx, dN2/dy: N2 = L2(2L2-1), dL2/dx = 0, dL2/dy = 1
        grad[:, 2, 0] = 0.0
        grad[:, 2, 1] = 4.0 * y - 1.0

        # N3 = 4*L0*L1, dN3/dx = 4(-L1 + L0) = 4(1-2x-y), dN3/dy = 4(-L1) = -4x
        grad[:, 3, 0] = 4.0 * (1.0 - 2.0 * x - y)
        grad[:, 3, 1] = -4.0 * x

        # N4 = 4*L1*L2, dN4/dx = 4*y, dN4/dy = 4*x
        grad[:, 4, 0] = 4.0 * y
        grad[:, 4, 1] = 4.0 * x

        # N5 = 4*L0*L2, dN5/dx = -4*y, dN5/dy = 4*(1-x-2y)
        grad[:, 5, 0] = -4.0 * y
        grad[:, 5, 1] = 4.0 * (1.0 - x - 2.0 * y)

        return grad

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return ((0, 1, 3), (1, 2, 4), (2, 0, 5))

    @property
    def face_element(self) -> ReferenceElement:
        from .line3 import Line3

        return Line3()

    @property
    def polynomial_degree(self) -> int:
        return 2
