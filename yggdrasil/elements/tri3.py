import numpy as np
from numpy.typing import NDArray

from ..domains.triangle import TriangleDomain
from .element import ReferenceElement


class Tri3(ReferenceElement):
    """3-node linear triangle on reference triangle (0,0)-(1,0)-(0,1)."""

    @property
    def domain(self) -> TriangleDomain:
        return TriangleDomain()

    @property
    def num_nodes(self) -> int:
        return 3

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        N0 = 1.0 - x - y
        N1 = x
        N2 = y
        return np.column_stack([N0, N1, N2])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        n = xi.shape[0]
        grad = np.empty((n, 3, 2))
        # dN0/dx = -1, dN0/dy = -1
        grad[:, 0, 0] = -1.0
        grad[:, 0, 1] = -1.0
        # dN1/dx = 1, dN1/dy = 0
        grad[:, 1, 0] = 1.0
        grad[:, 1, 1] = 0.0
        # dN2/dx = 0, dN2/dy = 1
        grad[:, 2, 0] = 0.0
        grad[:, 2, 1] = 1.0
        return grad
