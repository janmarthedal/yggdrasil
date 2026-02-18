from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..domains.hexahedron import HexahedronDomain
from .element import ReferenceElement


class Hex8(ReferenceElement):
    """8-node trilinear hexahedron on [0,1]^3.

    Node ordering:
      0: (0,0,0), 1: (1,0,0), 2: (1,1,0), 3: (0,1,0)  -- bottom face
      4: (0,0,1), 5: (1,0,1), 6: (1,1,1), 7: (0,1,1)  -- top face
    """

    @property
    def domain(self) -> HexahedronDomain:
        return HexahedronDomain()

    @property
    def num_nodes(self) -> int:
        return 8

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ])

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y, z = xi[:, 0], xi[:, 1], xi[:, 2]
        mx, my, mz = 1.0 - x, 1.0 - y, 1.0 - z

        return np.column_stack([
            mx * my * mz,  # 0
            x * my * mz,   # 1
            x * y * mz,    # 2
            mx * y * mz,   # 3
            mx * my * z,   # 4
            x * my * z,    # 5
            x * y * z,     # 6
            mx * y * z,    # 7
        ])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y, z = xi[:, 0], xi[:, 1], xi[:, 2]
        mx, my, mz = 1.0 - x, 1.0 - y, 1.0 - z
        n = xi.shape[0]

        grad = np.empty((n, 8, 3))

        # dN/dx, dN/dy, dN/dz for each node
        grad[:, 0, 0] = -my * mz
        grad[:, 0, 1] = -mx * mz
        grad[:, 0, 2] = -mx * my

        grad[:, 1, 0] = my * mz
        grad[:, 1, 1] = -x * mz
        grad[:, 1, 2] = -x * my

        grad[:, 2, 0] = y * mz
        grad[:, 2, 1] = x * mz
        grad[:, 2, 2] = -x * y

        grad[:, 3, 0] = -y * mz
        grad[:, 3, 1] = mx * mz
        grad[:, 3, 2] = -mx * y

        grad[:, 4, 0] = -my * z
        grad[:, 4, 1] = -mx * z
        grad[:, 4, 2] = mx * my

        grad[:, 5, 0] = my * z
        grad[:, 5, 1] = -x * z
        grad[:, 5, 2] = x * my

        grad[:, 6, 0] = y * z
        grad[:, 6, 1] = x * z
        grad[:, 6, 2] = x * y

        grad[:, 7, 0] = -y * z
        grad[:, 7, 1] = mx * z
        grad[:, 7, 2] = mx * y

        return grad

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return (
            (0, 3, 2, 1),
            (4, 5, 6, 7),
            (0, 1, 5, 4),
            (2, 3, 7, 6),
            (0, 4, 7, 3),
            (1, 2, 6, 5),
        )

    @property
    def face_element(self) -> ReferenceElement:
        from .quad4 import Quad4

        return Quad4()
