from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from ..refdomains.quadrilateral import QuadrilateralDomain
from .element import ReferenceElement


class Quad9(ReferenceElement):
    """9-node biquadratic quadrilateral on [0,1]^2.

    Node ordering:
      0: (0,0), 1: (1,0), 2: (1,1), 3: (0,1)       -- corners
      4: (0.5,0), 5: (1,0.5), 6: (0.5,1), 7: (0,0.5) -- edge midpoints
      8: (0.5,0.5)                                      -- center
    """

    @property
    def domain(self) -> QuadrilateralDomain:
        return QuadrilateralDomain()

    @property
    def num_nodes(self) -> int:
        return 9

    @property
    def node_coords(self) -> NDArray[np.float64]:
        return np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.0],
            [1.0, 0.5],
            [0.5, 1.0],
            [0.0, 0.5],
            [0.5, 0.5],
        ])

    @staticmethod
    def _lagrange_1d(
        t: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """1D quadratic Lagrange basis on [0,1] at nodes 0, 1, 0.5."""
        l0 = (1.0 - t) * (1.0 - 2.0 * t)  # node at 0
        l1 = t * (2.0 * t - 1.0)  # node at 1
        l2 = 4.0 * t * (1.0 - t)  # node at 0.5
        return l0, l1, l2

    @staticmethod
    def _lagrange_1d_deriv(
        t: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Derivatives of 1D quadratic Lagrange basis on [0,1]."""
        dl0 = 4.0 * t - 3.0
        dl1 = 4.0 * t - 1.0
        dl2 = 4.0 - 8.0 * t
        return dl0, dl1, dl2

    def shape_functions(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        lx0, lx1, lx2 = self._lagrange_1d(x)
        ly0, ly1, ly2 = self._lagrange_1d(y)

        # Tensor product, matching node ordering
        return np.column_stack([
            lx0 * ly0,  # node 0: (0,0)
            lx1 * ly0,  # node 1: (1,0)
            lx1 * ly1,  # node 2: (1,1)
            lx0 * ly1,  # node 3: (0,1)
            lx2 * ly0,  # node 4: (0.5,0)
            lx1 * ly2,  # node 5: (1,0.5)
            lx2 * ly1,  # node 6: (0.5,1)
            lx0 * ly2,  # node 7: (0,0.5)
            lx2 * ly2,  # node 8: (0.5,0.5)
        ])

    def shape_function_gradients(self, xi: NDArray[np.float64]) -> NDArray[np.float64]:
        x, y = xi[:, 0], xi[:, 1]
        n = xi.shape[0]

        lx0, lx1, lx2 = self._lagrange_1d(x)
        ly0, ly1, ly2 = self._lagrange_1d(y)
        dlx0, dlx1, dlx2 = self._lagrange_1d_deriv(x)
        dly0, dly1, dly2 = self._lagrange_1d_deriv(y)

        grad = np.empty((n, 9, 2))

        # (dN/dx, dN/dy) for each node
        nodes = [
            (dlx0, ly0, lx0, dly0),
            (dlx1, ly0, lx1, dly0),
            (dlx1, ly1, lx1, dly1),
            (dlx0, ly1, lx0, dly1),
            (dlx2, ly0, lx2, dly0),
            (dlx1, ly2, lx1, dly2),
            (dlx2, ly1, lx2, dly1),
            (dlx0, ly2, lx0, dly2),
            (dlx2, ly2, lx2, dly2),
        ]

        for i, (dlxi, lyj, lxi, dlyj) in enumerate(nodes):
            grad[:, i, 0] = dlxi * lyj
            grad[:, i, 1] = lxi * dlyj

        return grad

    @property
    def faces(self) -> tuple[tuple[int, ...], ...]:
        return ((0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7))

    @property
    def face_element(self) -> ReferenceElement:
        from .line3 import Line3

        return Line3()

    @property
    def polynomial_degree(self) -> int:
        return 2
