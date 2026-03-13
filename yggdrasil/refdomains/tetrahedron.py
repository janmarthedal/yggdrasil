import numpy as np
from numpy.typing import NDArray

from .refdomain import ReferenceDomain


class TetrahedronDomain(ReferenceDomain):
    """Reference tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).

    Volume = 1/6.
    """

    @property
    def topological_dimension(self) -> int:
        return 3

    def quadrature(
        self, order: int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Symmetric quadrature rules for the reference tetrahedron.

        Rules from Keast (1986) and other standard sources.
        """
        if order <= 1:
            return self._order1()
        elif order <= 2:
            return self._order2()
        elif order <= 3:
            return self._order3()
        else:
            raise ValueError(
                f"Tetrahedron quadrature not implemented for order {order} (max 3)"
            )

    def _order1(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """1-point rule, exact for order 1."""
        points = np.array([[0.25, 0.25, 0.25]])
        weights = np.array([1.0 / 6.0])
        return points, weights

    def _order2(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """4-point rule, exact for order 2 (Keast)."""
        a = 0.1381966011250105
        b = 0.5854101966249685
        points = np.array([
            [a, a, a],
            [b, a, a],
            [a, b, a],
            [a, a, b],
        ])
        weights = np.full(4, 1.0 / 24.0)
        return points, weights

    def _order3(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """5-point rule, exact for order 3 (Keast)."""
        points = np.array([
            [0.25, 0.25, 0.25],
            [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
            [0.5, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 0.5, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 0.5],
        ])
        weights = np.array([
            -2.0 / 15.0,
            3.0 / 40.0,
            3.0 / 40.0,
            3.0 / 40.0,
            3.0 / 40.0,
        ])
        return points, weights
