import numpy as np
from numpy.typing import NDArray

from .domain import ReferenceDomain


class TriangleDomain(ReferenceDomain):
    """Reference triangle with vertices (0,0), (1,0), (0,1)."""

    @property
    def topological_dimension(self) -> int:
        return 2

    def quadrature(
        self, order: int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Symmetric quadrature rules for the reference triangle.

        Area of reference triangle = 0.5.
        Rules from Dunavant (1985) and other standard sources.
        """
        if order <= 1:
            return self._order1()
        elif order <= 2:
            return self._order2()
        elif order <= 3:
            return self._order3()
        elif order <= 4:
            return self._order4()
        elif order <= 5:
            return self._order5()
        else:
            raise ValueError(
                f"Triangle quadrature not implemented for order {order} (max 5)"
            )

    def _order1(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """1-point rule, exact for order 1."""
        points = np.array([[1.0 / 3.0, 1.0 / 3.0]])
        weights = np.array([0.5])
        return points, weights

    def _order2(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """3-point rule, exact for order 2."""
        points = np.array([
            [1.0 / 6.0, 1.0 / 6.0],
            [2.0 / 3.0, 1.0 / 6.0],
            [1.0 / 6.0, 2.0 / 3.0],
        ])
        weights = np.full(3, 1.0 / 6.0)
        return points, weights

    def _order3(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """4-point rule, exact for order 3."""
        points = np.array([
            [1.0 / 3.0, 1.0 / 3.0],
            [1.0 / 5.0, 1.0 / 5.0],
            [3.0 / 5.0, 1.0 / 5.0],
            [1.0 / 5.0, 3.0 / 5.0],
        ])
        weights = np.array([
            -27.0 / 96.0,
            25.0 / 96.0,
            25.0 / 96.0,
            25.0 / 96.0,
        ])
        return points, weights

    def _order4(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """6-point rule, exact for order 4 (Dunavant)."""
        a1 = 0.445948490915965
        b1 = 0.108103018168070
        a2 = 0.091576213509771
        b2 = 0.816847572980459
        w1 = 0.111690794839006
        w2 = 0.054975871827661
        points = np.array([
            [a1, a1],
            [b1, a1],
            [a1, b1],
            [a2, a2],
            [b2, a2],
            [a2, b2],
        ])
        # Dunavant weights already sum to triangle area (0.5)
        weights = np.array([w1, w1, w1, w2, w2, w2])
        return points, weights

    def _order5(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """7-point rule, exact for order 5 (Dunavant)."""
        a0 = 1.0 / 3.0
        a1 = 0.470142064105115
        b1 = 0.059715871789770
        a2 = 0.101286507323456
        b2 = 0.797426985353087
        w0 = 0.225000000000000
        w1 = 0.132394152788506
        w2 = 0.125939180544827
        points = np.array([
            [a0, a0],
            [a1, a1],
            [b1, a1],
            [a1, b1],
            [a2, a2],
            [b2, a2],
            [a2, b2],
        ])
        # These weights sum to 1; multiply by area of reference triangle
        weights = np.array([w0, w1, w1, w1, w2, w2, w2]) * 0.5
        return points, weights
