import numpy as np
from numpy.typing import NDArray

from .refdomain import ReferenceDomain


class LineDomain(ReferenceDomain):
    """Reference line domain [0, 1]."""

    @property
    def topological_dimension(self) -> int:
        return 1

    def quadrature(
        self, order: int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Gauss-Legendre quadrature mapped from [-1,1] to [0,1].

        Uses ceil((order+1)/2) points to integrate polynomials of given degree exactly.
        """
        n_points = (order + 2) // 2  # ceil((order+1)/2)
        xi_ref, w_ref = np.polynomial.legendre.leggauss(n_points)
        # Map from [-1, 1] to [0, 1]: x = (xi + 1) / 2, dx = d(xi)/2
        points = ((xi_ref + 1.0) / 2.0).reshape(-1, 1)
        weights = w_ref / 2.0
        return points, weights
