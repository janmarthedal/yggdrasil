import numpy as np
from numpy.typing import NDArray

from .refdomain import ReferenceDomain


class PointDomain(ReferenceDomain):
    """Reference point domain (0-dimensional)."""

    @property
    def topological_dimension(self) -> int:
        return 0

    def quadrature(
        self, order: int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Single-point quadrature with weight 1."""
        points = np.empty((1, 0), dtype=np.float64)
        weights = np.array([1.0])
        return points, weights
