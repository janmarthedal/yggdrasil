import numpy as np
from numpy.typing import NDArray

from .refdomain import ReferenceDomain
from .line import LineDomain


class QuadrilateralDomain(ReferenceDomain):
    """Reference quadrilateral domain [0,1]^2."""

    @property
    def topological_dimension(self) -> int:
        return 2

    def quadrature(
        self, order: int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Tensor product of 1D Gauss-Legendre quadrature on [0,1]^2."""
        line = LineDomain()
        pts_1d, wts_1d = line.quadrature(order)
        pts_1d = pts_1d.ravel()

        # Tensor product
        px, py = np.meshgrid(pts_1d, pts_1d, indexing="ij")
        wx, wy = np.meshgrid(wts_1d, wts_1d, indexing="ij")

        points = np.column_stack([px.ravel(), py.ravel()])
        weights = (wx * wy).ravel()
        return points, weights
