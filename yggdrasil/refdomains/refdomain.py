from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray


class ReferenceDomain(ABC):
    @property
    @abstractmethod
    def topological_dimension(self) -> int: ...

    @abstractmethod
    def quadrature(
        self, order: int
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Return (points, weights) that integrate polynomials of given degree exactly.

        points: shape (num_quad, topological_dimension)
        weights: shape (num_quad,)
        """
        ...
