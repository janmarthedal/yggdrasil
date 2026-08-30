import numpy as np
from numpy.typing import NDArray


def grad_grad_form(N: NDArray[np.float64], grad_N: NDArray[np.float64]) -> NDArray[np.float64]:
    """Bilinear form ∫ ∇u·∇v dx — the stiffness integrand for the Laplacian."""
    return np.einsum("qia,qja->qij", grad_N, grad_N)


def mass_form(N: NDArray[np.float64], grad_N: NDArray[np.float64]) -> NDArray[np.float64]:
    """Bilinear form ∫ u v dx — the mass integrand."""
    return np.einsum("qi,qj->qij", N, N)
