import numpy as np


def grad_grad_form(N, grad_N):
    """Bilinear form ∫ ∇u·∇v dx — the stiffness integrand for the Laplacian."""
    return np.einsum("qia,qja->qij", grad_N, grad_N)
