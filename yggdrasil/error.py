from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray

from .mapping import compute_physical_gradients
from .mesh import Mesh


def l2_error(mesh: Mesh, uh: NDArray[np.float64], u_exact: Callable[[NDArray], NDArray], quadrature_order: int) -> float:
    """Return ||u_h - u_exact||_{L²} = sqrt(∫_Ω (u_h - u_exact)² dx)."""
    error_sq = 0.0
    for group in mesh.iter_element_groups():
        element = group.element
        xi, weights = element.domain.quadrature(quadrature_order)
        N = element.shape_functions(xi)  # (num_quad, npe)
        for e in range(group.num_elements):
            elem_nodes = group.connectivity[e]
            phys_coords = mesh.nodes[elem_nodes]  # (npe, spatial_dim)
            x_phys = N @ phys_coords  # (num_quad, spatial_dim)
            _, det_J = compute_physical_gradients(element, xi, phys_coords)
            jxw = weights * np.abs(det_J)  # (num_quad,)
            uh_q = N @ uh[elem_nodes]  # (num_quad,)
            diff = uh_q - u_exact(x_phys)  # (num_quad,)
            error_sq += np.dot(diff**2, jxw)
    return float(np.sqrt(error_sq))
