import numpy as np
from numpy.typing import NDArray

from .elements.element import ReferenceElement


def compute_jacobian(
    element: ReferenceElement,
    xi: NDArray[np.float64],
    physical_coords: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute the Jacobian of the mapping from reference to physical space.

    Parameters
    ----------
    element : ReferenceElement
    xi : (num_points, topo_dim) — evaluation points in reference space
    physical_coords : (num_nodes, spatial_dim) — physical node coordinates

    Returns
    -------
    J : (num_points, spatial_dim, topo_dim)
        J[q, i, j] = dx_i / d(xi_j)
    """
    # grad_N: (num_points, num_nodes, topo_dim)
    grad_N = element.shape_function_gradients(xi)
    # J = physical_coords^T @ grad_N  =>  einsum
    # physical_coords: (num_nodes, spatial_dim)
    # grad_N: (num_points, num_nodes, topo_dim)
    # result: (num_points, spatial_dim, topo_dim)
    return np.einsum("ni,qnj->qij", physical_coords, grad_N)


def compute_jacobian_det(
    element: ReferenceElement,
    xi: NDArray[np.float64],
    physical_coords: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute the determinant of the Jacobian (or its metric generalization).

    Parameters
    ----------
    element : ReferenceElement
    xi : (num_points, topo_dim) — evaluation points in reference space
    physical_coords : (num_nodes, spatial_dim) — physical node coordinates

    Returns
    -------
    det_J : (num_points,)
        Determinant of the Jacobian when topo_dim == spatial_dim,
        or sqrt(det(J^T J)) when topo_dim < spatial_dim.
    """
    J = compute_jacobian(element, xi, physical_coords)
    topo_dim = element.domain.topological_dimension
    spatial_dim = physical_coords.shape[1]

    if topo_dim == spatial_dim:
        return np.linalg.det(J)
    else:
        g = np.einsum("qij,qik->qjk", J, J)
        return np.sqrt(np.abs(np.linalg.det(g)))


def compute_physical_gradients(
    element: ReferenceElement,
    xi: NDArray[np.float64],
    physical_coords: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Map reference-space gradients to physical-space gradients.

    Parameters
    ----------
    element : ReferenceElement
    xi : (num_points, topo_dim)
    physical_coords : (num_nodes, spatial_dim)

    Returns
    -------
    grad_phys : (num_points, num_nodes, spatial_dim)
        Shape function gradients in physical coordinates.
    det_J : (num_points,)
        Determinant of the Jacobian (or its appropriate generalization).
    """
    J = compute_jacobian(element, xi, physical_coords)
    grad_N = element.shape_function_gradients(xi)

    topo_dim = element.domain.topological_dimension
    spatial_dim = physical_coords.shape[1]

    if topo_dim == spatial_dim:
        # Square Jacobian: invert directly
        J_inv = np.linalg.inv(J)  # (num_points, topo_dim, topo_dim)
        det_J = np.linalg.det(J)  # (num_points,)
        # grad_phys = grad_N @ J_inv
        grad_phys = np.einsum("qnj,qjk->qnk", grad_N, J_inv)
    else:
        # Non-square Jacobian (e.g., 2D surface in 3D): use pseudo-inverse
        # metric tensor g = J^T J
        g = np.einsum("qij,qik->qjk", J, J)
        det_g = np.linalg.det(g)
        det_J = np.sqrt(np.abs(det_g))
        g_inv = np.linalg.inv(g)
        # Tangent-plane gradient: grad_N @ g^{-1}
        grad_phys = np.einsum("qnj,qjk->qnk", grad_N, g_inv)

    return grad_phys, det_J
