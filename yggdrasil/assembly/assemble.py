from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp

from ..mapping import compute_physical_gradients
from ..mesh import Mesh


def assemble_bilinear_form(
    mesh: Mesh,
    bilinear_form: Callable[[NDArray, NDArray], NDArray],
    quadrature_order: int,
) -> sp.csr_matrix:
    """Assemble a global sparse matrix from a bilinear form.

    Parameters
    ----------
    mesh : Mesh
    bilinear_form : callable
        Signature: (N, grad_N) -> integrand
            N:      (num_quad, nodes_per_elem)
            grad_N: (num_quad, nodes_per_elem, spatial_dim)
            returns: (num_quad, nodes_per_elem, nodes_per_elem)
    quadrature_order : int
        Polynomial order for the quadrature rule.

    Returns
    -------
    K : scipy.sparse.csr_matrix of shape (num_nodes, num_nodes)
    """
    n_dofs = mesh.num_nodes
    rows_list = []
    cols_list = []
    vals_list = []

    for group in mesh.iter_element_groups():
        element = group.element
        xi, weights = element.domain.quadrature(quadrature_order)
        N = element.shape_functions(xi)  # (num_quad, nodes_per_elem)
        nodes_per_elem = element.num_nodes

        for e in range(group.num_elements):
            elem_nodes = group.connectivity[e]
            phys_coords = mesh.nodes[elem_nodes]  # (nodes_per_elem, spatial_dim)

            grad_N, det_J = compute_physical_gradients(element, xi, phys_coords)
            # grad_N: (num_quad, nodes_per_elem, spatial_dim)
            # det_J: (num_quad,)

            jxw = weights * np.abs(det_J)  # (num_quad,)

            integrand = bilinear_form(N, grad_N)  # (num_quad, npe, npe)
            Ke = np.einsum("qij,q->ij", integrand, jxw)  # (npe, npe)

            # Scatter into COO triplets
            local_rows = np.repeat(elem_nodes, nodes_per_elem)
            local_cols = np.tile(elem_nodes, nodes_per_elem)
            rows_list.append(local_rows)
            cols_list.append(local_cols)
            vals_list.append(Ke.ravel())

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    vals = np.concatenate(vals_list)

    K = sp.coo_matrix((vals, (rows, cols)), shape=(n_dofs, n_dofs))
    return K.tocsr()
