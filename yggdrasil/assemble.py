from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp

from .mapping import compute_physical_gradients
from .mesh import Mesh


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


def assemble_load_vector(
    mesh: Mesh,
    f: float | Callable[[NDArray], NDArray],
    quadrature_order: int,
) -> NDArray[np.float64]:
    """Assemble the load vector for a source term f.

    Parameters
    ----------
    mesh : Mesh
    f : float or callable
        Source term. Either a constant scalar or a callable with signature
        f(x) -> array where x has shape (num_quad, spatial_dim) and the
        return value has shape (num_quad,).
    quadrature_order : int
        Polynomial order for the quadrature rule.

    Returns
    -------
    b : ndarray of shape (num_nodes,)
    """
    b = np.zeros(mesh.num_nodes)

    for group in mesh.iter_element_groups():
        element = group.element
        xi, weights = element.domain.quadrature(quadrature_order)
        N = element.shape_functions(xi)

        for e in range(group.num_elements):
            elem_nodes = group.connectivity[e]
            phys_coords = mesh.nodes[elem_nodes]
            _, det_J = compute_physical_gradients(element, xi, phys_coords)
            jxw = weights * np.abs(det_J)

            if callable(f):
                x_phys = N @ phys_coords  # (num_quad, spatial_dim)
                f_vals = f(x_phys)  # (num_quad,)
                be = np.einsum("qi,q->i", N, f_vals * jxw)
            else:
                be = f * np.einsum("qi,q->i", N, jxw)
            b[elem_nodes] += be

    return b


def apply_dirichlet_bc(
    K: sp.spmatrix,
    b: NDArray[np.float64],
    bc_nodes: NDArray[np.intp],
    bc_val: float = 0.0,
) -> tuple[sp.csr_matrix, NDArray[np.float64]]:
    """Apply Dirichlet boundary conditions by zeroing rows/cols and setting diagonal to 1.

    Parameters
    ----------
    K : scipy.sparse matrix
        The global stiffness matrix.
    b : ndarray of shape (num_nodes,)
        The global load vector. Modified in place.
    bc_nodes : ndarray
        Indices of nodes where the Dirichlet BC is applied.
    bc_val : float
        The prescribed value at the boundary nodes.

    Returns
    -------
    K : scipy.sparse.csr_matrix
        The modified stiffness matrix.
    b : ndarray
        The modified load vector.
    """
    K = K.tolil()
    for node in bc_nodes:
        K[node, :] = 0
        K[:, node] = 0
        K[node, node] = 1.0
        b[node] = bc_val
    return K.tocsr(), b
