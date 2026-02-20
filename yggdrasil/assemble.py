from collections.abc import Callable

import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sp
import scipy.sparse.linalg

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
    bc_val: float | NDArray[np.float64] = 0.0,
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
    bc_val : float or ndarray of shape (len(bc_nodes),)
        The prescribed value(s) at the boundary nodes. Either a single scalar
        applied to all constrained nodes, or a per-node array.

    Returns
    -------
    K : scipy.sparse.csr_matrix
        The modified stiffness matrix.
    b : ndarray
        The modified load vector.
    """
    bc_vals = np.broadcast_to(bc_val, len(bc_nodes))
    b -= np.asarray(K.tocsc()[:, bc_nodes] @ bc_vals).ravel()
    K = K.tolil()
    for i, node in enumerate(bc_nodes):
        K[node, :] = 0
        K[:, node] = 0
        K[node, node] = 1.0
        b[node] = bc_vals[i]
    return K.tocsr(), b


def project_dirichlet_bc(
    boundary_mesh: Mesh,
    g: float | Callable[[NDArray], NDArray],
    quadrature_order: int,
) -> tuple[NDArray[np.intp], NDArray[np.float64]]:
    """L² projection of g onto the boundary FE trace space.

    Finds the function g_h in the trace space that minimises
    ``‖g − g_h‖_{L²(Γ)}``, by solving the boundary mass system
    ``M_Γ α = b_Γ`` where ``(M_Γ)_ij = ∫_Γ N_i N_j dS`` and
    ``(b_Γ)_i = ∫_Γ g N_i dS``.

    Parameters
    ----------
    boundary_mesh : Mesh
        A boundary sub-mesh as returned by ``extract_boundary`` or
        ``select_boundary_faces``. Must have
        ``point_data["original_node_index"]`` mapping local nodes to
        global DOF indices.
    g : float or callable
        The prescribed boundary function. Either a constant scalar or a
        callable with signature ``g(x) -> array`` where ``x`` has shape
        ``(num_quad, spatial_dim)`` and the return value has shape
        ``(num_quad,)``.
    quadrature_order : int
        Polynomial order for the quadrature rule used on the boundary.

    Returns
    -------
    bc_nodes : ndarray of shape (num_boundary_nodes,)
        Global node indices of the boundary nodes.
    bc_vals : ndarray of shape (num_boundary_nodes,)
        Projected DOF values minimising the L² error on the boundary.
        Suitable for passing directly to ``apply_dirichlet_bc``.
    """
    M = assemble_bilinear_form(
        boundary_mesh,
        lambda N, grad_N: np.einsum("qi,qj->qij", N, N),
        quadrature_order,
    )
    b = assemble_load_vector(boundary_mesh, g, quadrature_order)
    bc_vals = scipy.sparse.linalg.spsolve(M, b)
    bc_nodes = boundary_mesh.point_data["original_node_index"]
    return bc_nodes, bc_vals
