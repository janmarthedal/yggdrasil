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

            if isinstance(f, (int, float)):
                be = f * np.einsum("qi,q->i", N, jxw)
            else:
                x_phys = N @ phys_coords  # (num_quad, spatial_dim)
                f_vals = f(x_phys)  # (num_quad,)
                be = np.einsum("qi,q->i", N, f_vals * jxw)
            b[elem_nodes] += be

    return b


class CondensedSystem:
    """A Dirichlet-condensed linear system ready to solve.

    Attributes
    ----------
    K : sp.csr_matrix of shape (n_free, n_free)
        Reduced stiffness matrix for the free DOFs.
    b : ndarray of shape (n_free,)
        Reduced load vector with BC column contributions subtracted.
    """

    def __init__(
        self,
        K: sp.csr_matrix,
        b: NDArray[np.float64],
        free_nodes: NDArray[np.intp],
        bc_nodes: NDArray[np.intp],
        bc_vals: NDArray[np.float64],
        n_dofs: int,
    ):
        self.K = K
        self.b = b
        self._free_nodes = free_nodes
        self._bc_nodes = bc_nodes
        self._bc_vals = bc_vals
        self._n_dofs = n_dofs

    def reconstruct(self, u_f: NDArray[np.float64]) -> NDArray[np.float64]:
        """Assemble the full solution from the free-DOF solution.

        Parameters
        ----------
        u_f : ndarray of shape (n_free,)
            Solution of the reduced system (e.g. from ``spsolve(system.K, system.b)``).

        Returns
        -------
        u : ndarray of shape (n_dofs,)
            Full solution vector with free and constrained values placed at
            their correct global indices.
        """
        u = np.empty(self._n_dofs)
        u[self._free_nodes] = u_f
        u[self._bc_nodes] = self._bc_vals
        return u


def condense_dirichlet_bc(
    K: sp.spmatrix,
    b: NDArray[np.float64],
    bc_nodes: NDArray[np.intp],
    bc_val: float | NDArray[np.float64] = 0.0,
) -> CondensedSystem:
    """Condense Dirichlet BCs by eliminating constrained DOFs from the system.

    Partitions K and b into free (f) and constrained (c) blocks and returns
    the reduced system ``K_ff u_f = b_f - K_fc g``, where ``g`` are the
    prescribed boundary values.

    Parameters
    ----------
    K : scipy.sparse matrix of shape (n_dofs, n_dofs)
        The global stiffness matrix.
    b : ndarray of shape (n_dofs,)
        The global load vector.
    bc_nodes : ndarray
        Indices of nodes where Dirichlet BCs are applied.
    bc_val : float or ndarray of shape (len(bc_nodes),)
        Prescribed value(s) at the boundary nodes.

    Returns
    -------
    CondensedSystem
        Object with ``K`` (shape ``(n_free, n_free)``) and ``b`` (shape
        ``(n_free,)``), plus a ``reconstruct(u_f)`` method that assembles the
        full solution vector from the free-DOF solution.
    """
    bc_nodes = np.asarray(bc_nodes, dtype=np.intp)
    bc_vals = np.array(np.broadcast_to(bc_val, len(bc_nodes)), dtype=np.float64)

    n_dofs = K.shape[0]
    free_nodes = np.setdiff1d(np.arange(n_dofs, dtype=np.intp), bc_nodes)

    K_csc = sp.csc_matrix(K)
    K_ff = K_csc[free_nodes][:, free_nodes].tocsr()
    b_f = b[free_nodes] - np.asarray(K_csc[free_nodes][:, bc_nodes] @ bc_vals).ravel()

    return CondensedSystem(K_ff, b_f, free_nodes, bc_nodes, bc_vals, n_dofs)


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
        Must be at least ``2 * element.polynomial_degree`` for every element
        type in ``boundary_mesh`` so that the boundary mass matrix is
        integrated exactly. For example, use ``quadrature_order >= 2`` for
        linear boundary elements (Line2, Tri3) and ``>= 4`` for quadratic
        boundary elements (Line3, Tri6).

    Returns
    -------
    bc_nodes : ndarray of shape (num_boundary_nodes,)
        Global node indices of the boundary nodes.
    bc_vals : ndarray of shape (num_boundary_nodes,)
        Projected DOF values minimising the L² error on the boundary.
        Suitable for passing directly to ``apply_dirichlet_bc``.
    """
    for group in boundary_mesh.iter_element_groups():
        min_order = 2 * group.element.polynomial_degree
        assert quadrature_order >= min_order, (
            f"{type(group.element).__name__} elements have polynomial degree "
            f"{group.element.polynomial_degree}, so the boundary mass matrix "
            f"requires quadrature_order >= {min_order} "
            f"(got {quadrature_order})"
        )
    M = assemble_bilinear_form(
        boundary_mesh,
        lambda N, grad_N: np.einsum("qi,qj->qij", N, N),
        quadrature_order,
    )
    b = assemble_load_vector(boundary_mesh, g, quadrature_order)
    bc_vals = scipy.sparse.linalg.spsolve(M, b)
    bc_nodes = boundary_mesh.point_data["original_node_index"]
    return bc_nodes, bc_vals
