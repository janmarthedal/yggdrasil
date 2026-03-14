from collections.abc import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
from numpy.typing import NDArray

from .dof_map import DOFMap
from .forms import mass_form
from .mapping import compute_jacobian_det, compute_physical_gradients
from .mesh import Mesh


def assemble_bilinear_form(
    mesh: Mesh,
    bilinear_form: Callable[[NDArray, NDArray], NDArray],
    quadrature_order: int,
    dof_map: DOFMap | None = None,
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
    dof_map : DOFMap or None
        Mapping from mesh nodes to global DOF indices. Defaults to scalar
        (DOF index equals node index).

    Returns
    -------
    K : scipy.sparse.csr_matrix of shape (n_dofs, n_dofs)
    """
    if dof_map is None:
        dof_map = DOFMap(mesh)
    n_dofs = dof_map.n_dofs
    rows_list = []
    cols_list = []
    vals_list = []

    for group in mesh.iter_element_groups():
        element = group.element
        xi, weights = element.domain.quadrature(quadrature_order)
        N = element.shape_functions(xi)  # (num_quad, nodes_per_elem)

        for e in range(group.num_elements):
            elem_nodes = group.connectivity[e]
            phys_coords = mesh.nodes[elem_nodes]  # (nodes_per_elem, spatial_dim)

            grad_N, det_J = compute_physical_gradients(element, xi, phys_coords)
            # grad_N: (num_quad, nodes_per_elem, spatial_dim)
            # det_J: (num_quad,)

            jxw = weights * np.abs(det_J)  # (num_quad,)

            integrand = bilinear_form(N, grad_N)  # (num_quad, npe, npe)
            Ke = np.einsum("qij,q->ij", integrand, jxw)  # (npe, npe)

            elem_dofs = dof_map.element_dofs(elem_nodes)
            # Scatter into COO triplets
            local_rows = np.repeat(elem_dofs, len(elem_dofs))
            local_cols = np.tile(elem_dofs, len(elem_dofs))
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
    dof_map: DOFMap | None = None,
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
    dof_map : DOFMap or None
        Mapping from mesh nodes to global DOF indices. Defaults to scalar
        (DOF index equals node index).

    Returns
    -------
    b : ndarray of shape (n_dofs,)
    """
    if dof_map is None:
        dof_map = DOFMap(mesh)
    b = np.zeros(dof_map.n_dofs)

    for group in mesh.iter_element_groups():
        element = group.element
        xi, weights = element.domain.quadrature(quadrature_order)
        N = element.shape_functions(xi)

        for e in range(group.num_elements):
            elem_nodes = group.connectivity[e]
            phys_coords = mesh.nodes[elem_nodes]
            det_J = compute_jacobian_det(element, xi, phys_coords)
            jxw = weights * np.abs(det_J)

            if isinstance(f, (int, float)):
                be = f * np.einsum("qi,q->i", N, jxw)
            else:
                x_phys = N @ phys_coords  # (num_quad, spatial_dim)
                f_vals = f(x_phys)  # (num_quad,)
                be = np.einsum("qi,q->i", N, f_vals * jxw)

            elem_dofs = dof_map.element_dofs(elem_nodes)
            b[elem_dofs] += be

    return b


def assemble_neumann_bc(
    boundary_mesh: Mesh,
    g: float | Callable[[NDArray[np.float64]], NDArray[np.float64]],
    quadrature_order: int,
    n_dofs: int,
    dof_map: DOFMap | None = None,
) -> NDArray[np.float64]:
    """Assemble the Neumann BC contribution ∫_Γ g v dS into the global load vector.

    Add the result to ``b`` before calling ``condense_dirichlet_bc``.

    Parameters
    ----------
    boundary_mesh : Mesh
        Sub-mesh from ``extract_boundary`` or ``select_boundary_faces``.
        Must have ``point_data["original_node_index"]``.
    g : float or callable
        Prescribed flux (∂u/∂n = g on Γ_N).
        Callable signature: ``g(x) -> ndarray`` where x is (num_quad, spatial_dim).
    quadrature_order : int
        Quadrature order; use >= 2 * element.polynomial_degree for exactness.
    n_dofs : int
        Total DOFs in the global system (``mesh.num_nodes`` for scalar nodal FEM).
    dof_map : DOFMap or None
        When provided, ``n_dofs`` is ignored and replaced by ``dof_map.n_dofs``;
        boundary scatter uses ``dof_map.boundary_dofs``.

    Returns
    -------
    ndarray of shape (n_dofs,)
        Neumann contribution; scatter-add into the global load vector.
    """
    b_local = assemble_load_vector(boundary_mesh, g, quadrature_order)
    if dof_map is not None:
        n_dofs = dof_map.n_dofs
    b_global = np.zeros(n_dofs)
    original_idx = boundary_mesh.point_data["original_node_index"]
    if dof_map is not None:
        scatter_idx = dof_map.boundary_dofs(original_idx)
    else:
        scatter_idx = original_idx
    b_global[scatter_idx] += b_local
    return b_global


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
        free_dofs: NDArray[np.intp],
        bc_dofs: NDArray[np.intp],
        bc_vals: NDArray[np.float64],
        n_dofs: int,
    ):
        self.K = K
        self.b = b
        self._free_dofs = free_dofs
        self._bc_dofs = bc_dofs
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
            their correct global DOF indices.
        """
        u = np.empty(self._n_dofs)
        u[self._free_dofs] = u_f
        u[self._bc_dofs] = self._bc_vals
        return u


def condense_dirichlet_bc(
    K: sp.spmatrix,
    b: NDArray[np.float64],
    bc_dofs: NDArray[np.intp],
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
    bc_dofs : ndarray
        Global DOF indices (equals node indices for scalar problems) where
        Dirichlet BCs are applied.
    bc_val : float or ndarray of shape (len(bc_dofs),)
        Prescribed value(s) at the constrained DOFs.

    Returns
    -------
    CondensedSystem
        Object with ``K`` (shape ``(n_free, n_free)``) and ``b`` (shape
        ``(n_free,)``), plus a ``reconstruct(u_f)`` method that assembles the
        full solution vector from the free-DOF solution.
    """
    bc_dofs = np.asarray(bc_dofs, dtype=np.intp)
    bc_vals = np.array(np.broadcast_to(bc_val, len(bc_dofs)), dtype=np.float64)

    n_dofs = K.shape[0]
    free_dofs = np.setdiff1d(np.arange(n_dofs, dtype=np.intp), bc_dofs)

    K_csc = sp.csc_matrix(K)
    K_ff = K_csc[free_dofs][:, free_dofs].tocsr()
    b_f = b[free_dofs] - np.asarray(K_csc[free_dofs][:, bc_dofs] @ bc_vals).ravel()

    return CondensedSystem(K_ff, b_f, free_dofs, bc_dofs, bc_vals, n_dofs)


def project_dirichlet_bc(
    boundary_mesh: Mesh,
    g: float | Callable[[NDArray], NDArray],
    quadrature_order: int,
    dof_map: DOFMap | None = None,
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
        global DOF indices (equals node indices for scalar problems).
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
    dof_map : DOFMap or None
        When provided, ``dof_map.boundary_dofs`` is used to map boundary
        node indices to global DOF indices.

    Returns
    -------
    bc_dofs : ndarray of shape (num_boundary_nodes,)
        Global DOF indices (equals node indices for scalar problems) of the
        boundary nodes.
    bc_vals : ndarray of shape (num_boundary_nodes,)
        Projected DOF values minimising the L² error on the boundary.
        Suitable for passing directly to ``condense_dirichlet_bc``.
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
        mass_form,
        quadrature_order,
    )
    b = assemble_load_vector(boundary_mesh, g, quadrature_order)
    bc_vals = scipy.sparse.linalg.spsolve(M, b)
    original_idx = boundary_mesh.point_data["original_node_index"]
    if dof_map is not None:
        bc_dofs = dof_map.boundary_dofs(original_idx)
    else:
        bc_dofs = original_idx
    return bc_dofs, bc_vals
