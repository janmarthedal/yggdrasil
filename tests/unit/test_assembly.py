"""Tests for bilinear form assembly."""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg

from yggdrasil import ElementGroup, Mesh, assemble_bilinear_form, assemble_load_vector, condense_dirichlet_bc, grad_grad_form, project_dirichlet_bc
from yggdrasil.elements import Tri3


class TestSingleTriangle:
    def test_stiffness_right_triangle(self):
        """Stiffness matrix for the reference right triangle (0,0)-(1,0)-(0,1).

        For a right triangle with legs of length 1 and area = 0.5,
        the exact stiffness matrix is:
        K = 0.5 * [[2, -1, -1],
                    [-1, 1,  0],
                    [-1, 0,  1]]
        """
        nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        conn = np.array([[0, 1, 2]], dtype=np.intp)
        group = ElementGroup(element=Tri3(), connectivity=conn)
        mesh = Mesh(nodes, [group])

        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)

        K_expected = 0.5 * np.array([
            [2.0, -1.0, -1.0],
            [-1.0, 1.0, 0.0],
            [-1.0, 0.0, 1.0],
        ])
        np.testing.assert_allclose(K.toarray(), K_expected, atol=1e-14)

    def test_stiffness_scaled_triangle(self):
        """Stiffness matrix for a right triangle with legs of length 2.

        For this triangle, gradients scale by 1/2 and area scales by 4,
        so K scales by 4 * (1/4) = 1 relative to the unit triangle.
        The exact stiffness matrix is the same as the unit triangle.
        """
        nodes = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
        conn = np.array([[0, 1, 2]], dtype=np.intp)
        group = ElementGroup(element=Tri3(), connectivity=conn)
        mesh = Mesh(nodes, [group])

        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)

        # Same as unit triangle: grad scales as 1/s, area as s^2, net is independent of scale
        K_expected = 0.5 * np.array([
            [2.0, -1.0, -1.0],
            [-1.0, 1.0, 0.0],
            [-1.0, 0.0, 1.0],
        ])
        np.testing.assert_allclose(K.toarray(), K_expected, atol=1e-14)


class TestMultiElementMesh:
    def _make_two_tri_mesh(self):
        """Two-triangle mesh forming a unit square.

        3---2
        |T1/|
        | / |
        |/T0|
        0---1
        """
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ])
        conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
        group = ElementGroup(element=Tri3(), connectivity=conn)
        return Mesh(nodes, [group])

    def test_symmetry(self):
        """Assembled stiffness matrix should be symmetric."""
        mesh = self._make_two_tri_mesh()
        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
        diff = K - K.T
        assert sp.linalg.norm(diff) < 1e-14

    def test_row_sum_zero(self):
        """Row sums of the stiffness matrix should be zero (constant is in the null space)."""
        mesh = self._make_two_tri_mesh()
        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
        row_sums = np.array(K.sum(axis=1)).flatten()
        np.testing.assert_allclose(row_sums, 0.0, atol=1e-14)

    def test_shape(self):
        """Matrix should be n_nodes x n_nodes."""
        mesh = self._make_two_tri_mesh()
        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
        assert K.shape == (4, 4)


class TestLoadVector:
    """Tests for assemble_load_vector with scalar and callable source terms."""

    def _make_two_tri_mesh(self):
        """Two-triangle mesh forming a unit square.

        3---2
        |T1/|
        | / |
        |/T0|
        0---1
        """
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ])
        conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
        group = ElementGroup(element=Tri3(), connectivity=conn)
        return Mesh(nodes, [group])

    def test_scalar_float(self):
        """Scalar f_val should produce correct load vector."""
        mesh = self._make_two_tri_mesh()
        b = assemble_load_vector(mesh, 1.0, quadrature_order=1)
        # For a unit square with 2 triangles, each of area 0.5,
        # integral of 1 over the domain = 1.0, so sum(b) = 1.0
        np.testing.assert_allclose(b.sum(), 1.0, atol=1e-14)
        assert b.shape == (4,)

    def test_constant_callable_matches_scalar(self):
        """A callable returning ones should match scalar f=1."""
        mesh = self._make_two_tri_mesh()
        b_scalar = assemble_load_vector(mesh, 1.0, quadrature_order=1)
        b_callable = assemble_load_vector(
            mesh, lambda x: np.ones(x.shape[0]), quadrature_order=1
        )
        np.testing.assert_allclose(b_callable, b_scalar, atol=1e-14)

    def test_spatially_varying(self):
        """f(x) = x_0 on the unit square should integrate to 0.5."""
        mesh = self._make_two_tri_mesh()
        b = assemble_load_vector(mesh, lambda x: x[:, 0], quadrature_order=2)
        # integral of x over [0,1]^2 = 0.5
        np.testing.assert_allclose(b.sum(), 0.5, atol=1e-14)


class TestCondenseDirichletBC:
    """Tests for condense_dirichlet_bc and CondensedSystem.reconstruct."""

    def test_reduced_system_size(self):
        """Condensed K has shape (n_free, n_free) and b has shape (n_free,)."""
        K = sp.csr_matrix(np.array([
            [4.0, -1.0, -1.0],
            [-1.0, 4.0, -1.0],
            [-1.0, -1.0, 4.0],
        ]))
        b = np.array([1.0, 2.0, 3.0])
        bc_nodes = np.array([0], dtype=np.intp)

        system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

        assert system.K.shape == (2, 2)
        assert system.b.shape == (2,)

    def test_reconstruct_places_bc_values(self):
        """Reconstructed solution has the prescribed value at the BC node."""
        K = sp.csr_matrix(np.eye(3))
        b = np.zeros(3)
        bc_nodes = np.array([1], dtype=np.intp)

        system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=5.0)
        u_f = scipy.sparse.linalg.spsolve(system.K, system.b)
        u = system.reconstruct(u_f)

        assert u[1] == 5.0

    def test_multiple_bc_nodes_size(self):
        """Multiple BC nodes reduce the system size correctly."""
        K = sp.csr_matrix(np.array([
            [4.0, -1.0, -1.0, 0.0],
            [-1.0, 4.0, 0.0, -1.0],
            [-1.0, 0.0, 4.0, -1.0],
            [0.0, -1.0, -1.0, 4.0],
        ]))
        b = np.array([1.0, 2.0, 3.0, 4.0])
        bc_nodes = np.array([0, 3], dtype=np.intp)

        system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

        assert system.K.shape == (2, 2)
        assert system.b.shape == (2,)

    def test_rhs_column_contribution_subtracted(self):
        """b_f subtracts the K_fc * g column contribution from the free-DOF RHS."""
        K = sp.csr_matrix(np.array([
            [4.0, -1.0, -1.0],
            [-1.0, 4.0, -1.0],
            [-1.0, -1.0, 4.0],
        ]))
        b = np.array([0.0, 2.0, 3.0])
        bc_nodes = np.array([0], dtype=np.intp)

        system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=2.0)

        # b_f[0] = b[1] - K[1,0]*2 = 2.0 - (-1.0)*2.0 = 4.0
        # b_f[1] = b[2] - K[2,0]*2 = 3.0 - (-1.0)*2.0 = 5.0
        np.testing.assert_allclose(system.b[0], 4.0, atol=1e-14)
        np.testing.assert_allclose(system.b[1], 5.0, atol=1e-14)

    def test_poisson_condensed_solve(self):
        """Condensed solve of -∇²u=0 with node 2 free gives the correct discrete harmonic.

        Mesh: two triangles on [0,1]², nodes 0=(0,0), 1=(1,0), 2=(1,1), 3=(0,1).
        BCs: u[0]=0, u[1]=1, u[3]=0. Node 2 is free.

        K[2,2]=1, K[2,1]=-0.5, K[2,3]=-0.5 (off-diagonals from the two triangles).
        b_f = 0 - K[2,1]*1 - K[2,3]*0 = 0.5, so u[2] = 0.5.
        """
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ])
        conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
        group = ElementGroup(element=Tri3(), connectivity=conn)
        mesh = Mesh(nodes, [group])

        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
        b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

        bc_nodes = np.array([0, 1, 3], dtype=np.intp)
        bc_vals = np.array([0.0, 1.0, 0.0])
        system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=bc_vals)
        u_f = scipy.sparse.linalg.spsolve(system.K, system.b)
        u = system.reconstruct(u_f)

        np.testing.assert_allclose(u[2], 0.5, atol=1e-13)
        np.testing.assert_allclose(u[[0, 1, 3]], [0.0, 1.0, 0.0], atol=1e-13)


class TestProjectDirichletBC:
    """Tests for project_dirichlet_bc."""

    def _make_unit_square_mesh(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
        ])
        conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
        group = ElementGroup(element=Tri3(), connectivity=conn)
        return Mesh(nodes, [group])

    def test_constant_g_matches_nodal(self):
        """For constant g, L² projection should recover the constant at every node.

        Note: quadrature_order=2 is required because the boundary mass matrix
        for Line2 elements integrates products of linear shape functions
        (quadratic integrands), which 1-point midpoint rule cannot handle.
        """
        from yggdrasil import extract_boundary
        mesh = self._make_unit_square_mesh()
        bnd = extract_boundary(mesh)
        bc_nodes, bc_vals = project_dirichlet_bc(bnd, g=3.0, quadrature_order=2)
        np.testing.assert_allclose(bc_vals, 3.0, atol=1e-13)

    def test_linear_g_matches_nodal(self):
        """For g=x (affine), L² projection must agree with nodal sampling on linear elements."""
        from yggdrasil import extract_boundary
        mesh = self._make_unit_square_mesh()
        bnd = extract_boundary(mesh)
        bc_nodes, bc_vals = project_dirichlet_bc(bnd, g=lambda x: x[:, 0], quadrature_order=2)
        expected = mesh.nodes[bc_nodes, 0]
        np.testing.assert_allclose(bc_vals, expected, atol=1e-13)

    def test_returns_global_node_indices(self):
        """bc_nodes must be global indices into the original mesh."""
        from yggdrasil import extract_boundary
        mesh = self._make_unit_square_mesh()
        bnd = extract_boundary(mesh)
        bc_nodes, _ = project_dirichlet_bc(bnd, g=0.0, quadrature_order=2)
        assert set(bc_nodes).issubset(set(range(mesh.num_nodes)))

    def test_insufficient_quadrature_order_raises(self):
        """quadrature_order < 2*polynomial_degree must raise AssertionError."""
        import pytest
        from yggdrasil import extract_boundary
        mesh = self._make_unit_square_mesh()
        bnd = extract_boundary(mesh)
        with pytest.raises(AssertionError, match="quadrature_order >= 2"):
            project_dirichlet_bc(bnd, g=3.0, quadrature_order=1)

    def test_projection_solves_poisson(self):
        """project_dirichlet_bc + condense_dirichlet_bc should solve -∇²u=0, u=x correctly."""
        from yggdrasil import extract_boundary
        mesh = self._make_unit_square_mesh()

        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
        b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

        bnd = extract_boundary(mesh)
        bc_nodes, bc_vals = project_dirichlet_bc(bnd, g=lambda x: x[:, 0], quadrature_order=2)
        system = condense_dirichlet_bc(K, b, bc_nodes, bc_val=bc_vals)
        u = system.reconstruct(scipy.sparse.linalg.spsolve(system.K, system.b))

        np.testing.assert_allclose(u, mesh.nodes[:, 0], atol=1e-13)
