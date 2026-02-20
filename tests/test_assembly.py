"""Tests for bilinear form assembly."""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg

from yggdrasil import ElementGroup, Mesh, apply_dirichlet_bc, assemble_bilinear_form, assemble_load_vector
from yggdrasil.elements import Tri3


def stiffness_form(N, grad_N):
    """grad(u) · grad(v) integrand."""
    return np.einsum("qia,qja->qij", grad_N, grad_N)


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

        K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)

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

        K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)

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
        K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)
        diff = K - K.T
        assert sp.linalg.norm(diff) < 1e-14

    def test_row_sum_zero(self):
        """Row sums of the stiffness matrix should be zero (constant is in the null space)."""
        mesh = self._make_two_tri_mesh()
        K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)
        row_sums = np.array(K.sum(axis=1)).flatten()
        np.testing.assert_allclose(row_sums, 0.0, atol=1e-14)

    def test_shape(self):
        """Matrix should be n_nodes x n_nodes."""
        mesh = self._make_two_tri_mesh()
        K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)
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


class TestApplyDirichletBC:
    """Tests for apply_dirichlet_bc."""

    def test_zeroes_bc_rows_and_cols(self):
        """BC rows and columns should be zeroed with 1 on diagonal."""
        K = sp.csr_matrix(np.array([
            [4.0, -1.0, -1.0],
            [-1.0, 4.0, -1.0],
            [-1.0, -1.0, 4.0],
        ]))
        b = np.array([1.0, 2.0, 3.0])
        bc_nodes = np.array([0], dtype=np.intp)

        K_mod, b_mod = apply_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

        K_expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 4.0, -1.0],
            [0.0, -1.0, 4.0],
        ])
        np.testing.assert_allclose(K_mod.toarray(), K_expected)
        assert b_mod[0] == 0.0

    def test_nonzero_bc_value(self):
        """Prescribed nonzero value should appear in the load vector."""
        K = sp.csr_matrix(np.eye(3))
        b = np.array([1.0, 2.0, 3.0])
        bc_nodes = np.array([1], dtype=np.intp)

        _, b_mod = apply_dirichlet_bc(K, b, bc_nodes, bc_val=5.0)

        assert b_mod[1] == 5.0
        # Other entries unchanged
        assert b_mod[0] == 1.0
        assert b_mod[2] == 3.0

    def test_multiple_bc_nodes(self):
        """Multiple BC nodes should all be constrained."""
        K = sp.csr_matrix(np.array([
            [4.0, -1.0, -1.0, 0.0],
            [-1.0, 4.0, 0.0, -1.0],
            [-1.0, 0.0, 4.0, -1.0],
            [0.0, -1.0, -1.0, 4.0],
        ]))
        b = np.array([1.0, 2.0, 3.0, 4.0])
        bc_nodes = np.array([0, 3], dtype=np.intp)

        K_mod, b_mod = apply_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

        # BC rows/cols zeroed, diagonal = 1
        for node in bc_nodes:
            row = K_mod[node, :].toarray().ravel()
            col = K_mod[:, node].toarray().ravel()
            assert row[node] == 1.0
            assert np.count_nonzero(row) == 1
            assert col[node] == 1.0
            assert np.count_nonzero(col) == 1
            assert b_mod[node] == 0.0

        # Interior entries unchanged
        assert K_mod[1, 2] == 0.0
        assert K_mod[2, 1] == 0.0

    def test_nonzero_bc_subtracts_column_contribution(self):
        """Column contribution of constrained DOF should be subtracted from unconstrained rows."""
        # 3-node system: constrain node 0 to value 2.0
        # K[1,0] = -1, so b[1] should decrease by (-1)*2 = -2, i.e. 2 - (-2) = 4... wait:
        # b_new[1] = b[1] - K[1,0] * bc_val = 2.0 - (-1.0)*2.0 = 4.0
        K = sp.csr_matrix(np.array([
            [4.0, -1.0, -1.0],
            [-1.0, 4.0, -1.0],
            [-1.0, -1.0, 4.0],
        ]))
        b = np.array([0.0, 2.0, 3.0])
        bc_nodes = np.array([0], dtype=np.intp)

        _, b_mod = apply_dirichlet_bc(K, b, bc_nodes, bc_val=2.0)

        # b[1] -= K[1,0] * 2.0 = (-1.0) * 2.0 = -2.0, so b[1] = 2.0 - (-2.0) = 4.0
        np.testing.assert_allclose(b_mod[1], 4.0, atol=1e-14)
        # b[2] -= K[2,0] * 2.0 = (-1.0) * 2.0 = -2.0, so b[2] = 3.0 - (-2.0) = 5.0
        np.testing.assert_allclose(b_mod[2], 5.0, atol=1e-14)
        # constrained node gets the prescribed value
        np.testing.assert_allclose(b_mod[0], 2.0, atol=1e-14)

    def test_poisson_linear_solution(self):
        """Solve -∇²u = 0 on [0,1]² with u=x on boundary; exact solution is u=x.

        With correct column elimination, linear elements should reproduce the
        exact affine solution to machine precision.
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

        K = assemble_bilinear_form(mesh, stiffness_form, quadrature_order=1)
        b = assemble_load_vector(mesh, 0.0, quadrature_order=1)

        # All 4 nodes are on the boundary; prescribe u = x-coordinate
        bc_nodes = np.array([0, 1, 2, 3], dtype=np.intp)
        bc_vals = nodes[bc_nodes, 0]  # x-coordinates: [0, 1, 1, 0]

        # Apply each node individually to allow per-node values
        for node, val in zip(bc_nodes, bc_vals):
            K, b = apply_dirichlet_bc(K, b, np.array([node], dtype=np.intp), bc_val=val)

        u = scipy.sparse.linalg.spsolve(K, b)
        np.testing.assert_allclose(u, nodes[:, 0], atol=1e-13)
