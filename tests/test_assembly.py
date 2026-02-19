"""Tests for bilinear form assembly."""

import numpy as np
import scipy.sparse as sp

from yggdrasil import ElementGroup, Mesh, assemble_bilinear_form, assemble_load_vector
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
