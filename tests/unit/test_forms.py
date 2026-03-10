"""Tests for bilinear form functions in yggdrasil/forms.py."""

import numpy as np
import scipy.sparse as sp

from yggdrasil import ElementGroup, Mesh, assemble_bilinear_form, mass_form
from yggdrasil.elements import Tri3


def _make_unit_square_mesh():
    """Two-triangle mesh forming a unit square."""
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ])
    conn = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.intp)
    group = ElementGroup(element=Tri3(), connectivity=conn)
    return Mesh(nodes, [group])


class TestMassFormOutputShape:
    def test_output_shape(self):
        """mass_form returns shape (num_quad, npe, npe)."""
        # Tri3: 3 nodes per element, use 3 quadrature points
        num_quad, npe = 3, 3
        N = np.random.rand(num_quad, npe)
        grad_N = np.random.rand(num_quad, npe, 2)
        result = mass_form(N, grad_N)
        assert result.shape == (num_quad, npe, npe)

    def test_output_is_outer_product(self):
        """mass_form(N, grad_N)[q] == np.outer(N[q], N[q])."""
        num_quad, npe = 4, 3
        N = np.random.rand(num_quad, npe)
        grad_N = np.zeros((num_quad, npe, 2))
        result = mass_form(N, grad_N)
        for q in range(num_quad):
            np.testing.assert_allclose(result[q], np.outer(N[q], N[q]))


class TestMassMatrix:
    def test_symmetry(self):
        """Assembled mass matrix is symmetric."""
        mesh = _make_unit_square_mesh()
        M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)
        diff = M - M.T
        assert sp.linalg.norm(diff) < 1e-14

    def test_positive_definite(self):
        """All eigenvalues of the mass matrix are positive."""
        mesh = _make_unit_square_mesh()
        M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)
        eigvals = np.linalg.eigvalsh(M.toarray())
        assert np.all(eigvals > 0)

    def test_row_sum_equals_volume(self):
        """Sum of all mass matrix entries equals the domain volume (area = 1.0)."""
        mesh = _make_unit_square_mesh()
        M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)
        np.testing.assert_allclose(M.sum(), 1.0, atol=1e-14)

    def test_shape(self):
        """Mass matrix has shape (n_nodes, n_nodes)."""
        mesh = _make_unit_square_mesh()
        M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)
        assert M.shape == (4, 4)
