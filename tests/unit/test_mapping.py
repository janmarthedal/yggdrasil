"""Tests for physical-space mapping utilities."""

import numpy as np

from yggdrasil.elements import Hex8, Line2, Quad4, Tri3
from yggdrasil.mapping import (
    compute_jacobian,
    compute_mapped_quadrature,
    compute_physical_gradients,
)


class TestJacobian:
    def test_line_scaling(self):
        """A line element from x=2 to x=5 should have J = 3."""
        elem = Line2()
        xi = np.array([[0.25], [0.5], [0.75]])
        phys = np.array([[2.0], [5.0]])
        J = compute_jacobian(elem, xi, phys)
        np.testing.assert_allclose(J[:, 0, 0], 3.0, atol=1e-14)

    def test_tri_identity(self):
        """Triangle with vertices at reference positions: J = I."""
        elem = Tri3()
        xi = np.array([[0.25, 0.25]])
        phys = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        J = compute_jacobian(elem, xi, phys)
        np.testing.assert_allclose(J[0], np.eye(2), atol=1e-14)

    def test_tri_scaling(self):
        """Triangle scaled by 2: J = 2*I."""
        elem = Tri3()
        xi = np.array([[0.25, 0.25]])
        phys = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]])
        J = compute_jacobian(elem, xi, phys)
        np.testing.assert_allclose(J[0], 2.0 * np.eye(2), atol=1e-14)

    def test_quad_identity(self):
        """Quad with vertices at reference positions: J = I."""
        elem = Quad4()
        xi = np.array([[0.5, 0.5]])
        phys = elem.node_coords.copy()
        J = compute_jacobian(elem, xi, phys)
        np.testing.assert_allclose(J[0], np.eye(2), atol=1e-14)

    def test_hex_scaling(self):
        """Hex scaled by 3 in all directions."""
        elem = Hex8()
        xi = np.array([[0.5, 0.5, 0.5]])
        phys = elem.node_coords * 3.0
        J = compute_jacobian(elem, xi, phys)
        np.testing.assert_allclose(J[0], 3.0 * np.eye(3), atol=1e-14)


class TestPhysicalGradients:
    def test_line_gradient(self):
        """Physical gradient of linear element on [0, L] should be 1/L * reference."""
        elem = Line2()
        L = 4.0
        xi = np.array([[0.3]])
        phys = np.array([[0.0], [L]])
        grad_phys, det_J = compute_physical_gradients(elem, xi, phys)
        np.testing.assert_allclose(det_J, [L], atol=1e-14)
        # dN/dx_phys = dN/dxi * 1/L
        np.testing.assert_allclose(grad_phys[0, 0, 0], -1.0 / L, atol=1e-14)
        np.testing.assert_allclose(grad_phys[0, 1, 0], 1.0 / L, atol=1e-14)

    def test_tri_nonsymmetric_jacobian(self):
        """Physical gradients for a triangle with non-symmetric Jacobian.

        Triangle (0,0)-(1,0)-(1,1) has J = [[1,1],[0,1]], which is
        non-symmetric. This catches bugs where J_inv is accidentally
        transposed in the einsum.
        """
        elem = Tri3()
        xi = np.array([[0.25, 0.25]])
        phys = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]])
        grad_phys, det_J = compute_physical_gradients(elem, xi, phys)
        # J = [[1, 1], [0, 1]], det_J = 1
        np.testing.assert_allclose(det_J, [1.0], atol=1e-14)
        # J_inv = [[1, -1], [0, 1]]
        # grad_phys = ref_grad @ J_inv
        # ref_grad: N0=(-1,-1), N1=(1,0), N2=(0,1)
        # N0: (-1,-1) @ [[1,-1],[0,1]] = (-1, 0)
        # N1: (1, 0) @ [[1,-1],[0,1]] = (1, -1)
        # N2: (0, 1) @ [[1,-1],[0,1]] = (0, 1)
        expected = np.array([[[-1.0, 0.0], [1.0, -1.0], [0.0, 1.0]]])
        np.testing.assert_allclose(grad_phys, expected, atol=1e-14)

    def test_tri_scaled_gradient(self):
        """For scaled triangle, physical gradients = reference gradients / scale."""
        elem = Tri3()
        scale = 3.0
        xi = np.array([[0.25, 0.25]])
        phys = elem.node_coords * scale
        grad_phys, det_J = compute_physical_gradients(elem, xi, phys)
        # det_J = scale^2
        np.testing.assert_allclose(det_J, [scale**2], atol=1e-13)
        # Physical gradients = reference / scale
        ref_grad = elem.shape_function_gradients(xi)
        np.testing.assert_allclose(grad_phys, ref_grad / scale, atol=1e-13)


class TestMappedQuadrature:
    def test_line_area(self):
        """Integrate 1 over a physical line [0, 5]: should get 5."""
        elem = Line2()
        phys = np.array([[0.0], [5.0]])
        _, wts, det_J = compute_mapped_quadrature(elem, 1, phys)
        area = np.sum(wts * det_J)
        np.testing.assert_allclose(area, 5.0, atol=1e-14)

    def test_tri_area(self):
        """Area of triangle with vertices (0,0), (4,0), (0,3) = 6."""
        elem = Tri3()
        phys = np.array([[0.0, 0.0], [4.0, 0.0], [0.0, 3.0]])
        _, wts, det_J = compute_mapped_quadrature(elem, 1, phys)
        area = np.sum(wts * det_J)
        np.testing.assert_allclose(area, 6.0, atol=1e-13)

    def test_quad_area(self):
        """Area of [0,2] x [0,3] rectangle = 6."""
        elem = Quad4()
        phys = np.array([[0.0, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]])
        _, wts, det_J = compute_mapped_quadrature(elem, 1, phys)
        area = np.sum(wts * det_J)
        np.testing.assert_allclose(area, 6.0, atol=1e-13)

    def test_hex_volume(self):
        """Volume of [0,2]^3 = 8."""
        elem = Hex8()
        phys = elem.node_coords * 2.0
        _, wts, det_J = compute_mapped_quadrature(elem, 1, phys)
        vol = np.sum(wts * det_J)
        np.testing.assert_allclose(vol, 8.0, atol=1e-13)

    def test_physical_points_location(self):
        """Mapped quadrature points should lie inside the physical element."""
        elem = Tri3()
        phys = np.array([[1.0, 1.0], [3.0, 1.0], [1.0, 4.0]])
        phys_pts, _, _ = compute_mapped_quadrature(elem, 2, phys)
        # All x >= 1, y >= 1, and within the triangle
        assert np.all(phys_pts[:, 0] >= 1.0 - 1e-10)
        assert np.all(phys_pts[:, 1] >= 1.0 - 1e-10)
