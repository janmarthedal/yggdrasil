"""Tests for quadrature rules on all reference domains."""

import numpy as np
import pytest

from yggdrasil.refdomains import (
    HexahedronDomain,
    LineDomain,
    QuadrilateralDomain,
    TetrahedronDomain,
    TriangleDomain,
)


class TestLineDomainQuadrature:
    def test_dimension(self):
        d = LineDomain()
        assert d.topological_dimension == 1

    @pytest.mark.parametrize("order", [0, 1, 2, 3, 4, 5])
    def test_polynomial_exactness(self, order):
        """Integrate x^order over [0,1], exact result = 1/(order+1)."""
        d = LineDomain()
        pts, wts = d.quadrature(order)
        result = np.sum(wts * pts[:, 0] ** order)
        expected = 1.0 / (order + 1)
        np.testing.assert_allclose(result, expected, atol=1e-14)

    def test_weights_sum(self):
        """Weights should sum to domain measure = 1."""
        d = LineDomain()
        _, wts = d.quadrature(0)
        np.testing.assert_allclose(np.sum(wts), 1.0, atol=1e-14)


class TestTriangleDomainQuadrature:
    def test_dimension(self):
        d = TriangleDomain()
        assert d.topological_dimension == 2

    def test_constant_integration(self):
        """Integral of 1 over reference triangle = 0.5."""
        d = TriangleDomain()
        _, wts = d.quadrature(0)
        np.testing.assert_allclose(np.sum(wts), 0.5, atol=1e-14)

    @pytest.mark.parametrize("order", [1, 2, 3, 4, 5])
    def test_polynomial_exactness(self, order):
        """Integrate x^a * y^b with a+b <= order.

        For the reference triangle, integral of x^a * y^b =
        a! * b! / (a + b + 2)!
        """
        d = TriangleDomain()
        pts, wts = d.quadrature(order)

        for a in range(order + 1):
            for b in range(order + 1 - a):
                result = np.sum(wts * pts[:, 0] ** a * pts[:, 1] ** b)
                # Exact: a! * b! / (a + b + 2)!
                from math import factorial

                expected = factorial(a) * factorial(b) / factorial(a + b + 2)
                np.testing.assert_allclose(
                    result, expected, atol=1e-12, err_msg=f"Failed for x^{a}*y^{b}"
                )


class TestQuadrilateralDomainQuadrature:
    def test_dimension(self):
        d = QuadrilateralDomain()
        assert d.topological_dimension == 2

    def test_constant_integration(self):
        """Integral of 1 over [0,1]^2 = 1."""
        d = QuadrilateralDomain()
        _, wts = d.quadrature(0)
        np.testing.assert_allclose(np.sum(wts), 1.0, atol=1e-14)

    @pytest.mark.parametrize("order", [1, 2, 3, 4])
    def test_polynomial_exactness(self, order):
        """Integrate x^a * y^b over [0,1]^2 = 1/(a+1) * 1/(b+1)."""
        d = QuadrilateralDomain()
        pts, wts = d.quadrature(order)

        for a in range(order + 1):
            for b in range(order + 1):
                if a + b > order:
                    continue
                result = np.sum(wts * pts[:, 0] ** a * pts[:, 1] ** b)
                expected = 1.0 / (a + 1) / (b + 1)
                np.testing.assert_allclose(
                    result, expected, atol=1e-13, err_msg=f"Failed for x^{a}*y^{b}"
                )


class TestTetrahedronDomainQuadrature:
    def test_dimension(self):
        d = TetrahedronDomain()
        assert d.topological_dimension == 3

    def test_constant_integration(self):
        """Integral of 1 over reference tet = 1/6."""
        d = TetrahedronDomain()
        _, wts = d.quadrature(0)
        np.testing.assert_allclose(np.sum(wts), 1.0 / 6.0, atol=1e-14)

    @pytest.mark.parametrize("order", [1, 2, 3])
    def test_polynomial_exactness(self, order):
        """Integrate x^a * y^b * z^c over reference tet.

        Exact: a! * b! * c! / (a + b + c + 3)!
        """
        from math import factorial

        d = TetrahedronDomain()
        pts, wts = d.quadrature(order)

        for a in range(order + 1):
            for b in range(order + 1 - a):
                for c in range(order + 1 - a - b):
                    result = np.sum(
                        wts * pts[:, 0] ** a * pts[:, 1] ** b * pts[:, 2] ** c
                    )
                    expected = (
                        factorial(a) * factorial(b) * factorial(c)
                        / factorial(a + b + c + 3)
                    )
                    np.testing.assert_allclose(
                        result,
                        expected,
                        atol=1e-13,
                        err_msg=f"Failed for x^{a}*y^{b}*z^{c}",
                    )


class TestHexahedronDomainQuadrature:
    def test_dimension(self):
        d = HexahedronDomain()
        assert d.topological_dimension == 3

    def test_constant_integration(self):
        """Integral of 1 over [0,1]^3 = 1."""
        d = HexahedronDomain()
        _, wts = d.quadrature(0)
        np.testing.assert_allclose(np.sum(wts), 1.0, atol=1e-14)

    @pytest.mark.parametrize("order", [1, 2, 3])
    def test_polynomial_exactness(self, order):
        """Integrate x^a * y^b * z^c over [0,1]^3."""
        d = HexahedronDomain()
        pts, wts = d.quadrature(order)

        for a in range(order + 1):
            for b in range(order + 1):
                for c in range(order + 1):
                    if a + b + c > order:
                        continue
                    result = np.sum(
                        wts * pts[:, 0] ** a * pts[:, 1] ** b * pts[:, 2] ** c
                    )
                    expected = 1.0 / (a + 1) / (b + 1) / (c + 1)
                    np.testing.assert_allclose(
                        result,
                        expected,
                        atol=1e-13,
                        err_msg=f"Failed for x^{a}*y^{b}*z^{c}",
                    )
