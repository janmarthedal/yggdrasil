"""Tests for all reference elements: partition of unity, Kronecker delta,
gradient consistency."""

import numpy as np
import pytest

from yggdrasil.elements import Hex8, Line2, Line3, Quad4, Quad9, Tet4, Tri3, Tri6

ALL_ELEMENTS = [Line2, Line3, Tri3, Tri6, Quad4, Quad9, Tet4, Hex8]


def _random_interior_points(element, n=10, rng=None):
    """Generate random points inside the reference domain."""
    if rng is None:
        rng = np.random.default_rng(42)
    topo_dim = element.domain.topological_dimension

    if topo_dim == 1:
        return rng.uniform(0.05, 0.95, (n, 1))
    elif topo_dim == 2:
        if hasattr(element.domain, "__class__") and "Triangle" in type(
            element.domain
        ).__name__:
            # Generate points inside triangle using barycentric coordinates
            u = rng.uniform(0.05, 0.9, (n,))
            v = rng.uniform(0.05, 0.9, (n,))
            # Fold points outside triangle
            mask = u + v > 0.95
            u[mask] = 0.95 - u[mask]
            v[mask] = 0.95 - v[mask]
            return np.column_stack([u, v])
        else:
            return rng.uniform(0.05, 0.95, (n, 2))
    else:  # 3D
        if "Tetrahedron" in type(element.domain).__name__:
            pts = []
            while len(pts) < n:
                p = rng.uniform(0.05, 0.85, (3,))
                if p.sum() < 0.9:
                    pts.append(p)
            return np.array(pts)
        else:
            return rng.uniform(0.05, 0.95, (n, 3))


@pytest.mark.parametrize("ElementClass", ALL_ELEMENTS)
class TestPartitionOfUnity:
    def test_shape_functions_sum_to_one(self, ElementClass):
        """Shape functions must sum to 1 at any point in the domain."""
        elem = ElementClass()
        pts = _random_interior_points(elem)
        N = elem.shape_functions(pts)
        sums = np.sum(N, axis=1)
        np.testing.assert_allclose(sums, 1.0, atol=1e-14)

    def test_gradient_sums_to_zero(self, ElementClass):
        """Gradients of shape functions must sum to 0 (consequence of POU)."""
        elem = ElementClass()
        pts = _random_interior_points(elem)
        dN = elem.shape_function_gradients(pts)
        grad_sums = np.sum(dN, axis=1)
        np.testing.assert_allclose(grad_sums, 0.0, atol=1e-13)


@pytest.mark.parametrize("ElementClass", ALL_ELEMENTS)
class TestKroneckerDelta:
    def test_shape_at_nodes(self, ElementClass):
        """N_i(x_j) = delta_{ij}."""
        elem = ElementClass()
        nodes = elem.node_coords
        N = elem.shape_functions(nodes)
        np.testing.assert_allclose(N, np.eye(elem.num_nodes), atol=1e-14)


@pytest.mark.parametrize("ElementClass", ALL_ELEMENTS)
class TestGradientConsistency:
    def test_finite_difference(self, ElementClass):
        """Compare analytic gradients with centered finite differences."""
        elem = ElementClass()
        pts = _random_interior_points(elem, n=5)
        topo_dim = elem.domain.topological_dimension
        h = 1e-7

        analytic = elem.shape_function_gradients(pts)

        for d in range(topo_dim):
            pts_plus = pts.copy()
            pts_minus = pts.copy()
            pts_plus[:, d] += h
            pts_minus[:, d] -= h

            N_plus = elem.shape_functions(pts_plus)
            N_minus = elem.shape_functions(pts_minus)
            fd_grad = (N_plus - N_minus) / (2.0 * h)

            np.testing.assert_allclose(
                analytic[:, :, d],
                fd_grad,
                atol=1e-6,
                err_msg=f"Gradient mismatch for {ElementClass.__name__}, dim {d}",
            )
