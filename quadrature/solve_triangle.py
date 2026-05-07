#!/usr/bin/env python3
"""
Derive quadrature rules for the reference triangle from first principles.

Solves the moment equations for the symmetry-group parametrization
described in quadrature/CLAUDE.md. Points are Cartesian (x, y) on the
reference triangle with vertices (0,0), (1,0), (0,1); weights sum to 0.5.

Symmetry groups
---------------
n0 — centroid (1/3, 1/3): 1 point, 1 free parameter (weight only)
n1 — (a, a), (b, a), (a, b) where b = 1-2a: 3 points, 2 free parameters (a, weight)

Moment equations
----------------
For each independent monomial x^i y^j (i <= j, i+j <= order), the sum of
w * x^i * y^j over all quadrature points must equal the exact integral
  ∫∫_T x^i y^j dx dy = i! j! / (i+j+2)!

Monomials with i > j duplicate the i <= j equations for these symmetric groups,
so they are excluded to keep the system non-redundant.

Usage
-----
    uv run python quadrature/solve_triangle.py <order>   (order: 1-5)
"""

import math
import sys

import numpy as np
from scipy.optimize import least_squares


# (n0, n1): number of centroid and 3-point symmetry groups per order
ORDER_CONFIGS: dict[int, tuple[int, int]] = {
    1: (1, 0),
    2: (0, 1),
    3: (1, 1),
    4: (0, 2),
    5: (1, 2),
}


def exact_integral(i: int, j: int) -> float:
    """Exact value of ∫∫_T x^i y^j dx dy over the reference triangle."""
    return math.factorial(i) * math.factorial(j) / math.factorial(i + j + 2)


def n1_sum(a: float, i: int, j: int) -> float:
    """Sum of x^i y^j over the 3 points of an n1 group with parameter a."""
    b = 1.0 - 2.0 * a
    return a**i * a**j + b**i * a**j + a**i * b**j


def independent_monomials(order: int) -> list[tuple[int, int]]:
    """All (i, j) with i+j <= order and i <= j."""
    return [
        (i, total - i)
        for total in range(order + 1)
        for i in range(total + 1)
        if i <= total - i
    ]


def make_residuals(n0: int, n1: int, order: int):
    monomials = independent_monomials(order)
    exact = np.array([exact_integral(i, j) for i, j in monomials])
    # Centroid contributes (1/3)^(i+j) to each monomial — precomputed
    centroid_factor = np.array([(1.0 / 3.0) ** (i + j) for i, j in monomials])

    def residuals(x: np.ndarray) -> np.ndarray:
        computed = np.zeros(len(monomials))

        for k in range(n0):
            computed += x[k] * centroid_factor

        for k in range(n1):
            a = x[n0 + 2 * k]
            w = x[n0 + 2 * k + 1]
            computed += w * np.array([n1_sum(a, i, j) for i, j in monomials])

        return computed - exact

    return residuals


def initial_guess(n0: int, n1: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    w_per_point = 0.5 / (n0 + 3 * n1)
    x0 = []
    for _ in range(n0):
        x0.append(w_per_point)
    for k in range(n1):
        # Spread a values across (0.10, 0.45) so groups start distinguishable.
        # For n1=1 this gives a=0.10; for n1=2, a=0.10 and 0.45.
        a_base = 0.10 + k * (0.35 / max(n1 - 1, 1)) if n1 > 1 else 0.20
        a = float(np.clip(a_base + rng.uniform(-0.03, 0.03), 0.05, 0.48))
        x0.extend([a, 3.0 * w_per_point])
    return np.array(x0)


def solve(order: int) -> tuple[np.ndarray, np.ndarray]:
    n0, n1 = ORDER_CONFIGS[order]
    f = make_residuals(n0, n1, order)

    for seed in range(20):
        result = least_squares(
            f,
            initial_guess(n0, n1, seed=seed),
            method="lm",
            ftol=1e-14,
            xtol=1e-14,
            gtol=1e-14,
        )
        if float(np.max(np.abs(result.fun))) <= 1e-10:
            break
    else:
        max_res = float(np.max(np.abs(result.fun)))
        raise RuntimeError(f"Solver did not converge after 20 attempts: max residual = {max_res:.2e}")

    x = result.x
    points: list[list[float]] = []
    weights: list[float] = []

    for k in range(n0):
        points.append([1.0 / 3.0, 1.0 / 3.0])
        weights.append(float(x[k]))

    for k in range(n1):
        a = float(x[n0 + 2 * k])
        w = float(x[n0 + 2 * k + 1])
        b = 1.0 - 2.0 * a
        points += [[a, a], [b, a], [a, b]]
        weights += [w, w, w]

    return np.array(points), np.array(weights)


def verify(points: np.ndarray, weights: np.ndarray, order: int) -> bool:
    ok = True
    for total in range(order + 1):
        for i in range(total + 1):
            j = total - i
            computed = float(np.dot(weights, points[:, 0] ** i * points[:, 1] ** j))
            expected = exact_integral(i, j)
            err = abs(computed - expected)
            if err > 1e-10:
                print(f"  FAIL ({i},{j}): got {computed:.6e}  expected {expected:.6e}  err={err:.2e}")
                ok = False
    return ok


def main() -> None:
    valid = sorted(ORDER_CONFIGS)
    if len(sys.argv) != 2 or sys.argv[1] not in {str(k) for k in valid}:
        print(f"Usage: {sys.argv[0]} <order>   (order: {valid})")
        sys.exit(1)

    order = int(sys.argv[1])
    points, weights = solve(order)

    print(f"Order {order} triangle quadrature  ({len(points)} points)")
    print(f"  Sum of weights = {weights.sum():.15f}  (expected 0.5)")
    print()
    for i, (p, w) in enumerate(zip(points, weights)):
        print(f"  pt[{i}] = ({p[0]:.15f}, {p[1]:.15f})   w = {w:.15f}")
    print()

    ok = verify(points, weights, order)
    if ok:
        print(f"Exactness check passed: all monomials up to degree {order} correct.")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
