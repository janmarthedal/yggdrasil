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
n2 — all 6 permutations of barycentric (a, b, c) with c = 1-a-b:
     6 Cartesian points, 3 free parameters (a, b, weight)

Moment equations
----------------
For each independent monomial x^i y^j (i <= j, i+j <= order), the sum of
w * x^i * y^j over all quadrature points must equal the exact integral
  ∫∫_T x^i y^j dx dy = i! j! / (i+j+2)!

Monomials with i > j duplicate the i <= j equations for these symmetric groups,
so they are excluded to keep the system non-redundant.

Symmetry group configurations (n0, n1, n2) are from Dunavant (1985) Table IV.
Points or weights may be negative for some orders (the rules are still exact).

Usage
-----
    uv run python quadrature/solve_triangle.py <order>   (order: 1-20)
"""

import math
import sys
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


# (n0, n1, n2) from Dunavant (1985) Table IV — confirmed quadrature rule configurations
ORDER_CONFIGS: dict[int, tuple[int, int, int]] = {
    1:  (1, 0, 0),
    2:  (0, 1, 0),
    3:  (1, 1, 0),
    4:  (0, 2, 0),
    5:  (1, 2, 0),
    6:  (0, 2, 1),
    7:  (1, 2, 1),
    8:  (1, 3, 1),
    9:  (1, 4, 1),
    10: (1, 4, 2),
    11: (0, 5, 2),
    12: (0, 5, 3),
    13: (1, 6, 3),
    14: (0, 6, 4),
    15: (0, 6, 5),
    16: (1, 7, 5),
    17: (1, 8, 6),
    18: (1, 9, 7),
    19: (1, 8, 8),
    20: (1, 10, 8),
}


def exact_integral(i: int, j: int) -> float:
    """Exact value of ∫∫_T x^i y^j dx dy over the reference triangle."""
    return math.factorial(i) * math.factorial(j) / math.factorial(i + j + 2)



def independent_monomials(order: int) -> list[tuple[int, int]]:
    """All (i, j) with i+j <= order and i <= j."""
    return [
        (i, total - i)
        for total in range(order + 1)
        for i in range(total + 1)
        if i <= total - i
    ]


def make_residuals(n0: int, n1: int, n2: int, order: int):
    monomials = independent_monomials(order)
    exact = np.array([exact_integral(i, j) for i, j in monomials])
    is_ = np.array([i for i, j in monomials])
    js_ = np.array([j for i, j in monomials])
    centroid_factor = (1.0 / 3.0) ** (is_ + js_)

    def residuals(x: np.ndarray) -> np.ndarray:
        computed = np.zeros(len(monomials))
        idx = 0

        for _ in range(n0):
            computed += x[idx] * centroid_factor
            idx += 1

        for _ in range(n1):
            a_val, w = x[idx], x[idx + 1]
            b_val = 1.0 - 2.0 * a_val
            api, apj = a_val ** is_, a_val ** js_
            bpi, bpj = b_val ** is_, b_val ** js_
            computed += w * (api * apj + bpi * apj + api * bpj)
            idx += 2

        for _ in range(n2):
            a_val, b_val, w = x[idx], x[idx + 1], x[idx + 2]
            c_val = 1.0 - a_val - b_val
            api, apj = a_val ** is_, a_val ** js_
            bpi, bpj = b_val ** is_, b_val ** js_
            cpi, cpj = c_val ** is_, c_val ** js_
            computed += w * (bpi * cpj + cpi * bpj + api * cpj + cpi * apj + api * bpj + bpi * apj)
            idx += 3

        return computed - exact

    return residuals


def _polar_to_barycentric(r: float, angle_deg: float) -> tuple[float, float]:
    """Convert polar (r, angle in degrees) to barycentric (alpha, beta) for n2 initial guess.

    Uses the polar coordinate system of Dunavant (1985) Figure 2, centred on the
    triangle centroid; angle measured from the median alpha=0.
    """
    angle = math.radians(angle_deg)
    alpha = (1.0 - 2.0 * r * math.cos(angle)) / 3.0
    beta = (1.0 + r * math.cos(angle) - math.sqrt(3.0) * r * math.sin(angle)) / 3.0
    return alpha, beta


def initial_guess(n0: int, n1: int, n2: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    ng = n0 + 3 * n1 + 6 * n2
    w_init = 0.5 / ng  # equal per-point weight sums correctly to 0.5

    x0 = []

    for _ in range(n0):
        x0.append(w_init)

    # n1: alternate between two spread strategies across restarts.
    # Strategy A (even seeds): linear-a from 0.02 to 0.52 — works well when the
    # true a-values are packed in that range (most orders).
    # Strategy B (odd seeds): r-based from r=-0.90 to +0.90, mapping to
    # a=(1+r)/3 ∈ (0.03, 0.63) — covers groups with a>0.50 (outside triangle).
    use_linear_a = (seed % 2 == 0)
    for k in range(n1):
        if use_linear_a:
            a_base = 0.02 + k * (0.50 / max(n1 - 1, 1)) if n1 > 1 else 0.20
        else:
            r_base = -0.90 + k * (1.80 / max(n1 - 1, 1)) if n1 > 1 else 0.0
            a_base = (1.0 + r_base) / 3.0
        a = float(a_base + rng.uniform(-0.10, 0.10))
        x0.extend([a, w_init])

    # n2: polar coordinates. Angles span 30°–55° (Dunavant data lie in this range);
    # r spans 0.40–0.85. Observed true values: r ∈ [0.43, 0.93], θ ∈ [25°, 57°].
    for k in range(n2):
        r_base = 0.40 + k * (0.45 / max(n2 - 1, 1)) if n2 > 1 else 0.60
        ang_base = 30.0 + k * (25.0 / max(n2 - 1, 1)) if n2 > 1 else 45.0
        r = float(np.clip(r_base + rng.uniform(-0.10, 0.10), 0.20, 0.98))
        ang = float(np.clip(ang_base + rng.uniform(-12.0, 12.0), 2.0, 59.0))
        a, b = _polar_to_barycentric(r, ang)
        x0.extend([a, b, w_init])

    return np.array(x0)


@dataclass
class CentroidGroup:
    weight: float


@dataclass
class N1Group:
    a: float      # b = 1 - 2a is implied; a may be outside (0, 0.5) for some rules
    weight: float


@dataclass
class N2Group:
    """6-point group: all permutations of barycentric (a, b, c) with c = 1-a-b."""
    a: float
    b: float
    weight: float


def _lm(f, x0: np.ndarray, max_nfev: int = 5000) -> "least_squares result":
    return least_squares(f, x0, method="lm", ftol=1e-14, xtol=1e-14, gtol=1e-14,
                         max_nfev=max_nfev)


def solve(order: int) -> list[CentroidGroup | N1Group | N2Group]:
    n0, n1, n2 = ORDER_CONFIGS[order]
    f = make_residuals(n0, n1, n2, order)

    best = None
    for seed in range(100):
        result = _lm(f, initial_guess(n0, n1, n2, seed=seed))
        if best is None or np.max(np.abs(result.fun)) < np.max(np.abs(best.fun)):
            best = result
        if float(np.max(np.abs(best.fun))) <= 1e-9:
            break

    best_res = float(np.max(np.abs(best.fun)))
    if best_res > 1e-9:
        if best_res <= 1e-6:
            best = _lm(f, best.x, max_nfev=50_000)
        best_res = float(np.max(np.abs(best.fun)))
        if best_res > 1e-9:
            raise RuntimeError(
                f"Solver did not converge after 100 attempts + polish: max residual = {best_res:.2e}"
            )
    result = best

    x = result.x
    groups: list[CentroidGroup | N1Group | N2Group] = []
    idx = 0

    for _ in range(n0):
        groups.append(CentroidGroup(weight=float(x[idx])))
        idx += 1

    for _ in range(n1):
        groups.append(N1Group(a=float(x[idx]), weight=float(x[idx + 1])))
        idx += 2

    for _ in range(n2):
        groups.append(N2Group(a=float(x[idx]), b=float(x[idx + 1]), weight=float(x[idx + 2])))
        idx += 3

    return groups


def expand_groups(
    groups: list[CentroidGroup | N1Group | N2Group],
) -> tuple[np.ndarray, np.ndarray]:
    """Expand symmetry groups into flat (points, weights) arrays."""
    points: list[list[float]] = []
    weights: list[float] = []
    for g in groups:
        if isinstance(g, CentroidGroup):
            points.append([1.0 / 3.0, 1.0 / 3.0])
            weights.append(g.weight)
        elif isinstance(g, N1Group):
            b = 1.0 - 2.0 * g.a
            points += [[g.a, g.a], [b, g.a], [g.a, b]]
            weights += [g.weight, g.weight, g.weight]
        else:
            a, b, c = g.a, g.b, 1.0 - g.a - g.b
            points += [[b, c], [c, b], [a, c], [c, a], [a, b], [b, a]]
            weights += [g.weight] * 6
    return np.array(points), np.array(weights)


def verify(points: np.ndarray, weights: np.ndarray, order: int) -> bool:
    ok = True
    for total in range(order + 1):
        for i in range(total + 1):
            j = total - i
            computed = float(np.dot(weights, points[:, 0] ** i * points[:, 1] ** j))
            expected = exact_integral(i, j)
            err = abs(computed - expected)
            if err > 1e-9:
                print(f"  FAIL ({i},{j}): got {computed:.6e}  expected {expected:.6e}  err={err:.2e}")
                ok = False
    return ok


def main() -> None:
    valid = sorted(ORDER_CONFIGS)
    if len(sys.argv) != 2 or sys.argv[1] not in {str(k) for k in valid}:
        print(f"Usage: {sys.argv[0]} <order>   (order: {valid})")
        sys.exit(1)

    order = int(sys.argv[1])
    groups = solve(order)
    points, weights = expand_groups(groups)

    print(f"Order {order} triangle quadrature  ({len(points)} points, {len(groups)} groups)")
    print()

    print("Symmetry groups:")
    for g in groups:
        if isinstance(g, CentroidGroup):
            print(f"  n0  (1/3, 1/3)                                                            w = {g.weight:+.15f}")
        elif isinstance(g, N1Group):
            b = 1.0 - 2.0 * g.a
            print(f"  n1  a = {g.a:.15f}  b = {b:.15f}   w = {g.weight:+.15f}")
        else:
            c = 1.0 - g.a - g.b
            print(f"  n2  a = {g.a:.15f}  b = {g.b:.15f}  c = {c:.15f}   w = {g.weight:+.15f}")
    print()

    print("Quadrature points and weights:")
    print(f"  Sum of weights = {weights.sum():.15f}  (expected 0.5)")
    for i, (p, w) in enumerate(zip(points, weights)):
        print(f"  pt[{i:2d}] = ({p[0]:.15f}, {p[1]:.15f})   w = {w:+.15f}")
    print()

    ok = verify(points, weights, order)
    if ok:
        print(f"Exactness check passed: all monomials up to degree {order} correct.")
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()
