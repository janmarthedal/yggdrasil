#!/usr/bin/env python3
"""
Identify symmetric triangle quadrature rules using the Witherden-Vincent method.

Reference: Witherden & Vincent (2015), arXiv:1409.1865.

Triangle: vertices (0,0), (1,0), (0,1), area = 0.5.

Symmetry orbits in barycentric coordinates (λ₁+λ₂+λ₃=1):
  S₁       — centroid (1/3,1/3,1/3)    — 1 pt,  no free params
  S₂(α)    — Perm(α, α, 1-2α)          — 3 pts, 1 free param  (0 < α < 1/2)
  S₃(α,β)  — Perm(α, β, 1-α-β)         — 6 pts, 2 free params (α,β>0, α+β<1, α≠β)

Method:
  Points are parameterised by orbital parameters only. Given a candidate point set,
  weights are determined implicitly by solving the linear system Aω = b via SVD,
  where A[k,j] = Σ_{x in orbit j} ψₖ(x) and b = exact integrals of basis functions.
  The residual Aω - b is minimised with the Levenberg-Marquardt algorithm.

Usage:
  uv run python quadrature/witherden-vincent/solve_triangle.py <Np> <phi> [time_sec]
  Example: uv run python ... 7 5        # find φ=5 rule with 7 points
"""

import math
import sys
import time
from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares
from scipy.special import eval_jacobi


# ── Orthonormal Dubiner basis on T = {(x,y) : x,y≥0, x+y≤1} ─────────────────
#
# ψᵢⱼ(x,y) = √(2(2i+1)(i+j+1)) · Pᵢ^(0,0)(2r-1) · Pⱼ^(2i+1,0)(2y-1) · (1-y)^i
#   where r = x/(1-y)   (collapsed Duffy coordinate)
#
# Orthonormality: ∫_T ψᵢⱼ ψₘₙ dx dy = δᵢₘ δⱼₙ
# Only the constant mode integrates non-zero: ∫_T ψ₀₀ dx dy = √(1/2) = 1/√2

def _eval_basis(i: int, j: int, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    r = np.where(y < 1.0, x / np.maximum(1.0 - y, 1e-300), 0.0)
    norm = math.sqrt(2.0 * (2*i + 1) * (i + j + 1))
    return norm * eval_jacobi(i, 0, 0, 2*r - 1) * eval_jacobi(j, 2*i+1, 0, 2*y - 1) * (1.0 - y)**i


def _basis_list(phi: int) -> list[tuple[int, int]]:
    """All (i,j) with i,j≥0 and i+j≤phi, ordered by total degree."""
    return [(i, d - i) for d in range(phi + 1) for i in range(d + 1)]


# ── Symmetry orbit expansion ──────────────────────────────────────────────────
# Cartesian coords from barycentric (λ₁,λ₂,λ₃): (x,y) = (λ₂, λ₃)

def _expand_s1() -> np.ndarray:
    return np.array([[1/3, 1/3]])


def _expand_s2(alpha: float) -> np.ndarray:
    c = 1.0 - 2.0 * alpha
    return np.array([[alpha, c], [c, alpha], [alpha, alpha]])


def _expand_s3(alpha: float, beta: float) -> np.ndarray:
    gamma = 1.0 - alpha - beta
    return np.array([
        [beta, gamma], [gamma, beta],
        [alpha, gamma], [gamma, alpha],
        [alpha, beta], [beta, alpha],
    ])


def _clamp_s2(alpha: float, eps: float = 1e-8) -> float:
    return float(np.clip(alpha, eps, 0.5 - eps))


def _clamp_s3(alpha: float, beta: float, eps: float = 1e-8) -> tuple[float, float]:
    alpha = float(np.clip(alpha, eps, 1.0 - 2*eps))
    beta  = float(np.clip(beta,  eps, 1.0 - alpha - eps))
    return alpha, beta


# ── Decompositions ────────────────────────────────────────────────────────────

def _decompositions(Np: int):
    """Yield all (n0, n1, n2) with Np = n0·1 + n1·3 + n2·6, n0 ∈ {0,1}."""
    for n0 in [0, 1]:
        r0 = Np - n0
        if r0 < 0:
            continue
        for n1 in range(r0 // 3 + 1):
            r1 = r0 - 3 * n1
            if r1 % 6 == 0:
                yield (n0, n1, r1 // 6)


# ── Residual for Levenberg-Marquardt ──────────────────────────────────────────

_AREA = 0.5
_B0 = math.sqrt(_AREA)  # ∫_T ψ₀₀ dx dy


def _build_A(n0: int, n1: int, n2: int, params: np.ndarray,
             basis: list[tuple[int, int]]) -> tuple[np.ndarray, list[np.ndarray]]:
    """Build matrix A and collect orbit point arrays from orbital parameters."""
    idx = 0
    orbit_pts: list[np.ndarray] = []

    if n0:
        orbit_pts.append(_expand_s1())
    for _ in range(n1):
        orbit_pts.append(_expand_s2(_clamp_s2(params[idx])))
        idx += 1
    for _ in range(n2):
        alpha, beta = _clamp_s3(params[idx], params[idx + 1])
        orbit_pts.append(_expand_s3(alpha, beta))
        idx += 2

    nb = len(basis)
    A = np.empty((nb, len(orbit_pts)))
    for k, pts in enumerate(orbit_pts):
        xk, yk = pts[:, 0], pts[:, 1]
        for bi, (i, j) in enumerate(basis):
            A[bi, k] = np.sum(_eval_basis(i, j, xk, yk))

    return A, orbit_pts


def _make_residual(n0: int, n1: int, n2: int, phi: int):
    basis = _basis_list(phi)
    b = np.zeros(len(basis))
    b[0] = _B0

    def residual(params: np.ndarray) -> np.ndarray:
        A, _ = _build_A(n0, n1, n2, params, basis)
        omega, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
        return A @ omega - b

    return residual


# ── Seeding ───────────────────────────────────────────────────────────────────

def _seed(n1: int, n2: int, rng: np.random.Generator) -> np.ndarray:
    params = []
    for _ in range(n1):
        params.append(rng.uniform(0.05, 0.45))
    for _ in range(n2):
        lam = np.sort(rng.dirichlet([1.0, 1.0, 1.0]))
        # lam[0] < lam[1] < lam[2] (with probability 1)
        params.extend([float(lam[0]), float(lam[1])])
    return np.array(params)


# ── Rule dataclass ────────────────────────────────────────────────────────────

@dataclass
class Rule:
    n0: int
    n1: int
    n2: int
    orbital_params: np.ndarray  # free orbital parameters found by LMA
    points: np.ndarray          # shape (Np, 2) Cartesian
    weights: np.ndarray         # shape (Np,)  per-point weights
    residual: float             # max |Aω - b|
    phi: int
    is_pi: bool                 # positive weights and interior points


def _reconstruct(n0: int, n1: int, n2: int, params: np.ndarray, phi: int) -> Rule:
    basis = _basis_list(phi)
    b = np.zeros(len(basis))
    b[0] = _B0

    A, orbit_pts = _build_A(n0, n1, n2, params, basis)
    omega, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    residual = float(np.max(np.abs(A @ omega - b)))

    points = np.vstack(orbit_pts)
    weights = np.concatenate([np.full(len(pts), omega[k]) for k, pts in enumerate(orbit_pts)])

    eps = 1e-10
    interior = bool(np.all(points > -eps) and np.all(points.sum(axis=1) < 1 + eps))
    is_pi = bool(np.all(weights > -eps) and interior)

    return Rule(n0=n0, n1=n1, n2=n2, orbital_params=params,
                points=points, weights=weights, residual=residual,
                phi=phi, is_pi=is_pi)


# ── Verification ──────────────────────────────────────────────────────────────

def _exact_monomial(i: int, j: int) -> float:
    return math.factorial(i) * math.factorial(j) / math.factorial(i + j + 2)


def _verify(rule: Rule) -> bool:
    ok = True
    for d in range(rule.phi + 1):
        for i in range(d + 1):
            j = d - i
            computed = float(np.dot(rule.weights, rule.points[:, 0]**i * rule.points[:, 1]**j))
            expected = _exact_monomial(i, j)
            err = abs(computed - expected)
            if err > 1e-9:
                print(f"  FAIL monomial ({i},{j}): got {computed:.10e}, expected {expected:.10e}, err={err:.2e}")
                ok = False
    return ok


# ── Main solve loop ───────────────────────────────────────────────────────────

def solve(Np: int, phi: int, time_budget: float = 60.0) -> list[Rule]:
    decomps = list(_decompositions(Np))
    print(f"Nₚ={Np}, φ={phi}: {len(decomps)} decomposition(s) to try")

    rng = np.random.default_rng(0)
    found: list[Rule] = []

    # Split time budget evenly across decompositions that have free params.
    n_parameterised = sum(1 for n0, n1, n2 in decomps if n1 + 2*n2 > 0)
    per_decomp = time_budget / max(n_parameterised, 1)

    for n0, n1, n2 in decomps:
        n_params = n1 + 2 * n2
        label = f"(n₀={n0}, n₁={n1}, n₂={n2})"

        if n_params == 0:
            basis = _basis_list(phi)
            b = np.zeros(len(basis))
            b[0] = _B0
            A, _ = _build_A(n0, n1, n2, np.array([]), basis)
            omega, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            res = float(np.max(np.abs(A @ omega - b)))
            if res < 1e-10:
                rule = _reconstruct(n0, n1, n2, np.array([]), phi)
                print(f"  {label} found (no free params), residual={rule.residual:.2e}, PI={rule.is_pi}")
                found.append(rule)
            else:
                print(f"  {label} no free params, residual={res:.2e} — no rule")
            continue

        f = _make_residual(n0, n1, n2, phi)
        best_result = None
        best_res = np.inf
        attempts = 0
        t_decomp = time.time()

        while time.time() - t_decomp < per_decomp:
            x0 = _seed(n1, n2, rng)
            result = least_squares(f, x0, method="lm",
                                   ftol=1e-15, xtol=1e-15, gtol=1e-15,
                                   max_nfev=10_000)
            res = float(np.max(np.abs(result.fun)))
            attempts += 1

            if res < best_res:
                best_res = res
                best_result = result

            if best_res < 1e-10:
                break

        if best_res < 1e-10:
            rule = _reconstruct(n0, n1, n2, best_result.x, phi)
            print(f"  {label} found after {attempts} attempt(s), "
                  f"residual={rule.residual:.2e}, PI={rule.is_pi}")
            found.append(rule)
        else:
            why = "time" if time.time() - t_decomp >= per_decomp else "not converged"
            print(f"  {label} {why} after {attempts} attempt(s), best residual={best_res:.2e}")

    return found


# ── Output ────────────────────────────────────────────────────────────────────

def _print_rule(rule: Rule) -> None:
    Np = len(rule.points)
    print(f"\nRule: φ={rule.phi}, Nₚ={Np}, "
          f"(n₀={rule.n0}, n₁={rule.n1}, n₂={rule.n2}), "
          f"PI={rule.is_pi}, residual={rule.residual:.2e}")

    # Orbital parameters
    idx = 0
    orbit_idx = 0
    if rule.n0:
        print(f"  S₁  centroid (1/3, 1/3)   w = {rule.weights[orbit_idx]:+.15f}")
        orbit_idx += 1
    for _ in range(rule.n1):
        a = float(np.clip(rule.orbital_params[idx], 1e-8, 0.5 - 1e-8))
        c = 1.0 - 2 * a
        w = rule.weights[orbit_idx]
        print(f"  S₂  α={a:.15f}  (1-2α={c:.15f})   w={w:+.15f}")
        idx += 1
        orbit_idx += 3
    for _ in range(rule.n2):
        a, b = _clamp_s3(rule.orbital_params[idx], rule.orbital_params[idx + 1])
        c = 1.0 - a - b
        w = rule.weights[orbit_idx]
        print(f"  S₃  α={a:.15f}  β={b:.15f}  γ={c:.15f}   w={w:+.15f}")
        idx += 2
        orbit_idx += 6

    print(f"  Sum of weights = {rule.weights.sum():.15f}  (expected {_AREA})")

    print("  Quadrature points and weights:")
    for k, (p, w) in enumerate(zip(rule.points, rule.weights)):
        print(f"    pt[{k:2d}] = ({p[0]:.15f}, {p[1]:.15f})   w = {w:+.15f}")


def main() -> None:
    args = sys.argv[1:]
    if len(args) < 2 or len(args) > 3:
        print(f"Usage: {sys.argv[0]} <Np> <phi> [time_sec]")
        print("  Np       — number of quadrature points")
        print("  phi      — target quadrature strength (polynomial degree)")
        print("  time_sec — search time budget per run (default: 60)")
        sys.exit(1)

    Np   = int(args[0])
    phi  = int(args[1])
    tlim = float(args[2]) if len(args) == 3 else 60.0

    rules = solve(Np, phi, tlim)

    if not rules:
        print("\nNo rules found.")
        sys.exit(1)

    for rule in rules:
        _print_rule(rule)

        ok = _verify(rule)
        if ok:
            print(f"  Exactness check passed: all monomials up to degree {rule.phi} correct.")
        else:
            print("  Exactness check FAILED.")


if __name__ == "__main__":
    main()
