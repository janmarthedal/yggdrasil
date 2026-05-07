# Quadrature Sub-Project

## Goal

**Derive** quadrature weights and points from first principles by solving the moment equations for each reference element — not by copying tables from published sources. The initial focus is the triangle element. Output goes into `yggdrasil/refdomains/triangle.py` (and analogous files for other elements).

The published sources (Dunavant 1985, scikit-fem) serve only as post-hoc verification that the computed rules are correct. See `verification.md`.

## Reference elements in scope

Priority order matches main project's `CLAUDE.md`:
1. **Triangle** — vertices (0,0), (1,0), (0,1), area = 0.5
2. Line, Quadrilateral, Tetrahedron, Hexahedron (future)

---

## Triangle quadrature

### Conventions

- Points are Cartesian (x, y) coordinates on the reference triangle, stored as shape `(nqp, 2)`.
- Weights sum to the area of the reference triangle: `sum(weights) == 0.5`.

### Symmetry groups

Each quadrature rule is built from groups of symmetrically placed points:

| Group | Type | Points per group |
|-------|------|-----------------|
| n₀    | Centroid (1/3, 1/3) | 1 |
| n₁    | (a, a), (b, a), (a, b) — one free parameter | 3 |
| n₂    | all 6 permutations of (a, b, c) — two free parameters | 6 |

Total points: ng = n₀ + 3·n₁ + 6·n₂

For the rules currently in yggdrasil (orders 1–5), only n₀ and n₁ groups appear:

| Order (p) | ng | n₀ | n₁ | n₂ |
|-----------|----|----|----|----|
| 1         | 1  | 1  | 0  | 0  |
| 2         | 3  | 0  | 1  | 0  |
| 3         | 4  | 1  | 1  | 0  |
| 4         | 6  | 0  | 2  | 0  |
| 5         | 7  | 1  | 2  | 0  |

### Verification — exactness test

A quadrature rule of degree p is correct if and only if it integrates every monomial xⁱ yʲ with i+j ≤ p exactly. The exact value over the reference triangle is:

```
∫∫ xⁱ yʲ dx dy = i! j! / (i+j+2)!
```

See `verification.md` for deferred cross-checks against external sources.
