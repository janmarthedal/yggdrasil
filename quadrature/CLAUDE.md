# Quadrature Sub-Project

## Goal

Compute quadrature weights and points for the reference elements used by the main yggdrasil library. The initial focus is the triangle element. Output goes into `yggdrasil/refdomains/triangle.py` (and analogous files for other elements).

## Reference elements in scope

Priority order matches main project's `CLAUDE.md`:
1. **Triangle** — vertices (0,0), (1,0), (0,1), area = 0.5
2. Line, Quadrilateral, Tetrahedron, Hexahedron (future)

---

## Triangle quadrature

### Coordinate systems and conventions

Three different coordinate systems appear across the sources. Keep them straight:

| Source | Coordinates | Weights sum to |
|--------|-------------|----------------|
| Dunavant (1985) | Barycentric (α, β, γ), α+β+γ=1 | 1 |
| scikit-fem `get_quadrature_tri` | Cartesian (x, y), shape `(2, nqp)` | 0.5 |
| yggdrasil `TriangleDomain.quadrature` | Cartesian (x, y), shape `(nqp, 2)` | 0.5 |

### Mapping: Dunavant → yggdrasil

Dunavant's natural triangle has vertices labeled a=(1,0,0), 2=(0,1,0), 3=(0,0,1) in (α,β,γ).
Yggdrasil's reference triangle has vertices (0,0), (1,0), (0,1) in (x,y).

The consistent mapping is:

```
x = β,  y = γ,  α = 1 - x - y
```

Weight scaling (Dunavant weights sum to 1; yggdrasil weights must sum to area = 0.5):

```
w_yggdrasil = 0.5 × w_dunavant
```

### Symmetry groups (Dunavant Table II)

Each quadrature rule is built from groups of symmetrically placed points:

| Group | Type | Points per group | Parameters |
|-------|------|-----------------|------------|
| n₀    | Centroid (1/3, 1/3, 1/3) | 1 | weight w₀ only |
| n₁    | (a, b, b) and permutations | 3 | weight wᵢ, radius rᵢ |
| n₂    | (a, b, c) and all permutations | 6 | weight wᵢ, rᵢ, angle αᵢ |

Total points: ng = n₀ + 3·n₁ + 6·n₂

For the rules used in yggdrasil (orders 1–5), only n₀ and n₁ groups appear (no n₂):

| Order (p) | ng | n₀ | n₁ | n₂ |
|-----------|----|----|----|----|
| 1         | 1  | 1  | 0  | 0  |
| 2         | 3  | 0  | 1  | 0  |
| 3         | 4  | 1  | 1  | 0  |
| 4         | 6  | 0  | 2  | 0  |
| 5         | 7  | 1  | 2  | 0  |

### Dunavant Appendix II data for orders 1–5

Columns are (weight, α, β, γ). One row per symmetry group; multiply by group cardinality for actual points.

**p=1** (ng=1, n₀=1):
```
weight=1.0,  α=β=γ=1/3
```

**p=2** (ng=3, n₁=1):
```
weight=0.333333333333333,  α=0.666666666666667,  β=γ=0.166666666666667
```

**p=3** (ng=4, n₀=1, n₁=1):
```
weight=-0.562500000000000,  α=β=γ=0.333333333333333
weight= 0.520833333333333,  α=0.600000000000000,  β=γ=0.200000000000000
```

**p=4** (ng=6, n₁=2):
```
weight=0.223381589678011,  α=0.108103018168070,  β=γ=0.445948490915965
weight=0.109951743655322,  α=0.816847572980459,  β=γ=0.091576213509771
```

**p=5** (ng=7, n₀=1, n₁=2):
```
weight=0.225000000000000,  α=β=γ=0.333333333333333
weight=0.132394152788506,  α=0.059715871789770,  β=γ=0.470142064105115
weight=0.125939180544827,  α=0.797426985353087,  β=γ=0.101286507323456
```

### Deriving (x, y) points from a Dunavant n₁ group

For a group entry (α=a, β=γ=b):
- The 3 barycentric permutations are (a,b,b), (b,a,b), (b,b,a)
- Mapping x=β, y=γ gives Cartesian points: (b,b), (a,b), (b,a)

### Verification checklist

1. **Weights sum**: `sum(weights) == 0.5`
2. **Dunavant match**: each `w_yggdrasil == 0.5 × w_dunavant` (from Appendix II above)
3. **scikit-fem match**: values in `get_quadrature_tri` in `skfem/quadrature.py` use the same (x,y) convention and the same weights (sum to 0.5); points are stored as shape `(2, nqp)` vs yggdrasil's `(nqp, 2)`
4. **Exactness test**: integrate all monomials xⁱ yʲ with i+j ≤ p over the reference triangle using the quadrature rule and compare against the exact value `i! j! / (i+j+2)!` (times the area factor)

### Known issues in current `triangle.py`

- `_order4` has a misleading comment "Dunavant weights already sum to triangle area (0.5)" — the raw Dunavant weights actually sum to 1; the stored values are already halved, which is correct. The comment is wrong.
- `_order5` has a correct comment ("weights sum to 1; multiply by area") and correctly multiplies by 0.5.

---

## Key references

- **Dunavant (1985)**: D. A. Dunavant, "High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle", *Int. J. Numer. Methods Eng.*, 21, 1129–1148. PDF available locally at `quadrature/dunavant.pdf`. Appendix II (pages 12–14 of the PDF) contains the tabulated data.
- **scikit-fem**: `get_quadrature_tri` in `skfem/quadrature.py` (https://github.com/kinnala/scikit-fem). Uses the same reference triangle and same weight convention as yggdrasil (sum to 0.5). Points stored transposed: shape `(2, nqp)`.
