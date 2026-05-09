# Deferred Verification: Cross-checks Against External Sources

This document records how to verify triangle quadrature rules against Dunavant (1985) and scikit-fem. These are provenance checks — useful for diagnosing *where* a discrepancy came from if the primary exactness test fails, but not required for correctness.

---

## Coordinate system differences

Three different conventions appear across the sources:

| Source | Coordinates | Weights sum to |
|--------|-------------|----------------|
| Dunavant (1985) | Barycentric (α, β, γ), α+β+γ=1 | 1 |
| scikit-fem `get_quadrature_tri` | Cartesian (x, y), shape `(2, nqp)` | 0.5 |
| yggdrasil `TriangleDomain.quadrature` | Cartesian (x, y), shape `(nqp, 2)` | 0.5 |

---

## Dunavant (1985)

**Reference**: D. A. Dunavant, "High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle", *Int. J. Numer. Methods Eng.*, 21, 1129–1148. PDF at `quadrature/dunavant.pdf`; tabulated data in Appendix II (PDF pages 12–14).

### Coordinate mapping

Dunavant's natural triangle has vertices a=(1,0,0), 2=(0,1,0), 3=(0,0,1) in (α,β,γ). The mapping to yggdrasil's (x,y) is:

```
x = β,  y = γ,  α = 1 - x - y
```

### Weight scaling

Dunavant weights sum to 1; yggdrasil weights sum to 0.5 (area):

```
w_yggdrasil = 0.5 × w_dunavant
```

### Deriving (x, y) points from a Dunavant n₁ group

For a group entry (α=a, β=γ=b), the 3 barycentric permutations are (a,b,b), (b,a,b), (b,b,a). Mapping x=β, y=γ gives Cartesian points: (b,b), (a,b), (b,a).

### Appendix II data for orders 1–5

Columns are (weight, α, β, γ). One row per symmetry group.

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

### Known issue in `triangle.py`

`_order4` has a misleading comment "Dunavant weights already sum to triangle area (0.5)" — the raw Dunavant weights actually sum to 1; the stored values are already halved, which is correct. The comment is wrong.

---

## scikit-fem

**Reference**: `get_quadrature_tri` in `skfem/quadrature.py` (https://github.com/kinnala/scikit-fem).

Uses the same reference triangle (0,0), (1,0), (0,1) and the same weight convention as yggdrasil (sum to 0.5). Points are stored transposed: shape `(2, nqp)` with row 0 = x-coordinates, row 1 = y-coordinates.
