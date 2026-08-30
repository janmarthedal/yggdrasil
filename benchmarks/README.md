# Benchmark suite

Measures the cost of each stage of a finite-element solve across a range of mesh
resolutions, so optimization work (and regressions) can be tracked. Built on
[`pytest-benchmark`](https://pytest-benchmark.readthedocs.io/).

`benchmarks/` is **not** collected by `uv run pytest` (see `testpaths` in
`pyproject.toml`); it must be named explicitly.

## Running

```sh
# Core suite (numpy/scipy only)
uv run --group bench pytest benchmarks/

# Bigger meshes for profiling
uv run --group bench pytest benchmarks/ --bench-size large      # small | medium | large

# Narrow to one stage
uv run --group bench pytest benchmarks/test_assembly.py

# Save a baseline, then compare a later run against it
uv run --group bench pytest benchmarks/ --benchmark-save=baseline
uv run --group bench pytest benchmarks/ --benchmark-compare=baseline --benchmark-compare-fail=mean:10%

# Machine-readable output
uv run --group bench pytest benchmarks/ --benchmark-json=results.json
```

## scikit-fem comparison (opt-in)

`benchmarks/compare/` times the same problems in
[scikit-fem](https://scikit-fem.readthedocs.io/) where the feature sets overlap
(P1 Poisson assembly and solve on triangles/tetrahedra, consistent mass matrix).
Each test builds the *identical* discretisation in both libraries, asserts the
results agree, then benchmarks each side (`library` parameter = `yggdrasil` /
`skfem`). Without scikit-fem installed these tests are skipped.

```sh
uv run --group bench-compare pytest benchmarks/compare/
```

## What is measured

| File | Stage(s) |
|---|---|
| `test_meshgen.py` | `unit_square_tri_mesh`, `unit_cube_tet_mesh` |
| `test_assembly.py` | `assemble_bilinear_form` (stiffness, mass), `assemble_load_vector` |
| `test_boundary_bc.py` | `extract_boundary`, `tag_boundary_faces`, `condense_dirichlet_bc`, `assemble_neumann_bc`, `l2_project` |
| `test_solve.py` | `spsolve` on the condensed system; full Poisson (2D Dirichlet, 2D mixed, 3D) and Laplace eigenvalue solves |
| `test_error.py` | `l2_error` |

Each benchmark records `num_nodes`, `num_elements`, `num_free_dofs` and matrix
`nnz` in `extra_info` for throughput analysis.

## Scope

Problems and stages mirror `tests/system/` exactly. Only Tri3 (2D) and Tet4 (3D)
are exercised, since those are the only element types with structured mesh
generators today; extend `benchmarks/problems.py` as higher-order / quad / hex
generators land.
