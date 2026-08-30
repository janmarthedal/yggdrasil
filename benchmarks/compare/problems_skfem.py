"""scikit-fem builders for the problems that overlap with yggdrasil.

Each builder takes a yggdrasil ``Mesh`` and constructs the *same* discretisation
in scikit-fem (identical vertices and connectivity), so the two solvers can be
compared node-for-node as well as by timing.

All scikit-fem imports are deferred to call time; importing this module without
scikit-fem installed is harmless.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray

from yggdrasil.mesh import Mesh


@dataclass
class SkfemProblem:
    name: str
    assemble_lhs: Callable[[], Any]
    assemble_rhs: Callable[[], NDArray[np.float64]]
    solve_full: Callable[[], NDArray[np.float64]]


def _tri_basis(mesh: Mesh, intorder: int = 2):
    from skfem import Basis, ElementTriP1, MeshTri

    group = mesh.element_groups[0]
    m = MeshTri(mesh.nodes.T.copy(), group.connectivity.T.copy().astype(np.int64))
    return m, Basis(m, ElementTriP1(), intorder=intorder)


def _tet_basis(mesh: Mesh, intorder: int = 2):
    from skfem import Basis, ElementTetP1, MeshTet

    group = mesh.element_groups[0]
    m = MeshTet(mesh.nodes.T.copy(), group.connectivity.T.copy().astype(np.int64))
    return m, Basis(m, ElementTetP1(), intorder=intorder)


def poisson_2d_dirichlet(mesh: Mesh) -> SkfemProblem:
    """-Δu = 2π² sin(πx) sin(πy), u = 0 on ∂Ω."""

    def assemble_lhs():
        from skfem.models.poisson import laplace

        _, basis = _tri_basis(mesh)
        return laplace.assemble(basis)

    def assemble_rhs():
        from skfem import LinearForm

        _, basis = _tri_basis(mesh, intorder=5)

        @LinearForm
        def load(v, w):
            x = w.x
            return 2 * np.pi**2 * np.sin(np.pi * x[0]) * np.sin(np.pi * x[1]) * v

        return load.assemble(basis)

    def solve_full():
        from skfem import LinearForm, condense, solve
        from skfem.models.poisson import laplace

        _, basis = _tri_basis(mesh, intorder=5)
        K = laplace.assemble(basis)

        @LinearForm
        def load(v, w):
            x = w.x
            return 2 * np.pi**2 * np.sin(np.pi * x[0]) * np.sin(np.pi * x[1]) * v

        f = load.assemble(basis)
        return solve(*condense(K, f, D=basis.get_dofs()))

    return SkfemProblem("skfem:poisson_2d_dirichlet", assemble_lhs, assemble_rhs, solve_full)


def poisson_3d_dirichlet(mesh: Mesh) -> SkfemProblem:
    """-Δu = 1, u = 0 on ∂Ω."""

    def assemble_lhs():
        from skfem.models.poisson import laplace

        _, basis = _tet_basis(mesh)
        return laplace.assemble(basis)

    def assemble_rhs():
        from skfem.models.poisson import unit_load

        _, basis = _tet_basis(mesh)
        return unit_load.assemble(basis)

    def solve_full():
        from skfem import condense, solve
        from skfem.models.poisson import laplace, unit_load

        _, basis = _tet_basis(mesh)
        K = laplace.assemble(basis)
        f = unit_load.assemble(basis)
        return solve(*condense(K, f, D=basis.get_dofs()))

    return SkfemProblem("skfem:poisson_3d_dirichlet", assemble_lhs, assemble_rhs, solve_full)


def mass_2d(mesh: Mesh) -> SkfemProblem:
    """Consistent mass matrix on P1 triangles."""

    def assemble_lhs():
        from skfem.models.poisson import mass

        _, basis = _tri_basis(mesh)
        return mass.assemble(basis)

    def _unsupported():
        raise NotImplementedError

    return SkfemProblem("skfem:mass_2d", assemble_lhs, _unsupported, _unsupported)


def reorder_to_yggdrasil(skfem_mesh, values: NDArray[np.float64], ygg_mesh: Mesh) -> NDArray[np.float64]:
    """Map a scikit-fem nodal vector back to yggdrasil node ordering.

    ``MeshTri`` / ``MeshTet`` may permute vertices on construction; match by
    coordinate.
    """
    src = np.round(skfem_mesh.p.T, 9)
    dst = np.round(ygg_mesh.nodes, 9)
    src_order = np.lexsort(src.T[::-1])
    dst_order = np.lexsort(dst.T[::-1])
    out = np.empty_like(values)
    out[dst_order] = values[src_order]
    return out
