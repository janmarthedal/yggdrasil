"""Representative FEM problems and their solve stages, for benchmarking.

Each ``make_*`` function returns a :class:`Problem` whose attributes are the
individual stages of a finite-element solve (mesh generation, assembly, boundary
extraction, Dirichlet condensation, linear solve, error computation) exposed as
zero-argument callables plus the artifacts they produce.

The call sequences mirror ``tests/system/`` exactly, so a benchmark measures the
same code path the accuracy tests exercise.  This module has no pytest
dependency and can be imported and run on its own.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.typing import NDArray

from yggdrasil import (
    assemble_bilinear_form,
    assemble_load_vector,
    assemble_neumann_bc,
    condense_dirichlet_bc,
    extract_boundary,
    grad_grad_form,
    l2_error,
    l2_project,
    mass_form,
    select_boundary_faces,
    tag_boundary_faces,
    unit_cube_tet_mesh,
    unit_square_tri_mesh,
)
from yggdrasil.assemble import CondensedSystem
from yggdrasil.mesh import Mesh

_EPS = 1e-10


def _no_exact(_u: NDArray[np.float64]) -> float:
    return float("nan")


@dataclass
class Problem:
    """A single benchmark problem with its solve broken into stages.

    Attributes
    ----------
    name : str
        Human-readable identifier.
    mesh : Mesh
        The discretisation (built once at construction).
    build_mesh : callable
        Re-run mesh generation from scratch, returning a fresh ``Mesh``.
    assemble_lhs : callable
        Assemble the left-hand-side operator(s).  Returns the stiffness matrix,
        or a tuple ``(K, M)`` for the eigenvalue problem.
    assemble_rhs : callable
        Assemble the load vector (including any Neumann contribution).
    extract_bnd : callable
        Run ``extract_boundary`` (and any tagging) and return the boundary mesh.
    condense : callable
        Apply the Dirichlet condensation to the cached LHS/RHS, returning a
        ``CondensedSystem``.
    linear_solve : callable
        Solve the cached condensed system, returning the full solution vector.
    has_exact : bool
        Whether ``compute_error`` compares against a known closed-form solution
        (``False`` -> it returns ``nan``).
    compute_error : callable
        ``compute_error(u)`` -> L2 error against the exact solution.
    run_full : callable
        Execute the whole pipeline from scratch and return the full solution
        (or, for the eigenvalue problem, the eigenvalues).
    meta : dict
        Size information: ``num_nodes``, ``num_elements``, ``num_free_dofs``,
        ``nnz``.
    """

    name: str
    mesh: Mesh
    build_mesh: Callable[[], Mesh]
    assemble_lhs: Callable[[], Any]
    assemble_rhs: Callable[[], NDArray[np.float64]]
    extract_bnd: Callable[[], Mesh]
    condense: Callable[[], CondensedSystem]
    linear_solve: Callable[[], NDArray[np.float64]]
    run_full: Callable[[], Any]
    compute_error: Callable[[NDArray[np.float64]], float] = _no_exact
    has_exact: bool = True
    meta: dict[str, int] = field(default_factory=dict)


# --------------------------------------------------------------------------- #
# Exact solutions / source terms (copied from tests/system/*.py)
# --------------------------------------------------------------------------- #

def _sin_sin(x: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])


def _sin_sin_source(x: NDArray[np.float64]) -> NDArray[np.float64]:
    return 2 * np.pi**2 * np.sin(np.pi * x[:, 0]) * np.sin(np.pi * x[:, 1])


# --------------------------------------------------------------------------- #
# Problem builders
# --------------------------------------------------------------------------- #

def make_poisson_2d_dirichlet(n: int, quadrature_order: int = 1) -> Problem:
    """-Δu = f on the unit square, u = 0 on the boundary; u = sin(πx)sin(πy)."""

    def build_mesh() -> Mesh:
        return unit_square_tri_mesh(n)

    mesh = build_mesh()

    def assemble_lhs() -> sp.csr_matrix:
        return assemble_bilinear_form(mesh, grad_grad_form, quadrature_order)

    def assemble_rhs() -> NDArray[np.float64]:
        return assemble_load_vector(mesh, _sin_sin_source, quadrature_order=5)

    def extract_bnd() -> Mesh:
        return extract_boundary(mesh)

    K = assemble_lhs()
    b = assemble_rhs()
    bc_nodes = extract_bnd().point_data["original_node_index"]

    def condense() -> CondensedSystem:
        return condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    system = condense()

    def linear_solve() -> NDArray[np.float64]:
        return system.reconstruct(spla.spsolve(system.K, system.b))

    def compute_error(u: NDArray[np.float64]) -> float:
        return l2_error(mesh, u, _sin_sin, quadrature_order=5)

    def run_full() -> NDArray[np.float64]:
        m = build_mesh()
        kk = assemble_bilinear_form(m, grad_grad_form, quadrature_order)
        bb = assemble_load_vector(m, _sin_sin_source, quadrature_order=5)
        nodes = extract_boundary(m).point_data["original_node_index"]
        sys_ = condense_dirichlet_bc(kk, bb, nodes, bc_val=0.0)
        return sys_.reconstruct(spla.spsolve(sys_.K, sys_.b))

    return Problem(
        name=f"poisson_2d_dirichlet(n={n},q={quadrature_order})",
        mesh=mesh,
        build_mesh=build_mesh,
        assemble_lhs=assemble_lhs,
        assemble_rhs=assemble_rhs,
        extract_bnd=extract_bnd,
        condense=condense,
        linear_solve=linear_solve,
        run_full=run_full,
        compute_error=compute_error,
        meta=_meta(mesh, system, K),
    )


def make_poisson_2d_mixed_bc(n: int) -> Problem:
    """-Δu = 0 on the unit square, u = 0 on x=0, ∂u/∂n = 1 on x=1; u = x."""

    def build_mesh() -> Mesh:
        return unit_square_tri_mesh(n)

    mesh = build_mesh()

    def _tagged_boundary(m: Mesh) -> Mesh:
        bnd = extract_boundary(m)
        bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] < _EPS, tag=1)
        bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] > 1 - _EPS, tag=2)
        return bnd

    def assemble_lhs() -> sp.csr_matrix:
        return assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)

    def assemble_rhs() -> NDArray[np.float64]:
        b = assemble_load_vector(mesh, 0.0, quadrature_order=1)
        neumann = select_boundary_faces(_tagged_boundary(mesh), tag=2)
        b = b + assemble_neumann_bc(neumann, g=1.0, quadrature_order=1, dofs=mesh.num_nodes)
        return b

    def extract_bnd() -> Mesh:
        return _tagged_boundary(mesh)

    K = assemble_lhs()
    b = assemble_rhs()
    bc_nodes = select_boundary_faces(_tagged_boundary(mesh), tag=1).point_data["original_node_index"]

    def condense() -> CondensedSystem:
        return condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    system = condense()

    def linear_solve() -> NDArray[np.float64]:
        return system.reconstruct(spla.spsolve(system.K, system.b))

    def compute_error(u: NDArray[np.float64]) -> float:
        return l2_error(mesh, u, lambda x: x[:, 0], quadrature_order=2)

    def run_full() -> NDArray[np.float64]:
        m = build_mesh()
        kk = assemble_bilinear_form(m, grad_grad_form, quadrature_order=1)
        bb = assemble_load_vector(m, 0.0, quadrature_order=1)
        bnd = extract_boundary(m)
        bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] < _EPS, tag=1)
        bnd = tag_boundary_faces(bnd, lambda c: c[:, 0] > 1 - _EPS, tag=2)
        bb = bb + assemble_neumann_bc(
            select_boundary_faces(bnd, tag=2), g=1.0, quadrature_order=1, dofs=m.num_nodes
        )
        nodes = select_boundary_faces(bnd, tag=1).point_data["original_node_index"]
        sys_ = condense_dirichlet_bc(kk, bb, nodes, bc_val=0.0)
        return sys_.reconstruct(spla.spsolve(sys_.K, sys_.b))

    return Problem(
        name=f"poisson_2d_mixed_bc(n={n})",
        mesh=mesh,
        build_mesh=build_mesh,
        assemble_lhs=assemble_lhs,
        assemble_rhs=assemble_rhs,
        extract_bnd=extract_bnd,
        condense=condense,
        linear_solve=linear_solve,
        run_full=run_full,
        compute_error=compute_error,
        meta=_meta(mesh, system, K),
    )


def make_poisson_3d_dirichlet(n: int) -> Problem:
    """-Δu = 1 on the unit cube, u = 0 on ∂Ω (no closed form)."""

    def build_mesh() -> Mesh:
        return unit_cube_tet_mesh(n)

    mesh = build_mesh()

    def assemble_lhs() -> sp.csr_matrix:
        return assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)

    def assemble_rhs() -> NDArray[np.float64]:
        return assemble_load_vector(mesh, 1.0, quadrature_order=1)

    def extract_bnd() -> Mesh:
        return extract_boundary(mesh)

    K = assemble_lhs()
    b = assemble_rhs()
    bc_nodes = np.unique(extract_bnd().point_data["original_node_index"])

    def condense() -> CondensedSystem:
        return condense_dirichlet_bc(K, b, bc_nodes, bc_val=0.0)

    system = condense()

    def linear_solve() -> NDArray[np.float64]:
        return system.reconstruct(spla.spsolve(system.K, system.b))

    def run_full() -> NDArray[np.float64]:
        m = build_mesh()
        kk = assemble_bilinear_form(m, grad_grad_form, quadrature_order=1)
        bb = assemble_load_vector(m, 1.0, quadrature_order=1)
        nodes = np.unique(extract_boundary(m).point_data["original_node_index"])
        sys_ = condense_dirichlet_bc(kk, bb, nodes, bc_val=0.0)
        return sys_.reconstruct(spla.spsolve(sys_.K, sys_.b))

    return Problem(
        name=f"poisson_3d_dirichlet(n={n})",
        mesh=mesh,
        build_mesh=build_mesh,
        assemble_lhs=assemble_lhs,
        assemble_rhs=assemble_rhs,
        extract_bnd=extract_bnd,
        condense=condense,
        linear_solve=linear_solve,
        run_full=run_full,
        has_exact=False,
        meta=_meta(mesh, system, K),
    )


def make_mass_projection_2d(n: int) -> Problem:
    """L2 projection of sin(πx)sin(πy) onto P1 on the unit square."""

    def build_mesh() -> Mesh:
        return unit_square_tri_mesh(n)

    mesh = build_mesh()

    def assemble_lhs() -> sp.csr_matrix:
        return assemble_bilinear_form(mesh, mass_form, quadrature_order=2)

    def assemble_rhs() -> NDArray[np.float64]:
        return assemble_load_vector(mesh, _sin_sin, quadrature_order=2)

    def extract_bnd() -> Mesh:
        return extract_boundary(mesh)

    M = assemble_lhs()
    b = assemble_rhs()

    def condense() -> CondensedSystem:
        # No BCs for a plain projection; expose an identity "condensation".
        return condense_dirichlet_bc(M, b, np.empty(0, dtype=np.intp))

    system = condense()

    def linear_solve() -> NDArray[np.float64]:
        return spla.spsolve(M, b)

    def compute_error(u: NDArray[np.float64]) -> float:
        return l2_error(mesh, u, _sin_sin, quadrature_order=5)

    def run_full() -> NDArray[np.float64]:
        return l2_project(build_mesh(), _sin_sin, quadrature_order=2)

    return Problem(
        name=f"mass_projection_2d(n={n})",
        mesh=mesh,
        build_mesh=build_mesh,
        assemble_lhs=assemble_lhs,
        assemble_rhs=assemble_rhs,
        extract_bnd=extract_bnd,
        condense=condense,
        linear_solve=linear_solve,
        run_full=run_full,
        compute_error=compute_error,
        meta=_meta(mesh, system, M),
    )


def make_laplace_eigenvalue_2d(n: int, k: int = 6) -> Problem:
    """Smallest k eigenvalues of -Δ on the unit square with u = 0 on ∂Ω."""

    def build_mesh() -> Mesh:
        return unit_square_tri_mesh(n)

    mesh = build_mesh()

    def assemble_lhs() -> tuple[sp.csr_matrix, sp.csr_matrix]:
        K = assemble_bilinear_form(mesh, grad_grad_form, quadrature_order=1)
        M = assemble_bilinear_form(mesh, mass_form, quadrature_order=2)
        return K, M

    def assemble_rhs() -> NDArray[np.float64]:
        return np.zeros(mesh.num_nodes)

    def extract_bnd() -> Mesh:
        return extract_boundary(mesh)

    K, M = assemble_lhs()
    zero = assemble_rhs()
    bc_nodes = extract_bnd().point_data["original_node_index"]

    def condense() -> CondensedSystem:
        return condense_dirichlet_bc(K, zero, bc_nodes, bc_val=0.0)

    system_K = condense()
    system_M = condense_dirichlet_bc(M, zero, bc_nodes, bc_val=0.0)

    def linear_solve() -> NDArray[np.float64]:
        vals, _ = spla.eigsh(system_K.K, M=system_M.K, k=k, which="SM")
        return np.sort(vals)

    def run_full() -> NDArray[np.float64]:
        m = build_mesh()
        kk = assemble_bilinear_form(m, grad_grad_form, quadrature_order=1)
        mm = assemble_bilinear_form(m, mass_form, quadrature_order=2)
        z = np.zeros(m.num_nodes)
        nodes = extract_boundary(m).point_data["original_node_index"]
        sk = condense_dirichlet_bc(kk, z, nodes, bc_val=0.0)
        sm = condense_dirichlet_bc(mm, z, nodes, bc_val=0.0)
        vals, _ = spla.eigsh(sk.K, M=sm.K, k=k, which="SM")
        return np.sort(vals)

    return Problem(
        name=f"laplace_eigenvalue_2d(n={n})",
        mesh=mesh,
        build_mesh=build_mesh,
        assemble_lhs=assemble_lhs,
        assemble_rhs=assemble_rhs,
        extract_bnd=extract_bnd,
        condense=condense,
        linear_solve=linear_solve,
        run_full=run_full,
        has_exact=False,
        meta=_meta(mesh, system_K, K),
    )


def _meta(mesh: Mesh, system: CondensedSystem, operator: Any) -> dict[str, int]:
    return {
        "num_nodes": int(mesh.num_nodes),
        "num_elements": int(mesh.num_elements),
        "num_free_dofs": int(system.K.shape[0]),
        "nnz": int(operator.nnz),
    }


ALL_BUILDERS: dict[str, Callable[..., Problem]] = {
    "poisson_2d_dirichlet": make_poisson_2d_dirichlet,
    "poisson_2d_mixed_bc": make_poisson_2d_mixed_bc,
    "poisson_3d_dirichlet": make_poisson_3d_dirichlet,
    "mass_projection_2d": make_mass_projection_2d,
    "laplace_eigenvalue_2d": make_laplace_eigenvalue_2d,
}
