"""Direct-solver backends behind linalg.direct.factorise."""

import warnings

import numpy as np
import pytest

from aa540fem import geometry
from aa540fem.incompressible import FlowAssembler, FlowProblem
from aa540fem.linalg.direct import (
    available_direct_backends,
    direct_backend,
    estimate_direct_memory,
    factorise,
    petsc_available,
)
from aa540fem.linalg.dirichlet import eliminate
from aa540fem.linalg.solvers import LinearSolver


def saddle_point_system(elems=8, mu=1e-2):
    mesh = geometry(1.0, 1.0, elems, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    asm = FlowAssembler(FlowProblem(mesh, mu=mu, rho=1.0, stabilisation=True,
                                    bc={**walls, "top": (1.0, 0.0)}))
    fixed, vals = asm.dirichlet()
    U = np.zeros(asm.space.ndof)
    U[fixed] = vals
    U[:2 * asm.space.N] += 0.1 * np.sin(np.arange(2 * asm.space.N))
    U[fixed] = vals
    R, J = asm.steady_residual_jacobian(U, asm.body_load())
    elim = eliminate(J, fixed)
    return elim.K_bc, elim.apply_rhs(-R, np.zeros(fixed.size))


@pytest.mark.parametrize("backend", available_direct_backends())
def test_backends_solve_the_saddle_point_jacobian(backend):
    A, b = saddle_point_system()
    f = factorise(A, backend)
    assert f.backend == backend
    x = f.solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-10
    # repeated solves reuse the factorisation
    x2 = f.solve(2 * b)
    assert np.allclose(x2, 2 * x)
    assert f.nnz_factors is None or f.nnz_factors > A.nnz


def test_linear_solver_reports_backend():
    A, b = saddle_point_system(4)
    solver = LinearSolver(A, "direct", symmetric=False)
    x, info = solver.solve(b)
    assert info["backend"] in available_direct_backends()
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-10


def test_missing_backend_falls_back_with_a_warning(monkeypatch):
    monkeypatch.delenv("AA540FEM_DIRECT", raising=False)
    if petsc_available():
        pytest.skip("petsc4py is installed; nothing to fall back from")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        assert direct_backend("petsc") == "superlu"
    assert any("SuperLU" in str(x.message) for x in w)
    with pytest.raises(ValueError):
        direct_backend("pardiso")


def test_environment_selects_backend(monkeypatch):
    monkeypatch.setenv("AA540FEM_DIRECT", "superlu")
    assert direct_backend() == "superlu"
    A, b = saddle_point_system(4)
    assert factorise(A).backend == "superlu"


def test_memory_estimate_is_within_a_factor_of_the_factors():
    A, _ = saddle_point_system(20)
    f = factorise(A, "superlu")
    measured = 12 * f.nnz_factors
    estimate = estimate_direct_memory(A)
    assert measured / 4 < estimate < 4 * measured           # an order of magnitude
    assert estimate_direct_memory(A.nnz, n=A.shape[0]) == estimate
    with pytest.raises(ValueError):
        estimate_direct_memory(A.nnz)


def test_robust_mumps_retry_on_a_badly_scaled_projection_matrix():
    """The projection matrix [[M, -B^T], [B, 0]] of a mesh whose cells span many
    orders of magnitude overflows MUMPS's default workspace (INFOG(1) = -9);
    factorise() retries with more workspace and iterative refinement and stays
    on MUMPS, and the retry leaves the global PETSc options as it found them."""
    pytest.importorskip("petsc4py")
    from aa540fem.linalg.direct import PETScFactorisation, mumps_available

    if not mumps_available():
        pytest.skip("PETSc without MUMPS")
    from petsc4py import PETSc

    # a channel graded geometrically towards the wall: first cell 1e-7 of the height
    mesh = geometry(1.0, 1.0, 12, "quad9")
    eta = mesh.points[:, 1]
    mesh.points[:, 1] = (np.exp(16.0 * eta) - 1.0) / (np.exp(16.0) - 1.0)
    asm = FlowAssembler(FlowProblem(mesh, mu=1e-3, rho=1.0,
                                    bc={"left": (1.0, 0.0), "bottom": (0.0, 0.0),
                                        "top": (1.0, 0.0)}))
    fixed, vals = asm.dirichlet()
    A = eliminate(asm.pattern.matrix(asm.M_data - asm.BT_data + asm.B_data), fixed).K_bc
    b = np.random.default_rng(1).standard_normal(A.shape[0])
    robust = PETScFactorisation(A, robust=True)
    assert np.linalg.norm(A @ robust.solve(b) - b) / np.linalg.norm(b) < 1e-8
    opts = PETSc.Options()
    assert not opts.hasName("mat_mumps_icntl_10") and opts.getInt("mat_mumps_icntl_14") == 40
    with warnings.catch_warnings():
        warnings.simplefilter("error")                      # no fallback to SuperLU
        f = factorise(A, "petsc")
    assert f.backend == "petsc"
    assert np.linalg.norm(A @ f.solve(b) - b) / np.linalg.norm(b) < 1e-8
