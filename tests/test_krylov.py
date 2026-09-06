"""PETSc fieldsplit Krylov solver for the saddle-point Jacobians (skipped without petsc4py)."""

import numpy as np
import pytest

from aa540fem import geometry
from aa540fem.incompressible import FlowAssembler, FlowProblem, solve_flow
from aa540fem.linalg.direct import factorise, petsc_available
from aa540fem.linalg.dirichlet import eliminate

if not petsc_available():
    pytest.skip("petsc4py not available", allow_module_level=True)

from aa540fem.linalg.krylov import FieldSplitSolver, pressure_mass_matrix  # noqa: E402

SIDES = ("top", "right", "left", "bottom")


def poiseuille(elems=8):
    mesh = geometry(2.0, 1.0, elems, "quad9")
    inflow = lambda x, y: 4 * y * (1 - y)
    return FlowProblem(mesh, mu=1.0, rho=1.0, stabilisation=True,
                       bc={"left": (inflow, 0.0), "bottom": (0.0, 0.0), "top": (0.0, 0.0),
                           "right": "open"})


def cavity(elems=12, mu=0.01):
    mesh = geometry(1.0, 1.0, elems, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    return FlowProblem(mesh, mu=mu, rho=1.0, bc={**walls, "top": (1.0, 0.0)}, stabilisation=True)


def system(prob, perturb=0.0):
    asm = FlowAssembler(prob)
    fixed, vals = asm.dirichlet()
    U = np.zeros(asm.space.ndof)
    U[fixed] = vals
    if perturb:
        U[:2 * asm.space.N] += perturb * np.sin(np.arange(2 * asm.space.N))
        U[fixed] = vals
    R, J = asm.steady_residual_jacobian(U, asm.body_load())
    elim = eliminate(J, fixed)
    return asm, elim.K_bc, elim.apply_rhs(-R, np.zeros(fixed.size))


@pytest.mark.parametrize("case,velocity_pc", [("poiseuille", "ilu"), ("cavity", "ilu"),
                                              ("cavity", "gamg")])
def test_fieldsplit_matches_direct(case, velocity_pc):
    prob = poiseuille() if case == "poiseuille" else cavity()
    asm, A, b = system(prob, perturb=0.1)
    ref = factorise(A, "superlu").solve(b)
    solver = FieldSplitSolver(A, 2 * asm.space.N, pressure_mass_matrix(asm), rtol=1e-10,
                              velocity_pc=velocity_pc)
    x = solver.solve(b)
    assert np.linalg.norm(x - ref) / np.linalg.norm(ref) < 1e-6
    assert 0 < solver.iterations < 300


def test_solve_flow_with_fieldsplit_reproduces_direct():
    prob = cavity(10)
    direct = solve_flow(prob)
    krylov = solve_flow(prob, method="fieldsplit")
    assert krylov.info["converged"]
    assert np.allclose(krylov.U, direct.U, atol=1e-6)


def test_gpu_path_requires_cuda_build():
    from petsc4py import PETSc

    if not PETSc.Sys.hasExternalPackage("cuda"):
        pytest.skip("PETSc built without CUDA")
    asm, A, b = system(cavity(6))
    solver = FieldSplitSolver(A, 2 * asm.space.N, pressure_mass_matrix(asm), gpu=True)
    x = solver.solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-6
