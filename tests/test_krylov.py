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


def lumped_velocity_mass(asm):
    return np.asarray(abs(asm.M).sum(axis=1)).ravel()[:2 * asm.space.N]


@pytest.mark.parametrize("case,velocity_pc,schur", [("poiseuille", "ilu", "mass"),
                                                    ("poiseuille", "ilu", "lsc"),
                                                    ("cavity", "ilu", "lsc"),
                                                    ("cavity", "gamg", "lsc")])
def test_fieldsplit_matches_direct(case, velocity_pc, schur):
    prob = poiseuille() if case == "poiseuille" else cavity()
    asm, A, b = system(prob, perturb=0.1)
    ref = factorise(A, "superlu").solve(b)
    solver = FieldSplitSolver(A, 2 * asm.space.N, pressure_mass_matrix(asm), rtol=1e-10,
                              velocity_pc=velocity_pc, schur=schur,
                              Qdiag=lumped_velocity_mass(asm))
    x = solver.solve(b)
    assert np.linalg.norm(x - ref) / np.linalg.norm(ref) < 1e-6
    assert 0 < solver.iterations < 100


def converged_jacobian(prob, dt=None):
    """Eliminated Jacobian at the converged steady state (or the theta-scheme
    Jacobian with step ``dt`` there), a right-hand side and the assembler."""
    sol = solve_flow(prob)
    asm = FlowAssembler(prob)
    fixed, _ = asm.dirichlet()
    if dt is None:
        _, J = asm.steady_residual_jacobian(sol.U, asm.body_load())
    else:
        mt = asm.momentum_terms(sol.U, t=0.0, dt=dt, U_old=sol.U, jacobian=True)
        J = asm.pattern.matrix(asm.M_data / dt + 0.5 * (asm.K_data + mt.JN_data) + mt.JS_data
                               - asm.BT_data + asm.B_data)
    b = np.random.default_rng(1).standard_normal(J.shape[0])
    b[fixed] = 0.0
    return asm, eliminate(J, fixed).K_bc, b


@pytest.mark.parametrize("case,bound", [("cavity Re 100", 80), ("cavity Re 1000", 120),
                                        ("cylinder Re 20", 200), ("cylinder theta", 60)])
def test_lsc_iteration_counts_on_converged_jacobians(case, bound):
    """The LSC Schur preconditioner with ILU(0) on the velocity block converges
    on the stabilised Jacobians of the validation cases within recorded bounds
    (measured 35, 75, 133 and 28 iterations; the mass-matrix preconditioner
    does not converge on the last three within 300)."""
    import pathlib

    from aa540fem import read_mesh

    if case.startswith("cavity"):
        prob = cavity(24 if "100" in case else 32, 0.01 if "100" in case else 0.001)
        dt = None
    else:
        mesh = read_mesh(pathlib.Path(__file__).resolve().parent.parent / "examples" / "meshes"
                         / "cylinder_bl.msh")
        umax = 0.3 if "20" in case else 1.5
        inflow = lambda x, y: 4.0 * umax * y * (0.41 - y) / 0.41 ** 2
        prob = FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=True,
                           bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0),
                               "cylinder": (0.0, 0.0), "outlet": "open"})
        dt = None if "20" in case else 0.005
    asm, A, b = converged_jacobian(prob, dt)
    ref = factorise(A).solve(b)
    solver = FieldSplitSolver(A, 2 * asm.space.N, pressure_mass_matrix(asm), rtol=1e-8,
                              schur="lsc", Qdiag=lumped_velocity_mass(asm))
    x = solver.solve(b)
    assert np.linalg.norm(x - ref) / np.linalg.norm(ref) < 1e-6
    assert 0 < solver.iterations < bound


def test_lsc_passes_pinned_pressure_dof_through():
    """An enclosed flow pins one pressure dof; the LSC operator must not be singular."""
    from aa540fem.linalg.krylov import LSCPreconditioner

    asm, A, b = system(cavity(6))
    assert asm.problem.pins_pressure
    lsc = LSCPreconditioner(A, 2 * asm.space.N, Qdiag=lumped_velocity_mass(asm))
    assert lsc.pin.getArray().sum() == 1.0


def test_solve_flow_with_fieldsplit_reproduces_direct():
    prob = cavity(10)
    direct = solve_flow(prob)
    krylov = solve_flow(prob, method="fieldsplit")
    assert krylov.info["converged"]
    assert krylov.info["iterations"] <= direct.info["iterations"] + 2
    assert np.allclose(krylov.U, direct.U, atol=1e-6)


def test_theta_scheme_with_fieldsplit_reproduces_direct():
    from aa540fem.incompressible import solve_flow_transient

    prob = cavity(8)
    runs = {m: solve_flow_transient(prob, dt=0.05, t_end=0.2, scheme="theta", method=m,
                                    startup_steps=1, rtol=1e-8) for m in ("direct", "fieldsplit")}
    assert np.allclose(runs["fieldsplit"].snapshots[-1], runs["direct"].snapshots[-1], atol=1e-6)


def test_run_configuration_selects_the_krylov_path(monkeypatch):
    """``--linear petsc`` (``AA540FEM_LINEAR``) routes ``method="direct"`` to fieldsplit."""
    from aa540fem import hardware
    from aa540fem.incompressible import steady

    calls = []
    original = steady.fieldsplit_factory
    monkeypatch.setattr(steady, "fieldsplit_factory",
                        lambda asm, **kw: calls.append(kw) or original(asm, **kw))
    monkeypatch.setenv("AA540FEM_LINEAR", "petsc")
    hardware.configure()
    try:
        sol = solve_flow(cavity(6))
    finally:
        monkeypatch.delenv("AA540FEM_LINEAR")
        hardware.configure()
    assert sol.info["converged"] and len(calls) == 1 and calls[0]["gpu"] is False


def test_gpu_path_requires_cuda_build():
    from petsc4py import PETSc

    if not PETSc.Sys.hasExternalPackage("cuda"):
        pytest.skip("PETSc built without CUDA")
    asm, A, b = system(cavity(6))
    solver = FieldSplitSolver(A, 2 * asm.space.N, pressure_mass_matrix(asm), gpu=True)
    x = solver.solve(b)
    assert np.linalg.norm(A @ x - b) / np.linalg.norm(b) < 1e-6
