"""The residual-only evaluations of Newton and the continuation do not change the iterates."""

import numpy as np

from aa540fem import geometry
from aa540fem.incompressible import FlowAssembler, FlowProblem, solve_flow
from aa540fem.linalg.continuation import pseudo_transient
from aa540fem.linalg.newton import newton_iterate


def cavity(elems=10, mu=0.01):
    mesh = geometry(1.0, 1.0, elems, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    return FlowProblem(mesh, mu=mu, rho=1.0, bc={**walls, "top": (1.0, 0.0)},
                       stabilisation=True)


def counters(asm):
    calls = {"full": 0, "residual": 0}
    fixed, vals = asm.dirichlet()
    F = asm.body_load()

    def rj(U):
        calls["full"] += 1
        return asm.steady_residual_jacobian(U, F)

    def r(U):
        calls["residual"] += 1
        return asm.steady_residual(U, F)

    U = np.zeros(asm.space.ndof)
    U[fixed] = vals
    return rj, r, U, fixed, calls


def test_newton_iterates_identical_with_lazy_residual():
    asm = FlowAssembler(cavity())
    rj, r, U, fixed, calls = counters(asm)
    eager = newton_iterate(rj, U, fixed, rtol=1e-10)
    full_eager = calls["full"]
    calls["full"] = calls["residual"] = 0
    lazy = newton_iterate(rj, U, fixed, rtol=1e-10, residual=r)
    assert eager.converged and lazy.converged
    assert eager.iterations == lazy.iterations and eager.steps == lazy.steps
    assert np.array_equal(eager.T, lazy.T)
    assert np.array_equal(eager.residuals, lazy.residuals)
    assert calls["full"] <= full_eager                 # backtracked trials need no Jacobian
    # with a frozen Jacobian only the factorisations need one
    calls["full"] = calls["residual"] = 0
    frozen = newton_iterate(rj, U, fixed, rtol=1e-10, residual=r, frozen_jacobian=True)
    assert frozen.converged and calls["full"] < frozen.iterations
    assert calls["residual"] >= frozen.iterations


def test_pseudo_transient_identical_with_lazy_residual():
    asm = FlowAssembler(cavity(6))
    rj, r, U, fixed, calls = counters(asm)
    M = asm.M
    eager = pseudo_transient(rj, U, fixed, M, rtol=1e-9, dtau0=1.0, max_steps=30)
    full_eager = calls["full"]
    calls["full"] = calls["residual"] = 0
    lazy = pseudo_transient(rj, U, fixed, M, rtol=1e-9, dtau0=1.0, max_steps=30, residual=r)
    assert eager.converged and lazy.converged
    assert eager.iterations == lazy.iterations
    assert np.allclose(eager.T, lazy.T, atol=1e-12)
    assert np.allclose(eager.residuals, lazy.residuals, rtol=1e-8)
    assert calls["full"] <= full_eager


def test_solve_flow_uses_the_lazy_path_and_keeps_results():
    prob = cavity()
    sol = solve_flow(prob)
    assert sol.info["converged"]
    ref = solve_flow(prob, continuation="ptc")
    assert np.allclose(sol.U, ref.U, atol=1e-7)
