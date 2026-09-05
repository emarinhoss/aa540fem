"""Steady Navier-Stokes (or Stokes) by Newton's method, with pseudo-transient
continuation as the robust fallback at high Reynolds number.
"""

from __future__ import annotations

import warnings

import numpy as np

from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import FlowProblem
from aa540fem.incompressible.solution import FlowSolution
from aa540fem.linalg.dirichlet import DirichletEliminator
from aa540fem.linalg.newton import NewtonResult, newton_iterate
from aa540fem.linalg.solvers import LinearSolver

CONTINUATIONS = ("auto", "newton", "ptc")


def pseudo_transient(residual_jacobian, U, fixed, M, method="direct", rtol=1e-9, atol=1e-11,
                     dtau0=0.01, max_steps=200, dtau_max=1e6, inner_newton=4,
                     inner_rtol=1e-3, verbose=False) -> NewtonResult:
    """Pseudo-transient continuation: ``R(U) + (M / dtau)(U - U_k) = 0``.

    Each pseudo-time step is a backward-Euler step solved with up to
    ``inner_newton`` Newton iterations (to a relative tolerance
    ``inner_rtol`` of the step's initial residual); ``dtau`` grows by
    switched evolution relaxation ``dtau_{k+1} = dtau_k |R_k| / |R_{k+1}|``
    (bounded by factors 0.5 and 10) on the steady residual, so the iteration
    turns into plain Newton as the residual vanishes.  A step whose inner
    iteration diverges is retried with ``dtau / 4``.  ``M`` is the mass
    matrix (zero pressure block), so the continuity equation is never
    relaxed.
    """
    U = np.array(U, dtype=float, copy=True)
    fixed = np.asarray(fixed, dtype=int)
    zero = np.zeros(fixed.size)

    def evaluate(U):
        R, J = residual_jacobian(U)
        R = np.array(R, dtype=float, copy=True)
        R[fixed] = 0.0
        return R, J, np.linalg.norm(R)

    R, J, r = evaluate(U)
    result = NewtonResult(U, [r])
    target = max(atol, rtol * r)
    dtau = dtau0
    if verbose:
        print(f"  PTC 0: |R| = {r:.3e}")
    k = 0
    retries = 0
    while k < max_steps and r > target:
        # backward-Euler step from U with Newton on R(V) + M (V - U)/dtau
        V, RV, JV = U, R, J
        step_res = None
        ok = True
        for _ in range(inner_newton):
            G = RV + M @ ((V - U) / dtau)
            G[fixed] = 0.0
            g = np.linalg.norm(G)
            if step_res is None:
                step_res = g
            if g <= inner_rtol * step_res or g <= atol:
                break
            elim = DirichletEliminator((JV + M / dtau).tocsr(), fixed)
            delta, _ = LinearSolver(elim.K_bc, method, symmetric=False).solve(
                elim.apply_rhs(-G, zero))
            V = V + delta
            RV, JV, rV = evaluate(V)
            if not np.isfinite(rV) or rV > 1e3 * max(r, 1.0):
                ok = False
                break
        if not ok:
            retries += 1
            dtau *= 0.25
            if retries > 20:
                break
            if verbose:
                print(f"  PTC step rejected, dtau = {dtau:.3e}")
            continue
        rn = np.linalg.norm(RV)
        ratio = r / max(rn, 1e-300)
        dtau = min(dtau_max, dtau * min(10.0, max(0.5, ratio)))
        U, R, J, r = V, RV, JV, rn
        k += 1
        result.T = U
        result.residuals.append(r)
        result.steps.append(dtau)
        if verbose:
            print(f"  PTC {k}: |R| = {r:.3e}, dtau = {dtau:.3e}")
    result.converged = r <= target
    return result


def solve_flow(problem: FlowProblem, U0=None, method: str = "direct", verbose: bool = False,
               rtol: float = 1e-9, atol: float = 1e-11, max_newton: int = 30,
               damping: bool = True, stokes: bool = False, continuation: str = "auto",
               dtau0: float = 0.01, max_ptc: int = 200) -> FlowSolution:
    """Steady Navier-Stokes (or Stokes with ``stokes=True``).

    ``continuation``: ``"newton"`` (damped Newton from ``U0``, the first step
    from rest is the Stokes solution), ``"ptc"`` (pseudo-transient
    continuation, see :func:`pseudo_transient`) or ``"auto"`` (Newton, and
    PTC from the initial state if Newton does not converge).  Pass ``U0``
    (a previous ``FlowSolution.U``) for continuation in Reynolds number.
    ``method`` should be ``"direct"``: the saddle-point Jacobian is
    indefinite.
    """
    if continuation not in CONTINUATIONS:
        raise ValueError(f"Unknown continuation {continuation!r}; expected one of {CONTINUATIONS}")
    asm = FlowAssembler(problem)
    fixed, vals = asm.dirichlet()
    F = asm.body_load()
    U = np.zeros(asm.space.ndof) if U0 is None else np.array(U0, dtype=float, copy=True)
    U[fixed] = vals
    last = {}

    def residual_jacobian(U):
        if stokes:
            B = asm.Bx + asm.By
            R = asm.K @ U - F - B.T @ U + B @ U
            J = (asm.K - B.T + B).tocsr()
        else:
            R, J = asm.steady_residual_jacobian(U, F)
        last["J"] = J
        return R, J

    if verbose:
        print(f"Taylor-Hood: {asm.space.N} velocity nodes, {asm.space.Np} pressure nodes, "
              f"{asm.space.ndof} unknowns"
              + (", stabilised" if problem.stabilisation else ""))
    path = "newton"
    res = None
    if continuation in ("auto", "newton") or stokes:
        res = newton_iterate(residual_jacobian, U, fixed, method, rtol=rtol, atol=atol,
                             max_newton=max_newton, damping=damping, verbose=verbose,
                             abort_ratio=100.0 if continuation == "auto" else None)
    if not stokes and (continuation == "ptc" or (continuation == "auto" and not res.converged)):
        if verbose and res is not None:
            print("  Newton did not converge; switching to pseudo-transient continuation")
        path = "ptc" if res is None else "newton+ptc"
        res = pseudo_transient(residual_jacobian, U, fixed, asm.M, method, rtol, atol, dtau0,
                               max_ptc, verbose=verbose)
    if not res.converged:
        warnings.warn(f"steady solve did not converge in {res.iterations} iterations "
                      f"(|R| = {res.residuals[-1]:.2e})", stacklevel=2)
    info = {"iterations": res.iterations, "residuals": res.residuals,
            "converged": res.converged, "stokes": stokes, "continuation": path,
            "stabilisation": problem.stabilisation}
    if path != "newton":
        info["dtau"] = res.steps
    return FlowSolution(problem, asm.space, res.T, info, asm)
