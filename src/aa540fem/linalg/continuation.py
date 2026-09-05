"""Pseudo-transient continuation for steady nonlinear systems."""

from __future__ import annotations

import numpy as np

from aa540fem.linalg.dirichlet import DirichletEliminator
from aa540fem.linalg.newton import NewtonResult
from aa540fem.linalg.solvers import LinearSolver


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
    iteration diverges (its own residual growing tenfold, or non-finite) is
    retried with ``dtau / 4``.  ``M`` is the mass
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
        # backward-Euler step from U with Newton on G(V) = R(V) + M (V - U)/dtau
        V, RV, JV = U, R, J
        step_res = None
        g_prev = None
        ok = True
        for _ in range(inner_newton):
            G = RV + M @ ((V - U) / dtau)
            G[fixed] = 0.0
            g = np.linalg.norm(G)
            if not np.isfinite(g) or (g_prev is not None and g > 10.0 * g_prev):
                ok = False                      # the inner iteration is diverging
                break
            g_prev = g
            if step_res is None:
                step_res = g
            if g <= inner_rtol * step_res or g <= atol:
                break
            elim = DirichletEliminator((JV + M / dtau).tocsr(), fixed)
            delta, _ = LinearSolver(elim.K_bc, method, symmetric=False).solve(
                elim.apply_rhs(-G, zero))
            V = V + delta
            RV, JV, rV = evaluate(V)
        if ok and not np.isfinite(np.linalg.norm(RV)):
            ok = False
        if not ok:
            retries += 1
            dtau *= 0.25
            if retries > 40:
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
