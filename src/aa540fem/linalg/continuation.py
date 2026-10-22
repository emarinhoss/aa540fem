"""Pseudo-transient continuation for steady nonlinear systems."""

from __future__ import annotations

import numpy as np

from aa540fem.backends.pattern import add_matrices
from aa540fem.linalg.dirichlet import eliminate
from aa540fem.linalg.newton import NewtonResult
from aa540fem.linalg.solvers import LinearSolver


def pseudo_transient(residual_jacobian, U, fixed, M, method="direct", rtol=1e-9, atol=1e-11,
                     dtau0=0.01, max_steps=200, dtau_max=1e6, dtau_growth=1.0, inner_newton=10,
                     inner_rtol=1e-2, verbose=False, residual=None,
                     lagged=None) -> NewtonResult:
    """Pseudo-transient continuation: ``R(U) + (M / dtau)(U - U_k) = 0``.

    Each pseudo-time step is a backward-Euler step solved with up to
    ``inner_newton`` Newton iterations; it is accepted only when its own
    residual has dropped to ``inner_rtol`` of its initial value (otherwise
    the step is retried with ``dtau / 4``: an unconverged implicit step must
    never be accepted, its error accumulates).  ``dtau`` grows by switched
    evolution relaxation ``dtau_{k+1} = dtau_k |R_k| / |R_{k+1}|`` (bounded
    by factors 0.5 and 10) on the steady residual, so the iteration turns
    into plain Newton as the residual vanishes; ``dtau_growth`` is the least
    growth factor of a step that did not increase the residual (SER alone
    keeps ``dtau`` flat through a slow transient).  ``M`` is the mass
    matrix (zero pressure block), so the continuity equation is never
    relaxed.

    ``lagged``: optional callable ``(V, U_k, jacobian) -> (R, J)`` giving the
    residual with some coefficients frozen at the start ``U_k`` of the step
    (``J`` is ``None`` when ``jacobian`` is false).  The inner Newton then
    solves the step with those coefficients lagged, which removes their
    nonlinearity from the implicit step, while the step size control and the
    convergence test use the true residual.  At a steady state the lagged and
    true residuals coincide, so the converged solution is unchanged.
    """
    U = np.array(U, dtype=float, copy=True)
    fixed = np.asarray(fixed, dtype=int)
    zero = np.zeros(fixed.size)

    def evaluate(U, jacobian=True):
        if jacobian or residual is None:
            R, J = residual_jacobian(U)
        else:
            R, J = residual(U), None
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
        # backward-Euler step from U with Newton on G(V) = R(V) + M (V - U)/dtau;
        # the step is accepted only if the inner iteration actually converged
        if lagged is not None:
            RV, JV = lagged(U, U, True)
            RV = np.array(RV, dtype=float, copy=True)
            RV[fixed] = 0.0
            V = U
        else:
            V, RV, JV = U, R, J
        step_res = None
        ok = False
        for inner in range(inner_newton + 1):
            G = RV + M @ ((V - U) / dtau)
            G[fixed] = 0.0
            g = np.linalg.norm(G)
            if not np.isfinite(g):
                break
            if step_res is None:
                step_res = g
            if g <= inner_rtol * step_res or g <= atol:
                ok = True
                break
            if inner == inner_newton:
                break                           # out of inner iterations: reject
            if JV is None:
                if lagged is not None:
                    RV, JV = lagged(V, U, True)
                    RV = np.array(RV, dtype=float, copy=True)
                    RV[fixed] = 0.0
                else:
                    RV, JV, _ = evaluate(V)
            elim = eliminate(add_matrices(JV, M, 1.0 / dtau), fixed)
            solver = (method(elim.K_bc) if callable(method)
                      else LinearSolver(elim.K_bc, method, symmetric=False))
            delta, _ = solver.solve(elim.apply_rhs(-G, zero))
            alpha = 1.0
            for _ in range(4):                  # backtracking on the step residual
                Vn = V + alpha * delta
                if lagged is not None:
                    RVn, JVn = lagged(Vn, U, alpha == 1.0)
                    RVn = np.array(RVn, dtype=float, copy=True)
                    RVn[fixed] = 0.0
                    rVn = np.linalg.norm(RVn)
                else:
                    RVn, JVn, rVn = evaluate(Vn, jacobian=(alpha == 1.0))
                Gn = RVn + M @ ((Vn - U) / dtau)
                Gn[fixed] = 0.0
                if np.isfinite(rVn) and np.linalg.norm(Gn) < g:
                    break
                alpha *= 0.5
            V, RV, JV = Vn, RVn, JVn
            if not np.isfinite(rVn):
                break
        if not ok:
            retries += 1
            dtau *= 0.25
            if retries > 40:
                break
            if verbose:
                print(f"  PTC step rejected, dtau = {dtau:.3e}")
            continue
        if lagged is not None:                  # true residual of the accepted state
            RV, JV, _ = evaluate(V, jacobian=False)
        elif JV is None:                        # Jacobian of the accepted state
            RV, JV, _ = evaluate(V)
        rn = np.linalg.norm(RV)
        ratio = r / max(rn, 1e-300)
        dtau = min(dtau_max, dtau * min(10.0, max(dtau_growth if rn <= r else 0.5, ratio)))
        U, R, J, r = V, RV, JV, rn
        k += 1
        result.T = U
        result.residuals.append(r)
        result.steps.append(dtau)
        if verbose:
            print(f"  PTC {k}: |R| = {r:.3e}, dtau = {dtau:.3e}")
    result.converged = r <= target
    return result
