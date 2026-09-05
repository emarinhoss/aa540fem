"""Damped Newton (or Picard) iteration on a residual/Jacobian callable."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from aa540fem.linalg.dirichlet import DirichletEliminator
from aa540fem.linalg.solvers import LinearSolver


@dataclass
class NewtonResult:
    T: np.ndarray
    residuals: list = field(default_factory=list)   # ||R|| per iteration, incl. initial
    converged: bool = False
    steps: list = field(default_factory=list)       # damping factor of each update

    @property
    def iterations(self) -> int:
        return len(self.residuals) - 1


def newton_iterate(residual_jacobian, T, nodes, method: str = "direct", tol: float = 1e-10,
           maxiter=None, rtol: float = 1e-10, atol: float = 1e-12, max_newton: int = 25,
           damping: bool = True, verbose: bool = False,
           frozen_jacobian: bool = False, abort_ratio: float | None = None) -> NewtonResult:
    """Solve ``R(T) = 0`` with (damped) Newton iterations.

    Parameters
    ----------
    residual_jacobian : callable ``T -> (R, J)`` returning the residual
                        vector and its Jacobian (sparse); rows of ``nodes``
                        are ignored (Dirichlet nodes keep their values).
    T                 : initial guess, already satisfying the Dirichlet values.
    nodes             : Dirichlet node indices (no update there).
    method, tol, maxiter : linear solver settings (``cg`` is refused: the
                        Jacobian is not symmetric).
    rtol, atol        : stop when ``||R|| <= max(atol, rtol * ||R_0||)``.
    max_newton        : iteration limit.
    damping           : halve the update (up to 5 times) while the residual
                        does not decrease.
    frozen_jacobian   : reuse the factorised Jacobian for the following
                        iterations (modified Newton) as long as each one
                        reduces the residual by at least a factor 3;
                        otherwise it is refreshed.  The residual is always the
                        true one, so the converged answer is unchanged; useful
                        in time stepping where the Jacobian changes little
                        per step.
    abort_ratio       : stop early (not converged) as soon as the residual
                        exceeds this multiple of the initial residual, i.e.
                        Newton is diverging.
    """
    T = np.array(T, dtype=float, copy=True)
    nodes = np.asarray(nodes, dtype=int)

    def evaluate(T):
        R, J = residual_jacobian(T)
        R = np.array(R, dtype=float, copy=True)
        R[nodes] = 0.0
        return R, J, np.linalg.norm(R)

    R, J, r = evaluate(T)
    result = NewtonResult(T, [r])
    target = max(atol, rtol * r)
    if verbose:
        print(f"  Newton 0: |R| = {r:.3e}")

    solver = None
    for it in range(1, max_newton + 1):
        if r <= target:
            result.converged = True
            break
        if solver is None or not frozen_jacobian:
            elim = DirichletEliminator(J, nodes)
            solver = LinearSolver(elim.K_bc, method, tol, maxiter, symmetric=False)
        rhs = elim.apply_rhs(-R, np.zeros(nodes.size))
        delta, _ = solver.solve(rhs)

        alpha = 1.0
        for _ in range(6):
            Tn = T + alpha * delta
            Rn, Jn, rn = evaluate(Tn)
            if not damping or rn < r or alpha < 1.0 / 16:
                break
            alpha *= 0.5
        if frozen_jacobian and rn > r / 3.0:
            solver = None                       # poor contraction: refresh the Jacobian
        T, R, J, r = Tn, Rn, Jn, rn
        result.T = T
        result.residuals.append(r)
        result.steps.append(alpha)
        if verbose:
            print(f"  Newton {it}: |R| = {r:.3e}" + (f" (damping {alpha})" if alpha < 1 else ""))
        diverging = not np.isfinite(r) or r > (abort_ratio or np.inf) * result.residuals[0]
        if abort_ratio is not None and diverging:
            if verbose:
                print("  Newton diverging; aborted")
            break
    else:
        result.converged = r <= target
    return result
