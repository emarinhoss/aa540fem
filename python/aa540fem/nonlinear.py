"""Newton's method for temperature-dependent coefficients.

A coefficient callable with a parameter named ``T`` (see
:func:`aa540fem.util.call_coeff`) makes the problem nonlinear:

    u . grad T - div( kappa(x, y, T) grad T ) = f(x, y, T)

The residual ``R(T) = A(T) T - F(T)`` is driven to zero with Newton's
method; the Jacobian ``A + dA`` carries the derivatives of ``kappa`` and
``f`` with respect to ``T`` (central finite differences at the quadrature
points, see :func:`aa540fem.element.elem_operators`).  ``newton=False``
drops ``dA`` and gives Picard (fixed-point) iteration.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

from .boundary import DirichletEliminator
from .solver import (
    LinearSolver,
    Problem,
    Solution,
    assemble_operators,
    dirichlet_data,
    neumann_loads,
)
from .util import values_at


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
           frozen_jacobian: bool = False) -> NewtonResult:
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
    else:
        result.converged = r <= target
    return result


def solve_nonlinear(problem: Problem, verbose: bool = False, method: str = "direct",
                    tol: float = 1e-10, maxiter=None, T0=0.0, newton: bool = True,
                    rtol: float = 1e-10, atol: float = 1e-12, max_newton: int = 25,
                    damping: bool = True) -> Solution:
    """Steady solve of a problem with temperature-dependent coefficients.

    ``T0`` is the initial guess (constant or ``T0(x, y)``); the Dirichlet
    values are imposed on it.  ``newton=False`` uses Picard iteration.
    """
    mesh = problem.build_mesh()
    problem.validate(mesh)
    if not problem.has_dirichlet():
        raise ValueError("A pure Neumann problem is singular; fix T on at least one boundary")
    if verbose:
        print(f"Finished generating Grid: {mesh.n_nodes} nodes, {mesh.n_elems} elements.")

    nodes, vals = dirichlet_data(mesh, problem.bc_type, problem.bc_val)
    T = values_at(T0, mesh.x, mesh.y)
    T[nodes] = vals
    F_neu = neumann_loads(mesh, np.zeros(mesh.n_nodes), problem.bc_type, problem.bc_val)
    last = {}

    def residual_jacobian(T):
        ops = assemble_operators(mesh, problem, T=T)
        last["ops"] = ops
        R = ops.A @ T - ops.F - F_neu
        return R, (ops.J if newton else ops.A)

    if verbose:
        print("Solving Equations (Newton)..." if newton else "Solving Equations (Picard)...")
    res = newton_iterate(residual_jacobian, T, nodes, method, tol, maxiter, rtol, atol,
                         max_newton, damping, verbose)
    if not res.converged:
        warnings.warn(f"nonlinear iteration did not converge in {res.iterations} steps "
                      f"(|R| = {res.residuals[-1]:.2e})", stacklevel=2)

    ops = last["ops"]
    elim = DirichletEliminator(ops.J, nodes)
    info = {"method": method, "nonlinear": "newton" if newton else "picard",
            "iterations": res.iterations, "residuals": res.residuals,
            "converged": res.converged}
    return Solution(mesh, res.T, elim.K_bc, elim.apply_rhs(ops.F + F_neu, vals),
                    problem.material, info)
