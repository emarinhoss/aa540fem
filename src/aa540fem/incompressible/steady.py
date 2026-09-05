"""Steady Navier-Stokes (or Stokes) by Newton's method."""

from __future__ import annotations

import warnings

import numpy as np

from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import FlowProblem
from aa540fem.incompressible.solution import FlowSolution
from aa540fem.linalg.newton import newton_iterate


def solve_flow(problem: FlowProblem, U0=None, method: str = "direct", verbose: bool = False,
               rtol: float = 1e-9, atol: float = 1e-11, max_newton: int = 30,
               damping: bool = True, stokes: bool = False) -> FlowSolution:
    """Steady Navier-Stokes (or Stokes with ``stokes=True``) by Newton's method.

    The first Newton step from ``U0 = 0`` is the Stokes solution, which is
    the usual starting point; pass ``U0`` (a previous ``FlowSolution.U``)
    for continuation in Reynolds number.  ``method`` should be ``"direct"``:
    the saddle-point Jacobian is indefinite.
    """
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
              f"{asm.space.ndof} unknowns")
    res = newton_iterate(residual_jacobian, U, fixed, method, rtol=rtol, atol=atol,
                         max_newton=max_newton, damping=damping, verbose=verbose)
    if not res.converged:
        warnings.warn(f"Newton did not converge in {res.iterations} iterations "
                      f"(|R| = {res.residuals[-1]:.2e})", stacklevel=2)
    info = {"iterations": res.iterations, "residuals": res.residuals,
            "converged": res.converged, "stokes": stokes}
    return FlowSolution(problem, asm.space, res.T, info, asm)
