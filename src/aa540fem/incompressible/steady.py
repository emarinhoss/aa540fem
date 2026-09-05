"""Steady Navier-Stokes (or Stokes) by Newton's method, with pseudo-transient
continuation as the robust fallback at high Reynolds number.
"""

from __future__ import annotations

import warnings

import numpy as np
import scipy.sparse as sp

from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import FlowProblem
from aa540fem.incompressible.solution import FlowSolution
from aa540fem.linalg.continuation import pseudo_transient
from aa540fem.linalg.newton import newton_iterate

CONTINUATIONS = ("auto", "newton", "ptc")


def local_pseudo_time_scaling(mesh, nu, u_ref=1.0):
    """Nodal pseudo-time scale ``h / (u_ref + nu / h)`` with ``h`` the local
    mesh size: the local convective-diffusive time scale of each cell.  With
    it a pseudo-time step is a CFL-like number rather than a time, so the
    thin wall cells of a boundary-layer mesh and the coarse far field advance
    at their own pace (standard local time stepping)."""
    h = mesh.nodal_size()
    return h / (u_ref + nu / h)


def solve_flow(problem: FlowProblem, U0=None, method: str = "direct", verbose: bool = False,
               rtol: float = 1e-9, atol: float = 1e-11, max_newton: int = 30,
               damping: bool = True, stokes: bool = False, continuation: str = "auto",
               dtau0: float = 1.0, max_ptc: int = 200,
               local_timestep: bool = True) -> FlowSolution:
    """Steady Navier-Stokes (or Stokes with ``stokes=True``).

    ``continuation``: ``"newton"`` (damped Newton from ``U0``, the first step
    from rest is the Stokes solution), ``"ptc"`` (pseudo-transient
    continuation, see :func:`pseudo_transient`) or ``"auto"`` (Newton, and
    PTC from the initial state if Newton does not converge or stalls).  With
    ``local_timestep`` (default) the pseudo-time step is scaled by the local
    cell time scale (:func:`local_pseudo_time_scaling`) and ``dtau0`` is a
    CFL-like number; otherwise it is a global time.  Pass ``U0`` (a previous
    ``FlowSolution.U``) for continuation in Reynolds number.  ``method``
    should be ``"direct"``: the saddle-point Jacobian is indefinite.
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
                             abort_ratio=100.0 if continuation == "auto" else None,
                             stall_iterations=5 if continuation == "auto" else None)
    if not stokes and (continuation == "ptc" or (continuation == "auto" and not res.converged)):
        if verbose and res is not None:
            print("  Newton did not converge; switching to pseudo-transient continuation")
        path = "ptc" if res is None else "newton+ptc"
        if res is not None and res.residuals[-1] < res.residuals[0]:
            U = res.T                           # continue from Newton's best iterate
        M = asm.M
        if local_timestep:
            vel = fixed < 2 * asm.space.N
            u_ref = max(float(np.abs(vals[vel]).max()) if vel.any() else 0.0, 1e-3)
            scale = local_pseudo_time_scaling(asm.mesh, problem.mu / problem.rho, u_ref)
            inv = np.concatenate([1.0 / scale, 1.0 / scale, np.ones(asm.space.Np)])
            M = (sp.diags(inv) @ asm.M).tocsr()
        res = pseudo_transient(residual_jacobian, U, fixed, M, method, rtol, atol, dtau0,
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
