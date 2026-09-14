"""Steady Navier-Stokes (or Stokes) by Newton's method, with pseudo-transient
continuation as the robust fallback at high Reynolds number.
"""

from __future__ import annotations

import warnings

import numpy as np

from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import FlowProblem
from aa540fem.incompressible.solution import FlowSolution
from aa540fem.linalg.continuation import pseudo_transient
from aa540fem.linalg.dirichlet import eliminate
from aa540fem.linalg.newton import newton_iterate
from aa540fem.linalg.solvers import LinearSolver

CONTINUATIONS = ("auto", "newton", "ptc")


def local_pseudo_time_scaling(mesh, nu, u_ref):
    """Local pseudo-time scale of every node: ``h / u_ref`` with the longest
    edge ``h`` of the surrounding cells (``Mesh.nodal_size("max")``).

    The pseudo-time step of a node is ``dtau`` times this, so ``dtau`` is a
    CFL number and every cell advances at its own convective pace: with a
    global step the far field never moves, and with the thin dimension of
    stretched cells (the viscous scale) a boundary layer needs thousands of
    steps to convect the corrections of one cell to the next; the implicit
    steps do not need the viscous limit.  ``nu`` is accepted for interface
    stability and unused.
    """
    return mesh.nodal_size("max") / max(u_ref, 1e-12)


def project_divergence_free(asm: FlowAssembler, U, fixed, vals):
    """Closest discretely divergence-free velocity to ``U`` (in the mass norm).

    Solves the saddle-point system ``[[M, -B^T], [B, 0]] [w; q] = [M u; 0]``
    with the Dirichlet values imposed (the same projection the RK45 time
    integrator applies to its initial state) and returns ``U`` with the
    velocity replaced by ``w``; the pressure entries are kept.  The
    pseudo-transient continuation starts from this state: from a velocity
    that violates continuity, the first pseudo-time step needs a pressure
    jump of order ``1 / dtau`` to enforce it, which the stabilisation terms
    (quadratic in ``u`` and ``p``) turn into a residual that does not shrink
    with the step, so every step is rejected.
    """
    nv = asm.space.n_vel
    P = asm.pattern.matrix(asm.M_data - asm.BT_data + asm.B_data)   # [[M, -B^T], [B, 0]]
    elim = eliminate(P, fixed)
    b = np.zeros(asm.space.ndof)
    b[:nv] = (asm.M @ U)[:nv]
    sol, _ = LinearSolver(elim.K_bc, "direct", symmetric=False).solve(elim.apply_rhs(b, vals))
    V = np.array(U, dtype=float, copy=True)
    V[:nv] = sol[:nv]
    V[fixed] = vals
    return V


def solve_flow(problem: FlowProblem, U0=None, method: str = "direct", verbose: bool = False,
               rtol: float = 1e-9, atol: float = 1e-11, max_newton: int = 30,
               damping: bool = True, stokes: bool = False, continuation: str = "auto",
               dtau0: float = 1.0, max_ptc: int = 200, dtau_max: float = 1e6,
               local_timestep: bool = True) -> FlowSolution:
    """Steady Navier-Stokes (or Stokes with ``stokes=True``).

    ``continuation``: ``"newton"`` (damped Newton from ``U0``, the first step
    from rest is the Stokes solution), ``"ptc"`` (pseudo-transient
    continuation, see :func:`pseudo_transient`) or ``"auto"`` (Newton, and
    PTC from the initial state if Newton does not converge or stalls).  With
    ``local_timestep`` (default) the pseudo-time step is scaled by the local
    convective time scale (:func:`local_pseudo_time_scaling`) and ``dtau0``
    is a CFL number; otherwise it is a global time.  ``dtau_max`` caps the
    SER growth of the pseudo-time step: on stiff warm-started states the
    x10 growth overshoots the stability boundary and every overshoot costs
    a full rejected inner Newton.  Pass ``U0`` (a previous
    ``FlowSolution.U``) for continuation in Reynolds number.  ``method`` is
    ``"direct"`` (sparse LU; the saddle-point Jacobian is indefinite, so the
    Krylov methods of :class:`LinearSolver` do not apply) or
    ``"fieldsplit"`` (PETSc FGMRES with a Schur-complement block
    preconditioner, see :mod:`aa540fem.linalg.krylov`; the route to GPUs and
    to systems too large for a factorisation).
    """
    if continuation not in CONTINUATIONS:
        raise ValueError(f"Unknown continuation {continuation!r}; expected one of {CONTINUATIONS}")
    asm = FlowAssembler(problem)
    method = linear_method(method, asm)
    fixed, vals = asm.dirichlet()
    F = asm.body_load()
    U = np.zeros(asm.space.ndof) if U0 is None else np.array(U0, dtype=float, copy=True)
    U[fixed] = vals

    def residual_jacobian(U):
        if stokes:
            R = asm.K @ U - F - asm.BT @ U + asm.B @ U
            J = asm.pattern.matrix(asm.K_data - asm.BT_data + asm.B_data)
        else:
            R, J = asm.steady_residual_jacobian(U, F)
        return R, J

    def residual(U):
        if stokes:
            return asm.K @ U - F - asm.BT @ U + asm.B @ U
        return asm.steady_residual(U, F)

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
                             stall_iterations=5 if continuation == "auto" else None,
                             residual=residual)
    if not stokes and (continuation == "ptc" or (continuation == "auto" and not res.converged)):
        if verbose and res is not None:
            print("  Newton did not converge; switching to pseudo-transient continuation")
        path = "ptc" if res is None else "newton+ptc"
        if res is not None and res.residuals[-1] < 0.5 * res.residuals[0]:
            U = res.T       # continue from Newton's best iterate if it got somewhere;
            #                 heavily damped steps that barely reduced the residual
            #                 are a worse start than the initial state
        U = project_divergence_free(asm, U, fixed, vals)
        M = asm.M
        if local_timestep:
            vel = fixed < asm.space.n_vel
            u_ref = max(float(np.abs(vals[vel]).max()) if vel.any() else 0.0, 1e-3)
            scale = local_pseudo_time_scaling(asm.mesh, problem.mu / problem.rho, u_ref)
            inv = np.concatenate([1.0 / scale] * asm.space.dim + [np.ones(asm.space.Np)])
            M = asm.pattern.matrix(asm.M_data * inv[asm.pattern.rows])   # row scaling
        res = pseudo_transient(residual_jacobian, U, fixed, M, method, rtol, atol, dtau0,
                               max_ptc, dtau_max=dtau_max, verbose=verbose, residual=residual)
    if not res.converged:
        warnings.warn(f"steady solve did not converge in {res.iterations} iterations "
                      f"(|R| = {res.residuals[-1]:.2e})", stacklevel=2)
    info = {"iterations": res.iterations, "residuals": res.residuals,
            "converged": res.converged, "stokes": stokes, "continuation": path,
            "stabilisation": problem.stabilisation}
    if path != "newton":
        info["dtau"] = res.steps
    return FlowSolution(problem, asm.space, res.T, info, asm)


def linear_method(method, asm):
    """Resolve the ``method`` argument of the flow solvers.

    ``"direct"`` stays the sparse LU unless the run configuration asks for the
    PETSc Krylov path (``--linear petsc`` / ``petsc-cuda``,
    ``AA540FEM_LINEAR``); ``"fieldsplit"`` is that path explicitly; a
    callable ``A -> solver`` is passed through.
    """
    if callable(method):
        return method
    if method == "direct":
        try:
            from aa540fem.hardware import get_config

            backend = get_config().linear_backend
        except ImportError:
            backend = "scipy"
        if backend in ("petsc", "petsc-cuda"):
            return fieldsplit_factory(asm, gpu=backend == "petsc-cuda")
        return method
    if method == "fieldsplit":
        return fieldsplit_factory(asm)
    return method


def fieldsplit_factory(asm, rtol=1e-6, velocity_pc="ilu", gpu=None, schur="lsc", options=None):
    """``method`` callable for the Newton/PTC loops: the PETSc fieldsplit solver
    on the eliminated Jacobian (:class:`aa540fem.linalg.krylov.FieldSplitSolver`).

    ``schur``: ``"lsc"`` (least-squares commutator, default: robust for the
    convective and stabilised Jacobians) or ``"mass"`` (viscosity-scaled
    pressure mass matrix, adequate for Stokes-like problems);
    ``velocity_pc``: ``"ilu"`` or ``"gamg"``; ``rtol``: relative tolerance of
    the linear solves (the Newton residual is always evaluated exactly, so a
    moderate value only costs Newton iterations); ``options``: extra PETSc
    options (see the solver class).
    """
    from aa540fem.linalg.krylov import FieldSplitSolver, pressure_mass_matrix

    if gpu is None:
        try:
            from aa540fem.hardware import get_config

            gpu = get_config().linear_backend == "petsc-cuda"
        except ImportError:
            gpu = False
    Mp = pressure_mass_matrix(asm)
    n_vel = asm.space.n_vel
    Qdiag = np.asarray(abs(asm.M).sum(axis=1)).ravel()[:n_vel]     # lumped velocity mass

    class _Adapter:
        def __init__(self, A):
            self.solver = FieldSplitSolver(A, n_vel, Mp, rtol=rtol, velocity_pc=velocity_pc,
                                           gpu=gpu, schur=schur, Qdiag=Qdiag, options=options)

        def solve(self, b):
            x = self.solver.solve(b)
            return x, {"method": "fieldsplit", "iterations": self.solver.iterations,
                       "residual": self.solver.residual}

    return _Adapter
