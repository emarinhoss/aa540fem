"""Time integration of the Navier-Stokes system: adaptive RK45 (default,
projected explicit Runge-Kutta) or the theta-method with Newton per step.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse.linalg as spla

from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import SCHEMES, FlowProblem, values_at_pair
from aa540fem.incompressible.solution import FlowSolution, TransientFlowSolution
from aa540fem.linalg.dirichlet import eliminate
from aa540fem.linalg.newton import newton_iterate
from aa540fem.timestepping.rk import rk45


def _initial_state(asm, U0, fixed, vals):
    space, mesh = asm.space, asm.mesh
    if U0 is None:
        U = np.zeros(space.ndof)
    elif callable(U0):
        ux, uy = values_at_pair(U0, mesh.x, mesh.y)
        U = np.concatenate([np.broadcast_to(ux, mesh.x.shape), np.broadcast_to(uy, mesh.x.shape),
                            np.zeros(space.Np)]).astype(float)
    else:
        U = np.array(U0, dtype=float, copy=True)
    U[fixed] = vals
    return U


def solve_flow_transient(problem: FlowProblem, dt: float, t_end: float, theta: float = 0.5,
                         U0=None, method: str = "direct", store_every: int = 1,
                         verbose: bool = False, rtol: float = 1e-8, atol: float = 1e-10,
                         max_newton: int = 25, damping: bool = True, callback=None,
                         startup_steps: int = 0, frozen_jacobian: bool = True,
                         scheme: str = "rk45", output_interval: float | None = None,
                         dt_max: float | None = None,
                         adaptive: bool = True) -> TransientFlowSolution:
    """Time integration of the Navier-Stokes system.

    ``scheme="rk45"`` (default): explicit Dormand-Prince 5(4) on the velocity
    with a pressure projection at every stage (see
    :func:`_solve_flow_transient_rk45`); ``dt`` is the initial step, the
    error control uses ``rtol``/``atol`` (defaults 1e-8 / 1e-10 are the
    Newton tolerances of the theta scheme; the RK path uses 1e-4 / 1e-6
    unless you pass them explicitly), ``dt_max`` bounds the step,
    ``adaptive=False`` gives the fixed-step 5th-order method, and
    ``output_interval`` stores the fields exactly at multiples of that time.

    ``scheme="theta"``: theta-method with a Newton solve per step.
    Momentum: ``M (U - U_n)/dt + theta S(U) + (1 - theta) S(U_n) - B^T p = 0``
    with ``S(U) = K U + N(U) - F``; the pressure and the continuity equation
    are implicit.  ``startup_steps`` backward-Euler steps are taken first
    (Rannacher start-up), which damps the Crank-Nicolson ringing after an
    impulsive start.  With ``frozen_jacobian`` each step factorises its
    Jacobian once and iterates with it (modified Newton).

    ``U0`` may be a velocity/pressure vector or a callable
    ``(x, y) -> (ux, uy)`` for the initial velocity.  ``callback(step, t,
    solution)`` is called after every (accepted) step with a
    :class:`FlowSolution` (e.g. to log forces).
    """
    if scheme not in SCHEMES:
        raise ValueError(f"Unknown scheme {scheme!r}; expected one of {SCHEMES}")
    if scheme == "rk45":
        return _solve_flow_transient_rk45(problem, dt, t_end, U0, store_every, verbose,
                                          callback, rtol, atol, dt_max, adaptive,
                                          output_interval)
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must be in [0, 1]")
    nsteps = int(round(t_end / dt))
    if nsteps < 1 or abs(nsteps * dt - t_end) > 1e-8 * max(1.0, abs(t_end)):
        raise ValueError(f"dt = {dt} must divide t_end = {t_end}")

    asm = FlowAssembler(problem)
    space = asm.space
    B, BT, M = asm.B, asm.BT, asm.M
    time_dependent = problem.depends_on_time()

    fixed, vals = asm.dirichlet(0.0)
    U = _initial_state(asm, U0, fixed, vals)

    F_old = asm.body_load(0.0)
    S_old = asm.K @ U + asm.convection(U, jacobian=False) - F_old
    S_old[2 * space.N:] = 0.0

    times = [0.0]
    snapshots = [U.copy()]
    newton_iterations = []
    if verbose:
        print(f"Taylor-Hood: {space.ndof} unknowns; {nsteps} steps of dt = {dt} (theta = {theta})")

    for n in range(1, nsteps + 1):
        t = n * dt
        th = 1.0 if n <= startup_steps else theta
        if time_dependent:
            fixed, vals = asm.dirichlet(t)
        F_new = asm.body_load(t) if time_dependent else F_old
        U_old = U

        def evaluate(Un, jacobian, U_old=U_old, F_new=F_new, S_old=S_old, th=th, t=t):
            mt = asm.momentum_terms(Un, t=t, dt=dt, U_old=U_old, jacobian=jacobian)
            S = asm.K @ Un + mt.N - F_new
            S[2 * space.N:] = 0.0
            R = (M @ ((Un - U_old) / dt) + th * S + (1 - th) * S_old - BT @ Un + B @ Un
                 + mt.S)
            if not jacobian:
                return R
            J = asm.pattern.matrix(asm.M_data / dt + th * (asm.K_data + mt.JN_data)
                                   + mt.JS_data - asm.BT_data + asm.B_data)
            return R, J

        guess = U_old.copy()
        guess[fixed] = vals
        res = newton_iterate(lambda Un: evaluate(Un, True), guess, fixed, method, rtol=rtol,
                             atol=atol, max_newton=max_newton, damping=damping,
                             frozen_jacobian=frozen_jacobian,
                             residual=lambda Un: evaluate(Un, False))
        if not res.converged:
            raise RuntimeError(f"Newton did not converge at t = {t:.6g} "
                               f"(|R| = {res.residuals[-1]:.2e})")
        U = res.T
        S_old = asm.K @ U + asm.convection(U, jacobian=False) - F_new
        S_old[2 * space.N:] = 0.0
        F_old = F_new
        newton_iterations.append(res.iterations)
        if n % store_every == 0 or n == nsteps:
            times.append(t)
            snapshots.append(U.copy())
        if callback is not None:
            callback(n, t, FlowSolution(problem, space, U, {"step": n, "t": t}, asm))
        if verbose and (n % max(1, nsteps // 10) == 0 or n == nsteps):
            u, v, _ = space.split(U)
            print(f"  step {n}/{nsteps}, t = {t:.6g}, max |u| = {np.hypot(u, v).max():.6g}, "
                  f"{res.iterations} Newton iterations")

    info = {"scheme": "theta", "steps": nsteps, "dt": dt, "theta": theta,
            "startup_steps": startup_steps, "newton_iterations": newton_iterations}
    return TransientFlowSolution(problem, space, np.asarray(times), snapshots, info, asm)


def _solve_flow_transient_rk45(problem, dt, t_end, U0, store_every, verbose, callback,
                               rtol, atol, dt_max, adaptive, output_interval):
    """Projected explicit Runge-Kutta for the incompressible system.

    The state is the velocity.  Every stage derivative solves the constant
    saddle-point system ``[[M, -B^T], [B, 0]] [k; p] = [-(K u + N(u) - F); 0]``
    (Dirichlet velocity dofs carry the time derivative of their values), so
    every ``k`` is discretely divergence-free and so is every linear
    combination of stages: the continuity constraint is preserved exactly.
    The matrix is factorised once for the whole run.  The initial velocity
    is projected onto the divergence-free space in the same way.  The
    pressure stored with each snapshot is the stage pressure at that time
    (the FSAL stage of the accepted step).  With ``stabilisation`` the SUPG
    and grad-div terms use the pressure of the previous projection and the
    steady ``tau`` (PSPG is not applied: the projection keeps the continuity
    rows).
    """
    asm = FlowAssembler(problem)
    space = asm.space
    N = space.N
    P = asm.pattern.matrix(asm.M_data - asm.BT_data + asm.B_data)   # [[M, -B^T], [B, 0]]
    time_dependent = problem.depends_on_time()
    if rtol == 1e-8 and atol == 1e-10:          # theta-scheme Newton defaults: use RK defaults
        rtol, atol = 1e-4, 1e-6

    fixed, vals0 = asm.dirichlet(0.0)
    is_pressure = fixed >= 2 * N
    fixed_vel = fixed[~is_pressure]
    elim = eliminate(P, fixed)
    lu = spla.splu(elim.K_bc.tocsc())
    zero_rate = np.zeros(fixed.size)

    def fixed_values(t):
        if not time_dependent:
            return vals0, zero_rate
        v = asm.dirichlet(t)[1]
        r = asm.dirichlet(t, rate=True)[1]
        return v, r

    def project(b_vel, t, pin_rows_zero=False):
        """Solve the saddle-point system for ``(k, p)`` given the momentum rhs."""
        vals, rate = fixed_values(t)
        b = np.zeros(space.ndof)
        b[:2 * N] = b_vel[:2 * N]
        given = np.where(is_pressure, vals, rate)
        b = elim.apply_rhs(b, given)
        sol = lu.solve(b)
        return sol[:2 * N], sol[2 * N:]

    # initial velocity, projected onto the divergence-free space
    U = _initial_state(asm, U0, fixed, vals0)
    b = np.zeros(space.ndof)
    b[:2 * N] = (asm.M @ U)[:2 * N]
    sol = lu.solve(elim.apply_rhs(b, vals0))
    u = sol[:2 * N]
    F_const = asm.body_load(0.0)
    last = {"p": np.zeros(space.Np)}

    stabilised = problem.stabilisation

    def rhs(t, y):
        y = np.array(y, copy=True)
        vals, _ = fixed_values(t)
        y[fixed_vel] = vals[~is_pressure]
        U = np.concatenate([y, last["p"]])          # lagged pressure for the stabilisation
        F = asm.body_load(t) if time_dependent else F_const
        if stabilised:
            mt = asm.momentum_terms(U, t=t, pspg=False, jacobian=False)
            r = -(asm.K @ U + mt.N + mt.S - F)
        else:
            r = -(asm.K @ U + asm.convection(U, jacobian=False) - F)
        k, p = project(r, t)
        last["p"] = p
        return k

    def store_transform(t, y):
        return np.concatenate([y, last["p"]])

    def on_accept(step, t, y, h):
        if callback is not None:
            callback(step, t, FlowSolution(problem, space, np.concatenate([y, last["p"]]),
                                           {"step": step, "t": t, "dt": h}, asm))

    if verbose:
        print(f"Taylor-Hood: {space.ndof} unknowns; RK45 (Dormand-Prince), initial dt = {dt}, "
              f"rtol = {rtol}, atol = {atol}")
    free_vel = np.ones(2 * N, dtype=bool)
    free_vel[fixed_vel] = False
    res = rk45(rhs, u, 0.0, t_end, dt, rtol=rtol, atol=atol, dt_max=dt_max, adaptive=adaptive,
               output_interval=output_interval, store_every=store_every, error_mask=free_vel,
               on_accept=on_accept, store_transform=store_transform, verbose=verbose)
    # the stored initial state has no pressure yet: recover it from the first stage
    rhs(0.0, res.states[0][:2 * N])
    res.states[0][2 * N:] = last["p"]
    if verbose:
        print(f"  {res.info['steps']} steps ({res.info['rejected']} rejected), "
              f"dt in [{res.info['dt_min_used']:.3e}, {res.info['dt_max_used']:.3e}]")
    return TransientFlowSolution(problem, space, res.times, res.states, res.info, asm)
