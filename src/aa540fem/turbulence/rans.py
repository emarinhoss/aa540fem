"""Segregated coupling of the flow solver with the Spalart-Allmaras model.

Each outer iteration solves the SA equation for the frozen velocity field,
updates the eddy viscosity ``mu_t = rho nu_t`` in the momentum equations
(under-relaxed) and re-solves the steady flow from the previous state.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from aa540fem.incompressible.problem import OPEN, FlowProblem
from aa540fem.incompressible.solution import FlowSolution
from aa540fem.incompressible.steady import solve_flow
from aa540fem.turbulence.spalart_allmaras import SpalartAllmaras, SpalartAllmarasSolver


@dataclass
class RANSSolution:
    flow: FlowSolution
    nu_tilde: np.ndarray
    nu_t: np.ndarray
    distance: np.ndarray
    model: SpalartAllmaras
    history: list = field(default_factory=list)   # (outer, change of nu_t, wall force)
    converged: bool = False

    @property
    def mesh(self):
        return self.flow.mesh

    def save(self, path):
        """ParaView file with velocity, pressure, nu_tilde, nu_t/nu and wall distance."""
        from aa540fem.io.mesh_files import write_vtk

        f = self.flow
        return write_vtk(path, self.mesh, point_data={
            "velocity": f.velocity_3d(), "speed": f.speed, "p": f.p_nodal,
            "nu_tilde": self.nu_tilde,
            "nu_t_ratio": self.nu_t / self.model.nu, "wall_distance": self.distance})


def solve_rans(problem: FlowProblem, wall_tags, nu_tilde_inf=None, model=None,
               max_outer: int = 40, tol: float = 1e-3, relax: float = 0.7, verbose=False,
               force_tag=None, U0=None, viscosity_ramp=(1.0,),
               flow_options=None) -> RANSSolution:
    """Steady RANS solution with the Spalart-Allmaras model.

    Parameters
    ----------
    problem     : :class:`FlowProblem` (``stabilisation=True`` recommended);
                  ``mu``/``rho`` give the laminar viscosity.
    wall_tags   : boundary tags of no-slip walls (``nu_tilde = 0``, wall distance).
    nu_tilde_inf : freestream working variable, default ``3 nu`` (NASA TMR
                  recommendation, ``nu_t/nu`` about 0.21); imposed on every
                  boundary with a full velocity Dirichlet condition except the
                  walls.
    tol         : outer-iteration stop on the relative change of ``nu_t``.
    relax       : under-relaxation of ``nu_tilde`` between outer iterations.
    force_tag   : boundary whose force is reported in the history (default:
                  the first wall tag).
    viscosity_ramp : factors applied to the laminar viscosity in turn, e.g.
                  ``(100, 10, 1)``: the coupled problem is first converged
                  (loosely) at the higher viscosities and each stage starts
                  from the previous one, which is far more robust than
                  starting the wall-resolved high-Reynolds-number case from
                  rest.  The last factor must be 1.
    U0          : initial state for the first flow solve (a dof vector or a
                  callable ``(x, y[, z]) -> (ux, uy[, uz])``); a smooth boundary-layer
                  profile is a much better start than uniform flow on
                  wall-resolved meshes.
    flow_options : extra keyword arguments for :func:`solve_flow`
                  (e.g. ``dtau0`` of the pseudo-transient continuation).
    verbose     : print one line per outer iteration; ``verbose=2`` also
                  prints the Newton / pseudo-time history of the sub-solves.
    """
    wall_tags = [wall_tags] if isinstance(wall_tags, str) else list(wall_tags)
    if viscosity_ramp[-1] != 1.0:
        raise ValueError("viscosity_ramp must end with the factor 1")
    if len(viscosity_ramp) > 1:
        mu_target = problem.mu
        result = None
        for k, factor in enumerate(viscosity_ramp):
            problem.mu = mu_target * factor
            last = k == len(viscosity_ramp) - 1
            if verbose:
                print(f"RANS viscosity ramp: factor {factor:g}")
            result = _solve_rans(problem, wall_tags, nu_tilde_inf, None,
                                 max_outer if last else 6, tol if last else 10 * tol, relax,
                                 verbose, force_tag, U0 if result is None else result.flow.U,
                                 None if result is None else result.nu_tilde, flow_options)
        problem.mu = mu_target
        return result
    return _solve_rans(problem, wall_tags, nu_tilde_inf, model, max_outer, tol, relax, verbose,
                       force_tag, U0, None, flow_options)


def _initial_vector(problem, U0):
    if U0 is None or not callable(U0):
        return U0
    from aa540fem.incompressible.assembler import FlowAssembler
    from aa540fem.incompressible.transient import _initial_state

    asm = FlowAssembler(problem)
    fixed, vals = asm.dirichlet()
    return _initial_state(asm, U0, fixed, vals)


def _solve_rans(problem, wall_tags, nu_tilde_inf, model, max_outer, tol, relax, verbose,
                force_tag, U0, nu_tilde0, flow_options=None):
    flow_options = {"rtol": 1e-5, "atol": 1e-12, **(flow_options or {})}
    U0 = _initial_vector(problem, U0)
    mesh = problem.mesh
    nu = problem.mu / problem.rho
    model = model or SpalartAllmaras(nu)
    nt_inf = 3.0 * nu if nu_tilde_inf is None else nu_tilde_inf
    force_tag = force_tag or wall_tags[0]

    sa = SpalartAllmarasSolver(mesh, model, wall_tags, order=problem.order,
                               element_length=problem.element_length)
    fixed = {}
    for tag, spec in problem.bc.items():
        if spec == OPEN or tag in wall_tags or any(v is None for v in spec):
            continue
        fixed.update({int(n): nt_inf for n in mesh.bc_nodes[tag]})
    for tag in wall_tags:
        fixed.update({int(n): 0.0 for n in mesh.bc_nodes[tag]})
    fixed_nodes = np.array(sorted(fixed), dtype=int)
    fixed_vals = np.array([fixed[n] for n in fixed_nodes])

    nt = np.full(mesh.n_nodes, nt_inf) if nu_tilde0 is None else np.array(nu_tilde0, dtype=float)
    nt[fixed_nodes] = fixed_vals
    nu_t = model.eddy_viscosity(nt)
    problem.eddy_viscosity = problem.rho * nu_t
    inner = bool(verbose) and int(verbose) >= 2
    flow = solve_flow(problem, U0=U0, continuation="auto", verbose=inner, **flow_options)
    if verbose:
        print(f"RANS start: flow {flow.info['continuation']} in "
              f"{flow.info['iterations']} iterations")

    history = []
    converged = False
    for k in range(1, max_outer + 1):
        sa.set_velocity(*flow.velocity)
        res = sa.solve(nt, fixed_nodes, fixed_vals, rtol=1e-4, verbose=inner)
        nt_new = res.T
        nt = relax * nt_new + (1.0 - relax) * nt
        nu_t_new = model.eddy_viscosity(nt)
        change = np.linalg.norm(nu_t_new - nu_t) / max(np.linalg.norm(nu_t_new), 1e-300)
        nu_t = nu_t_new
        problem.eddy_viscosity = problem.rho * nu_t
        flow = solve_flow(problem, U0=flow.U, continuation="auto", verbose=inner, **flow_options)
        force = flow.forces(force_tag)
        history.append((k, change, force, res.iterations, res.converged, flow.info["iterations"]))
        if verbose:
            print(f"  outer {k}: |d nu_t| = {change:.2e}, SA {res.iterations} steps "
                  f"({'ok' if res.converged else 'NOT converged'}), flow {flow.info['iterations']} "
                  f"iterations, force on {force_tag} = ({force[0]:.5f}, {force[1]:.5f})")
        if change < tol and res.converged and flow.info["converged"]:
            converged = True
            break
    return RANSSolution(flow, nt, nu_t, sa.distance, model, history, converged)
