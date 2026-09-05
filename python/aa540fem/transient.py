"""Time integration of ``rho_c dT/dt + u . grad T - div(kappa grad T) = f``
with the theta-method (backward Euler for ``theta = 1``, Crank-Nicolson for
``theta = 0.5``).
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass, field

import numpy as np

from .boundary import DirichletEliminator
from .geometry import Mesh
from .solver import LinearSolver, Problem, assemble_operators, dirichlet_data, neumann_loads
from .util import values_at


@dataclass
class TransientSolution:
    mesh: Mesh
    times: np.ndarray            # stored times
    snapshots: list              # nodal fields at ``times``
    problem: Problem
    info: dict = field(default_factory=dict)

    @property
    def T(self) -> np.ndarray:
        """Final nodal field."""
        return self.snapshots[-1]

    def snapshot(self, t: float) -> np.ndarray:
        """Stored field closest to time ``t``."""
        return self.snapshots[int(np.argmin(np.abs(self.times - t)))]

    def save(self, path):
        """Write the final field to a ParaView file."""
        from .mesh_io import write_vtk

        return write_vtk(path, self.mesh, point_data={"T": self.T})

    def save_series(self, prefix):
        """Write one ``.vtu`` per stored step and a ``.pvd`` collection.

        ``prefix`` may include a directory; the files are
        ``<prefix>_0000.vtu`` ... and ``<prefix>.pvd`` (open the latter in
        ParaView to animate).  Returns the path of the ``.pvd`` file.
        """
        from .mesh_io import write_vtk

        prefix = pathlib.Path(prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        entries = []
        for i, (t, T) in enumerate(zip(self.times, self.snapshots)):
            name = f"{prefix.name}_{i:04d}.vtu"
            write_vtk(prefix.parent / name, self.mesh, point_data={"T": T})
            entries.append(
                f'    <DataSet timestep="{float(t)!r}" group="" part="0" file="{name}"/>')
        pvd = prefix.with_suffix(".pvd")
        pvd.write_text(
            '<?xml version="1.0"?>\n'
            '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
            "  <Collection>\n" + "\n".join(entries) + "\n  </Collection>\n</VTKFile>\n")
        return pvd


def solve_transient(problem: Problem, dt: float, t_end: float, theta: float = 1.0,
                    T0=0.0, method: str = "direct", tol: float = 1e-10, maxiter=None,
                    store_every: int = 1, verbose: bool = False) -> TransientSolution:
    """Integrate the transient problem from ``t = 0`` to ``t_end``.

    Parameters
    ----------
    problem     : :class:`Problem`; ``rho_c``, ``material``, ``velocity`` and
                  the boundary values may take a third argument ``t``.
    dt, t_end   : time step (must divide ``t_end``) and final time.
    theta       : 1 backward Euler (first order, unconditionally stable),
                  0.5 Crank-Nicolson (second order).
    T0          : initial field, constant or callable ``T0(x, y)``.
    method      : linear solver, see :class:`aa540fem.solver.LinearSolver`.
                  With time-independent coefficients the matrix is factorised
                  (or its preconditioner built) once.
    store_every : keep every n-th step (the initial and final fields are
                  always stored).
    """
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must be in [0, 1]")
    nsteps = int(round(t_end / dt))
    if nsteps < 1 or abs(nsteps * dt - t_end) > 1e-8 * max(1.0, abs(t_end)):
        raise ValueError(f"dt = {dt} must divide t_end = {t_end}")

    mesh = problem.build_mesh()
    problem.validate(mesh)
    if verbose:
        print(f"Finished generating Grid: {mesh.n_nodes} nodes, {mesh.n_elems} elements.")

    ops_time = problem.operators_depend_on_time()
    loads_time = problem.loads_depend_on_time()

    def system(t):
        ops = assemble_operators(mesh, problem, t=t, dt=dt)
        A = ops.A
        lhs = (ops.M / dt + theta * A).tocsr()
        rhs = (ops.M / dt - (1.0 - theta) * A).tocsr()
        return ops, lhs, rhs

    def load(ops, t):
        return neumann_loads(mesh, ops.F, problem.bc_type, problem.bc_val, t)

    ops, lhs, rhs = system(0.0)
    nodes, vals = dirichlet_data(mesh, problem.bc_type, problem.bc_val, 0.0)
    elim = DirichletEliminator(lhs, nodes)
    solver = LinearSolver(elim.K_bc, method, tol, maxiter, symmetric=problem.symmetric)

    T = values_at(T0, mesh.x, mesh.y)
    T[nodes] = vals
    F_old = load(ops, 0.0)

    times = [0.0]
    snapshots = [T.copy()]
    iterations = []
    if verbose:
        print(f"Time stepping: {nsteps} steps of dt = {dt} (theta = {theta})")

    for n in range(1, nsteps + 1):
        t = n * dt
        if ops_time:
            ops_new, lhs, rhs_new = system(t)
            elim = DirichletEliminator(lhs, nodes)
            solver = LinearSolver(elim.K_bc, method, tol, maxiter, symmetric=problem.symmetric)
            F_new = load(ops_new, t)
            # theta-method with time-varying operators: A^n on the explicit part
            b = rhs @ T + theta * F_new + (1.0 - theta) * F_old
            rhs = rhs_new
        else:
            F_new = load(ops, t) if loads_time else F_old
            b = rhs @ T + theta * F_new + (1.0 - theta) * F_old
        if loads_time:
            _, vals = dirichlet_data(mesh, problem.bc_type, problem.bc_val, t)
        b = elim.apply_rhs(b, vals)
        T, info = solver.solve(b)
        if "iterations" in info:
            iterations.append(info["iterations"])
        F_old = F_new
        if n % store_every == 0 or n == nsteps:
            times.append(t)
            snapshots.append(T.copy())
        if verbose and (n % max(1, nsteps // 10) == 0 or n == nsteps):
            print(f"  step {n}/{nsteps}, t = {t:.6g}, T in [{T.min():.6g}, {T.max():.6g}]")

    info = {"method": method, "steps": nsteps, "dt": dt, "theta": theta}
    if iterations:
        info["iterations"] = iterations
    return TransientSolution(mesh, np.asarray(times), snapshots, problem, info)
