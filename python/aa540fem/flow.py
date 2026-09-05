"""Incompressible Navier-Stokes with Taylor-Hood elements.

    rho ( du/dt + (u . grad) u ) - mu lap(u) + grad p = rho f
    div u = 0

Velocity is interpolated with the quadratic elements (``quad9`` or
``triangle6``), pressure with the linear element on their corner nodes
(Q2/Q1 or P2/P1), which is inf-sup stable, so no pressure stabilisation is
needed.  The viscous term uses the Laplacian form, whose natural boundary
condition ``mu du/dn - p n = 0`` is the usual "do-nothing" outflow.  The
nonlinear system is solved with Newton's method
(:func:`aa540fem.nonlinear.newton_iterate`) and time is integrated with the
theta-method, pressure and continuity being treated implicitly.

Unknown ordering: ``[u_x (all nodes), u_y (all nodes), p (corner nodes)]``.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from .element import physical_gradients
from .elements import PRESSURE_ELEMENT, get_element
from .geometry import Mesh
from .nonlinear import newton_iterate
from .quadrature import gauss_legendre_quad
from .util import accepts_time, values_at

OPEN = "open"


# ------------------------------------------------------------------ problem
@dataclass
class FlowProblem:
    """Definition of an incompressible flow case.

    Attributes
    ----------
    mesh        : :class:`Mesh` with ``quad9`` and/or ``triangle6`` blocks.
    mu, rho     : dynamic viscosity and density (constants).
    body_force  : optional ``(x, y[, t]) -> (fx, fy)`` acceleration.
    bc          : dict boundary tag -> ``(ux, uy)`` Dirichlet values (each a
                  constant, a callable ``(x, y[, t])`` or ``None`` for a free
                  component, e.g. ``(None, 0.0)`` for a symmetry line) or
                  ``"open"`` for the do-nothing outflow.  Unlisted tags are
                  open.  Where tags share a node the last listed wins.
    pin_pressure : fix the pressure at one node.  ``None`` pins it
                  automatically when no boundary is open (enclosed flow).
    pin_value   : value (constant or ``(x, y)``) of the pinned pressure.
    order       : quadrature order (``None``: element default).
    """

    mesh: Mesh
    mu: float = 1.0
    rho: float = 1.0
    body_force: object = None
    bc: dict = field(default_factory=dict)
    pin_pressure: bool | None = None
    pin_value: object = 0.0
    order: int | None = None

    def validate(self):
        for name in self.mesh.cells:
            if name not in PRESSURE_ELEMENT:
                raise ValueError(f"Taylor-Hood needs quadratic elements; got {name!r} "
                                 f"(supported: {sorted(PRESSURE_ELEMENT)})")
        unknown = [t for t in self.bc if t not in self.mesh.boundary]
        if unknown:
            raise ValueError(f"Unknown boundary tag(s) {unknown}; mesh has {self.mesh.tags}")
        for tag, spec in self.bc.items():
            if spec == OPEN:
                continue
            if not (isinstance(spec, (tuple, list)) and len(spec) == 2):
                raise ValueError(f"bc[{tag!r}] must be (ux, uy) or 'open'")

    @property
    def has_open_boundary(self) -> bool:
        return any(spec == OPEN for spec in self.bc.values()) or any(
            tag not in self.bc for tag in self.mesh.boundary)

    @property
    def pins_pressure(self) -> bool:
        return (not self.has_open_boundary) if self.pin_pressure is None else self.pin_pressure

    def depends_on_time(self) -> bool:
        if accepts_time(self.body_force):
            return True
        return any(accepts_time(v) for spec in self.bc.values() if spec != OPEN for v in spec)


# -------------------------------------------------------------------- space
class TaylorHoodSpace:
    """Degree-of-freedom layout ``[u_x, u_y, p]`` for a Taylor-Hood pair."""

    def __init__(self, mesh: Mesh):
        self.mesh = mesh
        self.N = mesh.n_nodes
        corners = np.concatenate([conn[:, :get_element(name).n_corners].ravel()
                                  for name, conn in mesh.cells.items()])
        self.pressure_nodes = np.unique(corners)
        self.Np = self.pressure_nodes.size
        self.p_index = -np.ones(self.N, dtype=int)
        self.p_index[self.pressure_nodes] = np.arange(self.Np)
        self.ndof = 2 * self.N + self.Np

    def dof_ux(self, nodes):
        return np.asarray(nodes, dtype=int)

    def dof_uy(self, nodes):
        return self.N + np.asarray(nodes, dtype=int)

    def dof_p(self, nodes):
        idx = self.p_index[np.asarray(nodes, dtype=int)]
        if (idx < 0).any():
            raise ValueError("pressure is only defined on corner nodes")
        return 2 * self.N + idx

    def split(self, U):
        """``(u_x, u_y, p)`` with ``p`` on the pressure nodes."""
        U = np.asarray(U)
        return U[:self.N], U[self.N:2 * self.N], U[2 * self.N:]

    def pressure_at_nodes(self, p):
        """Interpolate the linear pressure to every velocity node."""
        out = np.zeros(self.N)
        for name, conn in self.mesh.cells.items():
            el = get_element(name)
            pel = PRESSURE_ELEMENT[name]
            xi, eta = np.array(el.nodes, dtype=float).T
            psi, _, _ = pel.shape(xi, eta)                     # (n_nodes, n_corners)
            pe = p[self.p_index[conn[:, :pel.n_nodes]]]        # (ne, n_corners)
            out[conn] = pe @ psi.T
        return out


class _Block:
    """Precomputed quadrature data of one cell block."""

    def __init__(self, mesh, name, conn, order):
        self.name = name
        self.conn = conn
        self.el = get_element(name)
        self.pel = PRESSURE_ELEMENT[name]
        xi, eta, w = self.el.quadrature(self.el.full_order if order is None else order)
        self.phi, dxi, deta = self.el.shape(xi, eta)
        self.psi, _, _ = self.pel.shape(xi, eta)
        xe, ye = mesh.x[conn], mesh.y[conn]
        hs, self.dphi_dx, self.dphi_dy = physical_gradients(xe, ye, dxi, deta)
        self.wh = w[None, :] * hs
        self.X = xe @ self.phi.T
        self.Y = ye @ self.phi.T
        self.pconn = conn[:, :self.pel.n_nodes]


# -------------------------------------------------------------- assembly
class FlowAssembler:
    """Assembles residual and Jacobian of the Navier-Stokes system."""

    def __init__(self, problem: FlowProblem):
        problem.validate()
        self.problem = problem
        self.mesh = problem.mesh
        self.space = TaylorHoodSpace(self.mesh)
        self.blocks = [_Block(self.mesh, name, conn, problem.order)
                       for name, conn in self.mesh.cells.items()]
        self.K, self.M, self.Bx, self.By = self._linear_matrices()

    # -- helpers ------------------------------------------------------
    def _coo(self, entries):
        rows = np.concatenate([r.ravel() for r, _, _ in entries])
        cols = np.concatenate([c.ravel() for _, c, _ in entries])
        vals = np.concatenate([v.ravel() for _, _, v in entries])
        n = self.space.ndof
        return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    @staticmethod
    def _pair(rdofs, cdofs):
        """Row/column index arrays for batched element matrices (ne, a, b)."""
        return rdofs[:, :, None] + 0 * cdofs[:, None, :], cdofs[:, None, :] + 0 * rdofs[:, :, None]

    def _linear_matrices(self):
        sp_ = self.space
        mu, rho = self.problem.mu, self.problem.rho
        K, M, Bx, By = [], [], [], []
        for b in self.blocks:
            Ke = mu * (np.einsum("eq,eqi,eqj->eij", b.wh, b.dphi_dx, b.dphi_dx)
                       + np.einsum("eq,eqi,eqj->eij", b.wh, b.dphi_dy, b.dphi_dy))
            Me = rho * np.einsum("eq,qi,qj->eij", b.wh, b.phi, b.phi)
            Bxe = np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi_dx)
            Bye = np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi_dy)
            ux, uy, p = sp_.dof_ux(b.conn), sp_.dof_uy(b.conn), sp_.dof_p(b.pconn)
            for dof in (ux, uy):
                r, c = self._pair(dof, dof)
                K.append((r, c, Ke))
                M.append((r, c, Me))
            r, c = self._pair(p, ux)
            Bx.append((r, c, Bxe))
            r, c = self._pair(p, uy)
            By.append((r, c, Bye))
        return self._coo(K), self._coo(M), self._coo(Bx), self._coo(By)

    def body_load(self, t=0.0):
        """Load vector ``rho int phi_i f`` (zero without a body force)."""
        F = np.zeros(self.space.ndof)
        if self.problem.body_force is None:
            return F
        for b in self.blocks:
            fx, fy = (np.broadcast_to(np.asarray(v, dtype=float), b.X.shape)
                      for v in values_at_pair(self.problem.body_force, b.X, b.Y, t))
            np.add.at(F, self.space.dof_ux(b.conn).ravel(),
                      (self.problem.rho * np.einsum("eq,qi->ei", b.wh * fx, b.phi)).ravel())
            np.add.at(F, self.space.dof_uy(b.conn).ravel(),
                      (self.problem.rho * np.einsum("eq,qi->ei", b.wh * fy, b.phi)).ravel())
        return F

    def convection(self, U):
        """Convective residual ``N(u)`` and its Jacobian ``dN/dU`` (sparse)."""
        sp_ = self.space
        rho = self.problem.rho
        u, v, _ = sp_.split(U)
        N = np.zeros(sp_.ndof)
        J = []
        for b in self.blocks:
            ue, ve = u[b.conn], v[b.conn]
            uq, vq = ue @ b.phi.T, ve @ b.phi.T
            dudx = np.einsum("eqi,ei->eq", b.dphi_dx, ue)
            dudy = np.einsum("eqi,ei->eq", b.dphi_dy, ue)
            dvdx = np.einsum("eqi,ei->eq", b.dphi_dx, ve)
            dvdy = np.einsum("eqi,ei->eq", b.dphi_dy, ve)
            wr = b.wh * rho
            Nx = np.einsum("eq,qi->ei", wr * (uq * dudx + vq * dudy), b.phi)
            Ny = np.einsum("eq,qi->ei", wr * (uq * dvdx + vq * dvdy), b.phi)
            ux, uy = sp_.dof_ux(b.conn), sp_.dof_uy(b.conn)
            np.add.at(N, ux.ravel(), Nx.ravel())
            np.add.at(N, uy.ravel(), Ny.ravel())
            # Jacobian: (u . grad) delta_u  +  (delta_u . grad) u
            ugrad = uq[:, :, None] * b.dphi_dx + vq[:, :, None] * b.dphi_dy
            C = np.einsum("eq,qi,eqj->eij", wr, b.phi, ugrad)
            W = {key: np.einsum("eq,qi,qj->eij", wr * g, b.phi, b.phi)
                 for key, g in (("xx", dudx), ("xy", dudy), ("yx", dvdx), ("yy", dvdy))}
            for rd, cd, val in ((ux, ux, C + W["xx"]), (ux, uy, W["xy"]),
                                (uy, ux, W["yx"]), (uy, uy, C + W["yy"])):
                r, c = self._pair(rd, cd)
                J.append((r, c, val))
        return N, self._coo(J)

    # -- boundary conditions ------------------------------------------
    def dirichlet(self, t=0.0):
        """Fixed dofs and their values (velocity components and pressure pin)."""
        sp_ = self.space
        mesh = self.mesh
        fixed = {}
        for tag, spec in self.problem.bc.items():
            if spec == OPEN:
                continue
            nodes = mesh.bc_nodes[tag]
            for comp, val in enumerate(spec):
                if val is None:
                    continue
                dofs = sp_.dof_ux(nodes) if comp == 0 else sp_.dof_uy(nodes)
                vals = values_at(val, mesh.x[nodes], mesh.y[nodes], t)
                fixed.update(zip(dofs.tolist(), vals.tolist()))
        if self.problem.pins_pressure:
            node = sp_.pressure_nodes[0]
            fixed[int(sp_.dof_p([node])[0])] = float(
                values_at(self.problem.pin_value, mesh.x[[node]], mesh.y[[node]])[0])
        dofs = np.array(sorted(fixed), dtype=int)
        return dofs, np.array([fixed[d] for d in dofs])

    # -- steady operator ----------------------------------------------
    def steady_residual_jacobian(self, U, F, with_pressure=True):
        """``R = K u + N(u) - B^T p - F`` and ``B u`` with Jacobian."""
        N, dN = self.convection(U)
        A = self.K + dN
        B = self.Bx + self.By
        R = self.K @ U + N - F
        if with_pressure:
            R = R - B.T @ U + B @ U
            J = A - B.T + B
        else:
            J = A
        return R, J.tocsr()


def values_at_pair(fn, x, y, t=0.0):
    from .util import call_coeff

    return call_coeff(fn, x, y, t)


# ------------------------------------------------------------------ solution
@dataclass
class FlowSolution:
    problem: FlowProblem
    space: TaylorHoodSpace
    U: np.ndarray
    info: dict = field(default_factory=dict)
    assembler: FlowAssembler | None = field(default=None, repr=False)

    @property
    def mesh(self) -> Mesh:
        return self.problem.mesh

    @property
    def u(self) -> np.ndarray:
        return self.space.split(self.U)[0]

    @property
    def v(self) -> np.ndarray:
        return self.space.split(self.U)[1]

    @property
    def p(self) -> np.ndarray:
        """Pressure on the pressure (corner) nodes."""
        return self.space.split(self.U)[2]

    @property
    def p_nodal(self) -> np.ndarray:
        """Pressure interpolated to every node."""
        return self.space.pressure_at_nodes(self.p)

    @property
    def speed(self) -> np.ndarray:
        return np.hypot(self.u, self.v)

    def divergence_norm(self) -> float:
        """L2 norm of ``div u`` (a check of incompressibility)."""
        total = 0.0
        for b in self._assembler().blocks:
            div = (np.einsum("eqi,ei->eq", b.dphi_dx, self.u[b.conn])
                   + np.einsum("eqi,ei->eq", b.dphi_dy, self.v[b.conn]))
            total += np.sum(b.wh * div ** 2)
        return float(np.sqrt(total))

    def _assembler(self):
        if self.assembler is None:
            self.assembler = FlowAssembler(self.problem)
        return self.assembler

    def forces(self, tag: str, order: int = 3):
        """Force ``(Fx, Fy)`` exerted by the fluid on the boundary ``tag``.

        Integrates the traction of ``sigma = -p I + mu (grad u + grad u^T)``
        over the edges of the tag, using the parent elements' gradients.
        """
        return traction_forces(self, tag, order)

    def save(self, path):
        """Write velocity (3-component vector), speed and pressure for ParaView."""
        from .mesh_io import write_vtk

        vel = np.column_stack([self.u, self.v, np.zeros(self.space.N)])
        return write_vtk(path, self.mesh, point_data={"velocity": vel, "speed": self.speed,
                                                      "p": self.p_nodal})


def traction_forces(sol: FlowSolution, tag: str, order: int = 3):
    mesh = sol.mesh
    mu = sol.problem.mu
    if tag not in mesh.boundary:
        raise ValueError(f"Unknown boundary tag {tag!r}; mesh has {mesh.tags}")
    wanted = {tuple(sorted(e[:2])) for e in mesh.boundary[tag]}
    s, w1 = gauss_legendre_quad(order)
    p_nodal = sol.p_nodal
    Fx = Fy = 0.0
    for name, conn in mesh.cells.items():
        el = get_element(name)
        pel = PRESSURE_ELEMENT[name]
        ref = np.array(el.nodes, dtype=float)
        for face in el.faces:
            key = np.sort(conn[:, list(face[:2])], axis=1)
            sel = np.array([tuple(k) in wanted for k in key.tolist()])
            if not sel.any():
                continue
            ce = conn[sel]
            # quadrature points along the face in natural coordinates
            a, b_ = ref[face[0]], ref[face[1]]
            xi = 0.5 * (1 - s) * a[0] + 0.5 * (1 + s) * b_[0]
            eta = 0.5 * (1 - s) * a[1] + 0.5 * (1 + s) * b_[1]
            phi, dxi, deta = el.shape(xi, eta)
            psi, _, _ = pel.shape(xi, eta)
            xe, ye = mesh.x[ce], mesh.y[ce]
            _, dphi_dx, dphi_dy = physical_gradients(xe, ye, dxi, deta)
            ue, ve = sol.u[ce], sol.v[ce]
            dudx = np.einsum("eqi,ei->eq", dphi_dx, ue)
            dudy = np.einsum("eqi,ei->eq", dphi_dy, ue)
            dvdx = np.einsum("eqi,ei->eq", dphi_dx, ve)
            dvdy = np.einsum("eqi,ei->eq", dphi_dy, ve)
            pq = p_nodal[ce[:, :pel.n_nodes]] @ psi.T
            # tangent along the face (domain on the left) from the face nodes
            fn = ce[:, list(face)]
            if len(face) == 3:
                dl = np.column_stack([s - 0.5, s + 0.5, -2 * s])
            else:
                dl = np.column_stack([-0.5 * np.ones_like(s), 0.5 * np.ones_like(s)])
            tx = mesh.x[fn] @ dl.T
            ty = mesh.y[fn] @ dl.T
            ds = np.hypot(tx, ty)
            nx, ny = ty / ds, -tx / ds            # outward normal of the fluid domain
            sxx = -pq + 2 * mu * dudx
            syy = -pq + 2 * mu * dvdy
            sxy = mu * (dudy + dvdx)
            # force on the boundary object = - int sigma . n_fluid ds
            Fx -= np.sum(w1[None, :] * ds * (sxx * nx + sxy * ny))
            Fy -= np.sum(w1[None, :] * ds * (sxy * nx + syy * ny))
    return float(Fx), float(Fy)


# ------------------------------------------------------------------- solvers
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


@dataclass
class TransientFlowSolution:
    problem: FlowProblem
    space: TaylorHoodSpace
    times: np.ndarray
    snapshots: list
    info: dict = field(default_factory=dict)
    assembler: FlowAssembler | None = field(default=None, repr=False)

    @property
    def mesh(self) -> Mesh:
        return self.problem.mesh

    def at(self, index: int = -1) -> FlowSolution:
        """Stored step as a :class:`FlowSolution`."""
        return FlowSolution(self.problem, self.space, self.snapshots[index], self.info,
                            self.assembler)

    @property
    def final(self) -> FlowSolution:
        return self.at(-1)

    def save_series(self, prefix):
        """One ``.vtu`` per stored step plus a ``.pvd`` collection."""
        import pathlib

        from .mesh_io import write_vtk

        prefix = pathlib.Path(prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        entries = []
        for i, t in enumerate(self.times):
            sol = self.at(i)
            name = f"{prefix.name}_{i:04d}.vtu"
            vel = np.column_stack([sol.u, sol.v, np.zeros(self.space.N)])
            write_vtk(prefix.parent / name, self.mesh,
                      point_data={"velocity": vel, "speed": sol.speed, "p": sol.p_nodal})
            entries.append(
                f'    <DataSet timestep="{float(t)!r}" group="" part="0" file="{name}"/>')
        pvd = prefix.with_suffix(".pvd")
        pvd.write_text(
            '<?xml version="1.0"?>\n'
            '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
            "  <Collection>\n" + "\n".join(entries) + "\n  </Collection>\n</VTKFile>\n")
        return pvd


def solve_flow_transient(problem: FlowProblem, dt: float, t_end: float, theta: float = 0.5,
                         U0=None, method: str = "direct", store_every: int = 1,
                         verbose: bool = False, rtol: float = 1e-8, atol: float = 1e-10,
                         max_newton: int = 25, damping: bool = True) -> TransientFlowSolution:
    """Time integration with the theta-method and a Newton solve per step.

    Momentum: ``M (U - U_n)/dt + theta S(U) + (1 - theta) S(U_n) - B^T p = 0``
    with ``S(U) = K U + N(U) - F``; the pressure and the continuity equation
    are implicit.  ``U0`` may be a velocity/pressure vector or a callable
    ``(x, y) -> (ux, uy)`` for the initial velocity.
    """
    if not 0.0 <= theta <= 1.0:
        raise ValueError("theta must be in [0, 1]")
    nsteps = int(round(t_end / dt))
    if nsteps < 1 or abs(nsteps * dt - t_end) > 1e-8 * max(1.0, abs(t_end)):
        raise ValueError(f"dt = {dt} must divide t_end = {t_end}")

    asm = FlowAssembler(problem)
    space = asm.space
    mesh = problem.mesh
    B = (asm.Bx + asm.By).tocsr()
    M = asm.M
    time_dependent = problem.depends_on_time()

    fixed, vals = asm.dirichlet(0.0)
    if U0 is None:
        U = np.zeros(space.ndof)
    elif callable(U0):
        ux, uy = values_at_pair(U0, mesh.x, mesh.y)
        U = np.concatenate([np.broadcast_to(ux, mesh.x.shape), np.broadcast_to(uy, mesh.x.shape),
                            np.zeros(space.Np)]).astype(float)
    else:
        U = np.array(U0, dtype=float, copy=True)
    U[fixed] = vals

    F_old = asm.body_load(0.0)
    N_old, _ = asm.convection(U)
    S_old = asm.K @ U + N_old - F_old
    S_old[2 * space.N:] = 0.0

    times = [0.0]
    snapshots = [U.copy()]
    newton_iterations = []
    if verbose:
        print(f"Taylor-Hood: {space.ndof} unknowns; {nsteps} steps of dt = {dt} (theta = {theta})")

    for n in range(1, nsteps + 1):
        t = n * dt
        if time_dependent:
            fixed, vals = asm.dirichlet(t)
        F_new = asm.body_load(t) if time_dependent else F_old
        U_old = U

        def residual_jacobian(Un, U_old=U_old, F_new=F_new, S_old=S_old):
            N, dN = asm.convection(Un)
            S = asm.K @ Un + N - F_new
            S[2 * space.N:] = 0.0
            R = M @ ((Un - U_old) / dt) + theta * S + (1 - theta) * S_old - B.T @ Un + B @ Un
            J = (M / dt + theta * (asm.K + dN) - B.T + B).tocsr()
            return R, J

        guess = U_old.copy()
        guess[fixed] = vals
        res = newton_iterate(residual_jacobian, guess, fixed, method, rtol=rtol, atol=atol,
                             max_newton=max_newton, damping=damping)
        if not res.converged:
            raise RuntimeError(f"Newton did not converge at t = {t:.6g} "
                               f"(|R| = {res.residuals[-1]:.2e})")
        U = res.T
        N_new, _ = asm.convection(U)
        S_old = asm.K @ U + N_new - F_new
        S_old[2 * space.N:] = 0.0
        F_old = F_new
        newton_iterations.append(res.iterations)
        if n % store_every == 0 or n == nsteps:
            times.append(t)
            snapshots.append(U.copy())
        if verbose and (n % max(1, nsteps // 10) == 0 or n == nsteps):
            u, v, _ = space.split(U)
            print(f"  step {n}/{nsteps}, t = {t:.6g}, max |u| = {np.hypot(u, v).max():.6g}, "
                  f"{res.iterations} Newton iterations")

    info = {"steps": nsteps, "dt": dt, "theta": theta, "newton_iterations": newton_iterations}
    return TransientFlowSolution(problem, space, np.asarray(times), snapshots, info, asm)
