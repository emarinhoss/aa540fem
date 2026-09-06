"""Flow solutions: fields, forces on boundaries, ParaView output."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from aa540fem.core.elements import PRESSURE_ELEMENT, get_element
from aa540fem.core.mesh import Mesh
from aa540fem.core.quadrature import gauss_legendre_quad
from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import FlowProblem
from aa540fem.incompressible.space import TaylorHoodSpace
from aa540fem.transport.element import physical_gradients


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

    def wall_traction(self, tag: str, order: int = 3) -> dict:
        """Traction exerted by the fluid on the boundary ``tag`` at the edge
        quadrature points, sorted along ``x``.

        Returns a dict of flat arrays: ``x``, ``y``, ``tx``, ``ty`` (force per
        unit length on the wall), ``nx``, ``ny`` (outward normal of the fluid
        domain), ``p`` and ``weight`` (quadrature weight times edge length,
        so that ``sum(weight * tx)`` is the force).  The skin friction
        coefficient is ``2 tx / (rho U^2)`` on a wall aligned with ``x``.
        """
        return wall_traction(self, tag, order)

    def save(self, path):
        """Write velocity (3-component vector), speed and pressure for ParaView."""
        from aa540fem.io.mesh_files import write_vtk
        vel = np.column_stack([self.u, self.v, np.zeros(self.space.N)])
        return write_vtk(path, self.mesh, point_data={"velocity": vel, "speed": self.speed,
                                                      "p": self.p_nodal})


def wall_traction(sol: FlowSolution, tag: str, order: int = 3) -> dict:
    mesh = sol.mesh
    mu = sol.problem.mu
    if tag not in mesh.boundary:
        raise ValueError(f"Unknown boundary tag {tag!r}; mesh has {mesh.tags}")
    s, w1 = gauss_legendre_quad(order)
    p_nodal = sol.p_nodal
    out = {k: [] for k in ("x", "y", "tx", "ty", "nx", "ny", "p", "weight")}
    for name, conn in mesh.cells.items():
        el = get_element(name)
        pel = PRESSURE_ELEMENT[name]
        ref = np.array(el.nodes, dtype=float)
        for face in el.faces:
            sel = boundary_face_elements(mesh, tag, name, face)
            if sel.size == 0:
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
            # traction on the boundary object = - sigma . n_fluid
            out["tx"].append(-(sxx * nx + sxy * ny).ravel())
            out["ty"].append(-(sxy * nx + syy * ny).ravel())
            out["nx"].append(nx.ravel())
            out["ny"].append(ny.ravel())
            out["p"].append(pq.ravel())
            out["weight"].append((w1[None, :] * ds).ravel())
            out["x"].append((xe @ phi.T).ravel())
            out["y"].append((ye @ phi.T).ravel())
    if not out["x"]:
        return {k: np.zeros(0) for k in out}
    out = {k: np.concatenate(v) for k, v in out.items()}
    order_ = np.argsort(out["x"], kind="stable")
    return {k: v[order_] for k, v in out.items()}


def traction_forces(sol: FlowSolution, tag: str, order: int = 3):
    """Force ``(Fx, Fy)`` on the boundary ``tag``: the integral of :func:`wall_traction`."""
    tr = wall_traction(sol, tag, order)
    return float(np.sum(tr["weight"] * tr["tx"])), float(np.sum(tr["weight"] * tr["ty"]))


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
        from aa540fem.io.mesh_files import write_vtk
        from aa540fem.io.series import write_series

        def write_step(path, i):
            sol = self.at(i)
            vel = np.column_stack([sol.u, sol.v, np.zeros(self.space.N)])
            write_vtk(path, self.mesh,
                      point_data={"velocity": vel, "speed": sol.speed, "p": sol.p_nodal})

        return write_series(prefix, self.times, write_step)


def boundary_face_elements(mesh, tag, name, face):
    """Indices of the elements of block ``name`` whose local ``face`` lies on
    boundary ``tag`` (cached on the mesh: the force callbacks of the time
    integrators call this every step)."""
    cache = mesh.__dict__.setdefault("_face_cache", {})
    key = (tag, name, tuple(face))
    if key not in cache:
        edges = np.sort(mesh.boundary[tag][:, :2], axis=1)
        wanted = edges[:, 0] * mesh.n_nodes + edges[:, 1]
        conn = mesh.cells[name]
        pair = np.sort(conn[:, list(face[:2])], axis=1)
        keys = pair[:, 0] * mesh.n_nodes + pair[:, 1]
        cache[key] = np.nonzero(np.isin(keys, wanted))[0]
    return cache[key]

