"""Flow solutions: fields, forces on boundaries, ParaView output."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from aa540fem.core.elements import FACE_CORNERS, FACE_TYPES, PRESSURE_ELEMENT, get_element
from aa540fem.core.mesh import Mesh
from aa540fem.core.quadrature import gauss_legendre_quad, quadrature_points
from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import FlowProblem
from aa540fem.incompressible.space import TaylorHoodSpace
from aa540fem.transport.element import jacobian_nd, map_gradients_nd, physical_gradients

AXES = "xyz"


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
    def w(self) -> np.ndarray:
        """Third velocity component (3-D only)."""
        if self.space.dim < 3:
            raise AttributeError("a 2-D solution has no w component")
        return self.space.split(self.U)[2]

    @property
    def velocity(self) -> np.ndarray:
        """Velocity components as a ``(d, N)`` array."""
        return self.space.velocity(self.U)

    @property
    def p(self) -> np.ndarray:
        """Pressure on the pressure (corner) nodes."""
        return self.space.split(self.U)[-1]

    @property
    def p_nodal(self) -> np.ndarray:
        """Pressure interpolated to every node."""
        return self.space.pressure_at_nodes(self.p)

    @property
    def speed(self) -> np.ndarray:
        return np.linalg.norm(self.velocity, axis=0)

    def divergence_norm(self) -> float:
        """L2 norm of ``div u`` (a check of incompressibility)."""
        total = 0.0
        vel = self.velocity
        for b in self._assembler().blocks:
            div = sum(np.einsum("eqi,ei->eq", b.dphi[c], vel[c][b.conn]) for c in range(b.dim))
            total += np.sum(b.wh * div ** 2)
        return float(np.sqrt(total))

    def _assembler(self):
        if self.assembler is None:
            self.assembler = FlowAssembler(self.problem)
        return self.assembler

    def forces(self, tag: str, order: int = 3):
        """Force ``(Fx, Fy[, Fz])`` exerted by the fluid on the boundary ``tag``.

        Integrates the traction of ``sigma = -p I + mu (grad u + grad u^T)``
        over the faces of the tag, using the parent elements' gradients.
        """
        return traction_forces(self, tag, order)

    def wall_traction(self, tag: str, order: int = 3) -> dict:
        """Traction exerted by the fluid on the boundary ``tag`` at the face
        quadrature points, sorted along ``x``.

        Returns a dict of flat arrays: ``x``, ``y``[, ``z``], ``tx``, ``ty``
        [, ``tz``] (force per unit length / area on the wall), ``nx``, ``ny``
        [, ``nz``] (outward normal of the fluid domain), ``p`` and ``weight``
        (quadrature weight times edge length or face area, so that
        ``sum(weight * tx)`` is the force).  The skin friction coefficient is
        ``2 tx / (rho U^2)`` on a wall aligned with ``x``.  In 3-D ``order``
        is the quadrature order of the face element (``None``: its default).
        """
        return wall_traction(self, tag, order)

    def pressure_at(self, points) -> np.ndarray:
        """Pressure at arbitrary points ``(m, d)``.

        On simplex meshes (triangles, tetrahedra) the containing element is
        found from the corner (barycentric) coordinates and the linear
        pressure interpolated there, which is exact; on quadrilateral and
        hexahedral meshes the nearest pressure node is used.  Points outside
        the mesh take the value of the element they are closest to.
        """
        return pressure_at_points(self, np.atleast_2d(np.asarray(points, dtype=float)))

    def velocity_3d(self) -> np.ndarray:
        """Velocity as an ``(N, 3)`` array (zero third component in 2-D) for output."""
        vel = self.velocity.T
        return vel if vel.shape[1] == 3 else np.column_stack([vel, np.zeros(self.space.N)])

    def save(self, path):
        """Write velocity (3-component vector), speed and pressure for ParaView."""
        from aa540fem.io.mesh_files import write_vtk
        return write_vtk(path, self.mesh, point_data={"velocity": self.velocity_3d(),
                                                      "speed": self.speed, "p": self.p_nodal})


def wall_traction(sol: FlowSolution, tag: str, order: int = 3) -> dict:
    mesh = sol.mesh
    mu = sol.problem.mu
    if tag not in mesh.boundary:
        raise ValueError(f"Unknown boundary tag {tag!r}; mesh has {mesh.tags}")
    if mesh.dim == 3:
        return _wall_traction_3d(sol, tag, order)
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


def _wall_traction_3d(sol: FlowSolution, tag: str, order) -> dict:
    """3-D version of :func:`wall_traction`: quadrature on the boundary faces
    (triangles or quadrilaterals) mapped into the parent elements."""
    mesh = sol.mesh
    mu = sol.problem.mu
    vel = sol.velocity
    p_nodal = sol.p_nodal
    keys = [f"{a}" for a in AXES] + [f"t{a}" for a in AXES] + [f"n{a}" for a in AXES]
    out = {k: [] for k in keys + ["p", "weight"]}
    for name, conn in mesh.cells.items():
        el = get_element(name)
        pel = PRESSURE_ELEMENT[name]
        ref = np.array(el.nodes, dtype=float)
        for face, fel in el.face_elements():
            sel = boundary_face_elements(mesh, tag, name, face)
            if sel.size == 0:
                continue
            fcoords, fw = quadrature_points(fel.family, fel.full_order if order in (None, 3)
                                            else order)
            fphi, fdnat = fel.shape_at(fcoords)                 # face shape functions
            ce = conn[sel]
            # face quadrature points in the parent's natural coordinates (faces are
            # flat in reference space, so the face interpolation is exact)
            nat = tuple(fphi @ ref[list(face), k] for k in range(3))
            phi, dnat = el.shape_at(nat)
            psi = pel.shape(*nat)[0]
            xe = tuple(mesh.points[ce, k] for k in range(3))
            _, inverse = jacobian_nd(xe, dnat)
            dphi = map_gradients_nd(inverse, dnat)
            grad = [[np.einsum("eqi,ei->eq", dphi[k], vel[c][ce]) for k in range(3)]
                    for c in range(3)]
            pq = p_nodal[ce[:, :pel.n_nodes]] @ psi.T
            # tangents of the face from its own shape functions, normal = t1 x t2
            fn = ce[:, list(face)]
            t1 = [mesh.points[fn, k] @ fdnat[0].T for k in range(3)]
            t2 = [mesh.points[fn, k] @ fdnat[1].T for k in range(3)]
            nvec = [t1[1] * t2[2] - t1[2] * t2[1], t1[2] * t2[0] - t1[0] * t2[2],
                    t1[0] * t2[1] - t1[1] * t2[0]]
            dA = np.sqrt(sum(c * c for c in nvec))
            nrm = [c / dA for c in nvec]                        # outward normal of the fluid
            sigma = [[(-pq if c == k else 0.0) + mu * (grad[c][k] + grad[k][c]) for k in range(3)]
                     for c in range(3)]
            for c in range(3):
                out[f"t{AXES[c]}"].append(-sum(sigma[c][k] * nrm[k] for k in range(3)).ravel())
                out[f"n{AXES[c]}"].append(nrm[c].ravel())
                out[AXES[c]].append((xe[c] @ phi.T).ravel())
            out["p"].append(pq.ravel())
            out["weight"].append((fw[None, :] * dA).ravel())
    if not out["x"]:
        return {k: np.zeros(0) for k in out}
    out = {k: np.concatenate(v) for k, v in out.items()}
    order_ = np.argsort(out["x"], kind="stable")
    return {k: v[order_] for k, v in out.items()}


def pressure_at_points(sol: FlowSolution, points: np.ndarray) -> np.ndarray:
    mesh = sol.mesh
    d = mesh.dim
    p = sol.p
    out = np.full(points.shape[0], np.nan)
    best = np.full(points.shape[0], np.inf)
    for name, conn in mesh.cells.items():
        el = get_element(name)
        nc = el.n_corners
        if nc != d + 1:                                     # not a simplex: nearest node
            nodes = sol.space.pressure_nodes
            dist = np.linalg.norm(points[:, None, :] - mesh.points[nodes][None, :, :], axis=2)
            k = dist.argmin(axis=1)
            closer = dist[np.arange(points.shape[0]), k] < best
            out[closer] = p[k[closer]]
            best[closer] = dist[np.arange(points.shape[0]), k][closer]
            continue
        corners = mesh.points[conn[:, :nc]]                 # (ne, d+1, d)
        # barycentric coordinates of every point in every element (small meshes:
        # chunk the elements to bound the memory)
        T = np.transpose(corners[:, 1:, :] - corners[:, :1, :], (0, 2, 1))    # (ne, d, d)
        Tinv = np.linalg.inv(T)
        for q, x in enumerate(points):
            lam = np.einsum("eij,ej->ei", Tinv, x[None, :] - corners[:, 0, :])
            lam = np.column_stack([1.0 - lam.sum(axis=1), lam])
            violation = np.maximum(0.0, -lam).max(axis=1)      # 0 inside the element
            e = int(violation.argmin())
            if violation[e] < best[q]:
                best[q] = violation[e]
                pe = p[sol.space.p_index[conn[e, :nc]]]
                out[q] = float(np.clip(lam[e], 0.0, 1.0) @ pe / max(np.clip(lam[e], 0, 1).sum(),
                                                                       1e-300))
    return out


def traction_forces(sol: FlowSolution, tag: str, order: int = 3):
    """Force ``(Fx, Fy[, Fz])`` on the boundary ``tag``: the integral of
    :func:`wall_traction`."""
    tr = wall_traction(sol, tag, order)
    return tuple(float(np.sum(tr["weight"] * tr[f"t{AXES[c]}"])) for c in range(sol.space.dim))


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
            write_vtk(path, self.mesh, point_data={"velocity": sol.velocity_3d(),
                                                   "speed": sol.speed, "p": sol.p_nodal})

        return write_series(prefix, self.times, write_step)


def boundary_face_elements(mesh, tag, name, face):
    """Indices of the elements of block ``name`` whose local ``face`` lies on
    boundary ``tag`` (cached on the mesh: the force callbacks of the time
    integrators call this every step)."""
    cache = mesh.__dict__.setdefault("_face_cache", {})
    key = (tag, name, tuple(face))
    if key not in cache:
        el = get_element(name)
        ftype = el.face_types[el.faces.index(tuple(face))]
        nc = FACE_CORNERS[ftype]
        k = FACE_TYPES[mesh.dim][ftype]
        blocks = [b[:, :nc] for b in mesh.face_blocks(tag) if b.shape[1] == k]
        wanted = np.sort(np.vstack(blocks), axis=1) if blocks else np.zeros((0, nc), dtype=int)
        have = np.sort(mesh.cells[name][:, list(face[:nc])], axis=1)
        if wanted.size == 0:
            cache[key] = np.zeros(0, dtype=int)
        else:
            both = np.vstack([wanted, have])
            _, inverse = np.unique(both, axis=0, return_inverse=True)
            inverse = inverse.ravel()
            in_tag = np.zeros(inverse.max() + 1, dtype=bool)
            in_tag[inverse[:wanted.shape[0]]] = True
            cache[key] = np.nonzero(in_tag[inverse[wanted.shape[0]:]])[0]
    return cache[key]

