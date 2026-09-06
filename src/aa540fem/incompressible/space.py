"""Taylor-Hood degree-of-freedom layout and per-block quadrature data (2-D and 3-D)."""

from __future__ import annotations

import numpy as np

from aa540fem.core.elements import PRESSURE_ELEMENT, get_element
from aa540fem.core.mesh import Mesh
from aa540fem.core.quadrature import quadrature_points
from aa540fem.transport.element import (
    element_metric_nd,
    jacobian_nd,
    map_gradients_nd,
    physical_laplacian_nd,
)


class TaylorHoodSpace:
    """Degree-of-freedom layout ``[u_1, ..., u_d, p]`` for a Taylor-Hood pair:
    the ``d`` velocity components on all ``N`` nodes, then the pressure on
    the ``Np`` corner nodes."""

    def __init__(self, mesh: Mesh):
        self.mesh = mesh
        self.dim = mesh.dim
        self.N = mesh.n_nodes
        corners = np.concatenate([conn[:, :get_element(name).n_corners].ravel()
                                  for name, conn in mesh.cells.items()])
        self.pressure_nodes = np.unique(corners)
        self.Np = self.pressure_nodes.size
        self.p_index = -np.ones(self.N, dtype=int)
        self.p_index[self.pressure_nodes] = np.arange(self.Np)
        self.ndof = self.dim * self.N + self.Np

    @property
    def n_vel(self) -> int:
        """Number of velocity dofs (the first block)."""
        return self.dim * self.N

    def dof_u(self, component: int, nodes):
        """Dofs of velocity component ``component`` (0, 1[, 2]) at ``nodes``."""
        return component * self.N + np.asarray(nodes, dtype=int)

    def dof_ux(self, nodes):
        return self.dof_u(0, nodes)

    def dof_uy(self, nodes):
        return self.dof_u(1, nodes)

    def dof_uz(self, nodes):
        if self.dim < 3:
            raise ValueError("a 2-D space has no u_z")
        return self.dof_u(2, nodes)

    def dof_p(self, nodes):
        idx = self.p_index[np.asarray(nodes, dtype=int)]
        if (idx < 0).any():
            raise ValueError("pressure is only defined on corner nodes")
        return self.dim * self.N + idx

    def local_dofs(self, conn, pconn):
        """Element-local global dof numbers ``[u_1 nodes, ..., u_d nodes, p corners]``
        as an ``(n_elems, d n + n_corners)`` array, the layout of the element
        kernels."""
        return np.concatenate([self.dof_u(c, conn) for c in range(self.dim)]
                              + [self.dof_p(pconn)], axis=1)

    def split(self, U):
        """``(u_1, ..., u_d, p)`` with ``p`` on the pressure nodes."""
        U = np.asarray(U)
        return tuple(U[c * self.N:(c + 1) * self.N] for c in range(self.dim)) + (U[self.n_vel:],)

    def velocity(self, U):
        """The velocity components as a ``(d, N)`` array."""
        return np.asarray(U)[:self.n_vel].reshape(self.dim, self.N)

    def pressure_at_nodes(self, p):
        """Interpolate the linear pressure to every velocity node."""
        out = np.zeros(self.N)
        for name, conn in self.mesh.cells.items():
            el = get_element(name)
            pel = PRESSURE_ELEMENT[name]
            coords = np.array(el.nodes, dtype=float).T
            psi = pel.shape(*coords)[0]                       # (n_nodes, n_corners)
            pe = p[self.p_index[conn[:, :pel.n_nodes]]]        # (ne, n_corners)
            out[conn] = pe @ psi.T
        return out


class _Block:
    """Precomputed quadrature data of one cell block.

    ``dphi`` / ``dpsi`` are tuples of the ``d`` physical gradient arrays
    ``(ne, nq, n)`` of the velocity / pressure shape functions (also reachable
    as ``dphi_dx``, ``dphi_dy``[, ``dphi_dz``]); ``Gmat`` is the element metric
    tensor ``(ne, nq, d, d)`` and ``G`` its unique components
    ``(xx, xy, yy)`` or ``(xx, xy, xz, yy, yz, zz)``; ``Xq`` the quadrature
    point coordinates ``(x, y[, z])``.
    """

    def __init__(self, mesh, name, conn, order):
        self.name = name
        self.conn = conn
        self.el = get_element(name)
        self.pel = PRESSURE_ELEMENT[name]
        self.dim = d = self.el.dim
        coords, w = quadrature_points(self.el.family,
                                      self.el.full_order if order is None else order)
        self.phi, dnat = self.el.shape_at(coords)
        self.psi, dpnat = self.pel.shape_at(coords)
        xe = tuple(mesh.points[conn, k] for k in range(d))
        hs, inverse = jacobian_nd(xe, dnat)
        self.dphi = map_gradients_nd(inverse, dnat)
        self.dpsi = map_gradients_nd(inverse, dpnat)
        self.lap_phi = physical_laplacian_nd(inverse, self.el.hessian(*coords))
        self.Gmat = element_metric_nd(inverse, self.el.family)
        pairs = ((0, 0), (0, 1), (1, 1)) if d == 2 else \
            ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
        self.G = tuple(np.ascontiguousarray(self.Gmat[..., i, j]) for i, j in pairs)
        self.wh = w[None, :] * hs
        self.Xq = tuple(xe[k] @ self.phi.T for k in range(d))
        self.X, self.Y = self.Xq[0], self.Xq[1]
        if d == 3:
            self.Z = self.Xq[2]
        self.pconn = conn[:, :self.pel.n_nodes]
        self.n = conn.shape[1]                                     # velocity nodes per element
        self.nc = self.pconn.shape[1]                              # pressure (corner) nodes
        self.L = d * self.n + self.nc                              # local dofs
        self.ldof = None                                           # set by the assembler

    @property
    def dphi_dx(self):
        return self.dphi[0]

    @property
    def dphi_dy(self):
        return self.dphi[1]

    @property
    def dphi_dz(self):
        return self.dphi[2]

    @property
    def dpsi_dx(self):
        return self.dpsi[0]

    @property
    def dpsi_dy(self):
        return self.dpsi[1]

    @property
    def dpsi_dz(self):
        return self.dpsi[2]
