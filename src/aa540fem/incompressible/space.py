"""Taylor-Hood degree-of-freedom layout and per-block quadrature data."""

from __future__ import annotations

import numpy as np

from aa540fem.core.elements import PRESSURE_ELEMENT, get_element
from aa540fem.core.mesh import Mesh
from aa540fem.transport.element import jacobian, map_gradients, physical_laplacian


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
        self.psi, dpsi_dxi, dpsi_deta = self.pel.shape(xi, eta)
        xe, ye = mesh.x[conn], mesh.y[conn]
        hs, *inverse = jacobian(xe, ye, dxi, deta)
        self.dphi_dx, self.dphi_dy = map_gradients(inverse, dxi, deta)
        self.dpsi_dx, self.dpsi_dy = map_gradients(inverse, dpsi_dxi, dpsi_deta)
        self.lap_phi = physical_laplacian(inverse, *self.el.hessian(xi, eta))
        self.wh = w[None, :] * hs
        self.X = xe @ self.phi.T
        self.Y = ye @ self.phi.T
        self.pconn = conn[:, :self.pel.n_nodes]
