"""Residual and Jacobian assembly of the Navier-Stokes system."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from aa540fem.core.util import values_at, values_rate
from aa540fem.incompressible.problem import OPEN, FlowProblem, values_at_pair
from aa540fem.incompressible.space import TaylorHoodSpace, _Block


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

    def convection(self, U, jacobian: bool = True):
        """Convective residual ``N(u)`` and its Jacobian ``dN/dU`` (sparse).

        With ``jacobian=False`` only ``N(u)`` is returned.
        """
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
            if not jacobian:
                continue
            # Jacobian: (u . grad) delta_u  +  (delta_u . grad) u
            ugrad = uq[:, :, None] * b.dphi_dx + vq[:, :, None] * b.dphi_dy
            C = np.einsum("eq,qi,eqj->eij", wr, b.phi, ugrad)
            W = {key: np.einsum("eq,qi,qj->eij", wr * g, b.phi, b.phi)
                 for key, g in (("xx", dudx), ("xy", dudy), ("yx", dvdx), ("yy", dvdy))}
            for rd, cd, val in ((ux, ux, C + W["xx"]), (ux, uy, W["xy"]),
                                (uy, ux, W["yx"]), (uy, uy, C + W["yy"])):
                r, c = self._pair(rd, cd)
                J.append((r, c, val))
        if not jacobian:
            return N
        return N, self._coo(J)

    # -- boundary conditions ------------------------------------------
    def dirichlet(self, t=0.0, rate: bool = False):
        """Fixed dofs and their values (velocity components and pressure pin).

        With ``rate=True`` the time derivatives of the velocity values are
        returned (zero for the pressure pin and for constant values).
        """
        sp_ = self.space
        mesh = self.mesh
        evaluate = values_rate if rate else values_at
        fixed = {}
        for tag, spec in self.problem.bc.items():
            if spec == OPEN:
                continue
            nodes = mesh.bc_nodes[tag]
            for comp, val in enumerate(spec):
                if val is None:
                    continue
                dofs = sp_.dof_ux(nodes) if comp == 0 else sp_.dof_uy(nodes)
                vals = evaluate(val, mesh.x[nodes], mesh.y[nodes], t)
                fixed.update(zip(dofs.tolist(), vals.tolist()))
        if self.problem.pins_pressure:
            node = sp_.pressure_nodes[0]
            fixed[int(sp_.dof_p([node])[0])] = 0.0 if rate else float(
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
