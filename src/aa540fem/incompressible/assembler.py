"""Residual and Jacobian assembly of the Navier-Stokes system.

Galerkin terms: viscous ``K``, mass ``M`` and divergence ``B`` matrices
(assembled once) and the convective residual with its Jacobian.  Optional
residual-based stabilisation (``FlowProblem.stabilisation``): SUPG on the
momentum equations with the streamline weight ``tau (u . grad phi_i)``
applied to the full momentum residual

    R_m = rho (u . grad) u - mu lap u + grad p - rho f  [+ rho (u - u_n)/dt]

plus grad-div (LSIC) ``int gamma (div v)(div u)``, and optionally PSPG
``int tau grad q . R_m`` on the continuity equation.  The Laplacian of the
discrete velocity comes from the shape-function Hessians (exact on affine
elements), so the stabilisation is consistent: an exact polynomial solution
such as Poiseuille flow is reproduced exactly.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import scipy.sparse as sp

from aa540fem.core.util import values_at, values_rate
from aa540fem.incompressible.problem import OPEN, FlowProblem, values_at_pair
from aa540fem.incompressible.space import TaylorHoodSpace, _Block


class MomentumTerms(NamedTuple):
    """Nonlinear parts of the residual and their Jacobians (global, sparse)."""

    N: np.ndarray            # Galerkin convective residual
    J_N: sp.csr_matrix       # its Jacobian
    S: np.ndarray            # stabilisation residual (zero when off)
    J_S: sp.csr_matrix       # its Jacobian (empty when off)


class FlowAssembler:
    """Assembles residual and Jacobian of the Navier-Stokes system."""

    def __init__(self, problem: FlowProblem):
        problem.validate()
        self.problem = problem
        self.mesh = problem.mesh
        self.space = TaylorHoodSpace(self.mesh)
        self.blocks = [_Block(self.mesh, name, conn, problem.order)
                       for name, conn in self.mesh.cells.items()]
        self.mu_t = None
        if problem.eddy_viscosity is not None:
            self.mu_t = np.asarray(problem.eddy_viscosity, dtype=float)
            if self.mu_t.shape != (self.space.N,):
                raise ValueError("eddy_viscosity must be a nodal array")
        self.K, self.M, self.Bx, self.By = self._linear_matrices()

    def viscosity(self, b):
        """``mu_eff`` and its gradient at the quadrature points of block ``b``."""
        mu = self.problem.mu
        if self.mu_t is None:
            zero = np.zeros_like(b.X)
            return np.full_like(b.X, mu), zero, zero
        mt = self.mu_t[b.conn]
        return (mu + mt @ b.phi.T, np.einsum("eqi,ei->eq", b.dphi_dx, mt),
                np.einsum("eqi,ei->eq", b.dphi_dy, mt))

    # -- helpers ------------------------------------------------------
    def _coo(self, entries):
        n = self.space.ndof
        if not entries:
            return sp.csr_matrix((n, n))
        rows = np.concatenate([r.ravel() for r, _, _ in entries])
        cols = np.concatenate([c.ravel() for _, c, _ in entries])
        vals = np.concatenate([v.ravel() for _, _, v in entries])
        return sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()

    @staticmethod
    def _pair(rdofs, cdofs):
        """Row/column index arrays for batched element matrices (ne, a, b)."""
        return rdofs[:, :, None] + 0 * cdofs[:, None, :], cdofs[:, None, :] + 0 * rdofs[:, :, None]

    def _linear_matrices(self):
        sp_ = self.space
        rho = self.problem.rho
        K, M, Bx, By = [], [], [], []
        for b in self.blocks:
            mu_q, dmu_dx, dmu_dy = self.viscosity(b)
            wmu = b.wh * mu_q
            Ke = (np.einsum("eq,eqi,eqj->eij", wmu, b.dphi_dx, b.dphi_dx)
                  + np.einsum("eq,eqi,eqj->eij", wmu, b.dphi_dy, b.dphi_dy))
            Me = rho * np.einsum("eq,qi,qj->eij", b.wh, b.phi, b.phi)
            Bxe = np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi_dx)
            Bye = np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi_dy)
            ux, uy, p = sp_.dof_ux(b.conn), sp_.dof_uy(b.conn), sp_.dof_p(b.pconn)
            dofs = (ux, uy)
            for dof in dofs:
                r, c = self._pair(dof, dof)
                K.append((r, c, Ke))
                M.append((r, c, Me))
            if self.mu_t is not None:
                # variable viscosity term - grad(u)^T . grad(mu):
                # block (c, d) = - int phi_i (d_d mu) (d_c phi_j)
                dgrad = (b.dphi_dx, b.dphi_dy)
                dmu = (dmu_dx, dmu_dy)
                for c in range(2):
                    for d in range(2):
                        r, cc = self._pair(dofs[c], dofs[d])
                        K.append((r, cc, -np.einsum("eq,qi,eqj->eij", b.wh * dmu[d], b.phi,
                                                    dgrad[c])))
            r, c = self._pair(p, ux)
            Bx.append((r, c, Bxe))
            r, c = self._pair(p, uy)
            By.append((r, c, Bye))
        return self._coo(K), self._coo(M), self._coo(Bx), self._coo(By)

    def _body_force(self, b, t):
        if self.problem.body_force is None:
            return None
        return tuple(np.broadcast_to(np.asarray(v, dtype=float), b.X.shape)
                     for v in values_at_pair(self.problem.body_force, b.X, b.Y, t))

    def body_load(self, t=0.0):
        """Load vector ``rho int phi_i f`` (zero without a body force)."""
        F = np.zeros(self.space.ndof)
        if self.problem.body_force is None:
            return F
        for b in self.blocks:
            fx, fy = self._body_force(b, t)
            np.add.at(F, self.space.dof_ux(b.conn).ravel(),
                      (self.problem.rho * np.einsum("eq,qi->ei", b.wh * fx, b.phi)).ravel())
            np.add.at(F, self.space.dof_uy(b.conn).ravel(),
                      (self.problem.rho * np.einsum("eq,qi->ei", b.wh * fy, b.phi)).ravel())
        return F

    # -- nonlinear terms ----------------------------------------------
    def momentum_terms(self, U, t=0.0, dt=None, U_old=None, stabilise=None,
                       pspg=None, param_state=None) -> MomentumTerms:
        """Convective residual/Jacobian and, if enabled, the stabilisation terms.

        ``dt`` and ``U_old`` add ``rho (u - u_old)/dt`` to the momentum
        residual used by the stabilisation (theta scheme) and select the
        transient ``tau``.  ``stabilise``/``pspg`` override the problem
        settings (the RK45 path passes ``pspg=False``).  The stabilisation
        parameters ``tau`` and ``gamma`` are differentiated in the Jacobian
        (through the velocity magnitude and through the flow-direction
        dependence of the element length, which matters on stretched cells);
        with ``param_state`` given they are evaluated at that state and
        frozen instead.
        """
        prob = self.problem
        stab = prob.stabilisation if stabilise is None else stabilise
        pspg = prob.pspg if pspg is None else pspg
        sp_ = self.space
        rho = prob.rho
        u, v, p = sp_.split(U)
        if U_old is not None:
            u_old, v_old, _ = sp_.split(U_old)
        u_par, v_par, _ = sp_.split(U if param_state is None else param_state)

        N = np.zeros(sp_.ndof)
        S = np.zeros(sp_.ndof)
        JN, JS = [], []
        for b in self.blocks:
            ue, ve = u[b.conn], v[b.conn]
            uq, vq = ue @ b.phi.T, ve @ b.phi.T
            dudx = np.einsum("eqi,ei->eq", b.dphi_dx, ue)
            dudy = np.einsum("eqi,ei->eq", b.dphi_dy, ue)
            dvdx = np.einsum("eqi,ei->eq", b.dphi_dx, ve)
            dvdy = np.einsum("eqi,ei->eq", b.dphi_dy, ve)
            wh = b.wh
            phi = b.phi[None]                                        # (1, q, n)
            ugrad = uq[:, :, None] * b.dphi_dx + vq[:, :, None] * b.dphi_dy
            ux, uy, pd = sp_.dof_ux(b.conn), sp_.dof_uy(b.conn), sp_.dof_p(b.pconn)
            dofs = (ux, uy)

            # Galerkin convection: residual and Jacobian d/du [rho (u . grad) u]
            conv = rho * (uq * dudx + vq * dudy), rho * (uq * dvdx + vq * dvdy)
            np.add.at(N, ux.ravel(), np.einsum("eq,qi->ei", wh * conv[0], b.phi).ravel())
            np.add.at(N, uy.ravel(), np.einsum("eq,qi->ei", wh * conv[1], b.phi).ravel())
            dconv = ((rho * (ugrad + dudx[:, :, None] * phi), rho * dudy[:, :, None] * phi),
                     (rho * dvdx[:, :, None] * phi, rho * (ugrad + dvdy[:, :, None] * phi)))
            for c in range(2):
                for d in range(2):
                    r, cc = self._pair(dofs[c], dofs[d])
                    JN.append((r, cc, np.einsum("eq,qi,eqj->eij", wh, b.phi, dconv[c][d])))

            if not (stab or pspg):
                continue

            # full momentum residual at the quadrature points and its derivatives
            mu_q, dmu_dx, dmu_dy = self.viscosity(b)
            nu = mu_q / rho
            lap = np.einsum("eqi,ei->eq", b.lap_phi, ue), np.einsum("eqi,ei->eq", b.lap_phi, ve)
            grads = ((dudx, dudy), (dvdx, dvdy))
            # - div(mu grad u) - grad(u)^T grad(mu)
            #   = - mu lap u - grad mu . grad u - grad(u)^T grad mu
            visc = [-mu_q * lap[c] - (dmu_dx * grads[c][0] + dmu_dy * grads[c][1])
                    - (dmu_dx * grads[0][c] + dmu_dy * grads[1][c]) for c in range(2)]
            pe = p[sp_.p_index[b.pconn]]
            gradp = np.einsum("eqk,ek->eq", b.dpsi_dx, pe), np.einsum("eqk,ek->eq", b.dpsi_dy, pe)
            body = self._body_force(b, t)
            R = [conv[c] + visc[c] + gradp[c] - (rho * body[c] if body else 0.0)
                 for c in range(2)]
            dgradq = (b.dphi_dx, b.dphi_dy)
            dmu = (dmu_dx, dmu_dy)
            dR = [[dconv[c][d]
                   - (mu_q[:, :, None] * b.lap_phi
                      + dmu_dx[:, :, None] * b.dphi_dx + dmu_dy[:, :, None] * b.dphi_dy
                      if c == d else 0.0)
                   - dmu[d][:, :, None] * dgradq[c]          # d/du_{d,j} of -d_c u_d d_d mu
                   for d in range(2)] for c in range(2)]
            if dt is not None and U_old is not None:
                uo = u_old[b.conn] @ b.phi.T, v_old[b.conn] @ b.phi.T
                for c, (cur, old) in enumerate(((uq, uo[0]), (vq, uo[1]))):
                    R[c] = R[c] + rho * (cur - old) / dt
                    dR[c][c] = dR[c][c] + rho * phi / dt
            dpsi = (b.dpsi_dx, b.dpsi_dy)

            # stabilisation parameters (Tezduyar), from the parameter state
            upq, vpq = u_par[b.conn] @ b.phi.T, v_par[b.conn] @ b.phi.T
            umag = np.hypot(upq, vpq)
            moving = umag > 0
            safe = np.where(moving, umag, 1.0)
            sx = np.where(moving, upq / safe, 1.0)
            sy = np.where(moving, vpq / safe, 0.0)
            sgrad = sx[:, :, None] * b.dphi_dx + sy[:, :, None] * b.dphi_dy
            h = 2.0 / np.maximum(np.abs(sgrad).sum(axis=2), 1e-300)
            inv_dt2 = (2.0 / dt) ** 2 if dt is not None else 0.0
            tau = 1.0 / np.sqrt(inv_dt2 + (2.0 * umag / h) ** 2 + (4.0 * nu / h ** 2) ** 2)
            re_h = umag * h / (2.0 * nu)
            low = re_h < 3.0
            if prob.grad_div:
                gamma = 0.5 * h * umag * np.minimum(1.0, re_h / 3.0)
                dgamma_du = 0.5 * h * np.where(low, 2.0 * re_h / 3.0, 1.0)      # d gamma / d|u|
                dgamma_dh = np.where(low, h * umag ** 2 / (6.0 * nu), 0.5 * umag)  # d gamma / dh
            else:
                gamma = dgamma_du = dgamma_dh = 0.0 * umag
            dtau_du = -(4.0 / h ** 2) * umag * tau ** 3                            # d tau / d|u|
            dtau_dh = tau ** 3 * (4.0 * umag ** 2 / h ** 3 + 32.0 * nu ** 2 / h ** 5)  # d tau / dh
            # derivatives of the parameters with respect to the velocity components:
            # through |u| and through the flow-direction dependence of the element
            # length h(s), s = u/|u| (unless frozen at a separate state)
            follow = param_state is None
            dir_ = (np.where(moving, upq / safe, 0.0), np.where(moving, vpq / safe, 0.0))
            sgn = np.sign(sgrad)
            dh_ds = (-0.5 * h ** 2 * np.einsum("eqi,eqi->eq", sgn, b.dphi_dx),
                     -0.5 * h ** 2 * np.einsum("eqi,eqi->eq", sgn, b.dphi_dy))
            inv_u = np.where(moving, 1.0 / safe, 0.0)
            dh_du = ((dh_ds[0] * (1.0 - sx * sx) - dh_ds[1] * sy * sx) * inv_u,
                     (-dh_ds[0] * sx * sy + dh_ds[1] * (1.0 - sy * sy)) * inv_u)
            dtau_d = [dtau_du * dir_[d] + dtau_dh * dh_du[d] for d in range(2)]
            dgamma_d = [dgamma_du * dir_[d] + dgamma_dh * dh_du[d] for d in range(2)]
            dgrad = (b.dphi_dx, b.dphi_dy)

            if stab:
                w_i = tau[:, :, None] * ugrad                             # tau (u . grad phi_i)
                div = dudx + dvdy
                for c in range(2):
                    Sc = (np.einsum("eq,eqi->ei", wh * R[c], w_i)
                          + np.einsum("eq,eqi->ei", wh * gamma * div, dgrad[c]))
                    np.add.at(S, dofs[c].ravel(), Sc.ravel())
                    for d in range(2):
                        J = (np.einsum("eq,eqi,eqj->eij", wh, w_i, dR[c][d])
                             + np.einsum("eq,eqi,qj->eij", wh * tau * R[c], dgrad[d], b.phi)
                             + np.einsum("eq,eqi,eqj->eij", wh * gamma, dgrad[c], dgrad[d]))
                        if follow:
                            J = J + (np.einsum("eq,eqi,qj->eij", wh * R[c] * dtau_d[d],
                                               ugrad, b.phi)
                                     + np.einsum("eq,eqi,qj->eij", wh * div * dgamma_d[d],
                                                 dgrad[c], b.phi))
                        r, cc = self._pair(dofs[c], dofs[d])
                        JS.append((r, cc, J))
                    r, cc = self._pair(dofs[c], pd)
                    JS.append((r, cc, np.einsum("eq,eqi,eqk->eik", wh, w_i, dpsi[c])))

            if pspg:
                gradR = dpsi[0] * R[0][:, :, None] + dpsi[1] * R[1][:, :, None]
                Sp = np.einsum("eq,eqk->ek", wh * tau, gradR)
                np.add.at(S, pd.ravel(), Sp.ravel())
                for d in range(2):
                    J = (np.einsum("eq,eqk,eqj->ekj", wh * tau, dpsi[0], dR[0][d])
                         + np.einsum("eq,eqk,eqj->ekj", wh * tau, dpsi[1], dR[1][d]))
                    if follow:
                        J = J + np.einsum("eq,eqk,qj->ekj", wh * dtau_d[d], gradR, b.phi)
                    r, cc = self._pair(pd, dofs[d])
                    JS.append((r, cc, J))
                r, cc = self._pair(pd, pd)
                JS.append((r, cc, np.einsum("eq,eqk,eql->ekl", wh * tau, dpsi[0], dpsi[0])
                           + np.einsum("eq,eqk,eql->ekl", wh * tau, dpsi[1], dpsi[1])))
        return MomentumTerms(N, self._coo(JN), S, self._coo(JS))

    def convection(self, U, jacobian: bool = True):
        """Galerkin convective residual ``N(u)`` and its Jacobian (no stabilisation)."""
        mt = self.momentum_terms(U, stabilise=False, pspg=False)
        if not jacobian:
            return mt.N
        return mt.N, mt.J_N

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
    def steady_residual_jacobian(self, U, F, with_pressure=True, t=0.0):
        """``R = K u + N(u) + S(u) - B^T p - F`` and ``B u`` with the Jacobian."""
        mt = self.momentum_terms(U, t=t)
        A = self.K + mt.J_N + mt.J_S
        B = self.Bx + self.By
        R = self.K @ U + mt.N + mt.S - F
        if with_pressure:
            R = R - B.T @ U + B @ U
            J = A - B.T + B
        else:
            J = A
        return R, J.tocsr()
