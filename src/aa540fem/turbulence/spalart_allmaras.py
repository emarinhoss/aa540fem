"""Spalart-Allmaras one-equation turbulence model (negative variant, no trip).

The model transports the working variable ``nu_tilde``:

    d(nu_tilde)/dt + u . grad(nu_tilde)
        = c_b1 (1 - f_t2) S_tilde nu_tilde
          - [c_w1 f_w - c_b1/kappa^2 f_t2] (nu_tilde / d)^2
          + (1/sigma) [ div((nu + nu_tilde) grad nu_tilde) + c_b2 |grad nu_tilde|^2 ]

with the eddy viscosity ``nu_t = nu_tilde f_v1``.  Definitions and constants
follow the NASA Turbulence Modeling Resource description of "SA-neg" (the
negative-nu_tilde modification of Allmaras, Johnson & Spalart 2012, which
keeps the equation well posed when nu_tilde becomes negative) without the
f_t2 trip term ("SA-neg-noft2"); the original model is Spalart & Allmaras,
La Recherche Aerospatiale 1, 1994 (AIAA 92-0439).

Discretisation: the same quadratic elements and quadrature blocks as the
velocity, SUPG streamline upwinding with the convective velocity of the
flow solution, source terms and their derivatives evaluated at the
quadrature points (derivatives by central finite differences), Newton with
pseudo-transient continuation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from aa540fem.incompressible.space import _Block
from aa540fem.linalg.continuation import pseudo_transient
from aa540fem.linalg.newton import NewtonResult
from aa540fem.transport.element import supg_length

FD_STEP = 1e-6


@dataclass
class SpalartAllmaras:
    """Model functions and constants.  ``nu`` is the laminar kinematic viscosity."""

    nu: float
    cb1: float = 0.1355
    sigma: float = 2.0 / 3.0
    cb2: float = 0.622
    kappa: float = 0.41
    cw2: float = 0.3
    cw3: float = 2.0
    cv1: float = 7.1
    ct3: float = 1.2
    ct4: float = 0.5
    cv2: float = 0.7
    cv3: float = 0.9
    cn1: float = 16.0
    r_limit: float = 10.0
    ft2: bool = False          # trip term; off, as in the NASA TMR "SA-noft2" variant

    @property
    def cw1(self) -> float:
        return self.cb1 / self.kappa ** 2 + (1.0 + self.cb2) / self.sigma

    def fv1(self, chi):
        chi3 = chi ** 3
        return chi3 / (chi3 + self.cv1 ** 3)

    def eddy_viscosity(self, nu_tilde):
        """``nu_t = nu_tilde f_v1`` (zero where ``nu_tilde`` is negative)."""
        nt = np.asarray(nu_tilde, dtype=float)
        chi = nt / self.nu
        return np.where(nt > 0, nt * self.fv1(chi), 0.0)

    def diffusivity(self, nu_tilde):
        """``(nu + f_n nu_tilde) / sigma`` with ``f_n = 1`` for positive ``nu_tilde``."""
        nt = np.asarray(nu_tilde, dtype=float)
        chi = nt / self.nu
        chi3 = chi ** 3
        fn = np.where(nt >= 0, 1.0, (self.cn1 + chi3) / (self.cn1 - chi3))
        return (self.nu + fn * nt) / self.sigma

    def modified_vorticity(self, nu_tilde, omega, d):
        """``S_tilde`` with the SA-neg guard against negative values."""
        nt = np.asarray(nu_tilde, dtype=float)
        chi = nt / self.nu
        fv1 = self.fv1(chi)
        fv2 = 1.0 - chi / (1.0 + chi * fv1)
        sbar = nt * fv2 / (self.kappa ** 2 * d ** 2)
        low = sbar < -self.cv2 * omega
        safe = np.where(low, ((self.cv3 - 2 * self.cv2) * omega - sbar), 1.0)
        denom = np.where(safe == 0, 1e-300, safe)
        guarded = omega + omega * (self.cv2 ** 2 * omega + self.cv3 * sbar) / denom
        return np.where(low, guarded, omega + sbar)

    def source(self, nu_tilde, gx, gy, omega, d):
        """Production minus destruction plus the ``c_b2`` gradient term, at points."""
        nt = np.asarray(nu_tilde, dtype=float)
        d = np.maximum(d, 1e-12)
        chi = nt / self.nu
        grad2 = gx ** 2 + gy ** 2
        # positive branch
        s_tilde = np.maximum(self.modified_vorticity(nt, omega, d), 1e-300)
        r = np.minimum(nt / (s_tilde * self.kappa ** 2 * d ** 2), self.r_limit)
        g = r + self.cw2 * (r ** 6 - r)
        fw = g * ((1.0 + self.cw3 ** 6) / (g ** 6 + self.cw3 ** 6)) ** (1.0 / 6.0)
        ft2 = self.ct3 * np.exp(-self.ct4 * chi ** 2) if self.ft2 else 0.0
        prod = self.cb1 * (1.0 - ft2) * s_tilde * nt
        destr = (self.cw1 * fw - self.cb1 / self.kappa ** 2 * ft2) * (nt / d) ** 2
        positive = prod - destr
        # negative branch (Allmaras, Johnson & Spalart 2012)
        negative = self.cb1 * (1.0 - self.ct3) * omega * nt + self.cw1 * (nt / d) ** 2
        return np.where(nt >= 0, positive, negative) + self.cb2 / self.sigma * grad2


class SpalartAllmarasSolver:
    """Steady SA equation on the velocity mesh for a frozen velocity field."""

    def __init__(self, mesh, model: SpalartAllmaras, wall_tags, order=None):
        self.mesh = mesh
        self.model = model
        self.N = mesh.n_nodes
        self.blocks = [_Block(mesh, name, conn, order) for name, conn in mesh.cells.items()]
        self.distance = mesh.wall_distance(wall_tags)
        self.M = self._mass()
        self.u = np.zeros(self.N)
        self.v = np.zeros(self.N)

    def _mass(self):
        rows, cols, vals = [], [], []
        for b in self.blocks:
            Me = np.einsum("eq,qi,qj->eij", b.wh, b.phi, b.phi)
            n = b.conn.shape[1]
            rows.append(np.repeat(b.conn, n, axis=1).ravel())
            cols.append(np.tile(b.conn, (1, n)).ravel())
            vals.append(Me.ravel())
        return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.N, self.N)).tocsr()

    def set_velocity(self, u, v):
        self.u = np.asarray(u, dtype=float)
        self.v = np.asarray(v, dtype=float)

    def residual_jacobian(self, nu_tilde, supg: bool = True):
        """Weak residual of the steady equation and its Jacobian (CSR).

        Convection, diffusion and source are Galerkin terms; SUPG applies the
        streamline weight to the convective and source parts of the strong
        residual (the diffusion part is omitted).  Source and diffusivity
        derivatives are central finite differences at the quadrature points.
        """
        m = self.model
        nt = np.asarray(nu_tilde, dtype=float)
        R = np.zeros(self.N)
        rows, cols, vals = [], [], []
        for b in self.blocks:
            conn = b.conn
            nte = nt[conn]
            ntq = nte @ b.phi.T
            gx = np.einsum("eqi,ei->eq", b.dphi_dx, nte)
            gy = np.einsum("eqi,ei->eq", b.dphi_dy, nte)
            ue, ve = self.u[conn], self.v[conn]
            uq, vq = ue @ b.phi.T, ve @ b.phi.T
            omega = np.abs(np.einsum("eqi,ei->eq", b.dphi_dx, ve)
                           - np.einsum("eqi,ei->eq", b.dphi_dy, ue))
            dq = np.maximum(self.distance[conn] @ b.phi.T, 1e-12)
            wh = b.wh

            # diffusivity and its derivative
            D = m.diffusivity(ntq)
            h_n = FD_STEP * (m.nu + np.abs(ntq))
            dD = (m.diffusivity(ntq + h_n) - m.diffusivity(ntq - h_n)) / (2 * h_n)
            # source and its derivatives with respect to nu_tilde and its gradient
            s = m.source(ntq, gx, gy, omega, dq)
            s_n = (m.source(ntq + h_n, gx, gy, omega, dq)
                   - m.source(ntq - h_n, gx, gy, omega, dq)) / (2 * h_n)
            h_g = FD_STEP * (1.0 + np.abs(gx) + np.abs(gy))
            s_gx = (m.source(ntq, gx + h_g, gy, omega, dq)
                    - m.source(ntq, gx - h_g, gy, omega, dq)) / (2 * h_g)
            s_gy = (m.source(ntq, gx, gy + h_g, omega, dq)
                    - m.source(ntq, gx, gy - h_g, omega, dq)) / (2 * h_g)

            ugrad = uq[:, :, None] * b.dphi_dx + vq[:, :, None] * b.dphi_dy   # u . grad phi_j
            conv = uq * gx + vq * gy                                          # u . grad nu_tilde
            grad_i = b.dphi_dx * gx[:, :, None] + b.dphi_dy * gy[:, :, None]  # grad phi_i . grad nt
            dS = (s_n[:, :, None] * b.phi[None] + s_gx[:, :, None] * b.dphi_dx
                  + s_gy[:, :, None] * b.dphi_dy)                              # d s / d nt_j

            Re_ = (np.einsum("eq,qi->ei", wh * (conv - s), b.phi)
                   + np.einsum("eq,eqi->ei", wh * D, grad_i))
            Je = (np.einsum("eq,qi,eqj->eij", wh, b.phi, ugrad - dS)
                  + np.einsum("eq,eqi,eqj->eij", wh * D, b.dphi_dx, b.dphi_dx)
                  + np.einsum("eq,eqi,eqj->eij", wh * D, b.dphi_dy, b.dphi_dy)
                  + np.einsum("eq,eqi,qj->eij", wh * dD, grad_i, b.phi))
            if supg:
                umag = np.hypot(uq, vq)
                moving = umag > 0
                safe = np.where(moving, umag, 1.0)
                sx = np.where(moving, uq / safe, 1.0)
                sy = np.where(moving, vq / safe, 0.0)
                h = supg_length(sx, sy, b.dphi_dx, b.dphi_dy)
                tau = 1.0 / np.sqrt((2.0 * umag / h) ** 2 + (4.0 * D / h ** 2) ** 2)
                w_i = tau[:, :, None] * ugrad
                Re_ = Re_ + np.einsum("eq,eqi->ei", wh * (conv - s), w_i)
                Je = Je + np.einsum("eq,eqi,eqj->eij", wh, w_i, ugrad - dS)

            np.add.at(R, conn.ravel(), Re_.ravel())
            n = conn.shape[1]
            rows.append(np.repeat(conn, n, axis=1).ravel())
            cols.append(np.tile(conn, (1, n)).ravel())
            vals.append(Je.ravel())
        J = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.N, self.N)).tocsr()
        return R, J

    def solve(self, nu_tilde0, fixed, values, method="direct", rtol=1e-6, atol=1e-30,
              dtau0=0.01, max_steps=200, inner_newton=6, verbose=False) -> NewtonResult:
        """Steady solve by pseudo-transient continuation from ``nu_tilde0``.

        ``fixed``/``values``: Dirichlet nodes and values (walls: 0, inflow:
        the freestream level).
        """
        nt = np.array(nu_tilde0, dtype=float, copy=True)
        nt[fixed] = values
        return pseudo_transient(self.residual_jacobian, nt, fixed, self.M, method, rtol, atol,
                                dtau0, max_steps, inner_newton=inner_newton, verbose=verbose)
