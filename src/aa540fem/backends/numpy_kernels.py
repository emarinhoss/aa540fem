"""Reference (NumPy) element kernels of the Navier-Stokes assembly.

Every kernel works on one cell block (:class:`aa540fem.incompressible.space._Block`)
and returns element-local arrays in the layout ``[u_x nodes (n), u_y nodes
(n), p corner nodes (nc)]``, i.e. ``L = 2 n + nc`` local dofs, which the
assembler scatters with its :class:`~aa540fem.backends.pattern.ScatterPlan`.
The threaded numba kernels (:mod:`aa540fem.backends.numba_kernels`) produce
the same arrays; this module is the definition they are tested against.
"""

from __future__ import annotations

import numpy as np


# -- constant matrices ------------------------------------------------
def linear_local(b, mu_q, dmu_dx, dmu_dy, rho, variable_viscosity):
    """Element viscous ``K``, mass ``M``, divergence ``Bx``, ``By`` and the
    transposed divergence ``B^T`` (gradient) blocks, each ``(ne, L, L)``."""
    n, L = b.n, b.L
    ne = b.conn.shape[0]
    wmu = b.wh * mu_q
    Kuu = (np.einsum("eq,eqi,eqj->eij", wmu, b.dphi_dx, b.dphi_dx)
           + np.einsum("eq,eqi,eqj->eij", wmu, b.dphi_dy, b.dphi_dy))
    Muu = rho * np.einsum("eq,qi,qj->eij", b.wh, b.phi, b.phi)
    Bxe = np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi_dx)             # (ne, nc, n)
    Bye = np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi_dy)
    K = np.zeros((ne, L, L))
    M = np.zeros((ne, L, L))
    Bx = np.zeros((ne, L, L))
    By = np.zeros((ne, L, L))
    BT = np.zeros((ne, L, L))
    for c in range(2):
        K[:, c * n:(c + 1) * n, c * n:(c + 1) * n] = Kuu
        M[:, c * n:(c + 1) * n, c * n:(c + 1) * n] = Muu
    if variable_viscosity:
        # variable viscosity term - grad(u)^T . grad(mu):
        # block (c, d) = - int phi_i (d_d mu) (d_c phi_j)
        dgrad = (b.dphi_dx, b.dphi_dy)
        dmu = (dmu_dx, dmu_dy)
        for c in range(2):
            for d in range(2):
                K[:, c * n:(c + 1) * n, d * n:(d + 1) * n] -= np.einsum(
                    "eq,qi,eqj->eij", b.wh * dmu[d], b.phi, dgrad[c])
    Bx[:, 2 * n:, :n] = Bxe
    By[:, 2 * n:, n:2 * n] = Bye
    BT[:, :n, 2 * n:] = np.transpose(Bxe, (0, 2, 1))
    BT[:, n:2 * n, 2 * n:] = np.transpose(Bye, (0, 2, 1))
    return K, M, Bx, By, BT


def body_load_local(b, fx, fy, rho):
    """Element load ``rho int phi_i f`` as ``(ne, 2n)``."""
    return np.concatenate([rho * np.einsum("eq,qi->ei", b.wh * fx, b.phi),
                           rho * np.einsum("eq,qi->ei", b.wh * fy, b.phi)], axis=1)


# -- nonlinear terms --------------------------------------------------
def momentum_local(b, U, par, old, mu_q, dmu_dx, dmu_dy, body, options, backend="numpy"):
    """Element-local convective and stabilisation residuals and Jacobians.

    Returns ``(N_e (ne, 2n), S_e (ne, L) or None, JN_e (ne, 2n, 2n) or None,
    JS_e (ne, L, L) or None)``.  ``U``/``par``/``old`` are global dof vectors
    (the state, the state the stabilisation parameters are evaluated at,
    and the previous time level or ``None``); ``body`` is ``(fx, fy)`` at the
    quadrature points or ``None``; ``options`` holds ``rho, stab, pspg,
    grad_div, metric, follow, jacobian, inv_dt, inv_dt2``.
    """
    if backend == "numba":
        from aa540fem.backends.numba_kernels import momentum_local as _numba_momentum

        return _numba_momentum(b, U, par, old, mu_q, dmu_dx, dmu_dy, body, options)
    rho = options["rho"]
    stab, pspg, jac = options["stab"], options["pspg"], options["jacobian"]
    n, L = b.n, b.L
    ne = b.conn.shape[0]
    ux, uy = b.ldof[:, :n], b.ldof[:, n:2 * n]
    pd = b.ldof[:, 2 * n:]
    ue, ve = U[ux], U[uy]
    uq, vq = ue @ b.phi.T, ve @ b.phi.T
    dudx = np.einsum("eqi,ei->eq", b.dphi_dx, ue)
    dudy = np.einsum("eqi,ei->eq", b.dphi_dy, ue)
    dvdx = np.einsum("eqi,ei->eq", b.dphi_dx, ve)
    dvdy = np.einsum("eqi,ei->eq", b.dphi_dy, ve)
    wh = b.wh
    phi = b.phi[None]                                        # (1, q, n)
    ugrad = uq[:, :, None] * b.dphi_dx + vq[:, :, None] * b.dphi_dy

    # Galerkin convection: residual and Jacobian d/du [rho (u . grad) u]
    conv = rho * (uq * dudx + vq * dudy), rho * (uq * dvdx + vq * dvdy)
    N_e = np.concatenate([np.einsum("eq,qi->ei", wh * conv[0], b.phi),
                          np.einsum("eq,qi->ei", wh * conv[1], b.phi)], axis=1)
    dconv = ((rho * (ugrad + dudx[:, :, None] * phi), rho * dudy[:, :, None] * phi),
             (rho * dvdx[:, :, None] * phi, rho * (ugrad + dvdy[:, :, None] * phi)))
    JN_e = None
    if jac:
        JN_e = np.empty((ne, 2 * n, 2 * n))
        for c in range(2):
            for d in range(2):
                JN_e[:, c * n:(c + 1) * n, d * n:(d + 1) * n] = np.einsum(
                    "eq,qi,eqj->eij", wh, b.phi, dconv[c][d])
    if not (stab or pspg):
        return N_e, None, JN_e, (np.zeros((ne, L, L)) if jac else None)

    # full momentum residual at the quadrature points and its derivatives
    nu = mu_q / rho
    lap = np.einsum("eqi,ei->eq", b.lap_phi, ue), np.einsum("eqi,ei->eq", b.lap_phi, ve)
    grads = ((dudx, dudy), (dvdx, dvdy))
    # - div(mu grad u) - grad(u)^T grad(mu) = - mu lap u - grad mu . grad u - grad(u)^T grad mu
    visc = [-mu_q * lap[c] - (dmu_dx * grads[c][0] + dmu_dy * grads[c][1])
            - (dmu_dx * grads[0][c] + dmu_dy * grads[1][c]) for c in range(2)]
    pe = U[pd]
    gradp = np.einsum("eqk,ek->eq", b.dpsi_dx, pe), np.einsum("eqk,ek->eq", b.dpsi_dy, pe)
    R = [conv[c] + visc[c] + gradp[c] - (rho * body[c] if body is not None else 0.0)
         for c in range(2)]
    dgrad = (b.dphi_dx, b.dphi_dy)
    dmu = (dmu_dx, dmu_dy)
    if old is not None:
        uo = old[ux] @ b.phi.T, old[uy] @ b.phi.T
        R[0] = R[0] + rho * (uq - uo[0]) * options["inv_dt"]
        R[1] = R[1] + rho * (vq - uo[1]) * options["inv_dt"]
    dpsi = (b.dpsi_dx, b.dpsi_dy)

    # stabilisation parameters from the parameter state
    upq, vpq = par[ux] @ b.phi.T, par[uy] @ b.phi.T
    parameters = metric_parameters if options["metric"] else streamline_parameters
    tau, gamma, dtau_d, dgamma_d = parameters(b, upq, vpq, nu, options["inv_dt2"],
                                              options["grad_div"])
    follow = options["follow"]

    S_e = np.zeros((ne, L))
    if stab:
        w_i = tau[:, :, None] * ugrad                             # tau (u . grad phi_i)
        div = dudx + dvdy
        for c in range(2):
            S_e[:, c * n:(c + 1) * n] = (np.einsum("eq,eqi->ei", wh * R[c], w_i)
                                         + np.einsum("eq,eqi->ei", wh * gamma * div, dgrad[c]))
    gradR = None
    if pspg:
        gradR = dpsi[0] * R[0][:, :, None] + dpsi[1] * R[1][:, :, None]
        S_e[:, 2 * n:] = np.einsum("eq,eqk->ek", wh * tau, gradR)
    if not jac:
        return N_e, S_e, None, None

    dR = [[dconv[c][d]
           - (mu_q[:, :, None] * b.lap_phi
              + dmu_dx[:, :, None] * b.dphi_dx + dmu_dy[:, :, None] * b.dphi_dy
              if c == d else 0.0)
           - dmu[d][:, :, None] * dgrad[c]          # d/du_{d,j} of -d_c u_d d_d mu
           for d in range(2)] for c in range(2)]
    if old is not None:
        for c in range(2):
            dR[c][c] = dR[c][c] + rho * phi * options["inv_dt"]
    JS_e = np.zeros((ne, L, L))
    if stab:
        for c in range(2):
            rows = slice(c * n, (c + 1) * n)
            for d in range(2):
                J = (np.einsum("eq,eqi,eqj->eij", wh, w_i, dR[c][d])
                     + np.einsum("eq,eqi,qj->eij", wh * tau * R[c], dgrad[d], b.phi)
                     + np.einsum("eq,eqi,eqj->eij", wh * gamma, dgrad[c], dgrad[d]))
                if follow:
                    J = J + (np.einsum("eq,eqi,qj->eij", wh * R[c] * dtau_d[d], ugrad, b.phi)
                             + np.einsum("eq,eqi,qj->eij", wh * div * dgamma_d[d],
                                         dgrad[c], b.phi))
                JS_e[:, rows, d * n:(d + 1) * n] = J
            JS_e[:, rows, 2 * n:] = np.einsum("eq,eqi,eqk->eik", wh, w_i, dpsi[c])
    if pspg:
        for d in range(2):
            J = (np.einsum("eq,eqk,eqj->ekj", wh * tau, dpsi[0], dR[0][d])
                 + np.einsum("eq,eqk,eqj->ekj", wh * tau, dpsi[1], dR[1][d]))
            if follow:
                J = J + np.einsum("eq,eqk,qj->ekj", wh * dtau_d[d], gradR, b.phi)
            JS_e[:, 2 * n:, d * n:(d + 1) * n] = J
        JS_e[:, 2 * n:, 2 * n:] = (np.einsum("eq,eqk,eql->ekl", wh * tau, dpsi[0], dpsi[0])
                                   + np.einsum("eq,eqk,eql->ekl", wh * tau, dpsi[1], dpsi[1]))
    return N_e, S_e, JN_e, JS_e


# -- stabilisation parameters -----------------------------------------
def streamline_parameters(b, upq, vpq, nu, inv_dt2, grad_div):
    """Tezduyar's ``tau`` and ``gamma`` with the flow-direction element length.

    ``h = 2 / sum_i |s . grad phi_i|``, ``s = u / |u|``;
    ``tau = [(2/dt)^2 + (2|u|/h)^2 + (4 nu/h^2)^2]^(-1/2)``;
    ``gamma = h |u| / 2 min(1, Re_h / 3)``, ``Re_h = |u| h / (2 nu)``.
    Returns ``tau, gamma`` and their derivatives with respect to the two
    velocity components at the quadrature points (through ``|u|`` and through
    ``h(s)``), each ``(n_elems, nq)``.
    """
    umag = np.hypot(upq, vpq)
    moving = umag > 0
    safe = np.where(moving, umag, 1.0)
    sx = np.where(moving, upq / safe, 1.0)
    sy = np.where(moving, vpq / safe, 0.0)
    sgrad = sx[:, :, None] * b.dphi_dx + sy[:, :, None] * b.dphi_dy
    h = 2.0 / np.maximum(np.abs(sgrad).sum(axis=2), 1e-300)
    tau = 1.0 / np.sqrt(inv_dt2 + (2.0 * umag / h) ** 2 + (4.0 * nu / h ** 2) ** 2)
    re_h = umag * h / (2.0 * nu)
    low = re_h < 3.0
    if grad_div:
        gamma = 0.5 * h * umag * np.minimum(1.0, re_h / 3.0)
        dgamma_du = 0.5 * h * np.where(low, 2.0 * re_h / 3.0, 1.0)      # d gamma / d|u|
        dgamma_dh = np.where(low, h * umag ** 2 / (6.0 * nu), 0.5 * umag)  # d gamma / dh
    else:
        gamma = dgamma_du = dgamma_dh = 0.0 * umag
    dtau_du = -(4.0 / h ** 2) * umag * tau ** 3                            # d tau / d|u|
    dtau_dh = tau ** 3 * (4.0 * umag ** 2 / h ** 3 + 32.0 * nu ** 2 / h ** 5)  # d tau / dh
    dir_ = (np.where(moving, upq / safe, 0.0), np.where(moving, vpq / safe, 0.0))
    sgn = np.sign(sgrad)
    dh_ds = (-0.5 * h ** 2 * np.einsum("eqi,eqi->eq", sgn, b.dphi_dx),
             -0.5 * h ** 2 * np.einsum("eqi,eqi->eq", sgn, b.dphi_dy))
    inv_u = np.where(moving, 1.0 / safe, 0.0)
    dh_du = ((dh_ds[0] * (1.0 - sx * sx) - dh_ds[1] * sy * sx) * inv_u,
             (-dh_ds[0] * sx * sy + dh_ds[1] * (1.0 - sy * sy)) * inv_u)
    dtau_d = [dtau_du * dir_[d] + dtau_dh * dh_du[d] for d in range(2)]
    dgamma_d = [dgamma_du * dir_[d] + dgamma_dh * dh_du[d] for d in range(2)]
    return tau, gamma, dtau_d, dgamma_d


def metric_parameters(b, upq, vpq, nu, inv_dt2, grad_div):
    """``tau`` and ``gamma`` from the element metric tensor ``G`` (smooth in ``u``).

    ``tau = [(2/dt)^2 + u.G u + nu^2 G:G / 2]^(-1/2)`` (Shakib 1991, Bazilevs
    et al. 2007; ``u.G u = (2|u|/h_s)^2`` with the cell size ``h_s`` in the
    flow direction, and ``nu^2 G:G / 2 = (4 nu / h^2)^2`` on a square cell,
    the size of the smallest cell dimension on a stretched one);
    ``gamma = (h_s |u| / 2) min(1, Re_h / 3)`` with ``h_s |u| / 2 = |u|^2 / q``,
    ``q = sqrt(u.G u)``, ``Re_h = |u|^2 / (nu q)``.  Returns ``tau, gamma``
    and their derivatives with respect to the velocity components.
    """
    gxx, gxy, gyy = b.G
    gu = (gxx * upq + gxy * vpq, gxy * upq + gyy * vpq)                  # G u
    q2 = upq * gu[0] + vpq * gu[1]                                       # u . G u
    gg = gxx ** 2 + 2.0 * gxy ** 2 + gyy ** 2                            # G : G
    tau = 1.0 / np.sqrt(inv_dt2 + q2 + 0.5 * nu ** 2 * gg)
    dtau_d = [-tau ** 3 * gu[d] for d in range(2)]
    umag2 = upq ** 2 + vpq ** 2
    moving = q2 > 0
    q2s = np.where(moving, q2, 1.0)
    q = np.sqrt(q2s)
    if grad_div:
        hu2 = np.where(moving, umag2 / q, 0.0)                           # h_s |u| / 2
        re_h = hu2 / nu
        low = re_h < 3.0
        gamma = np.where(low, hu2 * re_h / 3.0, hu2)
        u_d = (upq, vpq)
        dgamma_d = [np.where(moving,
                             np.where(low,
                                      4.0 * umag2 * u_d[d] / (3.0 * nu * q2s)
                                      - 2.0 * umag2 ** 2 * gu[d] / (3.0 * nu * q2s ** 2),
                                      2.0 * u_d[d] / q - umag2 * gu[d] / q ** 3),
                             0.0) for d in range(2)]
    else:
        gamma = 0.0 * q2
        dgamma_d = [gamma, gamma]
    return tau, gamma, dtau_d, dgamma_d
