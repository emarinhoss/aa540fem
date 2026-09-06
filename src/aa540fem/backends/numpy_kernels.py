"""Reference (NumPy) element kernels of the Navier-Stokes assembly (2-D and 3-D).

Every kernel works on one cell block (:class:`aa540fem.incompressible.space._Block`)
of dimension ``d`` and returns element-local arrays in the layout
``[u_1 nodes (n), ..., u_d nodes (n), p corner nodes (nc)]``, i.e.
``L = d n + nc`` local dofs, which the assembler scatters with its
:class:`~aa540fem.backends.pattern.ScatterPlan`.  The threaded numba kernels
(:mod:`aa540fem.backends.numba_kernels`) produce the same arrays; this
module is the definition they are tested against.
"""

from __future__ import annotations

import numpy as np


# -- constant matrices ------------------------------------------------
def linear_local(b, mu_q, dmu, rho, variable_viscosity):
    """Element viscous ``K``, mass ``M``, the ``d`` divergence blocks ``B_k``
    (``q d_k u_k``, a list) and the transposed divergence ``B^T`` (gradient)
    block, each ``(ne, L, L)``.  ``dmu`` is the tuple of the ``d`` viscosity
    gradient components at the quadrature points."""
    d, n, L = b.dim, b.n, b.L
    ne = b.conn.shape[0]
    wmu = b.wh * mu_q
    Kuu = sum(np.einsum("eq,eqi,eqj->eij", wmu, b.dphi[k], b.dphi[k]) for k in range(d))
    Muu = rho * np.einsum("eq,qi,qj->eij", b.wh, b.phi, b.phi)
    Bk = [np.einsum("eq,qk,eqj->ekj", b.wh, b.psi, b.dphi[k]) for k in range(d)]   # (ne, nc, n)
    K = np.zeros((ne, L, L))
    M = np.zeros((ne, L, L))
    B = [np.zeros((ne, L, L)) for _ in range(d)]
    BT = np.zeros((ne, L, L))
    for c in range(d):
        K[:, c * n:(c + 1) * n, c * n:(c + 1) * n] = Kuu
        M[:, c * n:(c + 1) * n, c * n:(c + 1) * n] = Muu
    if variable_viscosity:
        # variable viscosity term - grad(u)^T . grad(mu):
        # block (c, k) = - int phi_i (d_k mu) (d_c phi_j)
        for c in range(d):
            for k in range(d):
                K[:, c * n:(c + 1) * n, k * n:(k + 1) * n] -= np.einsum(
                    "eq,qi,eqj->eij", b.wh * dmu[k], b.phi, b.dphi[c])
    for k in range(d):
        B[k][:, d * n:, k * n:(k + 1) * n] = Bk[k]
        BT[:, k * n:(k + 1) * n, d * n:] = np.transpose(Bk[k], (0, 2, 1))
    return K, M, B, BT


def body_load_local(b, f, rho):
    """Element load ``rho int phi_i f_c`` as ``(ne, d n)`` from the tuple ``f`` of
    the ``d`` force components at the quadrature points."""
    return np.concatenate([rho * np.einsum("eq,qi->ei", b.wh * f[c], b.phi)
                           for c in range(b.dim)], axis=1)


# -- nonlinear terms --------------------------------------------------
def momentum_local(b, U, par, old, mu_q, dmu, body, options, backend="numpy"):
    """Element-local convective and stabilisation residuals and Jacobians.

    Returns ``(N_e (ne, d n), S_e (ne, L) or None, JN_e (ne, d n, d n) or None,
    JS_e (ne, L, L) or None)``.  ``U``/``par``/``old`` are global dof vectors
    (the state, the state the stabilisation parameters are evaluated at,
    and the previous time level or ``None``); ``dmu`` the viscosity gradient
    components, ``body`` the ``d`` force components at the quadrature points
    or ``None``; ``options`` holds ``rho, stab, pspg, grad_div, metric,
    follow, jacobian, inv_dt, inv_dt2``.
    """
    if backend == "numba":
        from aa540fem.backends.numba_kernels import momentum_local as _numba_momentum

        return _numba_momentum(b, U, par, old, mu_q, dmu, body, options)
    rho = options["rho"]
    stab, pspg, jac = options["stab"], options["pspg"], options["jacobian"]
    d, n, L = b.dim, b.n, b.L
    ne = b.conn.shape[0]
    udof = [b.ldof[:, c * n:(c + 1) * n] for c in range(d)]
    pd = b.ldof[:, d * n:]
    ue = [U[udof[c]] for c in range(d)]
    uq = [ue[c] @ b.phi.T for c in range(d)]
    grad = [[np.einsum("eqi,ei->eq", b.dphi[k], ue[c]) for k in range(d)] for c in range(d)]
    wh = b.wh
    phi = b.phi[None]                                        # (1, q, n)
    dgrad = b.dphi
    ugrad = sum(uq[k][:, :, None] * dgrad[k] for k in range(d))     # u . grad phi_j

    # Galerkin convection: residual and Jacobian d/du [rho (u . grad) u]
    conv = [rho * sum(uq[k] * grad[c][k] for k in range(d)) for c in range(d)]
    N_e = np.concatenate([np.einsum("eq,qi->ei", wh * conv[c], b.phi) for c in range(d)], axis=1)
    dconv = [[rho * ((ugrad if c == k else 0.0) + grad[c][k][:, :, None] * phi)
              for k in range(d)] for c in range(d)]
    JN_e = None
    if jac:
        JN_e = np.empty((ne, d * n, d * n))
        for c in range(d):
            for k in range(d):
                JN_e[:, c * n:(c + 1) * n, k * n:(k + 1) * n] = np.einsum(
                    "eq,qi,eqj->eij", wh, b.phi, dconv[c][k])
    if not (stab or pspg):
        return N_e, None, JN_e, (np.zeros((ne, L, L)) if jac else None)

    # full momentum residual at the quadrature points and its derivatives
    nu = mu_q / rho
    lap = [np.einsum("eqi,ei->eq", b.lap_phi, ue[c]) for c in range(d)]
    # - div(mu grad u) - grad(u)^T grad(mu) = - mu lap u - grad mu . grad u - grad(u)^T grad mu
    visc = [-mu_q * lap[c] - sum(dmu[k] * grad[c][k] for k in range(d))
            - sum(dmu[k] * grad[k][c] for k in range(d)) for c in range(d)]
    pe = U[pd]
    dpsi = b.dpsi
    gradp = [np.einsum("eqk,ek->eq", dpsi[c], pe) for c in range(d)]
    R = [conv[c] + visc[c] + gradp[c] - (rho * body[c] if body is not None else 0.0)
         for c in range(d)]
    if old is not None:
        for c in range(d):
            uo = old[udof[c]] @ b.phi.T
            R[c] = R[c] + rho * (uq[c] - uo) * options["inv_dt"]

    # stabilisation parameters from the parameter state
    upq = [par[udof[c]] @ b.phi.T for c in range(d)]
    parameters = metric_parameters if options["metric"] else streamline_parameters
    tau, gamma, dtau_d, dgamma_d = parameters(b, upq, nu, options["inv_dt2"],
                                              options["grad_div"])
    follow = options["follow"]

    S_e = np.zeros((ne, L))
    w_i = tau[:, :, None] * ugrad                                 # tau (u . grad phi_i)
    div = sum(grad[c][c] for c in range(d))
    if stab:
        for c in range(d):
            S_e[:, c * n:(c + 1) * n] = (np.einsum("eq,eqi->ei", wh * R[c], w_i)
                                         + np.einsum("eq,eqi->ei", wh * gamma * div, dgrad[c]))
    gradR = None
    if pspg:
        gradR = sum(dpsi[c] * R[c][:, :, None] for c in range(d))
        S_e[:, d * n:] = np.einsum("eq,eqk->ek", wh * tau, gradR)
    if not jac:
        return N_e, S_e, None, None

    dvisc = mu_q[:, :, None] * b.lap_phi + sum(dmu[k][:, :, None] * dgrad[k] for k in range(d))
    dR = [[dconv[c][k] - (dvisc if c == k else 0.0)
           - dmu[k][:, :, None] * dgrad[c]          # d/du_{k,j} of -d_c u_k d_k mu
           for k in range(d)] for c in range(d)]
    if old is not None:
        for c in range(d):
            dR[c][c] = dR[c][c] + rho * phi * options["inv_dt"]
    JS_e = np.zeros((ne, L, L))
    if stab:
        for c in range(d):
            rows = slice(c * n, (c + 1) * n)
            for k in range(d):
                J = (np.einsum("eq,eqi,eqj->eij", wh, w_i, dR[c][k])
                     + np.einsum("eq,eqi,qj->eij", wh * tau * R[c], dgrad[k], b.phi)
                     + np.einsum("eq,eqi,eqj->eij", wh * gamma, dgrad[c], dgrad[k]))
                if follow:
                    J = J + (np.einsum("eq,eqi,qj->eij", wh * R[c] * dtau_d[k], ugrad, b.phi)
                             + np.einsum("eq,eqi,qj->eij", wh * div * dgamma_d[k],
                                         dgrad[c], b.phi))
                JS_e[:, rows, k * n:(k + 1) * n] = J
            JS_e[:, rows, d * n:] = np.einsum("eq,eqi,eqk->eik", wh, w_i, dpsi[c])
    if pspg:
        for k in range(d):
            J = sum(np.einsum("eq,eqk,eqj->ekj", wh * tau, dpsi[c], dR[c][k]) for c in range(d))
            if follow:
                J = J + np.einsum("eq,eqk,qj->ekj", wh * dtau_d[k], gradR, b.phi)
            JS_e[:, d * n:, k * n:(k + 1) * n] = J
        JS_e[:, d * n:, d * n:] = sum(np.einsum("eq,eqk,eql->ekl", wh * tau, dpsi[c], dpsi[c])
                                      for c in range(d))
    return N_e, S_e, JN_e, JS_e


# -- stabilisation parameters -----------------------------------------
def streamline_parameters(b, upq, nu, inv_dt2, grad_div):
    """Tezduyar's ``tau`` and ``gamma`` with the flow-direction element length.

    ``h = 2 / sum_i |s . grad phi_i|``, ``s = u / |u|``;
    ``tau = [(2/dt)^2 + (2|u|/h)^2 + (4 nu/h^2)^2]^(-1/2)``;
    ``gamma = h |u| / 2 min(1, Re_h / 3)``, ``Re_h = |u| h / (2 nu)``.
    ``upq`` is the list of the ``d`` velocity components at the quadrature
    points.  Returns ``tau, gamma`` and their derivatives with respect to the
    velocity components (through ``|u|`` and through ``h(s)``), each
    ``(n_elems, nq)``.
    """
    d = b.dim
    umag = np.sqrt(sum(u * u for u in upq))
    moving = umag > 0
    safe = np.where(moving, umag, 1.0)
    s = [np.where(moving, upq[k] / safe, 1.0 if k == 0 else 0.0) for k in range(d)]
    sgrad = sum(s[k][:, :, None] * b.dphi[k] for k in range(d))
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
    dir_ = [np.where(moving, upq[k] / safe, 0.0) for k in range(d)]
    sgn = np.sign(sgrad)
    dh_ds = [-0.5 * h ** 2 * np.einsum("eqi,eqi->eq", sgn, b.dphi[k]) for k in range(d)]
    inv_u = np.where(moving, 1.0 / safe, 0.0)
    dh_du = [sum(dh_ds[m] * ((1.0 if m == k else 0.0) - s[m] * s[k]) for m in range(d)) * inv_u
             for k in range(d)]
    dtau_d = [dtau_du * dir_[k] + dtau_dh * dh_du[k] for k in range(d)]
    dgamma_d = [dgamma_du * dir_[k] + dgamma_dh * dh_du[k] for k in range(d)]
    return tau, gamma, dtau_d, dgamma_d


def metric_parameters(b, upq, nu, inv_dt2, grad_div):
    """``tau`` and ``gamma`` from the element metric tensor ``G`` (smooth in ``u``).

    ``tau = [(2/dt)^2 + u.G u + nu^2 G:G / 2]^(-1/2)`` (Shakib 1991, Bazilevs
    et al. 2007; ``u.G u = (2|u|/h_s)^2`` with the cell size ``h_s`` in the
    flow direction, and ``nu^2 G:G / 2 = (4 nu / h^2)^2`` on a square cell,
    the size of the smallest cell dimension on a stretched one);
    ``gamma = (h_s |u| / 2) min(1, Re_h / 3)`` with ``h_s |u| / 2 = |u|^2 / q``,
    ``q = sqrt(u.G u)``, ``Re_h = |u|^2 / (nu q)``.  Returns ``tau, gamma``
    and their derivatives with respect to the velocity components.
    """
    d = b.dim
    G = b.Gmat
    gu = [sum(G[..., c, k] * upq[k] for k in range(d)) for c in range(d)]   # G u
    q2 = sum(upq[c] * gu[c] for c in range(d))                              # u . G u
    gg = np.einsum("...ij,...ij->...", G, G)                                # G : G
    tau = 1.0 / np.sqrt(inv_dt2 + q2 + 0.5 * nu ** 2 * gg)
    dtau_d = [-tau ** 3 * gu[c] for c in range(d)]
    umag2 = sum(u * u for u in upq)
    moving = q2 > 0
    q2s = np.where(moving, q2, 1.0)
    q = np.sqrt(q2s)
    if grad_div:
        hu2 = np.where(moving, umag2 / q, 0.0)                           # h_s |u| / 2
        re_h = hu2 / nu
        low = re_h < 3.0
        gamma = np.where(low, hu2 * re_h / 3.0, hu2)
        dgamma_d = [np.where(moving,
                             np.where(low,
                                      4.0 * umag2 * upq[c] / (3.0 * nu * q2s)
                                      - 2.0 * umag2 ** 2 * gu[c] / (3.0 * nu * q2s ** 2),
                                      2.0 * upq[c] / q - umag2 * gu[c] / q ** 3),
                             0.0) for c in range(d)]
    else:
        gamma = 0.0 * q2
        dgamma_d = [gamma] * d
    return tau, gamma, dtau_d, dgamma_d
