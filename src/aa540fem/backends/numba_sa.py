"""Threaded (numba) element kernel of the Spalart-Allmaras equation (2-D and 3-D).

Reproduces ``SpalartAllmarasSolver._local_numpy`` (residual ``(ne, n)`` and
Jacobian ``(ne, n, n)`` per element): Galerkin convection, diffusion and
source, SUPG on convection and source, source and diffusivity derivatives
by the same central finite differences at the quadrature points.  The model
functions are scalar ``@njit`` translations of :class:`SpalartAllmaras`;
its constants travel as a float array (see :func:`constants`).  Gradients
are passed stacked over the ``d`` directions (``dphi (d, ne, nq, n)``) and
the velocity as ``(d, ne, n)``.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

# positions in the constants array
(NU, CB1, SIGMA, CB2, KAPPA, CW1, CW2, CW3, CV1, CT3, CT4, CV2, CV3, CN1, R_LIMIT, FT2,
 FD_STEP) = range(17)


def constants(model, fd_step: float) -> np.ndarray:
    return np.array([model.nu, model.cb1, model.sigma, model.cb2, model.kappa, model.cw1,
                     model.cw2, model.cw3, model.cv1, model.ct3, model.ct4, model.cv2,
                     model.cv3, model.cn1, model.r_limit, 1.0 if model.ft2 else 0.0, fd_step])


@njit(cache=True, nogil=True)
def _diffusivity(nt, c):
    chi = nt / c[NU]
    chi3 = chi * chi * chi
    if nt >= 0.0:
        fn = 1.0
    else:
        fn = (c[CN1] + chi3) / (c[CN1] - chi3)
    return (c[NU] + fn * nt) / c[SIGMA]


@njit(cache=True, nogil=True)
def _source(nt, grad2, omega, d, c):
    """Source at a point from ``nu_tilde``, ``|grad nu_tilde|^2``, ``|omega|`` and ``d``."""
    if d < 1e-12:
        d = 1e-12
    nu = c[NU]
    chi = nt / nu
    chi3 = chi * chi * chi
    fv1 = chi3 / (chi3 + c[CV1] ** 3)
    kappa2 = c[KAPPA] * c[KAPPA]
    d2 = d * d
    # modified vorticity with the SA-neg guard
    fv2 = 1.0 - chi / (1.0 + chi * fv1)
    sbar = nt * fv2 / (kappa2 * d2)
    if sbar < -c[CV2] * omega:
        safe = (c[CV3] - 2.0 * c[CV2]) * omega - sbar
        if safe == 0.0:
            safe = 1e-300
        s_tilde = omega + omega * (c[CV2] * c[CV2] * omega + c[CV3] * sbar) / safe
    else:
        s_tilde = omega + sbar
    if s_tilde < 1e-300:
        s_tilde = 1e-300
    r = nt / (s_tilde * kappa2 * d2)
    if r > c[R_LIMIT]:
        r = c[R_LIMIT]
    elif r < -c[R_LIMIT]:
        r = -c[R_LIMIT]
    g = r + c[CW2] * (r ** 6 - r)
    cw36 = c[CW3] ** 6
    fw = g * ((1.0 + cw36) / (g ** 6 + cw36)) ** (1.0 / 6.0)
    ft2 = c[CT3] * np.exp(-c[CT4] * chi * chi) if c[FT2] != 0.0 else 0.0
    if nt >= 0.0:
        prod = c[CB1] * (1.0 - ft2) * s_tilde * nt
        destr = (c[CW1] * fw - c[CB1] / kappa2 * ft2) * (nt / d) ** 2
        s = prod - destr
    else:
        s = c[CB1] * (1.0 - c[CT3]) * omega * nt + c[CW1] * (nt / d) ** 2
    return s + c[CB2] / c[SIGMA] * grad2


@njit(cache=True, nogil=True)
def _grad2(g, d):
    s = 0.0
    for k in range(d):
        s += g[k] * g[k]
    return s


@njit(parallel=True, cache=True, nogil=True)
def sa_kernel(phi, dphi, psi, wh, G, nte, ue, de, c, supg, metric, R_e, J_e):
    d, ne, nq, n = dphi.shape
    nc = psi.shape[1]
    fd = c[FD_STEP]
    for e in prange(ne):
        ugrad = np.empty(n)
        grad_i = np.empty(n)
        dS = np.empty(n)
        g = np.empty(d)
        gp = np.empty(d)
        uq = np.empty(d)
        s_g = np.empty(d)
        vgrad = np.empty((d, d))                 # vgrad[c, k] = d u_c / d x_k
        for q in range(nq):
            w = wh[e, q]
            ntq = 0.0
            for k in range(d):
                g[k] = 0.0
                uq[k] = 0.0
                for m in range(d):
                    vgrad[k, m] = 0.0
            for i in range(n):
                ntq += phi[q, i] * nte[e, i]
                for k in range(d):
                    g[k] += dphi[k, e, q, i] * nte[e, i]
                    uq[k] += phi[q, i] * ue[k, e, i]
                    for m in range(d):
                        vgrad[k, m] += dphi[m, e, q, i] * ue[k, e, i]
            if d == 2:
                omega = abs(vgrad[1, 0] - vgrad[0, 1])
            else:
                cx = vgrad[2, 1] - vgrad[1, 2]
                cy = vgrad[0, 2] - vgrad[2, 0]
                cz = vgrad[1, 0] - vgrad[0, 1]
                omega = np.sqrt(cx * cx + cy * cy + cz * cz)
            dq = 0.0
            for k in range(nc):
                dq += psi[q, k] * de[e, k]
            if dq < 1e-12:
                dq = 1e-12
            # diffusivity, source and their finite-difference derivatives
            D = _diffusivity(ntq, c)
            h_n = fd * (c[NU] + abs(ntq))
            dD = (_diffusivity(ntq + h_n, c) - _diffusivity(ntq - h_n, c)) / (2.0 * h_n)
            grad2 = _grad2(g, d)
            s = _source(ntq, grad2, omega, dq, c)
            s_n = (_source(ntq + h_n, grad2, omega, dq, c)
                   - _source(ntq - h_n, grad2, omega, dq, c)) / (2.0 * h_n)
            h_g = fd
            for k in range(d):
                h_g += fd * abs(g[k])
            for k in range(d):
                for m in range(d):
                    gp[m] = g[m]
                gp[k] = g[k] + h_g
                plus = _source(ntq, _grad2(gp, d), omega, dq, c)
                gp[k] = g[k] - h_g
                minus = _source(ntq, _grad2(gp, d), omega, dq, c)
                s_g[k] = (plus - minus) / (2.0 * h_g)
            conv = 0.0
            for k in range(d):
                conv += uq[k] * g[k]
            for i in range(n):
                a = 0.0
                b = 0.0
                cdS = s_n * phi[q, i]
                for k in range(d):
                    a += uq[k] * dphi[k, e, q, i]
                    b += dphi[k, e, q, i] * g[k]
                    cdS += s_g[k] * dphi[k, e, q, i]
                ugrad[i] = a
                grad_i[i] = b
                dS[i] = cdS
            tau = 0.0
            if supg:
                if metric:
                    q2 = 0.0
                    gg = 0.0
                    for a_ in range(d):
                        for b_ in range(d):
                            q2 += uq[a_] * G[e, q, a_, b_] * uq[b_]
                            gg += G[e, q, a_, b_] * G[e, q, a_, b_]
                    tau = 1.0 / np.sqrt(q2 + 0.5 * D * D * gg)
                else:
                    umag2 = 0.0
                    for k in range(d):
                        umag2 += uq[k] * uq[k]
                    umag = np.sqrt(umag2)
                    ssum = 0.0
                    for i in range(n):
                        sg = 0.0
                        if umag > 0.0:
                            for k in range(d):
                                sg += uq[k] / umag * dphi[k, e, q, i]
                        else:
                            sg = dphi[0, e, q, i]
                        ssum += abs(sg)
                    if ssum < 1e-300:
                        ssum = 1e-300
                    h = 2.0 / ssum
                    tau = 1.0 / np.sqrt((2.0 * umag / h) ** 2 + (4.0 * D / (h * h)) ** 2)
            for i in range(n):
                wi = phi[q, i] + tau * ugrad[i]          # Galerkin + SUPG test function
                R_e[e, i] += w * ((conv - s) * wi + D * grad_i[i])
                for j in range(n):
                    diff = 0.0
                    for k in range(d):
                        diff += dphi[k, e, q, i] * dphi[k, e, q, j]
                    J_e[e, i, j] += w * (wi * (ugrad[j] - dS[j]) + D * diff
                                         + dD * grad_i[i] * phi[q, j])


def sa_local(b, nt, vel, distance, c, supg, metric):
    """Element residual ``(ne, n)`` and Jacobian ``(ne, n, n)`` with the numba kernel;
    ``vel`` is the ``(d, N)`` nodal velocity."""
    ne, n = b.conn.shape
    R_e = np.zeros((ne, n))
    J_e = np.zeros((ne, n, n))
    cache = b.__dict__.get("_numba_stacked")
    if cache is None:
        cache = (np.ascontiguousarray(np.stack(b.dphi)), np.ascontiguousarray(np.stack(b.dpsi)),
                 np.ascontiguousarray(b.Gmat))
        b.__dict__["_numba_stacked"] = cache
    dphi, _, G = cache
    sa_kernel(b.phi, dphi, b.psi, b.wh, G, np.ascontiguousarray(nt[b.conn]),
              np.ascontiguousarray(np.asarray(vel)[:, b.conn]),
              np.ascontiguousarray(distance[b.pconn]), c, bool(supg), bool(metric), R_e, J_e)
    return R_e, J_e
