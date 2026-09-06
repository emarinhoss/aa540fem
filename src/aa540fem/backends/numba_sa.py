"""Threaded (numba) element kernel of the Spalart-Allmaras equation.

Reproduces ``SpalartAllmarasSolver._local_numpy`` (residual ``(ne, n)`` and
Jacobian ``(ne, n, n)`` per element): Galerkin convection, diffusion and
source, SUPG on convection and source, source and diffusivity derivatives
by the same central finite differences at the quadrature points.  The model
functions are scalar ``@njit`` translations of :class:`SpalartAllmaras`;
its constants travel as a float array (see :func:`constants`).
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
def _source(nt, gx, gy, omega, d, c):
    if d < 1e-12:
        d = 1e-12
    nu = c[NU]
    chi = nt / nu
    chi3 = chi * chi * chi
    fv1 = chi3 / (chi3 + c[CV1] ** 3)
    grad2 = gx * gx + gy * gy
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


@njit(parallel=True, cache=True, nogil=True)
def sa_kernel(phi, dphi_dx, dphi_dy, psi, wh, gxx, gxy, gyy, nte, ue, ve, de, c,
              supg, metric, R_e, J_e):
    ne, nq, n = dphi_dx.shape
    nc = psi.shape[1]
    fd = c[FD_STEP]
    for e in prange(ne):
        ugrad = np.empty(n)
        grad_i = np.empty(n)
        dS = np.empty(n)
        for q in range(nq):
            w = wh[e, q]
            ntq = 0.0
            gx = 0.0
            gy = 0.0
            uq = 0.0
            vq = 0.0
            dudy = 0.0
            dvdx = 0.0
            for i in range(n):
                ntq += phi[q, i] * nte[e, i]
                gx += dphi_dx[e, q, i] * nte[e, i]
                gy += dphi_dy[e, q, i] * nte[e, i]
                uq += phi[q, i] * ue[e, i]
                vq += phi[q, i] * ve[e, i]
                dudy += dphi_dy[e, q, i] * ue[e, i]
                dvdx += dphi_dx[e, q, i] * ve[e, i]
            omega = abs(dvdx - dudy)
            dq = 0.0
            for k in range(nc):
                dq += psi[q, k] * de[e, k]
            if dq < 1e-12:
                dq = 1e-12
            # diffusivity, source and their finite-difference derivatives
            D = _diffusivity(ntq, c)
            h_n = fd * (c[NU] + abs(ntq))
            dD = (_diffusivity(ntq + h_n, c) - _diffusivity(ntq - h_n, c)) / (2.0 * h_n)
            s = _source(ntq, gx, gy, omega, dq, c)
            s_n = (_source(ntq + h_n, gx, gy, omega, dq, c)
                   - _source(ntq - h_n, gx, gy, omega, dq, c)) / (2.0 * h_n)
            h_g = fd * (1.0 + abs(gx) + abs(gy))
            s_gx = (_source(ntq, gx + h_g, gy, omega, dq, c)
                    - _source(ntq, gx - h_g, gy, omega, dq, c)) / (2.0 * h_g)
            s_gy = (_source(ntq, gx, gy + h_g, omega, dq, c)
                    - _source(ntq, gx, gy - h_g, omega, dq, c)) / (2.0 * h_g)
            conv = uq * gx + vq * gy
            for i in range(n):
                ugrad[i] = uq * dphi_dx[e, q, i] + vq * dphi_dy[e, q, i]
                grad_i[i] = dphi_dx[e, q, i] * gx + dphi_dy[e, q, i] * gy
                dS[i] = s_n * phi[q, i] + s_gx * dphi_dx[e, q, i] + s_gy * dphi_dy[e, q, i]
            tau = 0.0
            if supg:
                if metric:
                    q2 = gxx[e, q] * uq * uq + 2.0 * gxy[e, q] * uq * vq + gyy[e, q] * vq * vq
                    gg = gxx[e, q] ** 2 + 2.0 * gxy[e, q] ** 2 + gyy[e, q] ** 2
                    tau = 1.0 / np.sqrt(q2 + 0.5 * D * D * gg)
                else:
                    umag = np.hypot(uq, vq)
                    if umag > 0.0:
                        sx = uq / umag
                        sy = vq / umag
                    else:
                        sx = 1.0
                        sy = 0.0
                    ssum = 0.0
                    for i in range(n):
                        ssum += abs(sx * dphi_dx[e, q, i] + sy * dphi_dy[e, q, i])
                    if ssum < 1e-300:
                        ssum = 1e-300
                    h = 2.0 / ssum
                    tau = 1.0 / np.sqrt((2.0 * umag / h) ** 2 + (4.0 * D / (h * h)) ** 2)
            for i in range(n):
                wi = phi[q, i] + tau * ugrad[i]          # Galerkin + SUPG test function
                R_e[e, i] += w * ((conv - s) * wi + D * grad_i[i])
                for j in range(n):
                    J_e[e, i, j] += w * (wi * (ugrad[j] - dS[j])
                                         + D * (dphi_dx[e, q, i] * dphi_dx[e, q, j]
                                                + dphi_dy[e, q, i] * dphi_dy[e, q, j])
                                         + dD * grad_i[i] * phi[q, j])


def sa_local(b, nt, u, v, distance, c, supg, metric):
    """Element residual ``(ne, n)`` and Jacobian ``(ne, n, n)`` with the numba kernel."""
    ne, n = b.conn.shape
    R_e = np.zeros((ne, n))
    J_e = np.zeros((ne, n, n))
    gxx, gxy, gyy = b.G
    sa_kernel(b.phi, b.dphi_dx, b.dphi_dy, b.psi, b.wh, gxx, gxy, gyy,
              np.ascontiguousarray(nt[b.conn]), np.ascontiguousarray(u[b.conn]),
              np.ascontiguousarray(v[b.conn]), np.ascontiguousarray(distance[b.pconn]), c,
              bool(supg), bool(metric), R_e, J_e)
    return R_e, J_e
