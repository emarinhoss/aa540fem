"""Threaded (numba) element kernels of the Navier-Stokes assembly.

Same element-local outputs as :mod:`aa540fem.backends.numpy_kernels`
(``N_e (ne, 2n)``, ``S_e (ne, L)``, ``JN_e (ne, 2n, 2n)``, ``JS_e (ne, L, L)``
in the layout ``[u_x nodes, u_y nodes, p corner nodes]``), computed with one
explicit loop nest per element, the element loop in ``prange``.  Threads
write disjoint element slices only; the global scatter is the deterministic
:func:`segmented_reduce`, so results do not depend on the thread count.

The kernels compile on first use (a few seconds) and are cached on disk
(``NUMBA_CACHE_DIR`` or the package directory); :func:`warmup` compiles
them on a two-element mesh.
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange

# flag indices of the ``flags`` array
STAB, PSPG, GRAD_DIV, METRIC, FOLLOW, TRANSIENT, BODY, JACOBIAN = range(8)


# -- global scatter ---------------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def segmented_reduce(vals, perm, seg_start, out):
    """``out[k] = sum(vals[perm[seg_start[k]:seg_start[k + 1]]])`` in a fixed order."""
    for k in prange(out.shape[0]):
        s = 0.0
        for j in range(seg_start[k], seg_start[k + 1]):
            s += vals[perm[j]]
        out[k] = s


# -- stabilisation parameters at one quadrature point -----------------
@njit(cache=True, nogil=True)
def _params_metric(upq, vpq, nu, inv_dt2, grad_div, gxx, gxy, gyy):
    gux = gxx * upq + gxy * vpq
    guy = gxy * upq + gyy * vpq
    q2 = upq * gux + vpq * guy
    gg = gxx * gxx + 2.0 * gxy * gxy + gyy * gyy
    tau = 1.0 / np.sqrt(inv_dt2 + q2 + 0.5 * nu * nu * gg)
    tau3 = tau * tau * tau
    dtau_x = -tau3 * gux
    dtau_y = -tau3 * guy
    gamma = 0.0
    dgamma_x = 0.0
    dgamma_y = 0.0
    if grad_div and q2 > 0.0:
        umag2 = upq * upq + vpq * vpq
        q = np.sqrt(q2)
        hu2 = umag2 / q
        re_h = hu2 / nu
        if re_h < 3.0:
            gamma = hu2 * re_h / 3.0
            c1 = 4.0 * umag2 / (3.0 * nu * q2)
            c2 = 2.0 * umag2 * umag2 / (3.0 * nu * q2 * q2)
            dgamma_x = c1 * upq - c2 * gux
            dgamma_y = c1 * vpq - c2 * guy
        else:
            gamma = hu2
            q3 = q * q * q
            dgamma_x = 2.0 * upq / q - umag2 * gux / q3
            dgamma_y = 2.0 * vpq / q - umag2 * guy / q3
    return tau, gamma, dtau_x, dtau_y, dgamma_x, dgamma_y


@njit(cache=True, nogil=True)
def _params_streamline(upq, vpq, nu, inv_dt2, grad_div, dphi_dx_q, dphi_dy_q):
    n = dphi_dx_q.shape[0]
    umag = np.hypot(upq, vpq)
    moving = umag > 0.0
    if moving:
        sx = upq / umag
        sy = vpq / umag
        inv_u = 1.0 / umag
    else:
        sx = 1.0
        sy = 0.0
        inv_u = 0.0
    ssum = 0.0
    sgx = 0.0
    sgy = 0.0
    for i in range(n):
        sg = sx * dphi_dx_q[i] + sy * dphi_dy_q[i]
        ssum += abs(sg)
        if sg > 0.0:
            sgx += dphi_dx_q[i]
            sgy += dphi_dy_q[i]
        elif sg < 0.0:
            sgx -= dphi_dx_q[i]
            sgy -= dphi_dy_q[i]
    if ssum < 1e-300:
        ssum = 1e-300
    h = 2.0 / ssum
    h2 = h * h
    tau = 1.0 / np.sqrt(inv_dt2 + (2.0 * umag / h) ** 2 + (4.0 * nu / h2) ** 2)
    re_h = umag * h / (2.0 * nu)
    low = re_h < 3.0
    if grad_div:
        if low:
            gamma = 0.5 * h * umag * (re_h / 3.0)
            dgamma_du = 0.5 * h * (2.0 * re_h / 3.0)
            dgamma_dh = h * umag * umag / (6.0 * nu)
        else:
            gamma = 0.5 * h * umag
            dgamma_du = 0.5 * h
            dgamma_dh = 0.5 * umag
    else:
        gamma = 0.0
        dgamma_du = 0.0
        dgamma_dh = 0.0
    tau3 = tau * tau * tau
    dtau_du = -(4.0 / h2) * umag * tau3
    dtau_dh = tau3 * (4.0 * umag * umag / (h2 * h) + 32.0 * nu * nu / (h2 * h2 * h))
    dh_dsx = -0.5 * h2 * sgx
    dh_dsy = -0.5 * h2 * sgy
    dh_dux = (dh_dsx * (1.0 - sx * sx) - dh_dsy * sy * sx) * inv_u
    dh_duy = (-dh_dsx * sx * sy + dh_dsy * (1.0 - sy * sy)) * inv_u
    dir_x = sx * (1.0 if moving else 0.0)
    dir_y = sy * (1.0 if moving else 0.0)
    dtau_x = dtau_du * dir_x + dtau_dh * dh_dux
    dtau_y = dtau_du * dir_y + dtau_dh * dh_duy
    dgamma_x = dgamma_du * dir_x + dgamma_dh * dh_dux
    dgamma_y = dgamma_du * dir_y + dgamma_dh * dh_duy
    return tau, gamma, dtau_x, dtau_y, dgamma_x, dgamma_y


# -- momentum equations -----------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def momentum_kernel(phi, dphi_dx, dphi_dy, lap_phi, psi, dpsi_dx, dpsi_dy, wh, gxx, gxy, gyy,
                    ue, ve, pe, uoe, voe, upe, vpe, mu_q, dmu_dx, dmu_dy, fx, fy,
                    rho, inv_dt, inv_dt2, flags, N_e, S_e, JN_e, JS_e):
    ne, nq, n = dphi_dx.shape
    nc = psi.shape[1]
    stab = flags[STAB] != 0
    pspg = flags[PSPG] != 0
    grad_div = flags[GRAD_DIV] != 0
    metric = flags[METRIC] != 0
    follow = flags[FOLLOW] != 0
    transient = flags[TRANSIENT] != 0
    body = flags[BODY] != 0
    jac = flags[JACOBIAN] != 0
    stabilised = stab or pspg
    for e in prange(ne):
        ugrad = np.empty(n)
        dc00 = np.empty(n)
        dc01 = np.empty(n)
        dc10 = np.empty(n)
        dc11 = np.empty(n)
        dr00 = np.empty(n)
        dr01 = np.empty(n)
        dr10 = np.empty(n)
        dr11 = np.empty(n)
        for q in range(nq):
            w = wh[e, q]
            uq = 0.0
            vq = 0.0
            dudx = 0.0
            dudy = 0.0
            dvdx = 0.0
            dvdy = 0.0
            lapu = 0.0
            lapv = 0.0
            for i in range(n):
                uq += phi[q, i] * ue[e, i]
                vq += phi[q, i] * ve[e, i]
                dudx += dphi_dx[e, q, i] * ue[e, i]
                dudy += dphi_dy[e, q, i] * ue[e, i]
                dvdx += dphi_dx[e, q, i] * ve[e, i]
                dvdy += dphi_dy[e, q, i] * ve[e, i]
                lapu += lap_phi[e, q, i] * ue[e, i]
                lapv += lap_phi[e, q, i] * ve[e, i]
            for i in range(n):
                ugrad[i] = uq * dphi_dx[e, q, i] + vq * dphi_dy[e, q, i]
            conv0 = rho * (uq * dudx + vq * dudy)
            conv1 = rho * (uq * dvdx + vq * dvdy)
            for i in range(n):
                N_e[e, i] += w * conv0 * phi[q, i]
                N_e[e, n + i] += w * conv1 * phi[q, i]
            for j in range(n):
                dc00[j] = rho * (ugrad[j] + dudx * phi[q, j])
                dc01[j] = rho * dudy * phi[q, j]
                dc10[j] = rho * dvdx * phi[q, j]
                dc11[j] = rho * (ugrad[j] + dvdy * phi[q, j])
            if jac:
                for i in range(n):
                    wp = w * phi[q, i]
                    for j in range(n):
                        JN_e[e, i, j] += wp * dc00[j]
                        JN_e[e, i, n + j] += wp * dc01[j]
                        JN_e[e, n + i, j] += wp * dc10[j]
                        JN_e[e, n + i, n + j] += wp * dc11[j]
            if not stabilised:
                continue

            # strong momentum residual and its derivatives
            mu = mu_q[e, q]
            dmx = dmu_dx[e, q]
            dmy = dmu_dy[e, q]
            nu = mu / rho
            visc0 = -mu * lapu - (dmx * dudx + dmy * dudy) - (dmx * dudx + dmy * dvdx)
            visc1 = -mu * lapv - (dmx * dvdx + dmy * dvdy) - (dmx * dudy + dmy * dvdy)
            gpx = 0.0
            gpy = 0.0
            for k in range(nc):
                gpx += dpsi_dx[e, q, k] * pe[e, k]
                gpy += dpsi_dy[e, q, k] * pe[e, k]
            R0 = conv0 + visc0 + gpx
            R1 = conv1 + visc1 + gpy
            if body:
                R0 -= rho * fx[e, q]
                R1 -= rho * fy[e, q]
            if transient:
                uo = 0.0
                vo = 0.0
                for i in range(n):
                    uo += phi[q, i] * uoe[e, i]
                    vo += phi[q, i] * voe[e, i]
                R0 += rho * (uq - uo) * inv_dt
                R1 += rho * (vq - vo) * inv_dt
            # stabilisation parameters at the parameter state
            upq = 0.0
            vpq = 0.0
            for i in range(n):
                upq += phi[q, i] * upe[e, i]
                vpq += phi[q, i] * vpe[e, i]
            if metric:
                tau, gamma, dtau_x, dtau_y, dgamma_x, dgamma_y = _params_metric(
                    upq, vpq, nu, inv_dt2, grad_div, gxx[e, q], gxy[e, q], gyy[e, q])
            else:
                tau, gamma, dtau_x, dtau_y, dgamma_x, dgamma_y = _params_streamline(
                    upq, vpq, nu, inv_dt2, grad_div, dphi_dx[e, q], dphi_dy[e, q])
            div = dudx + dvdy
            if stab:
                for i in range(n):
                    wi = tau * ugrad[i]
                    S_e[e, i] += w * (R0 * wi + gamma * div * dphi_dx[e, q, i])
                    S_e[e, n + i] += w * (R1 * wi + gamma * div * dphi_dy[e, q, i])
            if pspg:
                for k in range(nc):
                    gradR = dpsi_dx[e, q, k] * R0 + dpsi_dy[e, q, k] * R1
                    S_e[e, 2 * n + k] += w * tau * gradR
            if not jac:
                continue
            for j in range(n):
                dvisc = (mu * lap_phi[e, q, j] + dmx * dphi_dx[e, q, j]
                         + dmy * dphi_dy[e, q, j])
                dr00[j] = dc00[j] - dvisc - dmx * dphi_dx[e, q, j]
                dr01[j] = dc01[j] - dmy * dphi_dx[e, q, j]
                dr10[j] = dc10[j] - dmx * dphi_dy[e, q, j]
                dr11[j] = dc11[j] - dvisc - dmy * dphi_dy[e, q, j]
                if transient:
                    dr00[j] += rho * phi[q, j] * inv_dt
                    dr11[j] += rho * phi[q, j] * inv_dt
            if stab:
                for i in range(n):
                    wi = tau * ugrad[i]
                    gx = dphi_dx[e, q, i]
                    gy = dphi_dy[e, q, i]
                    for j in range(n):
                        pj = phi[q, j]
                        gxj = dphi_dx[e, q, j]
                        gyj = dphi_dy[e, q, j]
                        # row u_x, i: tau (u.grad phi_i) dR_0/du_d + tau R_0 d_d phi_i phi_j
                        #             + gamma d_x phi_i d_d phi_j
                        j00 = wi * dr00[j] + tau * R0 * gx * pj + gamma * gx * gxj
                        j01 = wi * dr01[j] + tau * R0 * gy * pj + gamma * gx * gyj
                        # row u_y, i
                        j10 = wi * dr10[j] + tau * R1 * gx * pj + gamma * gy * gxj
                        j11 = wi * dr11[j] + tau * R1 * gy * pj + gamma * gy * gyj
                        if follow:
                            j00 += (R0 * dtau_x * ugrad[i] + div * dgamma_x * gx) * pj
                            j01 += (R0 * dtau_y * ugrad[i] + div * dgamma_y * gx) * pj
                            j10 += (R1 * dtau_x * ugrad[i] + div * dgamma_x * gy) * pj
                            j11 += (R1 * dtau_y * ugrad[i] + div * dgamma_y * gy) * pj
                        JS_e[e, i, j] += w * j00
                        JS_e[e, i, n + j] += w * j01
                        JS_e[e, n + i, j] += w * j10
                        JS_e[e, n + i, n + j] += w * j11
                    for k in range(nc):
                        JS_e[e, i, 2 * n + k] += w * wi * dpsi_dx[e, q, k]
                        JS_e[e, n + i, 2 * n + k] += w * wi * dpsi_dy[e, q, k]
            if pspg:
                for k in range(nc):
                    px = dpsi_dx[e, q, k]
                    py = dpsi_dy[e, q, k]
                    gradR = px * R0 + py * R1
                    for j in range(n):
                        pj = phi[q, j]
                        jx = tau * (px * dr00[j] + py * dr10[j])
                        jy = tau * (px * dr01[j] + py * dr11[j])
                        if follow:
                            jx += dtau_x * gradR * pj
                            jy += dtau_y * gradR * pj
                        JS_e[e, 2 * n + k, j] += w * jx
                        JS_e[e, 2 * n + k, n + j] += w * jy
                    for m in range(nc):
                        JS_e[e, 2 * n + k, 2 * n + m] += w * tau * (
                            px * dpsi_dx[e, q, m] + py * dpsi_dy[e, q, m])


def momentum_local(b, U, par, old, mu_q, dmu_dx, dmu_dy, body, options):
    """numba version of :func:`aa540fem.backends.numpy_kernels.momentum_local`."""
    n, L = b.n, b.L
    ne = b.conn.shape[0]
    ux, uy, pd = b.ldof[:, :n], b.ldof[:, n:2 * n], b.ldof[:, 2 * n:]
    stab, pspg, jac = options["stab"], options["pspg"], options["jacobian"]
    stabilised = stab or pspg
    ue = np.ascontiguousarray(U[ux])
    ve = np.ascontiguousarray(U[uy])
    pe = np.ascontiguousarray(U[pd])
    upe = ue if par is U else np.ascontiguousarray(par[ux])
    vpe = ve if par is U else np.ascontiguousarray(par[uy])
    transient = old is not None
    uoe = np.ascontiguousarray(old[ux]) if transient else ue
    voe = np.ascontiguousarray(old[uy]) if transient else ve
    if body is not None:
        fx = np.ascontiguousarray(body[0], dtype=float)
        fy = np.ascontiguousarray(body[1], dtype=float)
    else:
        fx = fy = mu_q
    flags = np.array([stab, pspg, options["grad_div"], options["metric"], options["follow"],
                      transient, body is not None, jac], dtype=np.int64)
    N_e = np.zeros((ne, 2 * n))
    S_e = np.zeros((ne, L)) if stabilised else np.zeros((1, L))
    JN_e = np.zeros((ne, 2 * n, 2 * n)) if jac else np.zeros((1, 2 * n, 2 * n))
    JS_e = np.zeros((ne, L, L)) if (jac and stabilised) else np.zeros((1, L, L))
    gxx, gxy, gyy = b.G
    momentum_kernel(b.phi, b.dphi_dx, b.dphi_dy, b.lap_phi, b.psi, b.dpsi_dx, b.dpsi_dy, b.wh,
                    gxx, gxy, gyy, ue, ve, pe, uoe, voe, upe, vpe,
                    np.ascontiguousarray(mu_q), np.ascontiguousarray(dmu_dx),
                    np.ascontiguousarray(dmu_dy), fx, fy, float(options["rho"]),
                    float(options["inv_dt"]), float(options["inv_dt2"]), flags,
                    N_e, S_e, JN_e, JS_e)
    if not stabilised:
        return N_e, None, (JN_e if jac else None), (np.zeros((ne, L, L)) if jac else None)
    if not jac:
        return N_e, S_e, None, None
    return N_e, S_e, JN_e, JS_e


def warmup(threads: int | None = None):
    """Compile the kernels on a two-element mesh (a few seconds once per machine)."""
    import numba

    from aa540fem.core.mesh import geometry
    from aa540fem.incompressible.assembler import FlowAssembler
    from aa540fem.incompressible.problem import FlowProblem

    if threads:
        numba.set_num_threads(threads)
    for elem in ("quad9", "triangle6"):
        mesh = geometry(1.0, 1.0, 1, elem)
        for stabilisation, pspg in ((False, False), (True, False), (True, True)):
            prob = FlowProblem(mesh, mu=0.1, rho=1.0, stabilisation=stabilisation, pspg=pspg,
                               bc={"left": (1.0, 0.0)}, eddy_viscosity=np.full(mesh.n_nodes, 0.01),
                               body_force=lambda x, y: (0.0 * x, 0.0 * y))
            asm = FlowAssembler(prob, backend="numba")
            U = np.linspace(0.0, 1.0, asm.space.ndof)
            asm.momentum_terms(U)
            asm.momentum_terms(U, dt=0.1, U_old=0.5 * U, param_state=U, jacobian=False)
    prob.element_length = "streamline"
    FlowAssembler(prob, backend="numba").momentum_terms(U)
