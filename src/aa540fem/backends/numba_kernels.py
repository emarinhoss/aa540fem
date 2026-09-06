"""Threaded (numba) element kernels of the Navier-Stokes assembly (2-D and 3-D).

Same element-local outputs as :mod:`aa540fem.backends.numpy_kernels`
(``N_e (ne, d n)``, ``S_e (ne, L)``, ``JN_e (ne, d n, d n)``, ``JS_e (ne, L, L)``
in the layout ``[u_1 nodes, ..., u_d nodes, p corner nodes]``), computed
with one explicit loop nest per element, the element loop in ``prange``.
Threads write disjoint element slices only; the global scatter is the
deterministic :func:`segmented_reduce`, so results do not depend on the
thread count.  Gradients are passed stacked over the ``d`` directions
(``dphi (d, ne, nq, n)``) and the velocity components as ``(d, ne, n)``.

The kernels compile on first use (a few seconds) and are cached on disk
(``NUMBA_CACHE_DIR`` or the package directory); :func:`warmup` compiles
them on small meshes.
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
def _params_metric(upq, nu, inv_dt2, grad_div, G, dtau, dgamma):
    """``tau, gamma`` from the metric tensor ``G (d, d)``; the derivatives with
    respect to the velocity components are written into ``dtau``, ``dgamma``."""
    d = upq.shape[0]
    gu = np.empty(d)
    q2 = 0.0
    gg = 0.0
    umag2 = 0.0
    for c in range(d):
        s = 0.0
        for k in range(d):
            s += G[c, k] * upq[k]
            gg += G[c, k] * G[c, k]
        gu[c] = s
        q2 += upq[c] * s
        umag2 += upq[c] * upq[c]
    tau = 1.0 / np.sqrt(inv_dt2 + q2 + 0.5 * nu * nu * gg)
    tau3 = tau * tau * tau
    for c in range(d):
        dtau[c] = -tau3 * gu[c]
        dgamma[c] = 0.0
    gamma = 0.0
    if grad_div and q2 > 0.0:
        q = np.sqrt(q2)
        hu2 = umag2 / q
        re_h = hu2 / nu
        if re_h < 3.0:
            gamma = hu2 * re_h / 3.0
            c1 = 4.0 * umag2 / (3.0 * nu * q2)
            c2 = 2.0 * umag2 * umag2 / (3.0 * nu * q2 * q2)
            for c in range(d):
                dgamma[c] = c1 * upq[c] - c2 * gu[c]
        else:
            gamma = hu2
            q3 = q * q * q
            for c in range(d):
                dgamma[c] = 2.0 * upq[c] / q - umag2 * gu[c] / q3
    return tau, gamma


@njit(cache=True, nogil=True)
def _params_streamline(upq, nu, inv_dt2, grad_div, dphi_q, dtau, dgamma):
    """Tezduyar's ``tau, gamma`` with the flow-direction length; ``dphi_q`` is
    ``(d, n)`` at the quadrature point."""
    d = upq.shape[0]
    n = dphi_q.shape[1]
    umag2 = 0.0
    for c in range(d):
        umag2 += upq[c] * upq[c]
    umag = np.sqrt(umag2)
    moving = umag > 0.0
    s = np.empty(d)
    if moving:
        for c in range(d):
            s[c] = upq[c] / umag
        inv_u = 1.0 / umag
    else:
        for c in range(d):
            s[c] = 0.0
        s[0] = 1.0
        inv_u = 0.0
    ssum = 0.0
    sg_sum = np.zeros(d)                       # sum_i sign(s . grad phi_i) grad phi_i
    for i in range(n):
        sg = 0.0
        for c in range(d):
            sg += s[c] * dphi_q[c, i]
        ssum += abs(sg)
        if sg > 0.0:
            for c in range(d):
                sg_sum[c] += dphi_q[c, i]
        elif sg < 0.0:
            for c in range(d):
                sg_sum[c] -= dphi_q[c, i]
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
    for k in range(d):
        dh_du = 0.0
        for m in range(d):
            dh_ds_m = -0.5 * h2 * sg_sum[m]
            delta = 1.0 if m == k else 0.0
            dh_du += dh_ds_m * (delta - s[m] * s[k])
        dh_du *= inv_u
        dir_k = s[k] if moving else 0.0
        dtau[k] = dtau_du * dir_k + dtau_dh * dh_du
        dgamma[k] = dgamma_du * dir_k + dgamma_dh * dh_du
    return tau, gamma


# -- momentum equations -----------------------------------------------
@njit(parallel=True, cache=True, nogil=True)
def momentum_kernel(phi, dphi, lap_phi, psi, dpsi, wh, G, ue, pe, uoe, upe, mu_q, dmu, f,
                    rho, inv_dt, inv_dt2, flags, N_e, S_e, JN_e, JS_e):
    d, ne, nq, n = dphi.shape
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
        uq = np.empty(d)
        upq = np.empty(d)
        grad = np.empty((d, d))                  # grad[c, k] = d u_c / d x_k
        lap = np.empty(d)
        conv = np.empty(d)
        R = np.empty(d)
        dtau = np.empty(d)
        dgamma = np.empty(d)
        ugrad = np.empty(n)
        dvisc = np.empty(n)
        dc = np.empty((d, d, n))
        dr = np.empty((d, d, n))
        gradR = np.empty(nc)
        for q in range(nq):
            w = wh[e, q]
            for c in range(d):
                uq[c] = 0.0
                lap[c] = 0.0
                for k in range(d):
                    grad[c, k] = 0.0
                for i in range(n):
                    uq[c] += phi[q, i] * ue[c, e, i]
                    lap[c] += lap_phi[e, q, i] * ue[c, e, i]
                    for k in range(d):
                        grad[c, k] += dphi[k, e, q, i] * ue[c, e, i]
            for i in range(n):
                s = 0.0
                for k in range(d):
                    s += uq[k] * dphi[k, e, q, i]
                ugrad[i] = s
            for c in range(d):
                s = 0.0
                for k in range(d):
                    s += uq[k] * grad[c, k]
                conv[c] = rho * s
                for i in range(n):
                    N_e[e, c * n + i] += w * conv[c] * phi[q, i]
            for c in range(d):
                for k in range(d):
                    for j in range(n):
                        dc[c, k, j] = rho * ((ugrad[j] if c == k else 0.0)
                                             + grad[c, k] * phi[q, j])
            if jac:
                for c in range(d):
                    for k in range(d):
                        for i in range(n):
                            wp = w * phi[q, i]
                            for j in range(n):
                                JN_e[e, c * n + i, k * n + j] += wp * dc[c, k, j]
            if not stabilised:
                continue

            # strong momentum residual and its derivatives
            mu = mu_q[e, q]
            nu = mu / rho
            for c in range(d):
                visc = -mu * lap[c]
                for k in range(d):
                    visc -= dmu[k, e, q] * grad[c, k] + dmu[k, e, q] * grad[k, c]
                gp = 0.0
                for m in range(nc):
                    gp += dpsi[c, e, q, m] * pe[e, m]
                R[c] = conv[c] + visc + gp
                if body:
                    R[c] -= rho * f[c, e, q]
                if transient:
                    uo = 0.0
                    for i in range(n):
                        uo += phi[q, i] * uoe[c, e, i]
                    R[c] += rho * (uq[c] - uo) * inv_dt
            # stabilisation parameters at the parameter state
            for c in range(d):
                s = 0.0
                for i in range(n):
                    s += phi[q, i] * upe[c, e, i]
                upq[c] = s
            if metric:
                tau, gamma = _params_metric(upq, nu, inv_dt2, grad_div, G[e, q], dtau, dgamma)
            else:
                tau, gamma = _params_streamline(upq, nu, inv_dt2, grad_div, dphi[:, e, q, :],
                                                dtau, dgamma)
            div = 0.0
            for c in range(d):
                div += grad[c, c]
            if stab:
                for i in range(n):
                    wi = tau * ugrad[i]
                    for c in range(d):
                        S_e[e, c * n + i] += w * (R[c] * wi + gamma * div * dphi[c, e, q, i])
            if pspg:
                for m in range(nc):
                    s = 0.0
                    for c in range(d):
                        s += dpsi[c, e, q, m] * R[c]
                    gradR[m] = s
                    S_e[e, d * n + m] += w * tau * s
            if not jac:
                continue
            for j in range(n):
                s = mu * lap_phi[e, q, j]
                for k in range(d):
                    s += dmu[k, e, q] * dphi[k, e, q, j]
                dvisc[j] = s
            for c in range(d):
                for k in range(d):
                    for j in range(n):
                        v = dc[c, k, j] - dmu[k, e, q] * dphi[c, e, q, j]
                        if c == k:
                            v -= dvisc[j]
                            if transient:
                                v += rho * phi[q, j] * inv_dt
                        dr[c, k, j] = v
            if stab:
                for i in range(n):
                    wi = tau * ugrad[i]
                    for c in range(d):
                        gci = dphi[c, e, q, i]
                        row = c * n + i
                        for k in range(d):
                            gki = dphi[k, e, q, i]
                            for j in range(n):
                                pj = phi[q, j]
                                val = (wi * dr[c, k, j] + tau * R[c] * gki * pj
                                       + gamma * gci * dphi[k, e, q, j])
                                if follow:
                                    val += (R[c] * dtau[k] * ugrad[i]
                                            + div * dgamma[k] * gci) * pj
                                JS_e[e, row, k * n + j] += w * val
                        for m in range(nc):
                            JS_e[e, row, d * n + m] += w * wi * dpsi[c, e, q, m]
            if pspg:
                for m in range(nc):
                    row = d * n + m
                    for k in range(d):
                        for j in range(n):
                            s = 0.0
                            for c in range(d):
                                s += dpsi[c, e, q, m] * dr[c, k, j]
                            val = tau * s
                            if follow:
                                val += dtau[k] * gradR[m] * phi[q, j]
                            JS_e[e, row, k * n + j] += w * val
                    for mm in range(nc):
                        s = 0.0
                        for c in range(d):
                            s += dpsi[c, e, q, m] * dpsi[c, e, q, mm]
                        JS_e[e, row, d * n + mm] += w * tau * s


def _stacked(b):
    """Gradient arrays stacked over the directions, cached on the block."""
    cache = b.__dict__.get("_numba_stacked")
    if cache is None:
        cache = (np.ascontiguousarray(np.stack(b.dphi)), np.ascontiguousarray(np.stack(b.dpsi)),
                 np.ascontiguousarray(b.Gmat))
        b.__dict__["_numba_stacked"] = cache
    return cache


def momentum_local(b, U, par, old, mu_q, dmu, body, options):
    """numba version of :func:`aa540fem.backends.numpy_kernels.momentum_local`."""
    d, n, L = b.dim, b.n, b.L
    ne = b.conn.shape[0]
    udof = b.ldof[:, :d * n].reshape(ne, d, n).transpose(1, 0, 2)          # (d, ne, n)
    pd = b.ldof[:, d * n:]
    stab, pspg, jac = options["stab"], options["pspg"], options["jacobian"]
    stabilised = stab or pspg
    ue = np.ascontiguousarray(U[udof])
    pe = np.ascontiguousarray(U[pd])
    upe = ue if par is U else np.ascontiguousarray(par[udof])
    transient = old is not None
    uoe = np.ascontiguousarray(old[udof]) if transient else ue
    dmu_s = np.ascontiguousarray(np.stack([np.asarray(g, dtype=float) for g in dmu]))
    if body is not None:
        f = np.ascontiguousarray(np.stack([np.asarray(c, dtype=float) for c in body]))
    else:
        f = dmu_s
    flags = np.array([stab, pspg, options["grad_div"], options["metric"], options["follow"],
                      transient, body is not None, jac], dtype=np.int64)
    N_e = np.zeros((ne, d * n))
    S_e = np.zeros((ne, L)) if stabilised else np.zeros((1, L))
    JN_e = np.zeros((ne, d * n, d * n)) if jac else np.zeros((1, d * n, d * n))
    JS_e = np.zeros((ne, L, L)) if (jac and stabilised) else np.zeros((1, L, L))
    dphi, dpsi, G = _stacked(b)
    momentum_kernel(b.phi, dphi, b.lap_phi, b.psi, dpsi, b.wh, G, ue, pe, uoe, upe,
                    np.ascontiguousarray(mu_q), dmu_s, f, float(options["rho"]),
                    float(options["inv_dt"]), float(options["inv_dt2"]), flags,
                    N_e, S_e, JN_e, JS_e)
    if not stabilised:
        return N_e, None, (JN_e if jac else None), (np.zeros((ne, L, L)) if jac else None)
    if not jac:
        return N_e, S_e, None, None
    return N_e, S_e, JN_e, JS_e


def warmup(threads: int | None = None):
    """Compile the kernels on small meshes (a few seconds once per machine)."""
    import numba

    from aa540fem.core.mesh import box, geometry
    from aa540fem.incompressible.assembler import FlowAssembler
    from aa540fem.incompressible.problem import FlowProblem

    if threads:
        numba.set_num_threads(threads)
    meshes = [geometry(1.0, 1.0, 1, e) for e in ("quad9", "triangle6")]
    meshes += [box(1.0, 1.0, 1.0, 1, e) for e in ("hexahedron27", "tetra10")]
    for mesh in meshes:
        d = mesh.dim
        zero = lambda *c: tuple(0.0 * c[0] for _ in range(d))   # noqa: E731
        for stabilisation, pspg in ((False, False), (True, False), (True, True)):
            prob = FlowProblem(mesh, mu=0.1, rho=1.0, stabilisation=stabilisation, pspg=pspg,
                               bc={"left": (1.0,) + (0.0,) * (d - 1)},
                               eddy_viscosity=np.full(mesh.n_nodes, 0.01), body_force=zero)
            asm = FlowAssembler(prob, backend="numba")
            U = np.linspace(0.0, 1.0, asm.space.ndof)
            asm.momentum_terms(U)
            asm.momentum_terms(U, dt=0.1, U_old=0.5 * U, param_state=U, jacobian=False)
        prob.element_length = "streamline"
        FlowAssembler(prob, backend="numba").momentum_terms(U)
