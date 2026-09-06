"""Element stiffness matrices and load vectors.  Port of ``elem_eqn.m``.

Unlike the MATLAB function, which handled one element at one quadrature
point, ``elem_eqn`` is vectorised over a batch of elements and sums over all
quadrature points at once.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from aa540fem.core.util import accepts_temperature, call_coeff, call_xyt
from aa540fem.transport.material import conductivity_and_forcing as _default_material


class ElementMatrices(NamedTuple):
    """Per-element matrices ``(n_elems, n, n)`` and load vectors ``(n_elems, n)``."""

    K: np.ndarray        # diffusion
    C: np.ndarray        # convection (+ SUPG)
    M: np.ndarray        # mass (+ SUPG)
    f: np.ndarray        # source (+ SUPG)
    dA: np.ndarray | None = None   # Newton terms d(K T - f)/dT beyond K, when T is given


FD_STEP = 1e-6


def jacobian(x, y, dphi_dxi, dphi_deta):
    """Isoparametric map derivatives at the quadrature points.

    ``x, y`` are ``(n_elems, n)`` nodal coordinates and ``dphi_*`` are
    ``(nq, n)``.  Returns ``(hs, dxi_dx, dxi_dy, deta_dx, deta_dy)``, each
    ``(n_elems, nq)``: the Jacobian determinant and the inverse Jacobian.
    """
    dx_dxi = x @ dphi_dxi.T
    dx_deta = x @ dphi_deta.T
    dy_dxi = y @ dphi_dxi.T
    dy_deta = y @ dphi_deta.T

    hs = dx_dxi * dy_deta - dx_deta * dy_dxi

    dxi_dx = dy_deta / hs
    dxi_dy = -dx_deta / hs
    deta_dx = -dy_dxi / hs
    deta_dy = dx_dxi / hs
    return hs, dxi_dx, dxi_dy, deta_dx, deta_dy


def map_gradients(inverse, dphi_dxi, dphi_deta):
    """Physical gradients ``(n_elems, nq, n)`` from an inverse Jacobian tuple.

    ``inverse`` is ``(dxi_dx, dxi_dy, deta_dx, deta_dy)`` as returned by
    :func:`jacobian` (without ``hs``); the natural derivatives may belong to
    a different element of the same reference geometry (e.g. the pressure
    element of a Taylor-Hood pair).
    """
    dxi_dx, dxi_dy, deta_dx, deta_dy = inverse
    dphi_dx = dxi_dx[:, :, None] * dphi_dxi[None] + deta_dx[:, :, None] * dphi_deta[None]
    dphi_dy = dxi_dy[:, :, None] * dphi_dxi[None] + deta_dy[:, :, None] * dphi_deta[None]
    return dphi_dx, dphi_dy


def physical_gradients(x, y, dphi_dxi, dphi_deta):
    """Shape-function gradients in physical coordinates.

    Returns ``(hs, dphi_dx, dphi_dy)`` with ``hs`` of shape ``(n_elems, nq)``
    and the gradients ``(n_elems, nq, n)``.
    """
    hs, dxi_dx, dxi_dy, deta_dx, deta_dy = jacobian(x, y, dphi_dxi, dphi_deta)
    dphi_dx, dphi_dy = map_gradients((dxi_dx, dxi_dy, deta_dx, deta_dy), dphi_dxi, dphi_deta)
    return hs, dphi_dx, dphi_dy


def physical_laplacian(inverse, d_xixi, d_xieta, d_etaeta):
    """Laplacian of the shape functions ``(n_elems, nq, n)`` from the natural Hessian.

    Uses the inverse Jacobian twice and neglects the curvature of the
    isoparametric map, which is exact on affine elements (straight-sided
    triangles, parallelogram quadrilaterals) and a small approximation on
    curved or distorted ones.
    """
    dxi_dx, dxi_dy, deta_dx, deta_dy = inverse
    a = dxi_dx ** 2 + dxi_dy ** 2
    b = dxi_dx * deta_dx + dxi_dy * deta_dy
    c = deta_dx ** 2 + deta_dy ** 2
    return (a[:, :, None] * d_xixi[None] + 2 * b[:, :, None] * d_xieta[None]
            + c[:, :, None] * d_etaeta[None])


def element_metric(inverse, family):
    """Element metric tensor ``G = J^-T T J^-1`` at the quadrature points.

    ``inverse`` is ``(dxi_dx, dxi_dy, deta_dx, deta_dy)`` from :func:`jacobian`.
    ``G`` measures lengths in units of the reference element: for a
    quadrilateral (reference ``[-1, 1]^2``) ``s . G s = (2 / h_s)^2`` with
    ``h_s`` the cell size in the unit direction ``s``, so a rectangle
    ``hx x hy`` gives ``G = diag(4/hx^2, 4/hy^2)``.  For triangles the
    constant ``T = 4 [[1, 1/2], [1/2, 1]]`` refers the metric to an
    equilateral reference triangle, so that ``G`` does not depend on which
    vertex carries the origin of the natural coordinates (an equilateral
    triangle of side ``h`` gives ``G = (4/h^2) I``, a right isosceles one
    the lengths of its legs along the legs).  Returns ``(Gxx, Gxy, Gyy)``,
    each ``(n_elems, nq)``.  Used by the metric form of the stabilisation
    parameters (Shakib 1991, Bazilevs et al. 2007), which is smooth in the
    velocity, unlike the flow-direction length of :func:`supg_length`.
    """
    dxi_dx, dxi_dy, deta_dx, deta_dy = inverse
    if family == "triangle":
        t11 = t22 = 4.0
        t12 = 2.0
    else:
        t11 = t22 = 1.0
        t12 = 0.0
    gxx = t11 * dxi_dx ** 2 + 2.0 * t12 * dxi_dx * deta_dx + t22 * deta_dx ** 2
    gxy = (t11 * dxi_dx * dxi_dy + t12 * (dxi_dx * deta_dy + dxi_dy * deta_dx)
           + t22 * deta_dx * deta_dy)
    gyy = t11 * dxi_dy ** 2 + 2.0 * t12 * dxi_dy * deta_dy + t22 * deta_dy ** 2
    return gxx, gxy, gyy


# -- d-dimensional mappings (the 2-D functions above are the d = 2 special case;
#    the arithmetic is the same, so 2-D results are unchanged) -------------
def jacobian_nd(coords, dphi_nat):
    """Isoparametric map in ``d`` dimensions.

    ``coords``: ``d`` arrays ``(n_elems, n)`` of nodal coordinates;
    ``dphi_nat``: ``d`` arrays ``(nq, n)`` of natural derivatives.  Returns
    ``(hs, inverse)`` with ``hs`` the Jacobian determinant ``(n_elems, nq)``
    and ``inverse[i][j] = d xi_i / d x_j`` (``d x d`` nested tuple of
    ``(n_elems, nq)`` arrays).
    """
    d = len(coords)
    J = [[np.asarray(coords[i]) @ np.asarray(dphi_nat[j]).T for j in range(d)] for i in range(d)]
    if d == 2:
        (a, b), (c, e) = J                              # dx/dxi dx/deta ; dy/dxi dy/deta
        hs = a * e - b * c
        inv = ((e / hs, -b / hs), (-c / hs, a / hs))    # dxi/dx dxi/dy ; deta/dx deta/dy
        return hs, inv
    if d == 3:
        (a, b, c), (e, f, g), (h, i, j) = J
        A = f * j - g * i
        B = -(e * j - g * h)
        C = e * i - f * h
        hs = a * A + b * B + c * C
        inv = ((A / hs, -(b * j - c * i) / hs, (b * g - c * f) / hs),
               (B / hs, (a * j - c * h) / hs, -(a * g - c * e) / hs),
               (C / hs, -(a * i - b * h) / hs, (a * f - b * e) / hs))
        return hs, inv
    raise ValueError("dimension must be 2 or 3")


def map_gradients_nd(inverse, dphi_nat):
    """Physical gradients, a tuple of ``d`` arrays ``(n_elems, nq, n)``."""
    d = len(inverse)
    return tuple(sum(inverse[i][j][:, :, None] * np.asarray(dphi_nat[i])[None]
                     for i in range(d)) for j in range(d))


_HESSIAN_PAIRS = {2: ((0, 0), (0, 1), (1, 1)),
                  3: ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))}


def physical_laplacian_nd(inverse, hessians):
    """Laplacian of the shape functions from the natural Hessian components
    (``(xx, xy, yy)`` or ``(xx, xy, xz, yy, yz, zz)``), neglecting the curvature
    of the map as :func:`physical_laplacian` does."""
    d = len(inverse)
    out = 0.0
    for (p, q), h in zip(_HESSIAN_PAIRS[d], hessians):
        coef = sum(inverse[p][k] * inverse[q][k] for k in range(d))
        out = out + (1.0 if p == q else 2.0) * coef[:, :, None] * np.asarray(h)[None]
    return out


def element_metric_nd(inverse, family):
    """Metric tensor ``G = J^-T T J^-1`` as a ``(n_elems, nq, d, d)`` array
    (see :func:`element_metric`; simplices use ``T`` of the equilateral
    reference simplex, ``T_ii = 4``, ``T_ij = 2``)."""
    d = len(inverse)
    if family in ("triangle", "tetra"):
        T = np.full((d, d), 2.0) + 2.0 * np.eye(d)
    else:
        T = np.eye(d)
    Ji = np.stack([np.stack(row, axis=-1) for row in inverse], axis=-2)   # (ne, nq, d, d)
    return np.einsum("...ki,kl,...lj->...ij", Ji, T, Ji)


def supg_length_nd(s, dphi):
    """Element length in the unit direction ``s`` (tuple of ``d`` arrays ``(n_elems, nq)``)
    from the gradient tuple ``dphi``: ``h = 2 / sum_i |s . grad phi_i|``."""
    sgrad = sum(s[k][:, :, None] * dphi[k] for k in range(len(s)))
    return 2.0 / np.maximum(np.abs(sgrad).sum(axis=2), 1e-300)


def supg_length(sx, sy, dphi_dx, dphi_dy):
    """Element length in the direction ``s = (sx, sy)`` (unit vectors, ``(n_elems, nq)``):
    ``h = 2 / sum_i |s . grad phi_i|`` (Tezduyar)."""
    sgrad = sx[:, :, None] * dphi_dx + sy[:, :, None] * dphi_dy
    return 2.0 / np.maximum(np.abs(sgrad).sum(axis=2), 1e-300)


def elem_eqn(x, y, phi, dphi_dxi, dphi_deta, w, material=None, t=0.0):
    """Compute element stiffness matrices ``Ke`` and load vectors ``fe``.

    Parameters
    ----------
    x, y        : ``(n_elems, n)`` nodal coordinates of each element
                  (a single element may be passed as ``(n,)`` arrays).
    phi         : ``(nq, n)`` shape functions at the quadrature points.
    dphi_dxi,
    dphi_deta   : ``(nq, n)`` derivatives w.r.t. the natural coordinates.
    w           : ``(nq,)`` quadrature weights.
    material    : callable ``(X, Y[, t]) -> (kxx, kxy, kyx, kyy, f)``; defaults to
                  :func:`aa540fem.conductivity_and_forcing.conductivity_and_forcing`.
    t           : time passed to ``material`` if it accepts one.

    Returns
    -------
    Ke : ``(n_elems, n, n)`` element stiffness matrices
    fe : ``(n_elems, n)`` element load vectors
    (with the leading dimension dropped if a single element was passed).
    """
    if material is None:
        material = _default_material

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    single = x.ndim == 1
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)
    phi = np.atleast_2d(phi)
    dphi_dxi = np.atleast_2d(dphi_dxi)
    dphi_deta = np.atleast_2d(dphi_deta)
    w = np.atleast_1d(np.asarray(w, dtype=float))

    # Physical coordinates of the quadrature points: (n_elems, nq)
    X = x @ phi.T
    Y = y @ phi.T

    hs, dphi_dx, dphi_dy = physical_gradients(x, y, dphi_dxi, dphi_deta)

    # Conductivity tensor and source at the quadrature points
    k11, k12, k21, k22, f = _broadcast(call_xyt(material, X, Y, t), X.shape)

    wh = w[None, :] * hs                 # (n_elems, nq)
    Ke = _diffusion(wh, k11, k12, k21, k22, dphi_dx, dphi_dy)
    fe = np.einsum("eq,qi->ei", wh * f, phi)

    if single:
        return Ke[0], fe[0]
    return Ke, fe


def _broadcast(values, shape):
    return tuple(np.broadcast_to(np.asarray(v, dtype=float), shape) for v in values)


def _diffusion(wh, k11, k12, k21, k22, dphi_dx, dphi_dy):
    def integrate(k, gi, gj):
        # sum_q w_q hs_q k_q gi[q, i] gj[q, j]
        return np.einsum("eq,eqi,eqj->eij", wh * k, gi, gj)

    return (integrate(k11, dphi_dx, dphi_dx)
            + integrate(k12, dphi_dx, dphi_dy)
            + integrate(k21, dphi_dy, dphi_dx)
            + integrate(k22, dphi_dy, dphi_dy))


def supg_tau(ux, uy, k11, k12, k21, k22, dphi_dx, dphi_dy, dt=None):
    """SUPG stabilisation parameter ``tau`` at the quadrature points.

    Steady form (``dt is None``): ``tau = h/(2|u|) (coth Pe - 1/Pe)`` with
    ``Pe = |u| h / (2 kappa_u)``, where ``kappa_u`` is the conductivity along
    the flow and ``h = 2 / sum_i |s . grad phi_i|`` the element length in the
    flow direction ``s = u/|u|`` (Tezduyar).  Transient form:
    ``tau = [ (2/dt)^2 + (2|u|/h)^2 + (4 kappa_u/h^2)^2 ]^(-1/2)``.
    ``tau`` is zero where the velocity vanishes.  All inputs are
    ``(n_elems, nq)`` except the gradients, ``(n_elems, nq, n)``.
    """
    umag = np.hypot(ux, uy)
    moving = umag > 0
    safe = np.where(moving, umag, 1.0)
    sx = ux / safe
    sy = uy / safe
    h = supg_length(sx, sy, dphi_dx, dphi_dy)
    kappa_u = (k11 * ux * ux + (k12 + k21) * ux * uy + k22 * uy * uy) / safe ** 2

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if dt is None:
            pe = np.where(kappa_u > 0, umag * h / (2.0 * kappa_u), np.inf)
            small = pe < 1e-3
            pe_safe = np.where(small, 1.0, pe)
            xi = np.where(small, pe / 3.0, 1.0 / np.tanh(pe_safe) - 1.0 / pe_safe)
            tau = h / (2.0 * safe) * xi
        else:
            tau = 1.0 / np.sqrt((2.0 / dt) ** 2 + (2.0 * umag / h) ** 2
                                + (4.0 * kappa_u / h ** 2) ** 2)
    return np.where(moving, tau, 0.0)


def elem_operators(x, y, phi, dphi_dxi, dphi_deta, w, material=None, velocity=None,
                   rho_c=1.0, supg=True, t=0.0, dt=None, T=None,
                   newton_terms: bool = True) -> ElementMatrices:
    """Element matrices of ``rho_c dT/dt + u.grad T - div(kappa grad T) = f``.

    Parameters are those of :func:`elem_eqn` plus

    velocity : callable ``(X, Y[, t]) -> (ux, uy)`` or ``None`` (no convection).
    rho_c    : heat capacity, constant or callable ``(X, Y[, t])``.
    supg     : add streamline-upwind Petrov-Galerkin stabilisation to the
               convection, mass and load terms (only with a velocity).  The
               diffusion part of the residual is omitted, which is exact for
               linear elements and a mild inconsistency for quadratic ones.
    dt       : time step, selects the transient form of ``tau``.
    T        : ``(n_elems, n)`` nodal temperatures.  Coefficients with a
               parameter named ``T`` are evaluated at the interpolated
               temperature and, for the material, the Newton terms
               ``dA = int grad phi_i . (dkappa/dT grad T) phi_j
               - int phi_i df/dT phi_j`` (plus the SUPG source part) are
               returned; the derivatives are central finite differences in
               ``T``.  A temperature-dependent ``velocity``/``rho_c`` is
               evaluated at ``T`` but not differentiated.
    newton_terms : set False to skip the ``dA`` computation (explicit time
               stepping only needs the residual).

    Returns an :class:`ElementMatrices` tuple ``(K, C, M, f, dA)``.
    """
    if material is None:
        material = _default_material

    x = np.atleast_2d(np.asarray(x, dtype=float))
    y = np.atleast_2d(np.asarray(y, dtype=float))
    phi = np.atleast_2d(phi)
    dphi_dxi = np.atleast_2d(dphi_dxi)
    dphi_deta = np.atleast_2d(dphi_deta)
    w = np.atleast_1d(np.asarray(w, dtype=float))

    X = x @ phi.T
    Y = y @ phi.T
    hs, dphi_dx, dphi_dy = physical_gradients(x, y, dphi_dxi, dphi_deta)

    Tq = None
    if T is not None:
        T = np.atleast_2d(np.asarray(T, dtype=float))
        Tq = T @ phi.T                                           # (n_elems, nq)

    k11, k12, k21, k22, f = _broadcast(call_coeff(material, X, Y, t, Tq), X.shape)
    rc = np.broadcast_to(np.asarray(call_coeff(rho_c, X, Y, t, Tq) if callable(rho_c)
                                    else rho_c, dtype=float), X.shape)
    wh = w[None, :] * hs

    Ke = _diffusion(wh, k11, k12, k21, k22, dphi_dx, dphi_dy)
    Me = np.einsum("eq,qi,qj->eij", wh * rc, phi, phi)
    fe = np.einsum("eq,qi->ei", wh * f, phi)
    Ce = np.zeros_like(Ke)
    wgt = None

    if velocity is not None:
        ux, uy = _broadcast(call_coeff(velocity, X, Y, t, Tq), X.shape)
        ugrad = ux[:, :, None] * dphi_dx + uy[:, :, None] * dphi_dy   # u . grad phi_j
        Ce = np.einsum("eq,qi,eqj->eij", wh, phi, ugrad)
        if supg:
            tau = supg_tau(ux, uy, k11, k12, k21, k22, dphi_dx, dphi_dy, dt)
            wgt = tau[:, :, None] * ugrad                             # tau u . grad phi_i
            Ce = Ce + np.einsum("eq,eqi,eqj->eij", wh, wgt, ugrad)
            Me = Me + np.einsum("eq,eqi,qj->eij", wh * rc, wgt, phi)
            fe = fe + np.einsum("eq,eqi->ei", wh * f, wgt)

    dAe = None
    if Tq is not None and newton_terms and accepts_temperature(material):
        dAe = _newton_terms(material, X, Y, t, Tq, T, wh, phi, dphi_dx, dphi_dy, wgt)
    return ElementMatrices(Ke, Ce, Me, fe, dAe)


def _newton_terms(material, X, Y, t, Tq, T, wh, phi, dphi_dx, dphi_dy, wgt):
    """Derivative of ``K(T) T - f(T)`` with respect to the nodal values beyond ``K``."""
    delta = FD_STEP * (1.0 + np.abs(Tq))
    plus = _broadcast(call_coeff(material, X, Y, t, Tq + delta), X.shape)
    minus = _broadcast(call_coeff(material, X, Y, t, Tq - delta), X.shape)
    dk11, dk12, dk21, dk22, df = ((a - b) / (2.0 * delta) for a, b in zip(plus, minus))

    gx = np.einsum("eqi,ei->eq", dphi_dx, T)          # grad T at the quadrature points
    gy = np.einsum("eqi,ei->eq", dphi_dy, T)
    hx = dk11 * gx + dk12 * gy                        # d(kappa grad T)/dT
    hy = dk21 * gx + dk22 * gy
    dA = (np.einsum("eq,eqi,qj->eij", wh, dphi_dx * hx[:, :, None], phi)
          + np.einsum("eq,eqi,qj->eij", wh, dphi_dy * hy[:, :, None], phi)
          - np.einsum("eq,qi,qj->eij", wh * df, phi, phi))
    if wgt is not None:
        dA = dA - np.einsum("eq,eqi,qj->eij", wh * df, wgt, phi)
    return dA
