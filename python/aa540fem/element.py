"""Element stiffness matrices and load vectors.  Port of ``elem_eqn.m``.

Unlike the MATLAB function, which handled one element at one quadrature
point, ``elem_eqn`` is vectorised over a batch of elements and sums over all
quadrature points at once.
"""

from __future__ import annotations

import numpy as np

from .conductivity_and_forcing import conductivity_and_forcing as _default_material


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


def physical_gradients(x, y, dphi_dxi, dphi_deta):
    """Shape-function gradients in physical coordinates.

    Returns ``(hs, dphi_dx, dphi_dy)`` with ``hs`` of shape ``(n_elems, nq)``
    and the gradients ``(n_elems, nq, n)``.
    """
    hs, dxi_dx, dxi_dy, deta_dx, deta_dy = jacobian(x, y, dphi_dxi, dphi_deta)
    dphi_dx = dxi_dx[:, :, None] * dphi_dxi[None] + deta_dx[:, :, None] * dphi_deta[None]
    dphi_dy = dxi_dy[:, :, None] * dphi_dxi[None] + deta_dy[:, :, None] * dphi_deta[None]
    return hs, dphi_dx, dphi_dy


def elem_eqn(x, y, phi, dphi_dxi, dphi_deta, w, material=None):
    """Compute element stiffness matrices ``Ke`` and load vectors ``fe``.

    Parameters
    ----------
    x, y        : ``(n_elems, n)`` nodal coordinates of each element
                  (a single element may be passed as ``(n,)`` arrays).
    phi         : ``(nq, n)`` shape functions at the quadrature points.
    dphi_dxi,
    dphi_deta   : ``(nq, n)`` derivatives w.r.t. the natural coordinates.
    w           : ``(nq,)`` quadrature weights.
    material    : callable ``(X, Y) -> (kxx, kxy, kyx, kyy, f)``; defaults to
                  :func:`aa540fem.conductivity_and_forcing.conductivity_and_forcing`.

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
    k11, k12, k21, k22, f = (np.broadcast_to(np.asarray(v, dtype=float), X.shape)
                             for v in material(X, Y))

    wh = w[None, :] * hs                 # (n_elems, nq)

    def integrate(k, gi, gj):
        # sum_q w_q hs_q k_q gi[q, i] gj[q, j]
        return np.einsum("eq,eqi,eqj->eij", wh * k, gi, gj)

    Ke = (integrate(k11, dphi_dx, dphi_dx)
          + integrate(k12, dphi_dx, dphi_dy)
          + integrate(k21, dphi_dy, dphi_dx)
          + integrate(k22, dphi_dy, dphi_dy))
    fe = np.einsum("eq,qi->ei", wh * f, phi)

    if single:
        return Ke[0], fe[0]
    return Ke, fe
