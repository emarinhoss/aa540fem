"""Lagrange interpolation (shape) functions on the reference elements.

Each function takes arrays of natural coordinates ``xi`` and ``eta`` of
length ``nq`` and returns three ``(nq, n)`` arrays: the shape functions and
their derivatives with respect to ``xi`` and ``eta``.  Local node ordering
follows the original MATLAB files.
"""

from __future__ import annotations

import numpy as np


def interpfunc_3(xi, eta):
    """3-node linear triangle.  Port of ``interpfunc_3.m``.

    Reference triangle with node 1 at (0, 0), node 2 at (0, 1) and node 3 at
    (1, 0).
    """
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    n = xi.size
    one = np.ones(n)
    zero = np.zeros(n)

    phi = np.column_stack([1.0 - xi - eta, eta, xi])
    dphi_dxi = np.column_stack([-one, zero, one])
    dphi_deta = np.column_stack([-one, one, zero])
    return phi, dphi_dxi, dphi_deta


def interpfunc_4(xi, eta):
    """4-node bilinear quadrilateral.  Port of ``interpfunc_4.m``.

    Counter-clockwise node ordering starting at (-1, -1).
    """
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    eta = np.atleast_1d(np.asarray(eta, dtype=float))

    phi = 0.25 * np.column_stack([
        (1 - xi) * (1 - eta),
        (1 + xi) * (1 - eta),
        (1 + xi) * (1 + eta),
        (1 - xi) * (1 + eta),
    ])
    dphi_dxi = 0.25 * np.column_stack([
        -(1 - eta),
        (1 - eta),
        (1 + eta),
        -(1 + eta),
    ])
    dphi_deta = 0.25 * np.column_stack([
        -(1 - xi),
        -(1 + xi),
        (1 + xi),
        (1 - xi),
    ])
    return phi, dphi_dxi, dphi_deta


def interpfunc_9(xi, eta):
    """9-node biquadratic (Lagrange) quadrilateral.  Port of ``interpfunc_9.m``.

    Row-major node ordering: nodes 1-3 along the bottom edge (eta = -1),
    4-6 along eta = 0 (node 5 is the centre), 7-9 along the top edge.
    """
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    eta = np.atleast_1d(np.asarray(eta, dtype=float))

    # 1-D quadratic Lagrange polynomials and derivatives on [-1, 1]
    def l1d(s):
        return np.column_stack([0.5 * (s * s - s), 1.0 - s * s, 0.5 * (s * s + s)])

    def dl1d(s):
        return np.column_stack([s - 0.5, -2.0 * s, s + 0.5])

    lx, ly = l1d(xi), l1d(eta)
    dlx, dly = dl1d(xi), dl1d(eta)

    # Tensor product with xi varying fastest (matches the MATLAB ordering)
    phi = (ly[:, :, None] * lx[:, None, :]).reshape(-1, 9)
    dphi_dxi = (ly[:, :, None] * dlx[:, None, :]).reshape(-1, 9)
    dphi_deta = (dly[:, :, None] * lx[:, None, :]).reshape(-1, 9)
    return phi, dphi_dxi, dphi_deta


_INTERP = {1: interpfunc_3, 2: interpfunc_4, 3: interpfunc_9}
NODES_PER_ELEMENT = {1: 3, 2: 4, 3: 9}


def interpfunc(elem_type: int, xi, eta):
    """Dispatch to the shape functions of ``elem_type`` (1, 2 or 3)."""
    try:
        return _INTERP[elem_type](xi, eta)
    except KeyError:
        raise ValueError(f"Unknown element type {elem_type!r}; expected 1, 2 or 3") from None
