"""Lagrange interpolation (shape) functions on the reference elements.

Each function takes arrays of natural coordinates ``xi`` and ``eta`` of
length ``nq`` and returns three ``(nq, n)`` arrays: the shape functions and
their derivatives with respect to ``xi`` and ``eta``.

Local node ordering follows the Gmsh / VTK convention (counter-clockwise
corners, then edge mid-nodes, then the centre node) so that meshes read from
files can be used directly.  This differs from the MATLAB code for the
3-node triangle (nodes 2 and 3 swapped) and the 9-node quadrilateral
(row-major there); see :mod:`aa540fem.elements`.
"""

from __future__ import annotations

import numpy as np


def _as_1d(xi, eta):
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    eta = np.atleast_1d(np.asarray(eta, dtype=float))
    return xi, eta


def interpfunc_3(xi, eta):
    """3-node linear triangle.

    Reference triangle with node 1 at (0, 0), node 2 at (1, 0) and node 3 at
    (0, 1) (counter-clockwise).
    """
    xi, eta = _as_1d(xi, eta)
    n = xi.size
    one = np.ones(n)
    zero = np.zeros(n)

    phi = np.column_stack([1.0 - xi - eta, xi, eta])
    dphi_dxi = np.column_stack([-one, one, zero])
    dphi_deta = np.column_stack([-one, zero, one])
    return phi, dphi_dxi, dphi_deta


def interpfunc_6(xi, eta):
    """6-node quadratic triangle (Gmsh ``triangle6`` ordering).

    Nodes 1-3 are the corners of :func:`interpfunc_3`; nodes 4, 5, 6 are the
    mid-points of edges 1-2, 2-3 and 3-1.
    """
    xi, eta = _as_1d(xi, eta)
    l1 = 1.0 - xi - eta
    l2 = xi
    l3 = eta

    phi = np.column_stack([
        l1 * (2 * l1 - 1),
        l2 * (2 * l2 - 1),
        l3 * (2 * l3 - 1),
        4 * l1 * l2,
        4 * l2 * l3,
        4 * l3 * l1,
    ])
    zero = np.zeros_like(xi)
    dphi_dxi = np.column_stack([
        -(4 * l1 - 1),
        4 * l2 - 1,
        zero,
        4 * (l1 - l2),
        4 * l3,
        -4 * l3,
    ])
    dphi_deta = np.column_stack([
        -(4 * l1 - 1),
        zero,
        4 * l3 - 1,
        -4 * l2,
        4 * l2,
        4 * (l1 - l3),
    ])
    return phi, dphi_dxi, dphi_deta


def interpfunc_4(xi, eta):
    """4-node bilinear quadrilateral.  Port of ``interpfunc_4.m``.

    Counter-clockwise node ordering starting at (-1, -1).
    """
    xi, eta = _as_1d(xi, eta)

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


def _lagrange_quadratic_1d(s):
    """Quadratic Lagrange polynomials on [-1, 1] at nodes -1, 0, 1 and derivatives."""
    lag = np.column_stack([0.5 * (s * s - s), 1.0 - s * s, 0.5 * (s * s + s)])
    dlag = np.column_stack([s - 0.5, -2.0 * s, s + 0.5])
    return lag, dlag


def interpfunc_9_rowmajor(xi, eta):
    """9-node biquadratic quadrilateral in the MATLAB (row-major) ordering.

    Nodes 1-3 along eta = -1, 4-6 along eta = 0 (node 5 is the centre), 7-9
    along eta = 1, with xi increasing within each row.
    """
    xi, eta = _as_1d(xi, eta)
    lx, dlx = _lagrange_quadratic_1d(xi)
    ly, dly = _lagrange_quadratic_1d(eta)

    phi = (ly[:, :, None] * lx[:, None, :]).reshape(-1, 9)
    dphi_dxi = (ly[:, :, None] * dlx[:, None, :]).reshape(-1, 9)
    dphi_deta = (dly[:, :, None] * lx[:, None, :]).reshape(-1, 9)
    return phi, dphi_dxi, dphi_deta


# Row-major index of each node in the Gmsh/VTK ``quad9`` ordering:
# corners (-1,-1) (1,-1) (1,1) (-1,1), mid-edges bottom right top left, centre.
QUAD9_FROM_ROWMAJOR = np.array([0, 2, 8, 6, 1, 5, 7, 3, 4])


def interpfunc_9(xi, eta):
    """9-node biquadratic quadrilateral (Gmsh ``quad9`` ordering).

    Corners counter-clockwise from (-1, -1), then the mid-edge nodes of the
    bottom, right, top and left edges, then the centre node.
    """
    phi, dphi_dxi, dphi_deta = interpfunc_9_rowmajor(xi, eta)
    p = QUAD9_FROM_ROWMAJOR
    return phi[:, p], dphi_dxi[:, p], dphi_deta[:, p]


NODES_PER_ELEMENT = {1: 3, 2: 4, 3: 9}


def interpfunc(elem_type, xi, eta):
    """Shape functions of an element given by legacy number (1, 2, 3) or name."""
    from aa540fem.core.elements import get_element
    return get_element(elem_type).shape(xi, eta)
