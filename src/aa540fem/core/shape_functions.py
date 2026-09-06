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


def hessian_3(xi, eta):
    """Second derivatives ``(d2/dxi2, d2/dxideta, d2/deta2)`` of the linear triangle: zero."""
    xi, eta = _as_1d(xi, eta)
    z = np.zeros((xi.size, 3))
    return z, z.copy(), z.copy()


def hessian_6(xi, eta):
    """Second derivatives of the 6-node triangle (constants)."""
    xi, eta = _as_1d(xi, eta)
    n = xi.size
    one = np.ones(n)
    zero = np.zeros(n)
    # L1 = 1 - xi - eta, L2 = xi, L3 = eta; phi = L(2L - 1), 4 L_a L_b
    d_xixi = np.column_stack([4 * one, 4 * one, zero, -8 * one, zero, zero])
    d_xieta = np.column_stack([4 * one, zero, zero, -4 * one, 4 * one, -4 * one])
    d_etaeta = np.column_stack([4 * one, zero, 4 * one, zero, zero, -8 * one])
    return d_xixi, d_xieta, d_etaeta


def hessian_4(xi, eta):
    """Second derivatives of the bilinear quadrilateral: only the mixed one is non-zero."""
    xi, eta = _as_1d(xi, eta)
    n = xi.size
    zero = np.zeros((n, 4))
    d_xieta = 0.25 * np.tile(np.array([1.0, -1.0, 1.0, -1.0]), (n, 1))
    return zero, d_xieta, zero.copy()


def hessian_9(xi, eta):
    """Second derivatives of the 9-node quadrilateral (Gmsh ordering)."""
    xi, eta = _as_1d(xi, eta)
    lx, dlx = _lagrange_quadratic_1d(xi)
    ly, dly = _lagrange_quadratic_1d(eta)
    ddl = np.tile(np.array([1.0, -2.0, 1.0]), (xi.size, 1))       # second derivative of each
    d_xixi = (ly[:, :, None] * ddl[:, None, :]).reshape(-1, 9)
    d_xieta = (dly[:, :, None] * dlx[:, None, :]).reshape(-1, 9)
    d_etaeta = (ddl[:, :, None] * lx[:, None, :]).reshape(-1, 9)
    pm = QUAD9_FROM_ROWMAJOR
    return d_xixi[:, pm], d_xieta[:, pm], d_etaeta[:, pm]


NODES_PER_ELEMENT = {1: 3, 2: 4, 3: 9}


def interpfunc(elem_type, xi, eta):
    """Shape functions of an element given by legacy number (1, 2, 3) or name."""
    from aa540fem.core.elements import get_element
    return get_element(elem_type).shape(xi, eta)


# -- three-dimensional elements (meshio / VTK local ordering) ------------
def _as_1d_3(xi, eta, zeta):
    return (np.atleast_1d(np.asarray(xi, dtype=float)),
            np.atleast_1d(np.asarray(eta, dtype=float)),
            np.atleast_1d(np.asarray(zeta, dtype=float)))


# barycentric gradients of the reference tetrahedron (0,0,0) (1,0,0) (0,1,0) (0,0,1)
_TET_DL = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
# mid-edge nodes 4..9 of ``tetra10`` (VTK ordering) lie on these corner pairs
TET10_EDGES = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))


def interpfunc_tet4(xi, eta, zeta):
    """4-node linear tetrahedron.  Returns ``(phi, dphi_dxi, dphi_deta, dphi_dzeta)``."""
    xi, eta, zeta = _as_1d_3(xi, eta, zeta)
    n = xi.size
    phi = np.column_stack([1.0 - xi - eta - zeta, xi, eta, zeta])
    d = [np.tile(_TET_DL[:, k], (n, 1)) for k in range(3)]
    return phi, d[0], d[1], d[2]


def interpfunc_tet10(xi, eta, zeta):
    """10-node quadratic tetrahedron (VTK / meshio ``tetra10`` ordering: corners,
    then the mid-points of the edges 01, 12, 02, 03, 13, 23)."""
    xi, eta, zeta = _as_1d_3(xi, eta, zeta)
    L = np.column_stack([1.0 - xi - eta - zeta, xi, eta, zeta])           # (nq, 4)
    phi = [L[:, i] * (2.0 * L[:, i] - 1.0) for i in range(4)]
    phi += [4.0 * L[:, a] * L[:, b] for a, b in TET10_EDGES]
    grads = []
    for k in range(3):
        g = [(4.0 * L[:, i] - 1.0) * _TET_DL[i, k] for i in range(4)]
        g += [4.0 * (L[:, a] * _TET_DL[b, k] + L[:, b] * _TET_DL[a, k]) for a, b in TET10_EDGES]
        grads.append(np.column_stack(g))
    return np.column_stack(phi), grads[0], grads[1], grads[2]


def hessian_tet4(xi, eta, zeta):
    xi, _, _ = _as_1d_3(xi, eta, zeta)
    z = np.zeros((xi.size, 4))
    return tuple(z.copy() for _ in range(6))


def hessian_tet10(xi, eta, zeta):
    """Second derivatives ``(xx, xy, xz, yy, yz, zz)`` of the 10-node tetrahedron (constants)."""
    xi, _, _ = _as_1d_3(xi, eta, zeta)
    n = xi.size
    out = []
    for (p, q) in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)):
        h = [4.0 * _TET_DL[i, p] * _TET_DL[i, q] for i in range(4)]
        h += [4.0 * (_TET_DL[a, p] * _TET_DL[b, q] + _TET_DL[b, p] * _TET_DL[a, q])
              for a, b in TET10_EDGES]
        out.append(np.tile(np.array(h), (n, 1)))
    return tuple(out)


# natural coordinates (xi, eta, zeta) of the 27 nodes of ``hexahedron27`` (VTK / meshio)
HEX27_NODES = np.array([
    (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
    (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
    (0, -1, -1), (1, 0, -1), (0, 1, -1), (-1, 0, -1),
    (0, -1, 1), (1, 0, 1), (0, 1, 1), (-1, 0, 1),
    (-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0),
    (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1),
    (0, 0, 0)], dtype=float)
HEX8_NODES = HEX27_NODES[:8]


def interpfunc_hex8(xi, eta, zeta):
    """8-node trilinear hexahedron (VTK ordering: bottom face counter-clockwise, then top)."""
    xi, eta, zeta = _as_1d_3(xi, eta, zeta)
    nodes = HEX8_NODES
    a, b, c = nodes[:, 0], nodes[:, 1], nodes[:, 2]
    fx, fy, fz = 1 + a * xi[:, None], 1 + b * eta[:, None], 1 + c * zeta[:, None]
    phi = 0.125 * fx * fy * fz
    return (phi, 0.125 * a * fy * fz, 0.125 * b * fx * fz, 0.125 * c * fx * fy)


def _hex27_factors(xi, eta, zeta):
    idx = (HEX27_NODES + 1).astype(int)                     # 0, 1, 2 <-> -1, 0, 1
    lx, dlx = _lagrange_quadratic_1d(xi)
    ly, dly = _lagrange_quadratic_1d(eta)
    lz, dlz = _lagrange_quadratic_1d(zeta)
    ddl = np.tile(np.array([1.0, -2.0, 1.0]), (xi.size, 1))
    return idx, (lx, dlx, ddl), (ly, dly, ddl), (lz, dlz, ddl)


def interpfunc_hex27(xi, eta, zeta):
    """27-node triquadratic hexahedron (VTK / meshio ``hexahedron27`` ordering:
    corners, 12 mid-edge nodes, 6 face centres, the centre)."""
    xi, eta, zeta = _as_1d_3(xi, eta, zeta)
    idx, (lx, dlx, _), (ly, dly, _), (lz, dlz, _) = _hex27_factors(xi, eta, zeta)
    i, j, k = idx[:, 0], idx[:, 1], idx[:, 2]
    phi = lx[:, i] * ly[:, j] * lz[:, k]
    return (phi, dlx[:, i] * ly[:, j] * lz[:, k], lx[:, i] * dly[:, j] * lz[:, k],
            lx[:, i] * ly[:, j] * dlz[:, k])


def hessian_hex8(xi, eta, zeta):
    """Second derivatives ``(xx, xy, xz, yy, yz, zz)`` of the trilinear hexahedron."""
    xi, eta, zeta = _as_1d_3(xi, eta, zeta)
    a, b, c = HEX8_NODES[:, 0], HEX8_NODES[:, 1], HEX8_NODES[:, 2]
    fx, fy, fz = 1 + a * xi[:, None], 1 + b * eta[:, None], 1 + c * zeta[:, None]
    zero = np.zeros((xi.size, 8))
    return (zero, 0.125 * a * b * fz, 0.125 * a * c * fy, zero.copy(), 0.125 * b * c * fx,
            zero.copy())


def hessian_hex27(xi, eta, zeta):
    """Second derivatives ``(xx, xy, xz, yy, yz, zz)`` of the triquadratic hexahedron."""
    xi, eta, zeta = _as_1d_3(xi, eta, zeta)
    idx, (lx, dlx, ddx), (ly, dly, ddy), (lz, dlz, ddz) = _hex27_factors(xi, eta, zeta)
    i, j, k = idx[:, 0], idx[:, 1], idx[:, 2]
    return (ddx[:, i] * ly[:, j] * lz[:, k], dlx[:, i] * dly[:, j] * lz[:, k],
            dlx[:, i] * ly[:, j] * dlz[:, k], lx[:, i] * ddy[:, j] * lz[:, k],
            lx[:, i] * dly[:, j] * dlz[:, k], lx[:, i] * ly[:, j] * ddz[:, k])
