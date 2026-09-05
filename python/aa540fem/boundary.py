"""Boundary conditions.  ``dirichlet`` ports ``dirichlet.m``; ``neumann`` fills
in the branch that was left unimplemented in ``main.m``.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .quadrature import gauss_legendre_quad


def _values_at(val, x, y):
    """Evaluate a boundary value that is either a constant or ``f(x, y)``."""
    if callable(val):
        return np.broadcast_to(np.asarray(val(x, y), dtype=float), np.shape(x)).copy()
    return np.full(np.shape(x), float(val))


def dirichlet(K, F, bc, val, x=None, y=None):
    """Impose ``T = val`` at the nodes ``bc``.  Port of ``dirichlet.m``.

    The corresponding rows and columns of ``K`` are zeroed, the diagonal set
    to one and ``F`` adjusted so that the remaining equations see the known
    values, exactly as the MATLAB routine does but with sparse matrices.

    Parameters
    ----------
    K   : ``(N, N)`` sparse stiffness matrix (any scipy.sparse format).
    F   : ``(N,)`` load vector.
    bc  : integer array of constrained node indices.
    val : constant, or callable ``val(x, y)`` (requires ``x`` and ``y``).

    Returns the modified ``(K, F)`` with ``K`` in CSR format.
    """
    bc = np.unique(np.asarray(bc, dtype=int))
    N = F.shape[0]
    if callable(val):
        if x is None or y is None:
            raise ValueError("nodal coordinates are needed for a callable Dirichlet value")
        vals = _values_at(val, np.asarray(x)[bc], np.asarray(y)[bc])
    else:
        vals = np.full(bc.size, float(val))

    K = sp.csr_matrix(K)
    F = np.array(F, dtype=float, copy=True)

    # Move the known values to the right-hand side ...
    F -= K[:, bc] @ vals
    # ... then zero the rows/columns and put a 1 on the diagonal.
    free = np.ones(N)
    free[bc] = 0.0
    D = sp.diags(free)
    fixed = sp.diags(1.0 - free)
    K = (D @ K @ D + fixed).tocsr()
    F[bc] = vals
    return K, F


def neumann(F, edges, val, x, y, order=3):
    """Add the flux ``q_n = n . (kappa grad T)`` prescribed on boundary edges.

    The weak form contributes ``F_i += int_Gamma q_n phi_i ds``, integrated
    with 1-D Gauss-Legendre quadrature along each (straight) edge.

    Parameters
    ----------
    F     : ``(N,)`` load vector (modified copy is returned).
    edges : ``(n_edges, m)`` node indices of each edge in Gmsh order
            ``(start, end)`` for ``m == 2`` or ``(start, end, mid)`` for
            ``m == 3``.
    val   : constant flux, or callable ``val(x, y)``.
    x, y  : nodal coordinate arrays.
    order : number of 1-D Gauss points per edge.
    """
    F = np.array(F, dtype=float, copy=True)
    edges = np.asarray(edges, dtype=int)
    if edges.size == 0:
        return F
    m = edges.shape[1]

    s, w = gauss_legendre_quad(order)
    if m == 2:
        phi = np.column_stack([0.5 * (1 - s), 0.5 * (1 + s)])
        dphi = np.column_stack([-0.5 * np.ones_like(s), 0.5 * np.ones_like(s)])
    elif m == 3:
        # nodes at s = -1 (start), +1 (end), 0 (mid)
        phi = np.column_stack([0.5 * (s * s - s), 0.5 * (s * s + s), 1 - s * s])
        dphi = np.column_stack([s - 0.5, s + 0.5, -2 * s])
    else:
        raise ValueError("edges must have 2 or 3 nodes")

    xe = np.asarray(x)[edges]                 # (n_edges, m)
    ye = np.asarray(y)[edges]
    X = xe @ phi.T                             # (n_edges, nq)
    Y = ye @ phi.T
    ds = np.hypot(xe @ dphi.T, ye @ dphi.T)    # |d(x,y)/ds|
    q = _values_at(val, X, Y)

    fe = np.einsum("eq,qi->ei", w[None, :] * ds * q, phi)
    np.add.at(F, edges.ravel(), fe.ravel())
    return F
