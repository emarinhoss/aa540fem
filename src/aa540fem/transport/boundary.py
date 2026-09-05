"""Boundary conditions of the scalar equation.  ``dirichlet`` ports
``dirichlet.m``; ``neumann`` fills in the branch that was left
unimplemented in ``main.m``.
"""

from __future__ import annotations

import numpy as np

from aa540fem.core.quadrature import gauss_legendre_quad
from aa540fem.core.util import values_at
from aa540fem.linalg.dirichlet import DirichletEliminator


def dirichlet(K, F, bc, val, x=None, y=None, t=0.0):
    """Impose ``T = val`` at the nodes ``bc``.  Port of ``dirichlet.m``.

    The corresponding rows and columns of ``K`` are zeroed, the diagonal set
    to one and ``F`` adjusted so that the remaining equations see the known
    values, exactly as the MATLAB routine does but with sparse matrices.

    Parameters
    ----------
    K   : ``(N, N)`` sparse stiffness matrix (any scipy.sparse format).
    F   : ``(N,)`` load vector.
    bc  : integer array of constrained node indices.
    val : constant, or callable ``val(x, y[, t])`` (requires ``x`` and ``y``).

    Returns the modified ``(K, F)`` with ``K`` in CSR format.
    """
    elim = DirichletEliminator(K, bc)
    if callable(val):
        if x is None or y is None:
            raise ValueError("nodal coordinates are needed for a callable Dirichlet value")
        vals = values_at(val, np.asarray(x)[elim.nodes], np.asarray(y)[elim.nodes], t)
    else:
        vals = np.full(elim.nodes.size, float(val))
    return elim.K_bc, elim.apply_rhs(F, vals)


def neumann(F, edges, val, x, y, order=3, t=0.0):
    """Add the flux ``q_n = n . (kappa grad T)`` prescribed on boundary edges.

    The weak form contributes ``F_i += int_Gamma q_n phi_i ds``, integrated
    with 1-D Gauss-Legendre quadrature along each (straight) edge.

    Parameters
    ----------
    F     : ``(N,)`` load vector (modified copy is returned).
    edges : ``(n_edges, m)`` node indices of each edge in Gmsh order
            ``(start, end)`` for ``m == 2`` or ``(start, end, mid)`` for
            ``m == 3``.
    val   : constant flux, or callable ``val(x, y[, t])``.
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
    q = values_at(val, X, Y, t)

    fe = np.einsum("eq,qi->ei", w[None, :] * ds * q, phi)
    np.add.at(F, edges.ravel(), fe.ravel())
    return F
