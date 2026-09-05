"""Elimination of prescribed nodal values from a linear system."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


class DirichletEliminator:
    """Symmetric elimination of prescribed nodal values, reusable across
    right-hand sides.

    Given the assembled matrix ``K`` and the constrained ``nodes``, the
    modified matrix ``K_bc`` has the corresponding rows and columns zeroed and
    a one on the diagonal (as in ``dirichlet.m``), and :meth:`apply_rhs`
    moves the known values to the right-hand side of a load vector.
    """

    def __init__(self, K, nodes):
        self.nodes = np.unique(np.asarray(nodes, dtype=int))
        K = sp.csr_matrix(K)
        self.N = K.shape[0]
        self.K_fixed = K[:, self.nodes].tocsc()
        free = np.ones(self.N)
        free[self.nodes] = 0.0
        D = sp.diags(free)
        self.K_bc = (D @ K @ D + sp.diags(1.0 - free)).tocsr()

    def apply_rhs(self, F, vals):
        """Return ``F`` adjusted for ``T[nodes] = vals`` (same order as ``nodes``)."""
        F = np.array(F, dtype=float, copy=True)
        vals = np.broadcast_to(np.asarray(vals, dtype=float), self.nodes.shape)
        F -= self.K_fixed @ vals
        F[self.nodes] = vals
        return F
