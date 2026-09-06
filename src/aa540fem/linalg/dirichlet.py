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


class PatternDirichlet:
    """The same elimination for a matrix on a fixed :class:`SparsityPattern`.

    The masks of the entries to zero and the positions of the diagonal are
    computed once per (pattern, nodes) and cached on the pattern, so a new
    matrix on the pattern is eliminated in O(nnz) vector operations instead
    of sparse triple products.  Same interface as :class:`DirichletEliminator`.
    """

    def __init__(self, K, nodes):
        pattern = K.pattern
        self.nodes = np.unique(np.asarray(nodes, dtype=int))
        self.N = pattern.n
        key = self.nodes.tobytes()
        cache = pattern.__dict__.setdefault("_dirichlet_cache", {})
        if key not in cache:
            fixed = np.zeros(pattern.n, dtype=bool)
            fixed[self.nodes] = True
            colfix = fixed[pattern.indices]
            cache[key] = (fixed[pattern.rows] | colfix, pattern.diagonal_positions()[self.nodes],
                          np.nonzero(colfix)[0])
        self.kill, self.diag, self.colfix = cache[key]
        data = np.array(K.data, dtype=float, copy=True)
        self._fixed_data = data[self.colfix]                      # original K[:, nodes] entries
        self._fixed_rows = pattern.rows[self.colfix]
        self._fixed_cols = pattern.indices[self.colfix]
        data[self.kill] = 0.0
        data[self.diag] = 1.0
        self.K_bc = pattern.matrix(data)

    def apply_rhs(self, F, vals):
        F = np.array(F, dtype=float, copy=True)
        vals = np.broadcast_to(np.asarray(vals, dtype=float), self.nodes.shape)
        v = np.zeros(self.N)
        v[self.nodes] = vals
        F -= np.bincount(self._fixed_rows, weights=self._fixed_data * v[self._fixed_cols],
                         minlength=self.N)
        F[self.nodes] = vals
        return F


def eliminate(K, nodes):
    """Eliminator for ``K``: the pattern-based one when ``K`` carries a pattern."""
    if getattr(K, "pattern", None) is not None:
        return PatternDirichlet(K, nodes)
    return DirichletEliminator(K, nodes)
