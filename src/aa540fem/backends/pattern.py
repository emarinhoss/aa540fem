"""Fixed sparsity patterns and element-to-global scatter plans.

The matrices of one discretisation (viscous, mass, divergence, convective
and stabilisation Jacobians, and every combination of them) share one
sparsity pattern: the union of the element dof blocks.  Building it once and
representing every matrix as a data vector on it removes the COO-to-CSR
conversion and the sparse additions from the assembly loop, and it gives the
element kernels a fixed place to write to.  The scatter of element-local
arrays into the data vector is a segmented sum over a permutation computed
once (:class:`ScatterPlan`), so the result does not depend on how many
threads perform it.  This split (local element arrays, then a scatter) is
what a threaded CPU kernel, a GPU kernel and a distributed matrix all
consume, which keeps the physics code independent of the backend.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


class PatternMatrix(sp.csr_matrix):
    """A CSR matrix that remembers the :class:`SparsityPattern` it lives on.

    Arithmetic returns plain SciPy matrices (as usual); the attribute is only
    used by the solvers to take the fast elimination path.
    """

    pattern: SparsityPattern | None = None


class SparsityPattern:
    """Canonical CSR structure (sorted indices, no duplicates) of a square matrix."""

    def __init__(self, n: int, indptr: np.ndarray, indices: np.ndarray):
        self.n = int(n)
        self.indptr = np.asarray(indptr, dtype=np.int64)
        self.indices = np.asarray(indices, dtype=np.int32)
        self.nnz = int(self.indices.size)
        self.rows = np.repeat(np.arange(self.n, dtype=np.int64), np.diff(self.indptr))
        # entry keys row * n + col are strictly increasing in CSR order
        self._keys = self.rows * self.n + self.indices
        self._diag = None

    @classmethod
    def from_element_dofs(cls, n: int, ldof_blocks) -> SparsityPattern:
        """Union of the ``ldof x ldof`` blocks of every element.

        ``ldof_blocks`` is a list of ``(n_elems, L)`` integer arrays, one per
        cell block (elements of different types have different ``L``).
        """
        rows, cols = [], []
        for ldof in ldof_blocks:
            ldof = np.asarray(ldof, dtype=np.int64)
            L = ldof.shape[1]
            rows.append(np.repeat(ldof, L, axis=1).ravel())
            cols.append(np.tile(ldof, (1, L)).ravel())
        rows = np.concatenate(rows) if rows else np.zeros(0, dtype=np.int64)
        cols = np.concatenate(cols) if cols else np.zeros(0, dtype=np.int64)
        A = sp.coo_matrix((np.ones(rows.size, dtype=np.int8), (rows, cols)), shape=(n, n)).tocsr()
        A.sum_duplicates()
        A.sort_indices()
        return cls(n, A.indptr, A.indices)

    # -- positions --------------------------------------------------------
    def scatter_map(self, rdof, cdof) -> np.ndarray:
        """Positions in the data vector of the entries ``(rdof[e, a], cdof[e, b])``.

        ``rdof`` is ``(n_elems, a)``, ``cdof`` ``(n_elems, b)``; the result is
        ``(n_elems, a, b)`` int64.  Every requested entry must be in the pattern.
        """
        rdof = np.asarray(rdof, dtype=np.int64)
        cdof = np.asarray(cdof, dtype=np.int64)
        keys = rdof[:, :, None] * self.n + cdof[:, None, :]
        pos = np.searchsorted(self._keys, keys)
        pos = np.minimum(pos, self.nnz - 1)
        if not np.array_equal(self._keys[pos], keys):
            raise ValueError("requested entries are not all in the sparsity pattern")
        return pos

    def diagonal_positions(self) -> np.ndarray:
        if self._diag is None:
            ar = np.arange(self.n, dtype=np.int64)
            self._diag = self.scatter_map(ar[:, None], ar[:, None])[:, 0, 0]
        return self._diag

    # -- matrices ---------------------------------------------------------
    def zeros(self) -> np.ndarray:
        return np.zeros(self.nnz)

    def matrix(self, data) -> PatternMatrix:
        """CSR matrix on this pattern sharing ``data`` (no copy, O(1))."""
        data = np.asarray(data, dtype=float)
        if data.shape != (self.nnz,):
            raise ValueError(f"data must have shape ({self.nnz},)")
        m = PatternMatrix((data, self.indices, self.indptr), shape=(self.n, self.n))
        m.pattern = self
        return m

    def data_of(self, A) -> np.ndarray:
        """Data vector of a SciPy matrix with (a subset of) this pattern."""
        A = sp.csr_matrix(A)
        A.sum_duplicates()
        A.sort_indices()
        rows = np.repeat(np.arange(self.n, dtype=np.int64), np.diff(A.indptr))
        pos = np.searchsorted(self._keys, rows * self.n + A.indices)
        pos = np.minimum(pos, self.nnz - 1)
        if not np.array_equal(self._keys[pos], rows * self.n + A.indices):
            raise ValueError("matrix has entries outside the pattern")
        data = self.zeros()
        np.add.at(data, pos, A.data)
        return data


class ScatterPlan:
    """Sum of element-local arrays into a data vector on a fixed pattern.

    ``maps`` is a list of position arrays ``(n_elems, a, b)`` from
    :meth:`SparsityPattern.scatter_map` (one per cell block or block pair);
    :meth:`assemble` takes the matching list of value arrays.  The reduction
    order is fixed (``perm``/``seg_start`` are computed once), so the result
    is bitwise reproducible for any number of threads.
    """

    def __init__(self, pattern: SparsityPattern, maps):
        self.pattern = pattern
        self.maps = [np.asarray(m, dtype=np.int64) for m in maps]
        self.sizes = [m.size for m in self.maps]
        flat = (np.concatenate([m.ravel() for m in self.maps]) if self.maps
                else np.zeros(0, dtype=np.int64))
        self.positions = flat
        self.perm = np.argsort(flat, kind="stable")
        self.seg_start = np.searchsorted(flat[self.perm], np.arange(pattern.nnz + 1))

    def assemble(self, values, out=None, backend: str = "numpy") -> np.ndarray:
        """Data vector of the sum of the element contributions.

        ``values`` lists arrays with the shapes of ``maps``; ``out`` (a data
        vector) is added to when given.
        """
        vals = (np.concatenate([np.ascontiguousarray(v, dtype=float).ravel() for v in values])
                if values else np.zeros(0))
        if vals.size != self.positions.size:
            raise ValueError("values do not match the scatter maps")
        if backend == "numba":
            from aa540fem.backends.numba_kernels import segmented_reduce

            data = np.zeros(self.pattern.nnz)
            segmented_reduce(vals, self.perm, self.seg_start, data)
        else:
            data = np.bincount(self.positions, weights=vals, minlength=self.pattern.nnz)
        if out is not None:
            out += data
            return out
        return data


def scatter_vector(n: int, dofs, values) -> np.ndarray:
    """Global vector from element-local values ``(n_elems, a)`` at ``dofs``."""
    return np.bincount(np.asarray(dofs).ravel(), weights=np.asarray(values, dtype=float).ravel(),
                       minlength=n)


def add_matrices(A, B, beta: float = 1.0):
    """``A + beta B`` on the shared pattern when both carry the same one,
    otherwise the ordinary sparse sum."""
    pa, pb = getattr(A, "pattern", None), getattr(B, "pattern", None)
    if pa is not None and pa is pb:
        return pa.matrix(A.data + beta * B.data)
    return (A + beta * B).tocsr()

