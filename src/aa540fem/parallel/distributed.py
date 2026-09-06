"""Distributed linear algebra on top of PETSc (``petsc4py`` + ``mpi4py``).

Stage A, :class:`DistributedFactorisation`: the matrix is replicated (every
rank assembled the same system); each rank inserts its own row range into
a distributed PETSc matrix, MUMPS factorises it across the ranks, and the
solution is scattered back to every rank, so the calling Newton loop runs
identically everywhere.  Under ``mpirun`` :func:`aa540fem.linalg.direct.factorise`
uses this class for the ``"petsc"`` backend.

Stage B, :func:`distributed_matrix`: each rank assembles only the elements
of its :class:`~aa540fem.parallel.partition.Partition` with the serial
element kernels and adds them into a distributed matrix with global dof
numbers; PETSc delivers the off-process entries to their owners.  This is
the data path of a fully distributed assembly; the prototype assembles the
constant matrices (viscous, mass, divergence) and is checked against the
serial ones in ``tests/test_mpi.py``.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from aa540fem.parallel.comm import world


def _petsc():
    from petsc4py import PETSc

    return PETSc


class DistributedFactorisation:
    """LU (MUMPS when available) of a replicated CSR matrix across the ranks of ``comm``."""

    backend = "petsc"

    def __init__(self, A, comm=None, solver_type: str | None = None):
        PETSc = _petsc()
        from aa540fem.linalg.direct import mumps_available

        comm = comm or world()
        A = sp.csr_matrix(A)
        self.n = A.shape[0]
        self.comm = comm
        self.solver_type = solver_type or ("mumps" if mumps_available() else "petsc")
        self.A = PETSc.Mat().create(comm=comm)
        self.A.setSizes(((None, self.n), (None, self.n)))
        self.A.setType("aij")
        self.A.setUp()
        r0, r1 = self.A.getOwnershipRange()
        self.rows = (r0, r1)
        local = A[r0:r1]
        self.A.setPreallocationCSR((local.indptr.astype(PETSc.IntType),
                                    local.indices.astype(PETSc.IntType), local.data))
        self.A.setValuesCSR(local.indptr.astype(PETSc.IntType),
                            local.indices.astype(PETSc.IntType), local.data)
        self.A.assemble()
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.A)
        self.ksp.setType("preonly")
        pc = self.ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverType(self.solver_type)
        if self.solver_type == "mumps":
            PETSc.Options()["mat_mumps_icntl_14"] = 40
        self.ksp.setFromOptions()
        pc.setUp()
        self._b = self.A.createVecRight()
        self._x = self.A.createVecRight()
        self._scatter, self._x_all = PETSc.Scatter.toAll(self._x)
        try:
            self.nnz_factors = int(pc.getFactorMatrix().getInfo().get("nz_used", 0)) or None
        except Exception:                                    # pragma: no cover
            self.nnz_factors = None

    def solve(self, b):
        b = np.asarray(b, dtype=float)
        r0, r1 = self.rows
        self._b.setArray(b[r0:r1])
        self.ksp.solve(self._b, self._x)
        self._scatter.scatter(self._x, self._x_all, mode=_petsc().ScatterMode.FORWARD)
        return self._x_all.getArray().copy()


def distributed_matrix(n, global_dofs, local_values, comm=None):
    """Distributed PETSc matrix assembled from one part's element contributions.

    ``n``: global number of dofs; ``global_dofs``: ``(n_local_elems, L)``
    global dof numbers of the part's elements (its element-local layout);
    ``local_values``: ``(n_local_elems, L, L)`` element matrices of the same
    elements.  Every rank calls this with its own part; PETSc routes the
    entries whose row another rank owns.  Returns the assembled ``Mat``.
    """
    PETSc = _petsc()
    comm = comm or world()
    A = PETSc.Mat().create(comm=comm)
    A.setSizes(((None, n), (None, n)))
    A.setType("aij")
    A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
    A.setUp()
    dofs = global_dofs.astype(PETSc.IntType)
    for e in range(dofs.shape[0]):
        A.setValues(dofs[e], dofs[e], local_values[e], addv=PETSc.InsertMode.ADD_VALUES)
    A.assemble()
    return A


def gather_csr(A, comm=None):
    """The rows of a distributed matrix collected on every rank as one SciPy CSR."""
    comm = comm or world()
    r0, r1 = A.getOwnershipRange()
    indptr, indices, data = A.getValuesCSR()
    n = A.getSize()[0]
    local = sp.csr_matrix((data, indices, indptr), shape=(r1 - r0, n))
    parts = comm.allgather(local)
    return sp.vstack(parts).tocsr()
