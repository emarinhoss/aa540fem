"""Sparse direct factorisations behind one interface.

``factorise(A)`` returns an object with ``solve(b)``; the backend is SciPy's
SuperLU (always available, single-threaded) or PETSc's LU with MUMPS when
``petsc4py`` is installed (multithreaded through OpenMP/BLAS, and distributed
across MPI ranks under ``mpirun``).  The saddle-point Jacobians of the flow
solver are indefinite with a zero pressure block; both backends pivot, and a
residual check after the first solve falls back to SuperLU when a backend
returns a poor solution.
"""

from __future__ import annotations

import importlib.util
import os
import warnings

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

DIRECT_BACKENDS = ("superlu", "petsc")
_petsc_state = {"checked": False, "ok": False, "mumps": False}


def petsc_available() -> bool:
    """``petsc4py`` is installed and loads (its PETSc build may be for another Python)."""
    if not _petsc_state["checked"]:
        _petsc_state["checked"] = True
        if importlib.util.find_spec("petsc4py") is not None:
            try:
                from petsc4py import PETSc

                _petsc_state["ok"] = True
                _petsc_state["mumps"] = bool(PETSc.Sys.hasExternalPackage("mumps"))
            except Exception:                                    # pragma: no cover
                _petsc_state["ok"] = False
    return _petsc_state["ok"]


def mumps_available() -> bool:
    return petsc_available() and _petsc_state["mumps"]


def available_direct_backends() -> list[str]:
    """Installed backends, preferred first."""
    return ["petsc", "superlu"] if petsc_available() else ["superlu"]


def direct_backend(requested: str | None = None) -> str:
    """Resolve the backend: argument > ``AA540FEM_DIRECT`` > run configuration > SuperLU."""
    name = requested or os.environ.get("AA540FEM_DIRECT") or None
    if name is None:
        try:
            from aa540fem.hardware import get_config

            name = get_config().direct_backend
        except ImportError:
            name = "auto"
    if name in (None, "auto"):
        return "petsc" if mumps_available() else "superlu"
    if name not in DIRECT_BACKENDS:
        raise ValueError(f"unknown direct backend {name!r}; expected one of {DIRECT_BACKENDS}")
    if name == "petsc" and not petsc_available():
        warnings.warn("petsc4py is not available; using SuperLU", stacklevel=2)
        return "superlu"
    return name


class SuperLUFactorisation:
    """SciPy ``splu`` (sequential SuperLU)."""

    backend = "superlu"

    def __init__(self, A):
        self.n = A.shape[0]
        self.lu = spla.splu(sp.csc_matrix(A))
        self.nnz_factors = int(self.lu.L.nnz + self.lu.U.nnz)

    def solve(self, b):
        return self.lu.solve(np.asarray(b, dtype=float))


class PETScFactorisation:
    """PETSc LU (MUMPS when the PETSc build has it).

    Sequential in a normal Python process; under ``mpirun`` with a replicated
    matrix the rows are split across the ranks and MUMPS factorises in
    parallel (see :mod:`aa540fem.parallel`).
    """

    backend = "petsc"

    def __init__(self, A, solver_type: str | None = None, threads: int | None = None):
        from petsc4py import PETSc

        A = sp.csr_matrix(A)
        self.n = A.shape[0]
        solver_type = solver_type or ("mumps" if mumps_available() else "petsc")
        self.solver_type = solver_type
        opts = PETSc.Options()
        if solver_type == "mumps":
            opts["mat_mumps_icntl_14"] = 40          # 40 % extra workspace for the pivoting
            if threads and threads > 1:
                opts["mat_mumps_use_omp_threads"] = int(threads)
        comm = PETSc.COMM_SELF
        self.A = PETSc.Mat().createAIJ(size=(self.n, self.n),
                                       csr=(A.indptr.astype(PETSc.IntType),
                                            A.indices.astype(PETSc.IntType), A.data), comm=comm)
        self.A.assemble()
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.A)
        self.ksp.setType("preonly")
        pc = self.ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverType(solver_type)
        self.ksp.setFromOptions()
        pc.setUp()
        self._x = self.A.createVecRight()
        self._b = self.A.createVecRight()
        try:
            info = pc.getFactorMatrix().getInfo()
            self.nnz_factors = int(info.get("nz_used", 0)) or None
        except Exception:                                        # pragma: no cover
            self.nnz_factors = None

    def solve(self, b):
        self._b.setArray(np.asarray(b, dtype=float))
        self.ksp.solve(self._b, self._x)
        return self._x.getArray().copy()


def factorise(A, backend: str | None = None, threads: int | None = None, check: bool = True):
    """Factorise the square sparse matrix ``A``.

    ``backend``: ``"superlu"``, ``"petsc"`` or ``None`` (resolved by
    :func:`direct_backend`).  With ``check`` the first backend's solution of
    a random right-hand side is verified and SuperLU is used instead when
    the relative residual exceeds 1e-8 (static pivoting on the zero pressure
    block of a saddle-point matrix can fail silently).
    """
    name = direct_backend(backend)
    if name == "petsc":
        f = PETScFactorisation(A, threads=threads)
        if check:
            rng = np.random.default_rng(0)
            b = rng.standard_normal(A.shape[0])
            x = f.solve(b)
            res = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
            if not np.isfinite(res) or res > 1e-8:
                warnings.warn(f"PETSc/{f.solver_type} factorisation inaccurate "
                              f"(relative residual {res:.1e}); using SuperLU", stacklevel=2)
                return SuperLUFactorisation(A)
        return f
    return SuperLUFactorisation(A)


# -- memory estimate ----------------------------------------------------
def estimate_direct_memory(A_or_nnz, n: int | None = None, fill: float | None = None) -> int:
    """Rough peak memory in bytes of a sparse LU of ``A`` (2D meshes).

    ``fill`` is the ratio ``nnz(L + U) / nnz(A)``; measured with SuperLU
    (COLAMD ordering) on the quadratic-element saddle-point Jacobians of this
    package it grows slowly with the size, 7-12 at 1e4 unknowns to 14-23 at
    2e4-5e4 (boundary-layer meshes fill more); ``1.8 ln(n)`` is used.  The
    estimate is 12 bytes per stored factor entry (value and index) plus the
    matrix; an order of magnitude, good enough to warn before a
    factorisation that cannot fit.
    """
    if isinstance(A_or_nnz, (int, np.integer)):
        nnz = int(A_or_nnz)
        if n is None:
            raise ValueError("n is needed when nnz is given")
    else:
        nnz = int(A_or_nnz.nnz)
        n = A_or_nnz.shape[0]
    if fill is None:
        fill = min(40.0, max(8.0, 1.8 * np.log(max(n, 2))))
    return int(12 * fill * nnz + 12 * nnz)
