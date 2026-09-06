"""Direct and preconditioned Krylov solvers for repeated solves with one matrix."""

from __future__ import annotations

import warnings

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from aa540fem.linalg.direct import factorise


class LinearSolver:
    """Direct or preconditioned Krylov solver for repeated solves with one matrix.

    ``method``: ``"direct"`` (sparse LU, factorised once; ``backend`` selects
    SuperLU or PETSc/MUMPS, see :mod:`aa540fem.linalg.direct`), ``"cg"``
    (symmetric systems only) or ``"gmres"``; the Krylov methods use pyamg
    smoothed aggregation as preconditioner when available, otherwise Jacobi
    (cg) or incomplete LU (gmres).
    """

    def __init__(self, A, method: str = "direct", tol: float = 1e-10, maxiter=None,
                 symmetric: bool = True, backend: str | None = None):
        self.method = method
        self.tol = tol
        self.maxiter = maxiter
        self.A = sp.csr_matrix(A)
        self.precond = None
        if method == "direct":
            self.lu = factorise(self.A, backend)
            self.backend = self.lu.backend
        elif method == "cg":
            if not symmetric:
                raise ValueError("cg needs a symmetric system; use method='gmres' "
                                 "for convection problems")
            self.M = self._amg("symmetric") or self._jacobi()
        elif method == "gmres":
            self.M = self._amg("nonsymmetric") or self._ilu()
        else:
            raise ValueError(f"Unknown method {method!r}; expected 'direct', 'cg' or 'gmres'")

    def _amg(self, symmetry):
        try:
            import pyamg
        except ImportError:
            return None
        ml = pyamg.smoothed_aggregation_solver(self.A, symmetry=symmetry)
        self.precond = "pyamg smoothed aggregation"
        return ml.aspreconditioner(cycle="V")

    def _jacobi(self):
        d = self.A.diagonal()
        self.precond = "Jacobi"
        return spla.LinearOperator(self.A.shape, matvec=lambda r: r / d)

    def _ilu(self):
        ilu = spla.spilu(self.A.tocsc(), drop_tol=1e-4, fill_factor=10)
        self.precond = "ILU"
        return spla.LinearOperator(self.A.shape, matvec=ilu.solve)

    def solve(self, F, verbose: bool = False):
        F = np.asarray(F, dtype=float)
        if self.method == "direct":
            return self.lu.solve(F), {"method": "direct", "backend": self.backend}

        count = [0]

        def callback(_):
            count[0] += 1

        if self.method == "cg":
            T, flag = spla.cg(self.A, F, M=self.M, rtol=self.tol, maxiter=self.maxiter,
                              callback=callback)
        else:
            T, flag = spla.gmres(self.A, F, M=self.M, rtol=self.tol, maxiter=self.maxiter,
                                 restart=50, callback=callback, callback_type="pr_norm")
        residual = np.linalg.norm(F - self.A @ T) / max(np.linalg.norm(F), 1e-300)
        if flag != 0:
            warnings.warn(f"{self.method} did not converge in {count[0]} iterations "
                          f"(relative residual {residual:.2e})", stacklevel=2)
        if verbose:
            print(f"{self.method} ({self.precond}): {count[0]} iterations, "
                  f"relative residual {residual:.2e}")
        return T, {"method": self.method, "preconditioner": self.precond,
                   "iterations": count[0], "residual": residual, "converged": flag == 0}
