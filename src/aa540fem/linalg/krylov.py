"""Iterative saddle-point solver: PETSc FGMRES with a block (fieldsplit) preconditioner.

The Navier-Stokes Jacobian has the block form ``[[F, -B^T], [B, C]]``
(velocity block ``F``, divergence ``B``, ``C`` zero without PSPG).  A
direct LU handles it up to some 1e5 unknowns; beyond that, and on GPUs,
the iterative route is flexible GMRES preconditioned by a block lower
triangular (Schur) factorisation: the velocity block is preconditioned by
algebraic multigrid (GAMG) or ILU, the pressure Schur complement
``B F^-1 B^T`` by its classical approximation, the pressure mass matrix
scaled by the effective viscosity (exact for Stokes, a robust choice at
moderate Reynolds numbers; Elman, Silvester & Wathen, *Finite Elements and
Fast Iterative Solvers*, ch. 9).  ``linear_backend="petsc-cuda"`` moves the
matrices and vectors to the GPU (``aijcusparse`` / ``cuda`` types) when the
PETSc build supports it.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from aa540fem.linalg.direct import petsc_available


def pressure_mass_matrix(asm, scaled: bool = True) -> sp.csr_matrix:
    """Pressure mass matrix on the corner nodes, optionally scaled by ``1 / nu_eff``."""
    Np = asm.space.Np
    rows, cols, vals = [], [], []
    for b in asm.blocks:
        mu_q, _, _ = asm.viscosity(b)
        weight = b.wh * (asm.problem.rho / mu_q if scaled else 1.0)
        Me = np.einsum("eq,qk,ql->ekl", weight, b.psi, b.psi)
        pd = asm.space.p_index[b.pconn]
        nc = pd.shape[1]
        rows.append(np.repeat(pd, nc, axis=1).ravel())
        cols.append(np.tile(pd, (1, nc)).ravel())
        vals.append(Me.ravel())
    return sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(Np, Np)).tocsr()


class FieldSplitSolver:
    """FGMRES + fieldsplit (Schur, lower) for the eliminated saddle-point Jacobian.

    ``A``: the CSR Jacobian after Dirichlet elimination; ``n_vel``: number of
    velocity dofs (the first block); ``Mp``: pressure mass matrix scaled by
    ``1 / nu`` (see :func:`pressure_mass_matrix`) used to precondition the
    Schur complement; ``velocity_pc``: ``"gamg"`` or ``"ilu"``; ``gpu``:
    use PETSc's CUDA types.
    """

    backend = "petsc-krylov"

    def __init__(self, A, n_vel, Mp, rtol=1e-8, maxiter=300, velocity_pc="ilu", gpu=False,
                 comm=None):
        if not petsc_available():
            raise ImportError("petsc4py is required for the fieldsplit solver")
        from petsc4py import PETSc

        A = sp.csr_matrix(A)
        self.n = A.shape[0]
        self.n_vel = int(n_vel)
        comm = comm or PETSc.COMM_SELF
        mat_type = "aijcusparse" if gpu else "aij"
        self.A = PETSc.Mat().createAIJ(size=A.shape, csr=(A.indptr.astype(PETSc.IntType),
                                                          A.indices.astype(PETSc.IntType), A.data),
                                       comm=comm)
        self.A.assemble()
        if gpu:
            self.A.convert(mat_type)
        Mp = sp.csr_matrix(Mp)
        self.Mp = PETSc.Mat().createAIJ(size=Mp.shape, csr=(Mp.indptr.astype(PETSc.IntType),
                                                            Mp.indices.astype(PETSc.IntType),
                                                            Mp.data), comm=comm)
        self.Mp.assemble()
        if gpu:
            self.Mp.convert(mat_type)

        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.A)
        self.ksp.setType("fgmres")
        self.ksp.setTolerances(rtol=rtol, max_it=maxiter)
        self.ksp.setGMRESRestart(50)
        pc = self.ksp.getPC()
        pc.setType("fieldsplit")
        is_u = PETSc.IS().createStride(self.n_vel, first=0, step=1, comm=comm)
        is_p = PETSc.IS().createStride(self.n - self.n_vel, first=self.n_vel, step=1, comm=comm)
        pc.setFieldSplitIS(("u", is_u), ("p", is_p))
        pc.setFieldSplitType(PETSc.PC.CompositeType.SCHUR)
        pc.setFieldSplitSchurFactType(PETSc.PC.SchurFactType.LOWER)
        pc.setFieldSplitSchurPreType(PETSc.PC.SchurPreType.USER, self.Mp)
        opts = PETSc.Options()
        prefix = self.ksp.getOptionsPrefix() or ""
        opts[f"{prefix}fieldsplit_u_ksp_type"] = "preonly"
        opts[f"{prefix}fieldsplit_u_pc_type"] = velocity_pc
        if velocity_pc == "gamg":
            opts[f"{prefix}fieldsplit_u_pc_gamg_type"] = "agg"
            opts[f"{prefix}fieldsplit_u_pc_gamg_agg_nsmooths"] = 0     # nonsymmetric block
            opts[f"{prefix}fieldsplit_u_mg_levels_ksp_type"] = "richardson"
            opts[f"{prefix}fieldsplit_u_mg_levels_pc_type"] = "sor"
        opts[f"{prefix}fieldsplit_p_ksp_type"] = "preonly"
        opts[f"{prefix}fieldsplit_p_pc_type"] = "jacobi"
        if gpu and velocity_pc == "ilu":                 # no ILU on the device
            opts[f"{prefix}fieldsplit_u_pc_type"] = "jacobi"
        self.ksp.setFromOptions()
        pc.setUp()
        self._x = self.A.createVecRight()
        self._b = self.A.createVecRight()
        if gpu:
            self._x.setType("cuda")
            self._b.setType("cuda")
        self.iterations = 0
        self.residual = None
        self.nnz_factors = None

    def solve(self, b):
        self._b.setArray(np.asarray(b, dtype=float))
        self.ksp.solve(self._b, self._x)
        self.iterations = self.ksp.getIterationNumber()
        self.residual = self.ksp.getResidualNorm()
        reason = self.ksp.getConvergedReason()
        if reason < 0:
            import warnings

            warnings.warn(f"fieldsplit FGMRES did not converge (reason {reason}) after "
                          f"{self.iterations} iterations", stacklevel=2)
        return self._x.getArray().copy()


def velocity_dof_count(asm, fixed) -> int:
    """Number of velocity dofs (the first block) of an assembler's system."""
    return 2 * asm.space.N
