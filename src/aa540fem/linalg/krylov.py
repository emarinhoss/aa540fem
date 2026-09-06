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


_instances = [0]                          # unique PETSc options prefix per solver


class LSCPreconditioner:
    """Least-squares commutator approximation of the inverse Schur complement.

    For the eliminated Jacobian ``[[F, -B^T], [B, C]]`` the Schur complement
    ``S = C + B F^-1 B^T`` is approximated by
    ``S^-1 ~ L^-1 (B Q^-1 F Q^-1 B^T) L^-1`` with ``L = B Q^-1 B^T`` and ``Q``
    the lumped velocity mass matrix (Elman, Howle, Shadid, Shuttleworth and
    Tuminaro 2006).  It is purely algebraic, so convection, grad-div and
    SUPG terms in ``F`` are accounted for automatically; ``L`` is a
    pressure-Laplacian-like matrix factorised once.  Pinned pressure dofs
    (identity rows of ``C``, zero rows of ``B``) are passed through.
    """

    def __init__(self, A, n_vel, Qdiag=None, comm=None, mat_type="aij", prefix="lsc_"):
        from petsc4py import PETSc

        A = sp.csr_matrix(A)
        F = A[:n_vel, :n_vel]
        B = A[n_vel:, :n_vel].tocsr()
        BT = A[:n_vel, n_vel:].tocsr() * -1.0          # A01 = -B^T
        if Qdiag is None:
            Qdiag = np.abs(F).sum(axis=1).A1           # fallback: row sums of |F|
        qinv = 1.0 / np.where(Qdiag > 0, Qdiag, 1.0)
        comm = comm or PETSc.COMM_SELF

        def mat(M):
            M = sp.csr_matrix(M)
            P = PETSc.Mat().createAIJ(size=M.shape, csr=(M.indptr.astype(PETSc.IntType),
                                                          M.indices.astype(PETSc.IntType),
                                                          M.data), comm=comm)
            P.assemble()
            if mat_type != "aij":
                P.convert(mat_type)
            return P

        F, B, BT = mat(F), mat(B), mat(BT)
        q = F.createVecRight()
        q.setArray(qinv)
        self._setup(F, B, BT, q, comm, prefix)

    @classmethod
    def from_petsc(cls, F, B, BT, qinv, comm=None, prefix="lsc_"):
        """Build from PETSc matrices (``BT`` already carrying the sign of ``B^T``)
        and the vector ``qinv = 1 / diag(Q)``; works on distributed matrices."""
        self = cls.__new__(cls)
        self._setup(F, B, BT, qinv, comm, prefix)
        return self

    def _setup(self, F, B, BT, qinv, comm, prefix):
        from petsc4py import PETSc

        from aa540fem.linalg.direct import mumps_available

        comm = comm or PETSc.COMM_SELF
        self.F, self.B, self.BT, self.qinv = F, B, BT, qinv
        Bq = B.duplicate(copy=True)
        Bq.diagonalScale(None, qinv)                   # B Q^-1
        L = Bq.matMult(BT)                             # B Q^-1 B^T
        L.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        diag = L.getDiagonal()
        pinned = diag.getArray() == 0.0                # pinned pressure dofs: zero rows of B
        self.pin = diag.duplicate()
        self.pin.setArray(pinned.astype(float))
        if comm.Get_size() > 1 or pinned.any():
            diag.axpy(1.0, self.pin)
            L.setDiagonal(diag, addv=PETSc.InsertMode.INSERT_VALUES)
            L.assemble()
        self.L = L
        self.ksp = PETSc.KSP().create(comm=comm)
        self.ksp.setOperators(self.L)
        self.ksp.setType("preonly")
        pc = self.ksp.getPC()
        pc.setType("lu")
        if comm.Get_size() > 1:
            if not mumps_available():
                raise RuntimeError("the distributed LSC preconditioner needs MUMPS")
            pc.setFactorSolverType("mumps")
        self.ksp.setOptionsPrefix(prefix)
        self.ksp.setFromOptions()
        self._t = self.L.createVecRight()
        self._u1 = self.F.createVecRight()
        self._u2 = self.F.createVecRight()
        self._p = self.L.createVecRight()

    def apply(self, pc, x, y):
        self.ksp.solve(x, self._t)                     # L^-1 x
        self.BT.mult(self._t, self._u1)                # B^T
        self._u1.pointwiseMult(self._u1, self.qinv)    # Q^-1
        self.F.mult(self._u1, self._u2)                # F
        self._u2.pointwiseMult(self._u2, self.qinv)    # Q^-1
        self.B.mult(self._u2, self._p)                 # B
        self.ksp.solve(self._p, y)                     # L^-1
        self._p.pointwiseMult(self.pin, x)             # pinned dofs: identity
        y.axpy(1.0, self._p)


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
                 comm=None, options=None, schur="lsc", restart=200, Qdiag=None):
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
        _instances[0] += 1
        self.ksp.setOptionsPrefix(f"fs{_instances[0]}_")    # options are global in PETSc
        self.ksp.setOperators(self.A)
        self.ksp.setType("fgmres")
        self.ksp.setTolerances(rtol=rtol, max_it=maxiter)
        self.ksp.setGMRESRestart(restart)
        pc = self.ksp.getPC()
        pc.setType("fieldsplit")
        is_u = PETSc.IS().createStride(self.n_vel, first=0, step=1, comm=comm)
        is_p = PETSc.IS().createStride(self.n - self.n_vel, first=self.n_vel, step=1, comm=comm)
        pc.setFieldSplitIS(("u", is_u), ("p", is_p))
        pc.setFieldSplitType(PETSc.PC.CompositeType.SCHUR)
        pc.setFieldSplitSchurFactType(PETSc.PC.SchurFactType.LOWER)
        opts = PETSc.Options()
        prefix = self.ksp.getOptionsPrefix() or ""
        if schur not in ("mass", "lsc"):
            raise ValueError(f"unknown Schur preconditioner {schur!r}; expected 'mass' or 'lsc'")
        self.schur = schur
        pc.setFieldSplitSchurPreType(PETSc.PC.SchurPreType.USER, self.Mp)
        opts[f"{prefix}fieldsplit_u_ksp_type"] = "preonly"
        opts[f"{prefix}fieldsplit_u_pc_type"] = velocity_pc
        if velocity_pc == "ilu" and "fieldsplit_u_pc_type" not in (options or {}):
            opts[f"{prefix}fieldsplit_u_pc_factor_levels"] = 1     # ILU(0) stalls at Re 20+
        if velocity_pc == "gamg":
            opts[f"{prefix}fieldsplit_u_pc_gamg_type"] = "agg"
            opts[f"{prefix}fieldsplit_u_pc_gamg_agg_nsmooths"] = 0     # nonsymmetric block
            opts[f"{prefix}fieldsplit_u_mg_levels_ksp_type"] = "richardson"
            opts[f"{prefix}fieldsplit_u_mg_levels_pc_type"] = "sor"
            opts[f"{prefix}fieldsplit_u_ksp_type"] = "richardson"      # two V-cycles
            opts[f"{prefix}fieldsplit_u_ksp_max_it"] = 2
        opts[f"{prefix}fieldsplit_p_ksp_type"] = "preonly"
        opts[f"{prefix}fieldsplit_p_pc_type"] = "lu" if schur == "mass" else "none"
        if gpu and velocity_pc == "ilu":                 # no ILU on the device
            opts[f"{prefix}fieldsplit_u_pc_type"] = "jacobi"
        for key, value in (options or {}).items():       # caller's PETSc options win
            opts[f"{prefix}{key}"] = value
        self.ksp.setFromOptions()
        pc.setUp()
        self.lsc = None
        if schur == "lsc":
            self.lsc = LSCPreconditioner(A, self.n_vel, Qdiag=Qdiag, comm=comm, mat_type=mat_type,
                                         prefix=f"{prefix}lsc_")
            _, ksp_p = pc.getFieldSplitSubKSP()
            pc_p = ksp_p.getPC()
            pc_p.setType(PETSc.PC.Type.PYTHON)
            pc_p.setPythonContext(self.lsc)
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
