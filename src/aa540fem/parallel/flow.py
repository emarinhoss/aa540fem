"""Domain-decomposed Navier-Stokes: every rank assembles its own elements.

:class:`DistributedFlowSystem` partitions the mesh (:mod:`partition`),
gives every rank a sub-mesh of its elements (owned nodes first, then the
one-layer halo of ghost nodes), assembles residual and Jacobian on that
sub-mesh with the ordinary :class:`~aa540fem.incompressible.FlowAssembler`
and its element kernels, and adds them into a distributed PETSc matrix and
vector through a local-to-global dof map, so PETSc delivers the halo rows
to their owners.  Global dofs are renumbered rank by rank
(``[u_x, u_y, p]`` of the nodes each rank owns), which makes PETSc's row
ranges coincide with the ownership.  The state lives in a distributed
vector; before every assembly it is scattered to the local (owned + ghost)
numbering, which is the halo exchange.  Dirichlet conditions are imposed
with ``MatZeroRowsColumns`` on the owned fixed rows, norms are PETSc
reductions, and the converged state is gathered to every rank in the
original numbering so that the serial post-processing (forces, output)
applies unchanged.

The Newton iteration (damped, optionally with a frozen Jacobian) and the
theta time stepping mirror :mod:`aa540fem.linalg.newton` and
:mod:`aa540fem.incompressible.transient`; the linear solves are distributed
MUMPS (``method="direct"``) or FGMRES with the fieldsplit/LSC preconditioner
of :mod:`aa540fem.linalg.krylov` on PETSc sub-matrices
(``method="fieldsplit"``).  In a serial process the same code runs on one
rank.
"""

from __future__ import annotations

import dataclasses
import time

import numpy as np

from aa540fem.core.mesh import Mesh
from aa540fem.incompressible.assembler import FlowAssembler, dirichlet_dofs
from aa540fem.incompressible.solution import FlowSolution, TransientFlowSolution
from aa540fem.incompressible.space import TaylorHoodSpace
from aa540fem.parallel.comm import world
from aa540fem.parallel.partition import partition_mesh


def _petsc():
    from petsc4py import PETSc

    return PETSc


class DistributedFlowSystem:
    """The distributed Navier-Stokes system of ``problem`` on ``comm``."""

    def __init__(self, problem, comm=None, backend=None, parts=None):
        PETSc = _petsc()
        self.PETSc = PETSc
        self.problem = problem
        self.comm = comm or world()
        self.rank, self.size = self.comm.Get_rank(), self.comm.Get_size()
        mesh = problem.mesh
        self.space = TaylorHoodSpace(mesh)                       # global layout, no pattern
        self.parts = partition_mesh(mesh, self.size, parts)     # identical on every rank
        part = self.parts[self.rank]
        self.part = part

        # -- rank-by-rank global numbering ----------------------------------
        dim = self.space.dim
        self.dim = dim
        is_p = self.space.p_index >= 0
        n_owned = np.array([p.owned.size for p in self.parts])
        np_owned = np.array([int(is_p[p.owned].sum()) for p in self.parts])
        offsets = np.concatenate([[0], np.cumsum(dim * n_owned + np_owned)])
        new_u = [np.empty(mesh.n_nodes, dtype=np.int64) for _ in range(dim)]
        new_p = -np.ones(mesh.n_nodes, dtype=np.int64)
        for r, p in enumerate(self.parts):
            for c in range(dim):
                new_u[c][p.owned] = offsets[r] + c * n_owned[r] + np.arange(n_owned[r])
            corners = p.owned[is_p[p.owned]]
            new_p[corners] = offsets[r] + dim * n_owned[r] + np.arange(np_owned[r])
        self.n = int(offsets[-1])
        self.r0, self.r1 = int(offsets[self.rank]), int(offsets[self.rank + 1])
        self.old2new = np.concatenate(new_u + [new_p[self.space.pressure_nodes]])
        assert self.old2new.size == self.space.ndof and (self.old2new >= 0).all()

        # -- the sub-mesh of this rank's elements ---------------------------
        nodes = np.concatenate([part.owned, part.ghost])
        self.nodes = nodes
        sub_mesh = Mesh(mesh.points[nodes], part.local_cells(mesh), {})
        self.sub_problem = dataclasses.replace(problem, mesh=sub_mesh, bc={}, pin_pressure=False,
                                               eddy_viscosity=None)
        if problem.eddy_viscosity is not None:
            self.sub_problem = dataclasses.replace(self.sub_problem,
                                                   eddy_viscosity=problem.eddy_viscosity[nodes])
        self.asm = FlowAssembler(self.sub_problem, backend)
        ls = self.asm.space
        self.l2g = np.concatenate([u[nodes] for u in new_u]
                                  + [new_p[nodes[ls.pressure_nodes]]]).astype(PETSc.IntType)
        assert (self.l2g >= 0).all()
        self.n_local = ls.ndof
        self.owned_rows = np.nonzero((self.l2g >= self.r0) & (self.l2g < self.r1))[0]

        # -- Dirichlet dofs this rank owns ----------------------------------
        fixed_old, values = dirichlet_dofs(problem, self.space)
        fixed_new = self.old2new[fixed_old]
        mine = (fixed_new >= self.r0) & (fixed_new < self.r1)
        self.fixed = fixed_new[mine].astype(PETSc.IntType)
        self.fixed_values = values[mine]
        self.time_dependent = problem.depends_on_time()

        # -- PETSc objects --------------------------------------------------
        self.lgmap = PETSc.LGMap().create(self.l2g, comm=self.comm)
        self.A = PETSc.Mat().create(comm=self.comm)
        self.A.setSizes(((self.r1 - self.r0, self.n), (self.r1 - self.r0, self.n)))
        self.A.setType("aij")
        d_nnz, o_nnz = self._preallocation()
        self.A.setPreallocationNNZ((d_nnz, o_nnz))
        self.A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        self.A.setLGMap(self.lgmap, self.lgmap)
        self.A.setUp()
        self.X = self.A.createVecRight()
        self.R = self.A.createVecRight()
        self.R.setLGMap(self.lgmap)
        self.delta = self.A.createVecRight()
        self.Ul = PETSc.Vec().createSeq(self.n_local, comm=PETSc.COMM_SELF)
        self.scatter = PETSc.Scatter().create(self.X, PETSc.IS().createGeneral(self.l2g,
                                                                              comm=self.comm),
                                              self.Ul, None)
        self.to_all, self.X_all = PETSc.Scatter.toAll(self.X)
        self.F_local = self.asm.body_load(0.0)
        self.ksp = None
        self.factorisations = 0
        self.linear_iterations = []

    # -- data movement -------------------------------------------------------
    def _preallocation(self):
        pattern = self.asm.pattern
        l2g = self.l2g
        d_nnz = np.ones(self.r1 - self.r0, dtype=np.int32)
        o_nnz = np.zeros(self.r1 - self.r0, dtype=np.int32)
        for row in self.owned_rows:
            cols = l2g[pattern.indices[pattern.indptr[row]:pattern.indptr[row + 1]]]
            inside = ((cols >= self.r0) & (cols < self.r1)).sum()
            i = l2g[row] - self.r0
            d_nnz[i] = max(1, inside)
            o_nnz[i] = cols.size - inside
        return d_nnz, o_nnz

    def local_state(self, X=None):
        """The state at this rank's owned and ghost dofs (the halo exchange)."""
        PETSc = self.PETSc
        self.scatter.scatter(self.X if X is None else X, self.Ul,
                             addv=PETSc.InsertMode.INSERT_VALUES, mode=PETSc.ScatterMode.FORWARD)
        return self.Ul.getArray().copy()

    def accumulate(self, local_values, vec):
        """Sum element-wise local contributions into the distributed vector ``vec``."""
        PETSc = self.PETSc
        vec.zeroEntries()
        vec.setValuesLocal(np.arange(self.n_local, dtype=PETSc.IntType), local_values,
                           addv=PETSc.InsertMode.ADD_VALUES)
        vec.assemblyBegin()
        vec.assemblyEnd()
        return vec

    def set_matrix(self, J_local):
        PETSc = self.PETSc
        self.A.zeroEntries()
        self.A.setValuesLocalCSR(J_local.indptr.astype(PETSc.IntType),
                                 J_local.indices.astype(PETSc.IntType), J_local.data,
                                 addv=PETSc.InsertMode.ADD_VALUES)
        self.A.assemble()
        self.A.zeroRowsColumns(self.fixed, diag=1.0)           # collective; delta is 0 there
        return self.A

    def gather(self, X=None):
        """The full state on every rank in the original dof numbering."""
        PETSc = self.PETSc
        self.to_all.scatter(self.X if X is None else X, self.X_all,
                            addv=PETSc.InsertMode.INSERT_VALUES, mode=PETSc.ScatterMode.FORWARD)
        return self.X_all.getArray()[self.old2new].copy()

    def set_state(self, U_old_numbering=None, t=0.0):
        """Initialise the distributed state from a full vector (original numbering)
        or zeros, and impose the Dirichlet values."""
        PETSc = self.PETSc
        x = self.X.getArray()
        if U_old_numbering is None:
            x[:] = 0.0
        else:
            rows = np.nonzero((self.old2new >= self.r0) & (self.old2new < self.r1))[0]
            x[self.old2new[rows] - self.r0] = np.asarray(U_old_numbering)[rows]
        if self.time_dependent:
            fixed_old, values = dirichlet_dofs(self.problem, self.space, t)
            fixed_new = self.old2new[fixed_old]
            mine = (fixed_new >= self.r0) & (fixed_new < self.r1)
            self.fixed, self.fixed_values = fixed_new[mine].astype(PETSc.IntType), values[mine]
        x[self.fixed - self.r0] = self.fixed_values

    # -- residuals ---------------------------------------------------------------
    def residual(self, evaluate, X=None, jacobian=False):
        """Distributed residual (and Jacobian) from the local ``evaluate``."""
        Ul = self.local_state(X)
        if jacobian:
            R_local, J_local = evaluate(Ul, True)
            self.set_matrix(J_local)
        else:
            R_local = evaluate(Ul, False)
        self.accumulate(R_local, self.R)
        self.R.setValues(self.fixed, np.zeros(self.fixed.size),
                         addv=self.PETSc.InsertMode.INSERT_VALUES)
        self.R.assemblyBegin()
        self.R.assemblyEnd()
        return self.R

    def steady_evaluate(self, t=0.0):
        asm, F = self.asm, self.F_local

        def evaluate(Ul, jacobian):
            if jacobian:
                return asm.steady_residual_jacobian(Ul, F, t=t)
            return asm.steady_residual(Ul, F, t=t)

        return evaluate

    # -- linear solver -----------------------------------------------------------
    def _make_ksp(self, method, rtol):
        PETSc = self.PETSc
        from aa540fem.linalg.direct import mumps_available

        ksp = PETSc.KSP().create(comm=self.comm)
        ksp.setOptionsPrefix(f"dist{id(self) % 100000}_")
        ksp.setOperators(self.A)
        opts = PETSc.Options()
        prefix = ksp.getOptionsPrefix()
        if method == "direct":
            if self.size > 1 and not mumps_available():
                raise RuntimeError("a distributed LU needs PETSc built with MUMPS; "
                                   "use method='fieldsplit'")
            ksp.setType("preonly")
            pc = ksp.getPC()
            pc.setType("lu")
            pc.setFactorSolverType("mumps" if mumps_available() else "petsc")
            opts["mat_mumps_icntl_14"] = 40
        elif method == "fieldsplit":
            n_own = self.dim * self.part.owned.size
            is_u = PETSc.IS().createStride(n_own, first=self.r0, step=1, comm=self.comm)
            is_p = PETSc.IS().createStride(self.r1 - self.r0 - n_own,
                                           first=self.r0 + n_own, step=1, comm=self.comm)
            ksp.setType("fgmres")
            ksp.setTolerances(rtol=rtol, max_it=300)
            ksp.setGMRESRestart(200)
            pc = ksp.getPC()
            pc.setType("fieldsplit")
            pc.setFieldSplitIS(("u", is_u), ("p", is_p))
            pc.setFieldSplitType(PETSc.PC.CompositeType.SCHUR)
            pc.setFieldSplitSchurFactType(PETSc.PC.SchurFactType.LOWER)
            pc.setFieldSplitSchurPreType(PETSc.PC.SchurPreType.SELF)
            opts[f"{prefix}fieldsplit_u_ksp_type"] = "preonly"
            opts[f"{prefix}fieldsplit_u_pc_type"] = "bjacobi"       # ILU(1) per rank
            opts[f"{prefix}fieldsplit_u_sub_pc_type"] = "ilu"
            opts[f"{prefix}fieldsplit_u_sub_pc_factor_levels"] = 1
            opts[f"{prefix}fieldsplit_p_ksp_type"] = "preonly"
            opts[f"{prefix}fieldsplit_p_pc_type"] = "none"
            self._is_u, self._is_p = is_u, is_p
        else:
            raise ValueError(f"unknown method {method!r}; expected 'direct' or 'fieldsplit'")
        ksp.setFromOptions()
        return ksp

    def _lumped_velocity_mass(self):
        """Lumped velocity mass on the owned velocity dofs (a global Vec restricted to ``is_u``)."""
        Qdiag = np.asarray(abs(self.asm.M).sum(axis=1)).ravel()
        Qdiag[self.asm.space.n_vel:] = 0.0
        q = self.A.createVecRight()
        q.setLGMap(self.lgmap)
        self.accumulate(Qdiag, q)
        return q

    def _setup_solver(self, method, rtol):
        """(Re)factorise or rebuild the preconditioner for the current ``A``."""
        from aa540fem.linalg.krylov import LSCPreconditioner

        PETSc = self.PETSc
        if self.ksp is None:
            self.ksp = self._make_ksp(method, rtol)
            self.method = method
        self.ksp.setOperators(self.A)
        if method == "fieldsplit":
            pc = self.ksp.getPC()
            pc.setUp()
            F = self.A.createSubMatrix(self._is_u, self._is_u)
            B = self.A.createSubMatrix(self._is_p, self._is_u)
            BT = self.A.createSubMatrix(self._is_u, self._is_p)
            BT.scale(-1.0)
            if not hasattr(self, "_qinv"):
                q = self._lumped_velocity_mass()
                qu = q.getSubVector(self._is_u).copy()
                arr = qu.getArray()
                arr[:] = 1.0 / np.where(arr > 0, arr, 1.0)
                self._qinv = qu
            lsc = LSCPreconditioner.from_petsc(F, B, BT, self._qinv, comm=self.comm,
                                               prefix=f"{self.ksp.getOptionsPrefix()}lsc_")
            _, ksp_p = pc.getFieldSplitSubKSP()
            pc_p = ksp_p.getPC()
            pc_p.setType(PETSc.PC.Type.PYTHON)
            pc_p.setPythonContext(lsc)
            self._lsc = lsc
        self.ksp.setUp()
        self.factorisations += 1

    # -- Newton ---------------------------------------------------------------------
    def newton(self, evaluate, method="direct", rtol=1e-9, atol=1e-11, max_newton=30,
               damping=True, frozen_jacobian=False, linear_rtol=1e-6, verbose=False,
               reuse_solver=False):
        """Damped Newton on the distributed state ``self.X`` (updated in place)."""
        X = self.X
        R = self.residual(evaluate, X, jacobian=False)
        r = R.norm()
        residuals = [r]
        target = max(atol, rtol * r)
        if verbose and self.rank == 0:
            print(f"  Newton 0: |R| = {r:.3e}")
        have_solver = reuse_solver and self.ksp is not None
        Xn = X.duplicate()
        converged = False
        for it in range(1, max_newton + 1):
            if r <= target:
                converged = True
                break
            if not have_solver or not frozen_jacobian:
                self.residual(evaluate, X, jacobian=True)
                self._setup_solver(method, linear_rtol)
                have_solver = True
            rhs = self.R.copy()
            rhs.scale(-1.0)
            self.ksp.solve(rhs, self.delta)
            self.linear_iterations.append(self.ksp.getIterationNumber())
            alpha = 1.0
            for _ in range(6):
                Xn.waxpy(alpha, self.delta, X)
                Rn = self.residual(evaluate, Xn, jacobian=False)
                rn = Rn.norm()
                if not damping or rn < r:
                    break
                alpha *= 0.5
            else:
                break
            if frozen_jacobian and rn > r / 3.0:
                have_solver = False
            Xn.copy(X)                                # X <- Xn
            r = rn
            residuals.append(r)
            if verbose and self.rank == 0:
                print(f"  Newton {it}: |R| = {r:.3e}"
                      + (f" (damping {alpha})" if alpha < 1 else ""))
        else:
            converged = r <= target
        return {"converged": converged, "iterations": len(residuals) - 1,
                "residuals": residuals}

    # -- drivers ------------------------------------------------------------------
    def solve_steady(self, U0=None, method="direct", rtol=1e-9, atol=1e-11, max_newton=30,
                     damping=True, verbose=False, linear_rtol=1e-6) -> FlowSolution:
        """Steady Navier-Stokes by distributed Newton; the returned solution holds the
        full state (original numbering) on every rank."""
        t0 = time.perf_counter()
        self.set_state(U0)
        info = self.newton(self.steady_evaluate(), method=method, rtol=rtol, atol=atol,
                           max_newton=max_newton, damping=damping, verbose=verbose,
                           linear_rtol=linear_rtol)
        info.update({"method": method, "ranks": self.size, "path": "distributed-newton",
                     "factorisations": self.factorisations,
                     "linear_iterations": list(self.linear_iterations),
                     "wall_time": time.perf_counter() - t0})
        return FlowSolution(self.problem, self.space, self.gather(), info)

    def solve_transient(self, dt, t_end, theta=0.5, U0=None, method="direct", rtol=1e-8,
                        atol=1e-10, max_newton=25, damping=True, startup_steps=0,
                        frozen_jacobian=True, reuse_jacobian=True, store_every=1,
                        callback=None, verbose=False,
                        linear_rtol=1e-6) -> TransientFlowSolution:
        """Theta scheme with a distributed Newton solve per step (see
        :func:`aa540fem.incompressible.solve_flow_transient`)."""
        nsteps = int(round(t_end / dt))
        if nsteps < 1 or abs(nsteps * dt - t_end) > 1e-8 * max(1.0, abs(t_end)):
            raise ValueError(f"dt = {dt} must divide t_end = {t_end}")
        asm = self.asm
        nv = asm.space.n_vel
        self.set_state(U0)
        Ul = self.local_state()
        F_old = self.F_local
        S_old = asm.K @ Ul + asm.convection(Ul, jacobian=False) - F_old
        S_old[nv:] = 0.0
        times, snapshots, newton_iterations = [0.0], [self.gather()], []
        X_old = self.X.duplicate()
        th_prev = None
        for n in range(1, nsteps + 1):
            t = n * dt
            th = 1.0 if n <= startup_steps else theta
            F_new = asm.body_load(t) if self.time_dependent else F_old
            self.X.copy(X_old)
            U_old = self.local_state(X_old)

            def evaluate(Un, jacobian, U_old=U_old, F_new=F_new, S_old=S_old, th=th, t=t):
                mt = asm.momentum_terms(Un, t=t, dt=dt, U_old=U_old, jacobian=jacobian)
                S = asm.K @ Un + mt.N - F_new
                S[nv:] = 0.0
                R = (asm.M @ ((Un - U_old) / dt) + th * S + (1 - th) * S_old - asm.BT @ Un
                     + asm.B @ Un + mt.S)
                if not jacobian:
                    return R
                J = asm.pattern.matrix(asm.M_data / dt + th * (asm.K_data + mt.JN_data)
                                       + mt.JS_data - asm.BT_data + asm.B_data)
                return R, J

            if self.time_dependent:
                self.set_state(self.gather(), t)
            reuse = reuse_jacobian and th == th_prev
            th_prev = th
            info = self.newton(evaluate, method=method, rtol=rtol, atol=atol,
                               max_newton=max_newton, damping=damping,
                               frozen_jacobian=frozen_jacobian, linear_rtol=linear_rtol,
                               reuse_solver=reuse, verbose=verbose)
            if not info["converged"]:
                raise RuntimeError(f"distributed Newton did not converge at t = {t:.6g}")
            Ul = self.local_state()
            S_old = asm.K @ Ul + asm.convection(Ul, jacobian=False) - F_new
            S_old[nv:] = 0.0
            F_old = F_new
            newton_iterations.append(info["iterations"])
            if n % store_every == 0 or n == nsteps:
                times.append(t)
                snapshots.append(self.gather())
            if callback is not None:
                callback(n, t, FlowSolution(self.problem, self.space, self.gather(),
                                            {"step": n, "t": t}))
        info = {"scheme": "theta", "steps": nsteps, "dt": dt, "theta": theta,
                "startup_steps": startup_steps, "newton_iterations": newton_iterations,
                "factorisations": self.factorisations, "ranks": self.size,
                "path": "distributed-newton", "method": method}
        return TransientFlowSolution(self.problem, self.space, np.asarray(times), snapshots, info)


def solve_flow_distributed(problem, comm=None, method="direct", backend=None, **kwargs):
    """Steady solve with :class:`DistributedFlowSystem` (see :meth:`solve_steady`)."""
    return DistributedFlowSystem(problem, comm=comm, backend=backend).solve_steady(method=method,
                                                                                    **kwargs)


def solve_flow_transient_distributed(problem, dt, t_end, comm=None, method="direct",
                                     backend=None, **kwargs):
    """Theta-scheme run with :class:`DistributedFlowSystem` (see :meth:`solve_transient`)."""
    system = DistributedFlowSystem(problem, comm=comm, backend=backend)
    return system.solve_transient(dt, t_end, method=method, **kwargs)
