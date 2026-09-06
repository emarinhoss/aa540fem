"""Residual and Jacobian assembly of the Navier-Stokes system.

Galerkin terms: viscous ``K``, mass ``M`` and divergence ``B`` matrices
(assembled once) and the convective residual with its Jacobian.  Optional
residual-based stabilisation (``FlowProblem.stabilisation``): SUPG on the
momentum equations with the streamline weight ``tau (u . grad phi_i)``
applied to the full momentum residual

    R_m = rho (u . grad) u - mu lap u + grad p - rho f  [+ rho (u - u_n)/dt]

plus grad-div (LSIC) ``int gamma (div v)(div u)``, and optionally PSPG
``int tau grad q . R_m`` on the continuity equation.  The Laplacian of the
discrete velocity comes from the shape-function Hessians (exact on affine
elements), so the stabilisation is consistent: an exact polynomial solution
such as Poiseuille flow is reproduced exactly.

All matrices live on one fixed sparsity pattern (:class:`SparsityPattern`,
the union of the element dof blocks); the element kernels
(:mod:`aa540fem.backends.numpy_kernels`, or the threaded numba versions)
return element-local arrays in the layout ``[u_x nodes, u_y nodes, p corner
nodes]`` that a :class:`ScatterPlan` sums into data vectors.  Jacobian
combinations are then vector operations on those data vectors.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import scipy.sparse as sp

from aa540fem.backends import ScatterPlan, SparsityPattern, assembly_backend, scatter_vector
from aa540fem.backends.numpy_kernels import (
    body_load_local,
    linear_local,
    momentum_local,
)
from aa540fem.core.util import values_at, values_rate
from aa540fem.incompressible.problem import OPEN, FlowProblem, values_at_pair
from aa540fem.incompressible.space import TaylorHoodSpace, _Block


class MomentumTerms(NamedTuple):
    """Nonlinear parts of the residual and their Jacobians.

    ``J_N``/``J_S`` are matrices on the assembler's sparsity pattern (``None``
    when assembled with ``jacobian=False``); ``JN_data``/``JS_data`` are the
    same Jacobians as data vectors on that pattern.
    """

    N: np.ndarray                    # Galerkin convective residual
    J_N: sp.csr_matrix | None        # its Jacobian
    S: np.ndarray                    # stabilisation residual (zero when off)
    J_S: sp.csr_matrix | None        # its Jacobian (zero matrix when off)
    JN_data: np.ndarray | None = None
    JS_data: np.ndarray | None = None


class FlowAssembler:
    """Assembles residual and Jacobian of the Navier-Stokes system.

    ``backend`` selects the element kernels (``"numpy"`` or ``"numba"``; the
    default comes from the run configuration, see :mod:`aa540fem.hardware`).
    """

    def __init__(self, problem: FlowProblem, backend: str | None = None):
        problem.validate()
        self.problem = problem
        self.mesh = problem.mesh
        self.space = TaylorHoodSpace(self.mesh)
        self.backend = assembly_backend(backend)
        self.blocks = [_Block(self.mesh, name, conn, problem.order)
                       for name, conn in self.mesh.cells.items()]
        for b in self.blocks:
            b.ldof = self.space.local_dofs(b.conn, b.pconn)
        self.mu_t = None
        if problem.eddy_viscosity is not None:
            self.mu_t = np.asarray(problem.eddy_viscosity, dtype=float)
            if self.mu_t.shape != (self.space.N,):
                raise ValueError("eddy_viscosity must be a nodal array")
        # one sparsity pattern for every matrix, and the element-to-pattern maps
        self.pattern = SparsityPattern.from_element_dofs(self.space.ndof,
                                                          [b.ldof for b in self.blocks])
        maps = [self.pattern.scatter_map(b.ldof, b.ldof) for b in self.blocks]
        self.plan = ScatterPlan(self.pattern, maps)                              # (ne, L, L)
        self.plan_uu = ScatterPlan(self.pattern, [m[:, :2 * b.n, :2 * b.n]
                                                  for m, b in zip(maps, self.blocks)])
        self._linear_matrices()

    # -- coefficients -------------------------------------------------
    def viscosity(self, b):
        """``mu_eff`` and its gradient at the quadrature points of block ``b``."""
        mu = self.problem.mu
        if self.mu_t is None:
            zero = np.zeros_like(b.X)
            return np.full_like(b.X, mu), zero, zero
        mt = self.mu_t[b.conn]
        return (mu + mt @ b.phi.T, np.einsum("eqi,ei->eq", b.dphi_dx, mt),
                np.einsum("eqi,ei->eq", b.dphi_dy, mt))

    def _body_force(self, b, t):
        if self.problem.body_force is None:
            return None
        return tuple(np.broadcast_to(np.asarray(v, dtype=float), b.X.shape)
                     for v in values_at_pair(self.problem.body_force, b.X, b.Y, t))

    # -- constant matrices --------------------------------------------
    def _linear_matrices(self):
        rho = self.problem.rho
        Ke, Me, Bxe, Bye, BTe = [], [], [], [], []
        for b in self.blocks:
            mu_q, dmu_dx, dmu_dy = self.viscosity(b)
            k, m, bx, by, bt = linear_local(b, mu_q, dmu_dx, dmu_dy, rho, self.mu_t is not None)
            Ke.append(k)
            Me.append(m)
            Bxe.append(bx)
            Bye.append(by)
            BTe.append(bt)
        self.K_data = self.plan.assemble(Ke, backend=self.backend)
        self.M_data = self.plan.assemble(Me, backend=self.backend)
        self.Bx_data = self.plan.assemble(Bxe, backend=self.backend)
        self.By_data = self.plan.assemble(Bye, backend=self.backend)
        self.B_data = self.Bx_data + self.By_data
        self.BT_data = self.plan.assemble(BTe, backend=self.backend)
        self.K = self.pattern.matrix(self.K_data)
        self.M = self.pattern.matrix(self.M_data)
        self.Bx = self.pattern.matrix(self.Bx_data)
        self.By = self.pattern.matrix(self.By_data)
        self.B = self.pattern.matrix(self.B_data)
        self.BT = self.pattern.matrix(self.BT_data)

    def body_load(self, t=0.0):
        """Load vector ``rho int phi_i f`` (zero without a body force)."""
        F = np.zeros(self.space.ndof)
        if self.problem.body_force is None:
            return F
        for b in self.blocks:
            fx, fy = self._body_force(b, t)
            F += scatter_vector(self.space.ndof, b.ldof[:, :2 * b.n],
                                body_load_local(b, fx, fy, self.problem.rho))
        return F

    # -- nonlinear terms ----------------------------------------------
    def momentum_terms(self, U, t=0.0, dt=None, U_old=None, stabilise=None,
                       pspg=None, param_state=None, jacobian: bool = True) -> MomentumTerms:
        """Convective residual/Jacobian and, if enabled, the stabilisation terms.

        ``dt`` and ``U_old`` add ``rho (u - u_old)/dt`` to the momentum
        residual used by the stabilisation (theta scheme) and select the
        transient ``tau``.  ``stabilise``/``pspg`` override the problem
        settings (the RK45 path passes ``pspg=False``).  The stabilisation
        parameters ``tau`` and ``gamma`` are differentiated in the Jacobian
        (through the velocity magnitude and through the flow-direction
        dependence of the element length, which matters on stretched cells);
        with ``param_state`` given they are evaluated at that state and
        frozen instead.  ``jacobian=False`` assembles the residuals only
        (time stepping and line searches need nothing else).
        """
        prob = self.problem
        stab = prob.stabilisation if stabilise is None else stabilise
        pspg = prob.pspg if pspg is None else pspg
        sp_ = self.space
        options = {"rho": prob.rho, "stab": bool(stab), "pspg": bool(pspg),
                   "grad_div": bool(prob.grad_div), "metric": prob.element_length == "metric",
                   "follow": param_state is None, "jacobian": bool(jacobian),
                   "inv_dt": 0.0 if dt is None else 1.0 / dt,
                   "inv_dt2": 0.0 if dt is None else (2.0 / dt) ** 2}
        U = np.asarray(U, dtype=float)
        par = U if param_state is None else np.asarray(param_state, dtype=float)
        old = None if (U_old is None or dt is None) else np.asarray(U_old, dtype=float)

        N = np.zeros(sp_.ndof)
        S = np.zeros(sp_.ndof)
        JN_e, JS_e = [], []
        for b in self.blocks:
            mu_q, dmu_dx, dmu_dy = self.viscosity(b)
            body = self._body_force(b, t) if (stab or pspg) else None
            Ne, Se, JNe, JSe = momentum_local(b, U, par, old, mu_q, dmu_dx, dmu_dy, body, options,
                                              self.backend)
            N += scatter_vector(sp_.ndof, b.ldof[:, :2 * b.n], Ne)
            if Se is not None:
                S += scatter_vector(sp_.ndof, b.ldof, Se)
            if jacobian:
                JN_e.append(JNe)
                JS_e.append(JSe)
        if not jacobian:
            return MomentumTerms(N, None, S, None)
        JN_data = self.plan_uu.assemble(JN_e, backend=self.backend)
        if stab or pspg:
            JS_data = self.plan.assemble(JS_e, backend=self.backend)
        else:
            JS_data = self.pattern.zeros()
        return MomentumTerms(N, self.pattern.matrix(JN_data), S, self.pattern.matrix(JS_data),
                             JN_data, JS_data)

    def convection(self, U, jacobian: bool = True):
        """Galerkin convective residual ``N(u)`` and its Jacobian (no stabilisation)."""
        mt = self.momentum_terms(U, stabilise=False, pspg=False, jacobian=jacobian)
        if not jacobian:
            return mt.N
        return mt.N, mt.J_N

    # -- boundary conditions ------------------------------------------
    def dirichlet(self, t=0.0, rate: bool = False):
        """Fixed dofs and their values (velocity components and pressure pin).

        With ``rate=True`` the time derivatives of the velocity values are
        returned (zero for the pressure pin and for constant values).
        """
        return dirichlet_dofs(self.problem, self.space, t, rate)

    # -- steady operator ----------------------------------------------
    def steady_residual(self, U, F, with_pressure=True, t=0.0):
        """``R = K u + N(u) + S(u) - B^T p - F`` and ``B u`` (no Jacobian)."""
        mt = self.momentum_terms(U, t=t, jacobian=False)
        R = self.K @ U + mt.N + mt.S - F
        if with_pressure:
            R = R - self.BT @ U + self.B @ U
        return R

    def steady_residual_jacobian(self, U, F, with_pressure=True, t=0.0):
        """``R = K u + N(u) + S(u) - B^T p - F`` and ``B u`` with the Jacobian."""
        mt = self.momentum_terms(U, t=t)
        R = self.K @ U + mt.N + mt.S - F
        J_data = self.K_data + mt.JN_data + mt.JS_data
        if with_pressure:
            R = R - self.BT @ U + self.B @ U
            J_data = J_data - self.BT_data + self.B_data
        return R, self.pattern.matrix(J_data)


def dirichlet_dofs(problem: FlowProblem, space: TaylorHoodSpace, t=0.0, rate: bool = False):
    """Fixed dofs of ``problem`` in the layout of ``space`` and their values
    (see :meth:`FlowAssembler.dirichlet`); needs no assembler, so a
    distributed solver can evaluate it from the mesh alone."""
    mesh = problem.mesh
    evaluate = values_rate if rate else values_at
    fixed = {}
    for tag, spec in problem.bc.items():
        if spec == OPEN:
            continue
        nodes = mesh.bc_nodes[tag]
        for comp, val in enumerate(spec):
            if val is None:
                continue
            dofs = space.dof_ux(nodes) if comp == 0 else space.dof_uy(nodes)
            vals = evaluate(val, mesh.x[nodes], mesh.y[nodes], t)
            fixed.update(zip(dofs.tolist(), vals.tolist()))
    if problem.pins_pressure:
        node = space.pressure_nodes[0]
        fixed[int(space.dof_p([node])[0])] = 0.0 if rate else float(
            values_at(problem.pin_value, mesh.x[[node]], mesh.y[[node]])[0])
    dofs = np.array(sorted(fixed), dtype=int)
    return dofs, np.array([fixed[d] for d in dofs])
