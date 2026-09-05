"""Steady solve of the scalar equation with temperature-dependent coefficients.

A coefficient callable with a parameter named ``T`` (see
:func:`aa540fem.core.util.call_coeff`) makes the problem nonlinear:

    u . grad T - div( kappa(x, y, T) grad T ) = f(x, y, T)

The residual ``R(T) = A(T) T - F(T)`` is driven to zero with Newton's
method (:func:`aa540fem.linalg.newton.newton_iterate`); the Jacobian
``A + dA`` carries the derivatives of ``kappa`` and ``f`` with respect to
``T`` (central finite differences at the quadrature points, see
:func:`aa540fem.transport.element.elem_operators`).  ``newton=False``
drops ``dA`` and gives Picard (fixed-point) iteration.
"""

from __future__ import annotations

import warnings

import numpy as np

from aa540fem.core.util import values_at
from aa540fem.linalg.dirichlet import DirichletEliminator
from aa540fem.linalg.newton import newton_iterate
from aa540fem.transport.problem import (
    Problem,
    Solution,
    assemble_operators,
    dirichlet_data,
    neumann_loads,
)


def solve_nonlinear(problem: Problem, verbose: bool = False, method: str = "direct",
                    tol: float = 1e-10, maxiter=None, T0=0.0, newton: bool = True,
                    rtol: float = 1e-10, atol: float = 1e-12, max_newton: int = 25,
                    damping: bool = True) -> Solution:
    """Steady solve of a problem with temperature-dependent coefficients.

    ``T0`` is the initial guess (constant or ``T0(x, y)``); the Dirichlet
    values are imposed on it.  ``newton=False`` uses Picard iteration.
    """
    mesh = problem.build_mesh()
    problem.validate(mesh)
    if not problem.has_dirichlet():
        raise ValueError("A pure Neumann problem is singular; fix T on at least one boundary")
    if verbose:
        print(f"Finished generating Grid: {mesh.n_nodes} nodes, {mesh.n_elems} elements.")

    nodes, vals = dirichlet_data(mesh, problem.bc_type, problem.bc_val)
    T = values_at(T0, mesh.x, mesh.y)
    T[nodes] = vals
    F_neu = neumann_loads(mesh, np.zeros(mesh.n_nodes), problem.bc_type, problem.bc_val)
    last = {}

    def residual_jacobian(T):
        ops = assemble_operators(mesh, problem, T=T)
        last["ops"] = ops
        R = ops.A @ T - ops.F - F_neu
        return R, (ops.J if newton else ops.A)

    if verbose:
        print("Solving Equations (Newton)..." if newton else "Solving Equations (Picard)...")
    res = newton_iterate(residual_jacobian, T, nodes, method, tol, maxiter, rtol, atol,
                         max_newton, damping, verbose)
    if not res.converged:
        warnings.warn(f"nonlinear iteration did not converge in {res.iterations} steps "
                      f"(|R| = {res.residuals[-1]:.2e})", stacklevel=2)

    ops = last["ops"]
    elim = DirichletEliminator(ops.J, nodes)
    info = {"method": method, "nonlinear": "newton" if newton else "picard",
            "iterations": res.iterations, "residuals": res.residuals,
            "converged": res.converged}
    return Solution(mesh, res.T, elim.K_bc, elim.apply_rhs(ops.F + F_neu, vals),
                    problem.material, info)
