"""Scalar transport: ``rho_c dT/dt + u . grad T - div(kappa grad T) = f``.

Steady (linear or Newton), transient (RK45 or theta), SUPG stabilisation,
post-processing.  The default material is in :mod:`aa540fem.transport.material`.
"""

from aa540fem.transport.boundary import dirichlet, neumann
from aa540fem.transport.element import ElementMatrices, elem_eqn, elem_operators, jacobian, supg_tau
from aa540fem.transport.material import conductivity_and_forcing
from aa540fem.transport.nonlinear import solve_nonlinear
from aa540fem.transport.postprocess import CellField, element_gradient, error_norms
from aa540fem.transport.problem import (
    DIRICHLET,
    NEUMANN,
    Operators,
    Problem,
    Solution,
    apply_boundary_conditions,
    assemble,
    assemble_operators,
    dirichlet_data,
    neumann_loads,
    solve,
)
from aa540fem.transport.transient import TransientSolution, solve_transient

__all__ = [
    "dirichlet", "neumann", "ElementMatrices", "elem_eqn", "elem_operators", "jacobian",
    "supg_tau", "conductivity_and_forcing", "solve_nonlinear", "CellField", "element_gradient",
    "error_norms", "DIRICHLET", "NEUMANN", "Operators", "Problem", "Solution",
    "apply_boundary_conditions", "assemble", "assemble_operators", "dirichlet_data",
    "neumann_loads", "solve", "TransientSolution", "solve_transient",
]
