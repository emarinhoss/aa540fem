"""Linear and nonlinear solvers."""

from aa540fem.linalg.dirichlet import DirichletEliminator
from aa540fem.linalg.newton import NewtonResult, newton_iterate
from aa540fem.linalg.solvers import LinearSolver

__all__ = ["DirichletEliminator", "NewtonResult", "newton_iterate", "LinearSolver"]
