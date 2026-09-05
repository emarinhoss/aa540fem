"""Finite element solvers for scalar transport and incompressible flow in 2-D.

Grown out of the AA 540 (University of Washington) anisotropic heat
conduction project.  Subpackages:

``core``            meshes, reference elements, shape functions, quadrature
``io``              mesh files (meshio) and ParaView output
``linalg``          Dirichlet elimination, linear solvers, Newton iteration
``timestepping``    explicit Runge-Kutta 45
``transport``       rho_c dT/dt + u . grad T - div(kappa grad T) = f
``incompressible``  Navier-Stokes with Taylor-Hood elements

The names below are re-exported for convenience.
"""

from aa540fem.core.elements import ELEMENTS, ReferenceElement, get_element
from aa540fem.core.mesh import Mesh, geometry
from aa540fem.core.quadrature import gauss_legendre_quad, gauss_trgl, quadrature_rule
from aa540fem.core.shape_functions import (
    interpfunc,
    interpfunc_3,
    interpfunc_4,
    interpfunc_6,
    interpfunc_9,
)
from aa540fem.incompressible import (
    FlowProblem,
    FlowSolution,
    TaylorHoodSpace,
    TransientFlowSolution,
    solve_flow,
    solve_flow_transient,
)
from aa540fem.linalg import DirichletEliminator, LinearSolver, NewtonResult, newton_iterate
from aa540fem.timestepping.rk import RKResult, rk45
from aa540fem.transport import (
    DIRICHLET,
    NEUMANN,
    CellField,
    ElementMatrices,
    Operators,
    Problem,
    Solution,
    TransientSolution,
    assemble,
    assemble_operators,
    dirichlet,
    elem_eqn,
    elem_operators,
    element_gradient,
    error_norms,
    jacobian,
    neumann,
    solve,
    solve_nonlinear,
    solve_transient,
    supg_tau,
)


def read_mesh(path, orient=True):
    """Read a mesh file with meshio; see :func:`aa540fem.io.mesh_files.read_mesh`."""
    from aa540fem.io.mesh_files import read_mesh as _read

    return _read(path, orient=orient)


def write_vtk(path, mesh, point_data=None, cell_data=None):
    """Write a ParaView file; see :func:`aa540fem.io.mesh_files.write_vtk`."""
    from aa540fem.io.mesh_files import write_vtk as _write

    return _write(path, mesh, point_data, cell_data)


__all__ = [
    "Mesh", "geometry", "read_mesh", "write_vtk",
    "ELEMENTS", "ReferenceElement", "get_element",
    "gauss_legendre_quad", "gauss_trgl", "quadrature_rule",
    "interpfunc", "interpfunc_3", "interpfunc_4", "interpfunc_6", "interpfunc_9",
    "ElementMatrices", "elem_eqn", "elem_operators", "jacobian", "supg_tau",
    "dirichlet", "neumann", "CellField", "element_gradient", "error_norms",
    "DIRICHLET", "NEUMANN", "DirichletEliminator", "LinearSolver", "Operators", "Problem",
    "Solution", "assemble", "assemble_operators", "solve",
    "TransientSolution", "solve_transient",
    "NewtonResult", "newton_iterate", "solve_nonlinear",
    "RKResult", "rk45",
    "FlowProblem", "FlowSolution", "TaylorHoodSpace", "TransientFlowSolution",
    "solve_flow", "solve_flow_transient",
]
