"""Finite element solver for 2-D steady anisotropic heat conduction.

Python port of the AA 540 (University of Washington) MATLAB project.

The equation solved is

    rho_c dT/dt + u . grad T - div( kappa(x, y) . grad T ) = f(x, y)

(steady and without convection this is the original
``div(kappa grad T) + f = 0``), with kappa a full 2x2 conductivity tensor
that may depend on T (Newton's method), an optional velocity field u with
SUPG stabilisation, and Dirichlet
(T = T0) or Neumann (n . kappa grad T = q_n) conditions on tagged boundaries.
Supported elements are 3- and 6-node triangles and 4- and 9-node
quadrilaterals, on the built-in structured rectangle or on meshes read from
files.
"""

from .boundary import dirichlet, neumann
from .element import ElementMatrices, elem_eqn, elem_operators, jacobian, supg_tau
from .elements import ELEMENTS, ReferenceElement, get_element
from .flow import (
    FlowProblem,
    FlowSolution,
    TaylorHoodSpace,
    TransientFlowSolution,
    solve_flow,
    solve_flow_transient,
)
from .geometry import Mesh, geometry
from .nonlinear import NewtonResult, newton_iterate, solve_nonlinear
from .postprocess import CellField, element_gradient, error_norms
from .quadrature import gauss_legendre_quad, gauss_trgl, quadrature_rule
from .rk import RKResult, rk45
from .shape_functions import interpfunc, interpfunc_3, interpfunc_4, interpfunc_6, interpfunc_9
from .solver import (
    DIRICHLET,
    NEUMANN,
    LinearSolver,
    Operators,
    Problem,
    Solution,
    assemble,
    assemble_operators,
    solve,
)
from .transient import TransientSolution, solve_transient


def read_mesh(path, orient=True):
    """Read a mesh file with meshio; see :func:`aa540fem.mesh_io.read_mesh`."""
    from .mesh_io import read_mesh as _read

    return _read(path, orient=orient)


def write_vtk(path, mesh, point_data=None, cell_data=None):
    """Write a ParaView file; see :func:`aa540fem.mesh_io.write_vtk`."""
    from .mesh_io import write_vtk as _write

    return _write(path, mesh, point_data, cell_data)


__all__ = [
    "Mesh", "geometry", "read_mesh", "write_vtk",
    "ELEMENTS", "ReferenceElement", "get_element",
    "gauss_legendre_quad", "gauss_trgl", "quadrature_rule",
    "interpfunc", "interpfunc_3", "interpfunc_4", "interpfunc_6", "interpfunc_9",
    "ElementMatrices", "elem_eqn", "elem_operators", "jacobian", "supg_tau",
    "dirichlet", "neumann",
    "CellField", "element_gradient", "error_norms",
    "DIRICHLET", "NEUMANN", "LinearSolver", "Operators", "Problem", "Solution",
    "assemble", "assemble_operators", "solve",
    "TransientSolution", "solve_transient",
    "NewtonResult", "newton_iterate", "solve_nonlinear",
    "RKResult", "rk45",
    "FlowProblem", "FlowSolution", "TaylorHoodSpace", "TransientFlowSolution",
    "solve_flow", "solve_flow_transient",
]
