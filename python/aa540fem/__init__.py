"""Finite element solver for 2-D steady anisotropic heat conduction.

Python port of the AA 540 (University of Washington) MATLAB project.

The equation solved on a rectangular domain [0, a] x [0, b] is

    div( kappa(x, y) . grad T ) + f(x, y) = 0

with kappa a full 2x2 conductivity tensor, and Dirichlet (T = T0) or
Neumann (n . kappa grad T = q_n) boundary conditions on each of the four
sides.  Supported elements are 3-node linear triangles, 4-node bilinear
quadrilaterals and 9-node biquadratic quadrilaterals.
"""

from .geometry import Mesh, geometry
from .quadrature import gauss_legendre_quad, gauss_trgl, quadrature_rule
from .shape_functions import interpfunc_3, interpfunc_4, interpfunc_9, interpfunc
from .element import elem_eqn
from .boundary import dirichlet, neumann
from .solver import Problem, Solution, solve

__all__ = [
    "Mesh",
    "geometry",
    "gauss_legendre_quad",
    "gauss_trgl",
    "quadrature_rule",
    "interpfunc_3",
    "interpfunc_4",
    "interpfunc_9",
    "interpfunc",
    "elem_eqn",
    "dirichlet",
    "neumann",
    "Problem",
    "Solution",
    "solve",
]
