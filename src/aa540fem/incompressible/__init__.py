"""Incompressible Navier-Stokes with Taylor-Hood elements.

    rho ( du/dt + (u . grad) u ) - mu lap(u) + grad p = rho f,   div u = 0

Velocity on the quadratic elements (``quad9`` / ``triangle6``), pressure on
their corner nodes (Q2/Q1, P2/P1), Newton for the steady problem, adaptive
RK45 (projected) or the theta-method in time, forces from the traction.
"""

from aa540fem.incompressible.assembler import FlowAssembler
from aa540fem.incompressible.problem import OPEN, SCHEMES, FlowProblem
from aa540fem.incompressible.solution import FlowSolution, TransientFlowSolution, traction_forces
from aa540fem.incompressible.space import TaylorHoodSpace
from aa540fem.incompressible.steady import solve_flow
from aa540fem.incompressible.transient import solve_flow_transient

__all__ = [
    "FlowAssembler", "OPEN", "SCHEMES", "FlowProblem", "FlowSolution", "TransientFlowSolution",
    "traction_forces", "TaylorHoodSpace", "solve_flow", "solve_flow_transient",
]
