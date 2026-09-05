"""Assembly, boundary conditions and linear solve.  Port of ``main.m``'s
computational part, so that the driver script only holds the user inputs.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .boundary import dirichlet, neumann
from .element import elem_eqn
from .geometry import SIDES, Mesh, geometry
from .quadrature import MIN_ORDER, default_order, quadrature_rule
from .shape_functions import interpfunc

DIRICHLET = 0
NEUMANN = 1


@dataclass
class Problem:
    """User inputs, mirroring the top of ``main.m``.

    Attributes
    ----------
    a, b      : domain size in x and y.
    elems     : number of cells along each direction (equal in x and y).
    elem_type : 1 --> 3-node linear, 2 --> 4-node linear, 3 --> 9-node quadratic.
    bc_type   : dict side -> 0 (Dirichlet) or 1 (Neumann) for
                ``"top"``, ``"right"``, ``"left"``, ``"bottom"``.
    bc_val    : dict side -> value (constant or callable ``val(x, y)``):
                the temperature for Dirichlet sides, the normal flux
                ``n . (kappa grad T)`` for Neumann sides.
    order     : quadrature order.  Triangles: 1, 3, 4, 6, 7, 9, 12, 13
                (Gauss points on the triangle); quadrilaterals: number of
                Gauss-Legendre points per direction.  ``None`` selects the
                lowest order that fully integrates the element (1, 2 and 3
                for element types 1, 2 and 3); a lower order gives a
                rank-deficient, usually singular, system.
    material  : optional callable ``(x, y) -> (kxx, kxy, kyx, kyy, f)``
                overriding ``conductivity_and_forcing``.
    """

    a: float = 4.0
    b: float = 6.0
    elems: int = 50
    elem_type: int = 3
    bc_type: dict = field(default_factory=lambda: {
        "top": DIRICHLET, "right": NEUMANN, "left": NEUMANN, "bottom": DIRICHLET})
    bc_val: dict = field(default_factory=lambda: {
        "top": 100.0, "right": 0.0, "left": 0.0, "bottom": 0.0})
    order: int | None = None
    material: object = None


@dataclass
class Solution:
    mesh: Mesh
    T: np.ndarray                # nodal temperatures, flat
    K: sp.csr_matrix             # stiffness matrix after boundary conditions
    F: np.ndarray                # load vector after boundary conditions

    @property
    def grid(self) -> np.ndarray:
        """Temperatures on the ``(ny, nx)`` grid of ``mesh.X``/``mesh.Y``."""
        return self.mesh.to_grid(self.T)


def assemble(mesh: Mesh, order: int | None = None, material=None, verbose: bool = False):
    """Assemble the global stiffness matrix ``K`` (CSR) and load vector ``F``.

    ``order=None`` uses :func:`aa540fem.quadrature.default_order`.
    """
    if order is None:
        order = default_order(mesh.elem_type)
    elif order < MIN_ORDER[mesh.elem_type]:
        warnings.warn(
            f"quadrature order {order} under-integrates element type {mesh.elem_type} "
            f"(minimum {MIN_ORDER[mesh.elem_type]}); the system may be singular",
            stacklevel=2)
    xi, eta, w = quadrature_rule(mesh.elem_type, order)
    phi, dphi_dxi, dphi_deta = interpfunc(mesh.elem_type, xi, eta)

    xe = mesh.x[mesh.conn]
    ye = mesh.y[mesh.conn]
    Ke, fe = elem_eqn(xe, ye, phi, dphi_dxi, dphi_deta, w, material=material)

    n = mesh.nodes_per_element
    rows = np.repeat(mesh.conn, n, axis=1).ravel()   # conn[e, i] for each j
    cols = np.tile(mesh.conn, (1, n)).ravel()        # conn[e, j] for each i
    N = mesh.n_nodes
    K = sp.coo_matrix((Ke.ravel(), (rows, cols)), shape=(N, N)).tocsr()
    F = np.zeros(N)
    np.add.at(F, mesh.conn.ravel(), fe.ravel())

    if verbose:
        print("Finished Assembling Global Stiffness Matrix.")
    return K, F


def apply_boundary_conditions(mesh: Mesh, K, F, bc_type: dict, bc_val: dict, verbose=False):
    """Apply Neumann fluxes first, then Dirichlet values (which win at corners)."""
    for side in SIDES:
        if bc_type[side] == NEUMANN:
            F = neumann(F, mesh.bc_edges[side], bc_val[side], mesh.x, mesh.y)
    for side in SIDES:
        if bc_type[side] == DIRICHLET:
            K, F = dirichlet(K, F, mesh.bc_nodes[side], bc_val[side], mesh.x, mesh.y)
        elif bc_type[side] != NEUMANN:
            raise ValueError(f"Unknown boundary condition type {bc_type[side]!r} on {side}")
    if verbose:
        print("Boundary Conditions Applied.")
    return K, F


def solve(problem: Problem, verbose: bool = False) -> Solution:
    """Mesh, assemble, apply boundary conditions and solve for the temperatures."""
    mesh = geometry(problem.a, problem.b, problem.elems, problem.elem_type)
    if verbose:
        print("Finished generating Grid.")

    K, F = assemble(mesh, problem.order, problem.material, verbose)
    K, F = apply_boundary_conditions(mesh, K, F, problem.bc_type, problem.bc_val, verbose)

    if all(t == NEUMANN for t in problem.bc_type.values()):
        raise ValueError("A pure Neumann problem is singular; fix T on at least one side")

    if verbose:
        print("Solving Equations ...")
    T = spla.spsolve(K.tocsc(), F)
    return Solution(mesh, T, K, F)
