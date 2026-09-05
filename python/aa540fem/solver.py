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
from .elements import get_element
from .geometry import Mesh, geometry

DIRICHLET = 0
NEUMANN = 1


@dataclass
class Problem:
    """User inputs, mirroring the top of ``main.m``.

    Attributes
    ----------
    a, b      : domain size in x and y (structured rectangle only).
    elems     : number of cells along each direction (structured only).
    elem_type : 1 --> 3-node linear, 2 --> 4-node linear, 3 --> 9-node
                quadratic; element names (``"triangle6"`` ...) also work.
    mesh      : an explicit :class:`Mesh` (e.g. from
                :func:`aa540fem.mesh_io.read_mesh`).  When given, ``a``,
                ``b``, ``elems`` and ``elem_type`` are ignored.
    bc_type   : dict boundary tag -> 0 (Dirichlet) or 1 (Neumann).  Tags of
                the mesh that are not listed get the natural (zero-flux)
                condition.
    bc_val    : dict boundary tag -> value (constant or callable ``val(x, y)``):
                the temperature for Dirichlet tags, the normal flux
                ``n . (kappa grad T)`` for Neumann tags.
    order     : quadrature order.  Triangles: 1, 3, 4, 6, 7, 9, 12, 13
                (Gauss points on the triangle); quadrilaterals: number of
                Gauss-Legendre points per direction.  ``None`` selects the
                lowest order that fully integrates each element type; a
                lower order gives a rank-deficient, usually singular, system.
    material  : optional callable ``(x, y) -> (kxx, kxy, kyx, kyy, f)``
                overriding ``conductivity_and_forcing``.
    """

    a: float = 4.0
    b: float = 6.0
    elems: int = 50
    elem_type: object = 3
    mesh: Mesh | None = None
    bc_type: dict = field(default_factory=lambda: {
        "top": DIRICHLET, "right": NEUMANN, "left": NEUMANN, "bottom": DIRICHLET})
    bc_val: dict = field(default_factory=lambda: {
        "top": 100.0, "right": 0.0, "left": 0.0, "bottom": 0.0})
    order: int | None = None
    material: object = None

    def build_mesh(self) -> Mesh:
        if self.mesh is not None:
            return self.mesh
        return geometry(self.a, self.b, self.elems, self.elem_type)


@dataclass
class Solution:
    mesh: Mesh
    T: np.ndarray                # nodal temperatures, flat
    K: sp.csr_matrix             # stiffness matrix after boundary conditions
    F: np.ndarray                # load vector after boundary conditions
    material: object = None      # material used (None = default module)
    info: dict = field(default_factory=dict)   # solver statistics

    @property
    def grid(self) -> np.ndarray:
        """Temperatures on the ``(ny, nx)`` grid of ``mesh.X``/``mesh.Y``."""
        return self.mesh.to_grid(self.T)

    def flux(self):
        """Element-centroid gradient and flux, see :func:`aa540fem.postprocess.element_gradient`."""
        from .postprocess import element_gradient

        return element_gradient(self.mesh, self.T, self.material)

    def save(self, path, cell_fields: bool = True):
        """Write the solution to a ParaView file (``.vtu``, ``.vtk``, ...)."""
        from .mesh_io import write_vtk

        cell_data = None
        if cell_fields:
            cf = self.flux()
            cell_data = {"gradT": cf.grad, "flux": cf.flux}
        return write_vtk(path, self.mesh, point_data={"T": self.T}, cell_data=cell_data)


def assemble(mesh: Mesh, order: int | None = None, material=None, verbose: bool = False):
    """Assemble the global stiffness matrix ``K`` (CSR) and load vector ``F``.

    Every cell block of the mesh is integrated with its own element; ``order``
    applies to all blocks and ``None`` uses each element's ``min_order``.
    """
    N = mesh.n_nodes
    rows, cols, vals = [], [], []
    F = np.zeros(N)

    for name, conn in mesh.cells.items():
        el = get_element(name)
        o = el.min_order if order is None else order
        if o < el.min_order:
            warnings.warn(
                f"quadrature order {o} under-integrates {el.name} elements "
                f"(minimum {el.min_order}); the system may be singular",
                stacklevel=2)
        xi, eta, w = el.quadrature(o)
        phi, dphi_dxi, dphi_deta = el.shape(xi, eta)

        Ke, fe = elem_eqn(mesh.x[conn], mesh.y[conn], phi, dphi_dxi, dphi_deta, w,
                          material=material)

        n = el.n_nodes
        rows.append(np.repeat(conn, n, axis=1).ravel())   # conn[e, i] for each j
        cols.append(np.tile(conn, (1, n)).ravel())        # conn[e, j] for each i
        vals.append(Ke.ravel())
        np.add.at(F, conn.ravel(), fe.ravel())

    K = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(N, N)).tocsr()

    if verbose:
        print("Finished Assembling Global Stiffness Matrix.")
    return K, F


def apply_boundary_conditions(mesh: Mesh, K, F, bc_type: dict, bc_val: dict, verbose=False):
    """Apply Neumann fluxes first, then Dirichlet values (which win at corners).

    Tags in ``bc_type`` must exist in ``mesh.boundary``; mesh tags that are
    not mentioned are left natural (zero flux).
    """
    unknown = [t for t in bc_type if t not in mesh.boundary]
    if unknown:
        raise ValueError(f"Unknown boundary tag(s) {unknown}; mesh has {mesh.tags}")
    for tag, kind in bc_type.items():
        if kind not in (DIRICHLET, NEUMANN):
            raise ValueError(f"Unknown boundary condition type {kind!r} on {tag!r}")
        if tag not in bc_val:
            raise ValueError(f"No boundary value given for tag {tag!r}")

    for tag, kind in bc_type.items():
        if kind == NEUMANN:
            F = neumann(F, mesh.boundary[tag], bc_val[tag], mesh.x, mesh.y)
    for tag, kind in bc_type.items():
        if kind == DIRICHLET:
            K, F = dirichlet(K, F, mesh.bc_nodes[tag], bc_val[tag], mesh.x, mesh.y)
    if verbose:
        print("Boundary Conditions Applied.")
    return K, F


def _solve_cg(K, F, tol, maxiter, verbose):
    """Conjugate gradients, AMG-preconditioned when pyamg is available."""
    K = K.tocsr()
    try:
        import pyamg

        ml = pyamg.smoothed_aggregation_solver(K)
        M = ml.aspreconditioner(cycle="V")
        precond = "pyamg smoothed aggregation"
    except ImportError:
        d = K.diagonal()
        M = spla.LinearOperator(K.shape, matvec=lambda r: r / d)
        precond = "Jacobi"

    count = [0]

    def callback(_):
        count[0] += 1

    T, flag = spla.cg(K, F, M=M, rtol=tol, maxiter=maxiter, callback=callback)
    residual = np.linalg.norm(F - K @ T) / max(np.linalg.norm(F), 1e-300)
    if flag != 0:
        warnings.warn(f"CG did not converge in {count[0]} iterations (relative residual "
                      f"{residual:.2e})", stacklevel=3)
    if verbose:
        print(f"CG ({precond}): {count[0]} iterations, relative residual {residual:.2e}")
    return T, {"method": "cg", "preconditioner": precond, "iterations": count[0],
               "residual": residual, "converged": flag == 0}


def solve(problem: Problem, verbose: bool = False, method: str = "direct",
          tol: float = 1e-10, maxiter: int | None = None) -> Solution:
    """Mesh, assemble, apply boundary conditions and solve for the temperatures.

    ``method`` is ``"direct"`` (sparse LU) or ``"cg"`` (conjugate gradients with
    algebraic multigrid preconditioning when pyamg is installed).
    """
    mesh = problem.build_mesh()
    if verbose:
        print(f"Finished generating Grid: {mesh.n_nodes} nodes, {mesh.n_elems} elements.")

    K, F = assemble(mesh, problem.order, problem.material, verbose)
    K, F = apply_boundary_conditions(mesh, K, F, problem.bc_type, problem.bc_val, verbose)

    if not any(t == DIRICHLET for t in problem.bc_type.values()):
        raise ValueError("A pure Neumann problem is singular; fix T on at least one boundary")

    if verbose:
        print("Solving Equations ...")
    if method == "direct":
        T = spla.spsolve(K.tocsc(), F)
        info = {"method": "direct"}
    elif method == "cg":
        T, info = _solve_cg(K, F, tol, maxiter, verbose)
    else:
        raise ValueError(f"Unknown method {method!r}; expected 'direct' or 'cg'")
    return Solution(mesh, T, K, F, problem.material, info)
