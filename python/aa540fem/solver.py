"""Assembly, boundary conditions and linear solve.  Port of ``main.m``'s
computational part, so that the driver script only holds the user inputs.

The steady problem solved here is

    u . grad T - div( kappa grad T ) = f

which reduces to the original ``div(kappa grad T) + f = 0`` without a
velocity.  :mod:`aa540fem.transient` adds the ``rho_c dT/dt`` term.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .boundary import DirichletEliminator, neumann
from .element import elem_operators
from .elements import get_element
from .geometry import Mesh, geometry
from .util import accepts_time, values_at

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
    bc_val    : dict boundary tag -> value (constant or callable
                ``val(x, y)`` / ``val(x, y, t)``): the temperature for
                Dirichlet tags, the normal flux ``n . (kappa grad T)`` for
                Neumann tags.  Where Dirichlet tags share a node the last
                listed tag wins.
    order     : quadrature order.  Triangles: 1, 3, 4, 6, 7, 9, 12, 13
                (Gauss points on the triangle); quadrilaterals: number of
                Gauss-Legendre points per direction.  ``None`` selects each
                element's ``full_order`` (exact stiffness, mass and
                convection matrices for constant coefficients); an order
                below ``min_order`` gives a rank-deficient system.
    material  : optional callable ``(x, y[, t]) -> (kxx, kxy, kyx, kyy, f)``
                overriding ``conductivity_and_forcing``.
    velocity  : optional callable ``(x, y[, t]) -> (ux, uy)`` adding the
                convection term ``u . grad T``.
    rho_c     : volumetric heat capacity (constant or callable), used by the
                transient solver only.
    supg      : streamline-upwind Petrov-Galerkin stabilisation of the
                convection term (ignored without a velocity).
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
    velocity: object = None
    rho_c: object = 1.0
    supg: bool = True

    def build_mesh(self) -> Mesh:
        if self.mesh is not None:
            return self.mesh
        return geometry(self.a, self.b, self.elems, self.elem_type)

    @property
    def symmetric(self) -> bool:
        """True if the assembled system is symmetric (no convection)."""
        return self.velocity is None

    def has_dirichlet(self) -> bool:
        return any(kind == DIRICHLET for kind in self.bc_type.values())

    def validate(self, mesh: Mesh):
        unknown = [t for t in self.bc_type if t not in mesh.boundary]
        if unknown:
            raise ValueError(f"Unknown boundary tag(s) {unknown}; mesh has {mesh.tags}")
        for tag, kind in self.bc_type.items():
            if kind not in (DIRICHLET, NEUMANN):
                raise ValueError(f"Unknown boundary condition type {kind!r} on {tag!r}")
            if tag not in self.bc_val:
                raise ValueError(f"No boundary value given for tag {tag!r}")

    def operators_depend_on_time(self) -> bool:
        return any(accepts_time(fn) for fn in (self.material, self.velocity, self.rho_c))

    def loads_depend_on_time(self) -> bool:
        return self.operators_depend_on_time() or any(
            accepts_time(v) for v in self.bc_val.values())


@dataclass
class Operators:
    """Assembled global matrices (CSR) and load vector."""

    K: sp.csr_matrix     # diffusion
    C: sp.csr_matrix     # convection (+ SUPG)
    M: sp.csr_matrix     # mass (+ SUPG)
    F: np.ndarray        # source (+ SUPG)

    @property
    def A(self) -> sp.csr_matrix:
        """Steady operator ``K + C``."""
        return (self.K + self.C).tocsr()


@dataclass
class Solution:
    mesh: Mesh
    T: np.ndarray                # nodal temperatures, flat
    K: sp.csr_matrix             # system matrix after boundary conditions
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


# ------------------------------------------------------------------ assembly
def assemble_operators(mesh: Mesh, problem: Problem, t: float = 0.0, dt=None,
                       verbose: bool = False) -> Operators:
    """Assemble diffusion, convection and mass matrices and the load vector.

    Every cell block of the mesh is integrated with its own element;
    ``problem.order`` applies to all blocks and ``None`` uses each element's
    ``full_order``.  ``t`` is passed to time-dependent coefficients and
    ``dt`` selects the transient SUPG parameter.
    """
    N = mesh.n_nodes
    rows, cols = [], []
    kv, cv, mv = [], [], []
    F = np.zeros(N)

    for name, conn in mesh.cells.items():
        el = get_element(name)
        o = el.full_order if problem.order is None else problem.order
        if o < el.min_order:
            warnings.warn(
                f"quadrature order {o} under-integrates {el.name} elements "
                f"(minimum {el.min_order}); the system may be singular",
                stacklevel=2)
        xi, eta, w = el.quadrature(o)
        phi, dphi_dxi, dphi_deta = el.shape(xi, eta)

        Ke, Ce, Me, fe = elem_operators(
            mesh.x[conn], mesh.y[conn], phi, dphi_dxi, dphi_deta, w,
            material=problem.material, velocity=problem.velocity, rho_c=problem.rho_c,
            supg=problem.supg, t=t, dt=dt)

        n = el.n_nodes
        rows.append(np.repeat(conn, n, axis=1).ravel())   # conn[e, i] for each j
        cols.append(np.tile(conn, (1, n)).ravel())        # conn[e, j] for each i
        kv.append(Ke.ravel())
        cv.append(Ce.ravel())
        mv.append(Me.ravel())
        np.add.at(F, conn.ravel(), fe.ravel())

    rows = np.concatenate(rows)
    cols = np.concatenate(cols)

    def mat(vals):
        m = sp.coo_matrix((np.concatenate(vals), (rows, cols)), shape=(N, N)).tocsr()
        m.eliminate_zeros()
        return m

    ops = Operators(mat(kv), mat(cv), mat(mv), F)
    if verbose:
        print("Finished Assembling Global Stiffness Matrix.")
    return ops


def assemble(mesh: Mesh, order: int | None = None, material=None, verbose: bool = False):
    """Assemble the diffusion matrix ``K`` (CSR) and load vector ``F`` only."""
    ops = assemble_operators(mesh, Problem(mesh=mesh, order=order, material=material),
                             verbose=verbose)
    return ops.K, ops.F


# ------------------------------------------------------- boundary conditions
def neumann_loads(mesh: Mesh, F, bc_type: dict, bc_val: dict, t: float = 0.0):
    """Add the prescribed fluxes of all Neumann tags to ``F``."""
    for tag, kind in bc_type.items():
        if kind == NEUMANN:
            F = neumann(F, mesh.boundary[tag], bc_val[tag], mesh.x, mesh.y, t=t)
    return F


def dirichlet_data(mesh: Mesh, bc_type: dict, bc_val: dict, t: float = 0.0):
    """Constrained nodes and their values; the last listed tag wins at shared nodes."""
    values = {}
    for tag, kind in bc_type.items():
        if kind == DIRICHLET:
            nodes = mesh.bc_nodes[tag]
            vals = values_at(bc_val[tag], mesh.x[nodes], mesh.y[nodes], t)
            values.update(zip(nodes.tolist(), vals.tolist()))
    nodes = np.array(sorted(values), dtype=int)
    return nodes, np.array([values[n] for n in nodes])


def apply_boundary_conditions(mesh: Mesh, K, F, bc_type: dict, bc_val: dict,
                              verbose: bool = False, t: float = 0.0):
    """Apply Neumann fluxes and then eliminate the Dirichlet values.

    Tags in ``bc_type`` must exist in ``mesh.boundary``; mesh tags that are
    not mentioned are left natural (zero flux).
    """
    Problem(mesh=mesh, bc_type=bc_type, bc_val=bc_val).validate(mesh)
    F = neumann_loads(mesh, F, bc_type, bc_val, t)
    nodes, vals = dirichlet_data(mesh, bc_type, bc_val, t)
    elim = DirichletEliminator(K, nodes)
    if verbose:
        print("Boundary Conditions Applied.")
    return elim.K_bc, elim.apply_rhs(F, vals)


# --------------------------------------------------------------- linear solve
class LinearSolver:
    """Direct or preconditioned Krylov solver for repeated solves with one matrix.

    ``method``: ``"direct"`` (sparse LU, factorised once), ``"cg"`` (symmetric
    systems only) or ``"gmres"``; the Krylov methods use pyamg smoothed
    aggregation as preconditioner when available, otherwise Jacobi (cg) or
    incomplete LU (gmres).
    """

    def __init__(self, A, method: str = "direct", tol: float = 1e-10, maxiter=None,
                 symmetric: bool = True):
        self.method = method
        self.tol = tol
        self.maxiter = maxiter
        self.A = sp.csr_matrix(A)
        self.precond = None
        if method == "direct":
            self.lu = spla.splu(self.A.tocsc())
        elif method == "cg":
            if not symmetric:
                raise ValueError("cg needs a symmetric system; use method='gmres' "
                                 "for convection problems")
            self.M = self._amg("symmetric") or self._jacobi()
        elif method == "gmres":
            self.M = self._amg("nonsymmetric") or self._ilu()
        else:
            raise ValueError(f"Unknown method {method!r}; expected 'direct', 'cg' or 'gmres'")

    def _amg(self, symmetry):
        try:
            import pyamg
        except ImportError:
            return None
        ml = pyamg.smoothed_aggregation_solver(self.A, symmetry=symmetry)
        self.precond = "pyamg smoothed aggregation"
        return ml.aspreconditioner(cycle="V")

    def _jacobi(self):
        d = self.A.diagonal()
        self.precond = "Jacobi"
        return spla.LinearOperator(self.A.shape, matvec=lambda r: r / d)

    def _ilu(self):
        ilu = spla.spilu(self.A.tocsc(), drop_tol=1e-4, fill_factor=10)
        self.precond = "ILU"
        return spla.LinearOperator(self.A.shape, matvec=ilu.solve)

    def solve(self, F, verbose: bool = False):
        F = np.asarray(F, dtype=float)
        if self.method == "direct":
            return self.lu.solve(F), {"method": "direct"}

        count = [0]

        def callback(_):
            count[0] += 1

        if self.method == "cg":
            T, flag = spla.cg(self.A, F, M=self.M, rtol=self.tol, maxiter=self.maxiter,
                              callback=callback)
        else:
            T, flag = spla.gmres(self.A, F, M=self.M, rtol=self.tol, maxiter=self.maxiter,
                                 restart=50, callback=callback, callback_type="pr_norm")
        residual = np.linalg.norm(F - self.A @ T) / max(np.linalg.norm(F), 1e-300)
        if flag != 0:
            warnings.warn(f"{self.method} did not converge in {count[0]} iterations "
                          f"(relative residual {residual:.2e})", stacklevel=2)
        if verbose:
            print(f"{self.method} ({self.precond}): {count[0]} iterations, "
                  f"relative residual {residual:.2e}")
        return T, {"method": self.method, "preconditioner": self.precond,
                   "iterations": count[0], "residual": residual, "converged": flag == 0}


def solve(problem: Problem, verbose: bool = False, method: str = "direct",
          tol: float = 1e-10, maxiter: int | None = None) -> Solution:
    """Mesh, assemble, apply boundary conditions and solve the steady problem.

    ``method`` is ``"direct"`` (sparse LU), ``"cg"`` (conjugate gradients,
    symmetric problems only) or ``"gmres"``; both Krylov methods use pyamg
    preconditioning when it is installed.
    """
    mesh = problem.build_mesh()
    problem.validate(mesh)
    if not problem.has_dirichlet():
        raise ValueError("A pure Neumann problem is singular; fix T on at least one boundary")
    if verbose:
        print(f"Finished generating Grid: {mesh.n_nodes} nodes, {mesh.n_elems} elements.")

    ops = assemble_operators(mesh, problem, verbose=verbose)
    K, F = apply_boundary_conditions(mesh, ops.A, ops.F, problem.bc_type, problem.bc_val,
                                     verbose)

    if verbose:
        print("Solving Equations ...")
    solver = LinearSolver(K, method, tol, maxiter, symmetric=problem.symmetric)
    T, info = solver.solve(F, verbose)
    return Solution(mesh, T, K, F, problem.material, info)
