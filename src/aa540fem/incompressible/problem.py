"""Definition of an incompressible flow case."""

from __future__ import annotations

from dataclasses import dataclass, field

from aa540fem.core.elements import PRESSURE_ELEMENT
from aa540fem.core.mesh import Mesh
from aa540fem.core.util import accepts_time, call_coeff_nd

OPEN = "open"
SCHEMES = ("rk45", "theta")
ELEMENT_LENGTHS = ("streamline", "metric")


@dataclass
class FlowProblem:
    """Definition of an incompressible flow case.

    Attributes
    ----------
    mesh        : :class:`Mesh` with ``quad9`` and/or ``triangle6`` blocks.
    mu, rho     : dynamic viscosity and density (constants).
    body_force  : optional ``(x, y[, t]) -> (fx, fy)`` acceleration.
    bc          : dict boundary tag -> ``(ux, uy)`` Dirichlet values (each a
                  constant, a callable ``(x, y[, t])`` or ``None`` for a free
                  component, e.g. ``(None, 0.0)`` for a symmetry line) or
                  ``"open"`` for the do-nothing outflow.  Unlisted tags are
                  open.  Where tags share a node the last listed wins.
    pin_pressure : fix the pressure at one node.  ``None`` pins it
                  automatically when no boundary is open (enclosed flow).
    pin_value   : value (constant or ``(x, y)``) of the pinned pressure.
    order       : quadrature order (``None``: element default).
    stabilisation : residual-based SUPG stabilisation of the momentum
                  equations plus grad-div (LSIC) stabilisation; consistent
                  (an exact solution stays exact), needed at high cell
                  Reynolds numbers.  Off by default (plain Galerkin).
    grad_div    : include the grad-div (LSIC) term in the stabilisation
                  (``gamma = h |u| / 2 min(1, Re_h / 3)``); it helps the steady
                  Newton iteration at high Reynolds number but adds numerical
                  dissipation to unsteady wakes, so switch it off for
                  vortex-shedding runs.
    pspg        : add PSPG pressure stabilisation to the continuity equation
                  (not needed for the inf-sup stable Taylor-Hood pair; ignored
                  by the RK45 time integrator).
    eddy_viscosity : optional nodal array of the turbulent viscosity
                  ``mu_t`` (set by the RANS coupling); the momentum equations
                  then use ``mu + mu_t(x)`` with the extra term
                  ``- grad(u)^T . grad(mu_t)`` of the variable-viscosity
                  stress divergence.
    element_length : how the stabilisation parameters measure the cell:
                  ``"metric"`` (default) uses the element metric tensor ``G``
                  (``tau = [(2/dt)^2 + u.G u + nu^2 G:G / 2]^(-1/2)``), which
                  is smooth in the velocity and lets Newton converge on
                  strongly stretched boundary-layer cells; ``"streamline"``
                  uses Tezduyar's element length in the flow direction,
                  ``h = 2 / sum_i |s . grad phi_i|``, which jumps between the
                  cell length and its height with the slightest rotation of
                  the velocity on such cells and stalls Newton there.
    """

    mesh: Mesh
    mu: float = 1.0
    rho: float = 1.0
    body_force: object = None
    bc: dict = field(default_factory=dict)
    pin_pressure: bool | None = None
    pin_value: object = 0.0
    order: int | None = None
    stabilisation: bool = False
    grad_div: bool = True
    pspg: bool = False
    eddy_viscosity: object = None
    element_length: str = "metric"

    def validate(self):
        if self.element_length not in ELEMENT_LENGTHS:
            raise ValueError(f"element_length must be one of {ELEMENT_LENGTHS}, "
                             f"got {self.element_length!r}")
        for name in self.mesh.cells:
            if name not in PRESSURE_ELEMENT:
                raise ValueError(f"Taylor-Hood needs quadratic elements; got {name!r} "
                                 f"(supported: {sorted(PRESSURE_ELEMENT)})")
        unknown = [t for t in self.bc if t not in self.mesh.boundary]
        if unknown:
            raise ValueError(f"Unknown boundary tag(s) {unknown}; mesh has {self.mesh.tags}")
        for tag, spec in self.bc.items():
            if spec == OPEN:
                continue
            if not (isinstance(spec, (tuple, list)) and len(spec) == self.mesh.dim):
                raise ValueError(f"bc[{tag!r}] must have {self.mesh.dim} velocity components "
                                 "or be 'open'")

    @property
    def has_open_boundary(self) -> bool:
        return any(spec == OPEN for spec in self.bc.values()) or any(
            tag not in self.bc for tag in self.mesh.boundary)

    @property
    def pins_pressure(self) -> bool:
        return (not self.has_open_boundary) if self.pin_pressure is None else self.pin_pressure

    def depends_on_time(self) -> bool:
        d = self.mesh.dim
        if accepts_time(self.body_force, d):
            return True
        return any(accepts_time(v, d) for spec in self.bc.values() if spec != OPEN for v in spec)


def values_at_pair(fn, coords, t=0.0):
    """Evaluate a ``(x, y[, z][, t]) -> (a, b[, c])`` callable at the coordinate tuple."""
    return call_coeff_nd(fn, coords, t)
