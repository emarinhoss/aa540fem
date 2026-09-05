"""Definition of an incompressible flow case."""

from __future__ import annotations

from dataclasses import dataclass, field

from aa540fem.core.elements import PRESSURE_ELEMENT
from aa540fem.core.mesh import Mesh
from aa540fem.core.util import accepts_time, call_coeff

OPEN = "open"
SCHEMES = ("rk45", "theta")


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
    pspg        : add PSPG pressure stabilisation to the continuity equation
                  (not needed for the inf-sup stable Taylor-Hood pair; ignored
                  by the RK45 time integrator).
    eddy_viscosity : optional nodal array of the turbulent viscosity
                  ``mu_t`` (set by the RANS coupling); the momentum equations
                  then use ``mu + mu_t(x)`` with the extra term
                  ``- grad(u)^T . grad(mu_t)`` of the variable-viscosity
                  stress divergence.
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
    pspg: bool = False
    eddy_viscosity: object = None

    def validate(self):
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
            if not (isinstance(spec, (tuple, list)) and len(spec) == 2):
                raise ValueError(f"bc[{tag!r}] must be (ux, uy) or 'open'")

    @property
    def has_open_boundary(self) -> bool:
        return any(spec == OPEN for spec in self.bc.values()) or any(
            tag not in self.bc for tag in self.mesh.boundary)

    @property
    def pins_pressure(self) -> bool:
        return (not self.has_open_boundary) if self.pin_pressure is None else self.pin_pressure

    def depends_on_time(self) -> bool:
        if accepts_time(self.body_force):
            return True
        return any(accepts_time(v) for spec in self.bc.values() if spec != OPEN for v in spec)


def values_at_pair(fn, x, y, t=0.0):
    """Evaluate a ``(x, y[, t]) -> (a, b)`` callable."""
    return call_coeff(fn, x, y, t)
