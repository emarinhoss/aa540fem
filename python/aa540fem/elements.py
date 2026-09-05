"""Reference-element registry.

A :class:`ReferenceElement` bundles everything the assembly needs to know
about one element type: node count, shape functions, quadrature, the local
node lists of its edges, and how to reverse its orientation.  Elements are
named after the meshio / Gmsh cell types (``triangle``, ``triangle6``,
``quad``, ``quad9``); the MATLAB element numbers 1, 2 and 3 are accepted as
aliases for ``triangle``, ``quad`` and ``quad9``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .quadrature import quadrature_rule
from .shape_functions import interpfunc_3, interpfunc_4, interpfunc_6, interpfunc_9


@dataclass(frozen=True)
class ReferenceElement:
    """Description of a 2-D reference element.

    Attributes
    ----------
    name       : meshio cell type.
    family     : ``"triangle"`` or ``"quad"`` (selects the quadrature rule).
    n_nodes    : nodes per element.
    min_order  : lowest quadrature order that fully integrates the stiffness
                 matrix for constant conductivity.
    faces      : local node indices of each edge, ordered ``(start, end[, mid])``
                 with the element on the left when walking from start to end.
    face_type  : cell type of an edge, ``"line"`` or ``"line3"``.
    reverse    : node permutation that flips the element orientation.
    centroid   : natural coordinates of the element centroid.
    legacy_id  : MATLAB element number, if any.
    """

    name: str
    family: str
    n_nodes: int
    min_order: int
    faces: tuple
    face_type: str
    reverse: tuple
    centroid: tuple
    shape_fn: Callable = field(repr=False, compare=False)
    legacy_id: int | None = None

    def shape(self, xi, eta):
        """Shape functions and natural derivatives, each ``(nq, n_nodes)``."""
        return self.shape_fn(xi, eta)

    def quadrature(self, order: int | None = None):
        """Quadrature points and weights; ``None`` selects ``min_order``."""
        if order is None:
            order = self.min_order
        return quadrature_rule(self.family, order)

    @property
    def default_order(self) -> int:
        return self.min_order

    @property
    def nodes_per_face(self) -> int:
        return len(self.faces[0])

    @property
    def n_corners(self) -> int:
        return 3 if self.family == "triangle" else 4


_third = 1.0 / 3.0

TRIANGLE = ReferenceElement(
    name="triangle", family="triangle", n_nodes=3, min_order=1,
    faces=((0, 1), (1, 2), (2, 0)), face_type="line",
    reverse=(0, 2, 1), centroid=(_third, _third),
    shape_fn=interpfunc_3, legacy_id=1,
)

TRIANGLE6 = ReferenceElement(
    name="triangle6", family="triangle", n_nodes=6, min_order=3,
    faces=((0, 1, 3), (1, 2, 4), (2, 0, 5)), face_type="line3",
    reverse=(0, 2, 1, 5, 4, 3), centroid=(_third, _third),
    shape_fn=interpfunc_6,
)

QUAD = ReferenceElement(
    name="quad", family="quad", n_nodes=4, min_order=2,
    faces=((0, 1), (1, 2), (2, 3), (3, 0)), face_type="line",
    reverse=(0, 3, 2, 1), centroid=(0.0, 0.0),
    shape_fn=interpfunc_4, legacy_id=2,
)

QUAD9 = ReferenceElement(
    name="quad9", family="quad", n_nodes=9, min_order=3,
    faces=((0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7)), face_type="line3",
    reverse=(0, 3, 2, 1, 7, 6, 5, 4, 8), centroid=(0.0, 0.0),
    shape_fn=interpfunc_9, legacy_id=3,
)

ELEMENTS = {e.name: e for e in (TRIANGLE, TRIANGLE6, QUAD, QUAD9)}
LEGACY_ELEMENTS = {e.legacy_id: e for e in ELEMENTS.values() if e.legacy_id is not None}
FACE_ELEMENTS = {"line": 2, "line3": 3}


def get_element(key) -> ReferenceElement:
    """Look up an element by name, MATLAB number (1, 2, 3) or instance."""
    if isinstance(key, ReferenceElement):
        return key
    if isinstance(key, (int, np.integer)) and not isinstance(key, bool):
        try:
            return LEGACY_ELEMENTS[int(key)]
        except KeyError:
            raise ValueError(f"Unknown element type {key!r}; expected 1, 2 or 3") from None
    try:
        return ELEMENTS[key]
    except (KeyError, TypeError):
        raise ValueError(
            f"Unknown element type {key!r}; expected one of {sorted(ELEMENTS)}") from None
