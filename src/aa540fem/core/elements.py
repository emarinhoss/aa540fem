"""Reference-element registry.

A :class:`ReferenceElement` bundles everything the assembly needs to know
about one element type: dimension, node count, shape functions, quadrature,
the local node lists of its faces and edges, and how to reverse its
orientation.  Elements are named after the meshio cell types (``triangle``,
``triangle6``, ``quad``, ``quad9`` in 2D; ``tetra``, ``tetra10``,
``hexahedron``, ``hexahedron27``, ``wedge``, ``wedge18`` in 3D, in meshio's VTK
local ordering); the
MATLAB element numbers 1, 2 and 3 are accepted as aliases for ``triangle``,
``quad`` and ``quad9``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from aa540fem.core.quadrature import quadrature_rule
from aa540fem.core.shape_functions import (
    HEX27_NODES,
    TET10_EDGES,
    WEDGE6_NODES,
    WEDGE18_NODES,
    hessian_3,
    hessian_4,
    hessian_6,
    hessian_9,
    hessian_hex8,
    hessian_hex27,
    hessian_tet4,
    hessian_tet10,
    hessian_wedge6,
    hessian_wedge18,
    interpfunc_3,
    interpfunc_4,
    interpfunc_6,
    interpfunc_9,
    interpfunc_hex8,
    interpfunc_hex27,
    interpfunc_tet4,
    interpfunc_tet10,
    interpfunc_wedge6,
    interpfunc_wedge18,
)


@dataclass(frozen=True)
class ReferenceElement:
    """Description of a reference element (2-D or 3-D).

    Attributes
    ----------
    name       : meshio cell type.
    family     : ``"triangle"``, ``"quad"``, ``"tetra"``, ``"hexahedron"`` or
                 ``"wedge"`` (selects the quadrature rule).
    n_nodes    : nodes per element.
    min_order  : lowest quadrature order that fully integrates the stiffness
                 matrix for constant conductivity.
    full_order : quadrature order that also integrates the mass and
                 convection matrices exactly (the default for assembly).
    faces      : local node indices of each face: in 2D an edge ordered
                 ``(start, end[, mid])`` with the element on the left when
                 walking from start to end; in 3D a triangle or quadrilateral
                 in the face element's own ordering with the outward normal.
    face_type  : cell type of a face (``"line"``, ``"line3"``, ``"triangle6"``,
                 ``"quad9"``, ...); prisms have two face types, see ``face_types``.
    face_types : cell type of every face (defaults to ``face_type`` for all).
    reverse    : node permutation that flips the element orientation.
    centroid   : natural coordinates of the element centroid.
    nodes      : natural coordinates of the nodes.
    legacy_id  : MATLAB element number, if any.
    dim        : 2 or 3.
    n_corners  : corner (vertex) nodes, the first ``n_corners`` local nodes.
    edges      : corner pairs of the element edges (cell size measures).
    """

    name: str
    family: str
    n_nodes: int
    min_order: int
    full_order: int
    faces: tuple
    face_type: str
    reverse: tuple
    centroid: tuple
    nodes: tuple
    shape_fn: Callable = field(repr=False, compare=False)
    hessian_fn: Callable = field(repr=False, compare=False)
    legacy_id: int | None = None
    dim: int = 2
    n_corners: int = 0
    edges: tuple = ()
    face_types: tuple = ()

    def __post_init__(self):
        if not self.face_types:
            object.__setattr__(self, "face_types", (self.face_type,) * len(self.faces))

    def shape(self, *coords):
        """Shape functions and natural derivatives, each ``(nq, n_nodes)``:
        ``(phi, dphi_dxi, dphi_deta[, dphi_dzeta])``."""
        return self.shape_fn(*coords)

    def hessian(self, *coords):
        """Second natural derivatives, each ``(nq, n)``: ``(xx, xy, yy)`` in 2D,
        ``(xx, xy, xz, yy, yz, zz)`` in 3D."""
        return self.hessian_fn(*coords)

    def shape_at(self, coords):
        """``(phi, (dphi_1, ..., dphi_d))`` for a tuple of natural-coordinate arrays."""
        out = self.shape_fn(*coords)
        return out[0], tuple(out[1:])

    def quadrature(self, order: int | None = None):
        """Quadrature points and weights; ``None`` selects ``full_order``."""
        if order is None:
            order = self.full_order
        return quadrature_rule(self.family, order)

    @property
    def default_order(self) -> int:
        return self.full_order

    @property
    def nodes_per_face(self) -> int:
        return len(self.faces[0])

    def face_elements(self):
        """``(face, face element)`` pairs, the face element being the reference
        element of the face's own type."""
        return [(f, get_element(t)) for f, t in zip(self.faces, self.face_types)]

    @property
    def n_faces(self) -> int:
        return len(self.faces)


def _cyclic_edges(n):
    return tuple((i, (i + 1) % n) for i in range(n))


_third = 1.0 / 3.0

TRIANGLE = ReferenceElement(
    name="triangle", family="triangle", n_nodes=3, min_order=1, full_order=3,
    faces=((0, 1), (1, 2), (2, 0)), face_type="line",
    reverse=(0, 2, 1), centroid=(_third, _third),
    nodes=((0, 0), (1, 0), (0, 1)),
    shape_fn=interpfunc_3, hessian_fn=hessian_3, legacy_id=1,
    n_corners=3, edges=_cyclic_edges(3),
)

TRIANGLE6 = ReferenceElement(
    name="triangle6", family="triangle", n_nodes=6, min_order=3, full_order=6,
    faces=((0, 1, 3), (1, 2, 4), (2, 0, 5)), face_type="line3",
    reverse=(0, 2, 1, 5, 4, 3), centroid=(_third, _third),
    nodes=((0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)),
    shape_fn=interpfunc_6, hessian_fn=hessian_6, n_corners=3, edges=_cyclic_edges(3),
)

QUAD = ReferenceElement(
    name="quad", family="quad", n_nodes=4, min_order=2, full_order=2,
    faces=((0, 1), (1, 2), (2, 3), (3, 0)), face_type="line",
    reverse=(0, 3, 2, 1), centroid=(0.0, 0.0),
    nodes=((-1, -1), (1, -1), (1, 1), (-1, 1)),
    shape_fn=interpfunc_4, hessian_fn=hessian_4, legacy_id=2,
    n_corners=4, edges=_cyclic_edges(4),
)

QUAD9 = ReferenceElement(
    name="quad9", family="quad", n_nodes=9, min_order=3, full_order=3,
    faces=((0, 1, 4), (1, 2, 5), (2, 3, 6), (3, 0, 7)), face_type="line3",
    reverse=(0, 3, 2, 1, 7, 6, 5, 4, 8), centroid=(0.0, 0.0),
    nodes=((-1, -1), (1, -1), (1, 1), (-1, 1), (0, -1), (1, 0), (0, 1), (-1, 0), (0, 0)),
    shape_fn=interpfunc_9, hessian_fn=hessian_9, legacy_id=3,
    n_corners=4, edges=_cyclic_edges(4),
)

# -- 3-D elements (meshio / VTK ordering) ----------------------------------
_TET_EDGE_MID = {frozenset(e): 4 + i for i, e in enumerate(TET10_EDGES)}
_TET_FACES = ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2))          # outward normals


def _tet10_face(corners):
    a, b, c = corners
    return (a, b, c, _TET_EDGE_MID[frozenset((a, b))], _TET_EDGE_MID[frozenset((b, c))],
            _TET_EDGE_MID[frozenset((c, a))])


_TET_NODES = tuple(map(tuple, [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]))
_TET10_NODES = _TET_NODES + tuple(tuple(0.5 * (np.array(_TET_NODES[a]) + np.array(_TET_NODES[b])))
                                  for a, b in TET10_EDGES)

TETRA = ReferenceElement(
    name="tetra", family="tetra", n_nodes=4, min_order=1, full_order=2,
    faces=_TET_FACES, face_type="triangle", reverse=(0, 2, 1, 3),
    centroid=(0.25, 0.25, 0.25), nodes=_TET_NODES,
    shape_fn=interpfunc_tet4, hessian_fn=hessian_tet4, dim=3, n_corners=4, edges=TET10_EDGES,
)

TETRA10 = ReferenceElement(
    name="tetra10", family="tetra", n_nodes=10, min_order=2, full_order=4,
    faces=tuple(_tet10_face(f) for f in _TET_FACES), face_type="triangle6",
    reverse=(0, 2, 1, 3, 6, 5, 4, 7, 9, 8), centroid=(0.25, 0.25, 0.25), nodes=_TET10_NODES,
    shape_fn=interpfunc_tet10, hessian_fn=hessian_tet10, dim=3, n_corners=4, edges=TET10_EDGES,
)

_HEX_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
              (0, 4), (1, 5), (2, 6), (3, 7))
_HEX_EDGE_MID = {frozenset(e): 8 + i for i, e in enumerate(_HEX_EDGES)}
# outward faces (corner quads) and the index of their centre node
_HEX_FACES = (((0, 3, 2, 1), 24), ((0, 1, 5, 4), 22), ((1, 2, 6, 5), 21),
              ((2, 3, 7, 6), 23), ((0, 4, 7, 3), 20), ((4, 5, 6, 7), 25))


def _hex27_face(corners, centre):
    mids = tuple(_HEX_EDGE_MID[frozenset((corners[i], corners[(i + 1) % 4]))] for i in range(4))
    return corners + mids + (centre,)


HEXAHEDRON = ReferenceElement(
    name="hexahedron", family="hexahedron", n_nodes=8, min_order=2, full_order=2,
    faces=tuple(f for f, _ in _HEX_FACES), face_type="quad",
    reverse=(1, 0, 3, 2, 5, 4, 7, 6), centroid=(0.0, 0.0, 0.0),
    nodes=tuple(map(tuple, HEX27_NODES[:8])),
    shape_fn=interpfunc_hex8, hessian_fn=hessian_hex8, dim=3, n_corners=8, edges=_HEX_EDGES,
)

HEXAHEDRON27 = ReferenceElement(
    name="hexahedron27", family="hexahedron", n_nodes=27, min_order=3, full_order=3,
    faces=tuple(_hex27_face(f, c) for f, c in _HEX_FACES), face_type="quad9",
    reverse=(1, 0, 3, 2, 5, 4, 7, 6, 8, 11, 10, 9, 12, 15, 14, 13, 17, 16, 19, 18,
             21, 20, 22, 23, 24, 25, 26),
    centroid=(0.0, 0.0, 0.0), nodes=tuple(map(tuple, HEX27_NODES)),
    shape_fn=interpfunc_hex27, hessian_fn=hessian_hex27, dim=3, n_corners=8, edges=_HEX_EDGES,
)

# prisms: bottom and top triangles (outward normals -z, +z) and three quadrilaterals
_WEDGE_EDGES = ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5))
_WEDGE_EDGE_MID = {frozenset(e): 6 + i for i, e in enumerate(_WEDGE_EDGES)}
_WEDGE_TRI_FACES = ((0, 2, 1), (3, 4, 5))
_WEDGE_QUAD_FACES = (((0, 1, 4, 3), 15), ((1, 2, 5, 4), 16), ((2, 0, 3, 5), 17))


def _wedge18_faces():
    faces = []
    for c in _WEDGE_TRI_FACES:
        faces.append(c + tuple(_WEDGE_EDGE_MID[frozenset((c[i], c[(i + 1) % 3]))]
                               for i in range(3)))
    for c, centre in _WEDGE_QUAD_FACES:
        faces.append(c + tuple(_WEDGE_EDGE_MID[frozenset((c[i], c[(i + 1) % 4]))]
                               for i in range(4)) + (centre,))
    return tuple(faces)


WEDGE = ReferenceElement(
    name="wedge", family="wedge", n_nodes=6, min_order=1, full_order=3,
    faces=_WEDGE_TRI_FACES + tuple(c for c, _ in _WEDGE_QUAD_FACES), face_type="triangle",
    face_types=("triangle", "triangle", "quad", "quad", "quad"),
    reverse=(0, 2, 1, 3, 5, 4), centroid=(_third, _third, 0.0),
    nodes=tuple(map(tuple, WEDGE6_NODES)),
    shape_fn=interpfunc_wedge6, hessian_fn=hessian_wedge6, dim=3, n_corners=6,
    edges=_WEDGE_EDGES,
)

WEDGE18 = ReferenceElement(
    name="wedge18", family="wedge", n_nodes=18, min_order=3, full_order=6,
    faces=_wedge18_faces(), face_type="triangle6",
    face_types=("triangle6", "triangle6", "quad9", "quad9", "quad9"),
    reverse=(0, 2, 1, 3, 5, 4, 8, 7, 6, 11, 10, 9, 12, 14, 13, 17, 16, 15),
    centroid=(_third, _third, 0.0), nodes=tuple(map(tuple, WEDGE18_NODES)),
    shape_fn=interpfunc_wedge18, hessian_fn=hessian_wedge18, dim=3, n_corners=6,
    edges=_WEDGE_EDGES,
)

ELEMENTS = {e.name: e for e in (TRIANGLE, TRIANGLE6, QUAD, QUAD9, TETRA, TETRA10, HEXAHEDRON,
                                HEXAHEDRON27, WEDGE, WEDGE18)}

# Taylor-Hood pairs: quadratic velocity element -> linear pressure element on
# its corner nodes (the first n_corners local nodes).
PRESSURE_ELEMENT = {"triangle6": TRIANGLE, "quad9": QUAD, "tetra10": TETRA,
                    "hexahedron27": HEXAHEDRON, "wedge18": WEDGE}
LEGACY_ELEMENTS = {e.legacy_id: e for e in ELEMENTS.values() if e.legacy_id is not None}
# boundary cell types (nodes per face) of 2-D and 3-D meshes
FACE_TYPES = {2: {"line": 2, "line3": 3},
              3: {"triangle": 3, "triangle6": 6, "quad": 4, "quad9": 9}}
FACE_ELEMENTS = FACE_TYPES[2]
FACE_CORNERS = {"line": 2, "line3": 2, "triangle": 3, "triangle6": 3, "quad": 4, "quad9": 4}
# face type of a 3-D boundary face from its node count
FACE_TYPE_BY_NODES = {n: t for t, n in FACE_TYPES[3].items()}


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
