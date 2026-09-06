"""Mesh data structure and structured rectangle generator.

:class:`Mesh` holds nodes, one connectivity block per element type and named
boundary edge sets, so that both the structured rectangle of the original
MATLAB code (:func:`geometry`) and unstructured meshes read from files
(:mod:`aa540fem.mesh_io`) go through the same solver.

Node numbering of the structured mesh is row-major with ``x`` varying
fastest: node ``k`` sits at ``(xx[k % nx], yy[k // nx])``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aa540fem.core.elements import ELEMENTS, FACE_CORNERS, FACE_ELEMENTS, FACE_TYPES, get_element

SIDES = ("top", "right", "left", "bottom")
BOX_SIDES = ("left", "right", "bottom", "top", "front", "back")     # x=0, x=a, y=0, y=b, z=0, z=c


@dataclass
class Mesh:
    """Nodes, connectivity and boundary information of a 2-D or 3-D mesh.

    Attributes
    ----------
    points   : ``(n_nodes, d)`` nodal coordinates, ``d`` = 2 or 3.
    cells    : dict cell type -> ``(n_elems, nodes_per_element)`` 0-based
               connectivity in meshio/VTK local ordering (see
               :mod:`aa540fem.core.elements`).
    boundary : dict tag -> ``(n_faces, k)`` boundary faces: in 2D edges
               ordered ``(start, end[, mid])`` with the domain on the left,
               in 3D triangles or quadrilaterals (``k`` = 3, 6, 4 or 9) in the
               face element's ordering with the outward normal.
    X, Y, nx, ny : grid form of the coordinates, only for structured 2-D meshes.
    """

    points: np.ndarray
    cells: dict
    boundary: dict
    X: np.ndarray | None = None
    Y: np.ndarray | None = None
    nx: int | None = None
    ny: int | None = None

    def __post_init__(self):
        self.points = np.asarray(self.points, dtype=float)
        if self.points.ndim != 2 or self.points.shape[1] not in (2, 3):
            raise ValueError("points must be an (n, 2) or (n, 3) array")
        self.cells = {name: np.asarray(c, dtype=int) for name, c in self.cells.items()}
        dims = set()
        for name, conn in self.cells.items():
            el = get_element(name)
            if conn.ndim != 2 or conn.shape[1] != el.n_nodes:
                raise ValueError(f"{name} connectivity must be (n, {el.n_nodes})")
            dims.add(el.dim)
        if len(dims) > 1:
            raise ValueError("cannot mix 2-D and 3-D cell blocks")
        if dims and dims != {self.points.shape[1]}:
            raise ValueError(f"{dims.pop()}-D cells need points with that many coordinates")
        empty = np.zeros((0, 2 if self.dim == 2 else 3), dtype=int)
        self.boundary = {tag: np.asarray(e, dtype=int).reshape(-1, np.asarray(e).shape[-1])
                         if np.asarray(e).size else empty
                         for tag, e in self.boundary.items()}

    # ------------------------------------------------------------ basic info
    @property
    def x(self) -> np.ndarray:
        return self.points[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.points[:, 1]

    @property
    def z(self) -> np.ndarray:
        if self.dim < 3:
            raise AttributeError("a 2-D mesh has no z coordinate")
        return self.points[:, 2]

    @property
    def dim(self) -> int:
        """Spatial dimension (2 or 3)."""
        return self.points.shape[1]

    @property
    def coords(self) -> tuple:
        """The coordinate columns ``(x, y[, z])``."""
        return tuple(self.points[:, k] for k in range(self.dim))

    @property
    def n_nodes(self) -> int:
        return self.points.shape[0]

    @property
    def n_elems(self) -> int:
        return sum(c.shape[0] for c in self.cells.values())

    @property
    def elements(self) -> dict:
        """dict cell type -> :class:`ReferenceElement` for each block."""
        return {name: get_element(name) for name in self.cells}

    @property
    def tags(self) -> list:
        return list(self.boundary)

    @property
    def is_structured(self) -> bool:
        return self.X is not None

    # ------------------------------------------- single-block compatibility
    def _single(self):
        if len(self.cells) != 1:
            raise ValueError("mesh has several cell blocks; use mesh.cells")
        return next(iter(self.cells.items()))

    @property
    def conn(self) -> np.ndarray:
        """Connectivity of the (single) cell block."""
        return self._single()[1]

    @property
    def elem_type(self):
        """Element name of the (single) cell block."""
        return self._single()[0]

    @property
    def nodes_per_element(self) -> int:
        return get_element(self.elem_type).n_nodes

    @property
    def bc_edges(self) -> dict:
        return self.boundary

    @property
    def bc_nodes(self) -> dict:
        """dict tag -> sorted unique node indices on that boundary (cached; the
        cache is rebuilt when the set of tags or their edge counts change)."""
        stamp = tuple((tag, edges.shape) for tag, edges in self.boundary.items())
        cache = self.__dict__.get("_bc_nodes_cache")
        if cache is None or cache[0] != stamp:
            cache = (stamp, {tag: np.unique(edges) for tag, edges in self.boundary.items()})
            self.__dict__["_bc_nodes_cache"] = cache
        return cache[1]

    # ---------------------------------------------------------------- utils
    def to_grid(self, values) -> np.ndarray:
        """Reshape a nodal vector to the ``(ny, nx)`` grid of ``X``/``Y``."""
        if not self.is_structured:
            raise ValueError("to_grid is only available for structured meshes")
        return np.asarray(values).reshape(self.ny, self.nx)

    def jacobian_at_centroids(self) -> dict:
        """dict cell type -> Jacobian determinant of each element at its centroid."""
        out = {}
        for name, conn in self.cells.items():
            el = get_element(name)
            _, dnat = el.shape_at(tuple([c] for c in el.centroid))
            d = self.dim
            J = np.stack([np.stack([self.points[conn, i] @ dnat[j].T for j in range(d)], -1)
                          for i in range(d)], -2)                  # (ne, 1, d, d)
            out[name] = np.linalg.det(J).ravel()
        return out

    def check_orientation(self) -> int:
        """Reverse elements with a negative Jacobian so all are positively
        oriented (counter-clockwise in 2D).

        Returns the number of elements that were flipped.
        """
        flipped = 0
        for name, det in self.jacobian_at_centroids().items():
            neg = det < 0
            if neg.any():
                el = get_element(name)
                conn = self.cells[name]
                conn[neg] = conn[neg][:, list(el.reverse)]
                flipped += int(neg.sum())
        return flipped

    def boundary_faces(self) -> np.ndarray:
        """Faces that belong to exactly one element, i.e. the outer boundary.

        Returned as an ``(n_faces, k)`` array in the same ordering as
        ``boundary`` entries (domain on the left in 2D, outward normal in 3D).
        Mixed face types (triangles and quadrilaterals of a 3-D mesh with both
        tetrahedra and hexahedra) are returned as a dict face type -> array.
        """
        quadratic = {get_element(n).face_type in ("line3", "triangle6", "quad9")
                     for n in self.cells}
        if len(quadratic) > 1:
            raise ValueError("cannot mix linear and quadratic cell blocks")
        by_type = {}
        for name, conn in self.cells.items():
            el = get_element(name)
            for f in el.faces:
                by_type.setdefault(el.face_type, []).append(conn[:, list(f)])
        out = {}
        for ftype, parts in by_type.items():
            faces = np.vstack(parts)
            key = np.sort(faces[:, :FACE_CORNERS[ftype]], axis=1)
            _, idx, counts = np.unique(key, axis=0, return_index=True, return_counts=True)
            out[ftype] = faces[np.sort(idx[counts == 1])]
        return out if len(out) > 1 else next(iter(out.values()))

    def triangulation(self) -> np.ndarray:
        """Corner-node triangles covering every element, for plotting (2-D)."""
        if self.dim != 2:
            raise ValueError("triangulation is only available for 2-D meshes")
        tris = []
        for name, conn in self.cells.items():
            el = get_element(name)
            if el.family == "triangle":
                tris.append(conn[:, :3])
            else:
                tris.append(conn[:, [0, 1, 2]])
                tris.append(conn[:, [0, 2, 3]])
        return np.vstack(tris)

    def nodal_size(self, reduce: str = "min") -> np.ndarray:
        """Corner-edge length of the elements touching each node.

        ``reduce="min"``: the shortest edge, i.e. the thin dimension of
        stretched boundary-layer cells (the viscous length scale);
        ``reduce="max"``: the longest edge, their streamwise dimension (the
        convective length scale, which sets the local pseudo-time step of an
        implicit continuation: a step limited by the thin dimension needs
        thousands of steps to convect anything along a boundary layer).
        """
        if reduce not in ("min", "max"):
            raise ValueError("reduce must be 'min' or 'max'")
        pick = np.minimum if reduce == "min" else np.maximum
        h = np.full(self.n_nodes, np.inf if reduce == "min" else 0.0)
        for name, conn in self.cells.items():
            el = get_element(name)
            lengths = np.full(conn.shape[0], np.inf if reduce == "min" else 0.0)
            for i, j in el.edges:
                a, b = conn[:, i], conn[:, j]
                lengths = pick(lengths, np.linalg.norm(self.points[a] - self.points[b], axis=1))
            for j in range(conn.shape[1]):
                pick.at(h, conn[:, j], lengths)
        return h

    def wall_distance(self, tags) -> np.ndarray:
        """Nodal distance to the boundary edges of ``tags`` (a name or a list)."""
        from aa540fem.core.wall_distance import wall_distance

        return wall_distance(self, tags)

    def centroids(self) -> np.ndarray:
        """Physical centroid of every element, concatenated in block order."""
        out = []
        for name, conn in self.cells.items():
            el = get_element(name)
            phi = el.shape(*[[c] for c in el.centroid])[0]
            out.append(np.column_stack([self.points[conn, k] @ phi[0] for k in range(self.dim)]))
        return np.vstack(out)


def geometry(a: float, b: float, elems: int, elem_type=2) -> Mesh:
    """Build a structured mesh of ``elems x elems`` cells on ``[0, a] x [0, b]``.

    Port of ``geometry.m``.

    Parameters
    ----------
    a, b      : domain size in the x and y directions.
    elems     : number of cells along each direction.  For triangles every
                cell is split into two elements.
    elem_type : 1 / ``"triangle"``   --> 3-node linear triangles
                2 / ``"quad"``       --> 4-node bilinear quadrilaterals
                3 / ``"quad9"``      --> 9-node biquadratic quadrilaterals
                ``"triangle6"``      --> 6-node quadratic triangles

    The boundary tags are ``"top"``, ``"right"``, ``"left"`` and ``"bottom"``.
    """
    elems = int(elems)
    if elems < 1:
        raise ValueError("elems must be >= 1")
    el = get_element(elem_type)

    p = 2 if el.face_type == "line3" else 1     # nodes per cell edge minus one
    nx = ny = p * elems + 1

    xx = np.linspace(0.0, a, nx)
    yy = np.linspace(0.0, b, ny)
    X, Y = np.meshgrid(xx, yy)          # shape (ny, nx), X[i, j] = xx[j]
    points = np.column_stack([X.ravel(), Y.ravel()])

    # Global index of the lower-left node of each cell
    ex, ey = np.meshgrid(np.arange(elems), np.arange(elems))
    n0 = (p * ey * nx + p * ex).ravel()

    if el.name == "triangle":
        lower = np.column_stack([n0, n0 + 1, n0 + nx])
        upper = np.column_stack([n0 + nx + 1, n0 + nx, n0 + 1])
        conn = np.empty((2 * n0.size, 3), dtype=int)
        conn[0::2] = lower
        conn[1::2] = upper
    elif el.name == "triangle6":
        # corner offsets: c0 = origin, c1 = +2 (right), c2 = +2nx (up), c3 = far corner
        c1, c2, c3 = 2, 2 * nx, 2 * nx + 2
        lower = np.column_stack([n0, n0 + c1, n0 + c2,
                                 n0 + 1, n0 + nx + 1, n0 + nx])
        upper = np.column_stack([n0 + c3, n0 + c2, n0 + c1,
                                 n0 + 2 * nx + 1, n0 + nx + 1, n0 + nx + 2])
        conn = np.empty((2 * n0.size, 6), dtype=int)
        conn[0::2] = lower
        conn[1::2] = upper
    elif el.name == "quad":
        conn = np.column_stack([n0, n0 + 1, n0 + nx + 1, n0 + nx])
    else:  # quad9, Gmsh ordering
        offsets = np.array([0, 2, 2 * nx + 2, 2 * nx,          # corners
                            1, nx + 2, 2 * nx + 1, nx,          # bottom right top left
                            nx + 1])                            # centre
        conn = n0[:, None] + offsets[None, :]

    ix = np.arange(nx)
    iy = np.arange(ny)
    # Nodes along each side, traversed with the domain on the left
    side_nodes = {
        "bottom": ix,
        "right": iy * nx + (nx - 1),
        "top": ((ny - 1) * nx + ix)[::-1],
        "left": (iy * nx)[::-1],
    }

    def edges(nodes):
        starts = np.arange(0, nodes.size - 1, p)
        if p == 1:
            return np.column_stack([nodes[starts], nodes[starts + 1]])
        return np.column_stack([nodes[starts], nodes[starts + 2], nodes[starts + 1]])

    boundary = {side: edges(side_nodes[side]) for side in SIDES}
    return Mesh(points, {el.name: conn}, boundary, X, Y, nx, ny)


def box(a: float, b: float, c: float, elems, elem_type="hexahedron27") -> Mesh:
    """Structured mesh of a box ``[0, a] x [0, b] x [0, c]`` in 27-node
    hexahedra (``"hexahedron27"``) or 10-node tetrahedra (``"tetra10"``, six
    per cell), ``elems`` cells per direction (an int or a triple).

    Boundary tags: ``left``/``right`` (x = 0 / a), ``bottom``/``top``
    (y = 0 / b), ``front``/``back`` (z = 0 / c).
    """
    el = get_element(elem_type)
    if el.name not in ("hexahedron27", "tetra10"):
        raise ValueError("box supports 'hexahedron27' and 'tetra10'")
    ex, ey, ez = (elems, elems, elems) if np.isscalar(elems) else tuple(int(e) for e in elems)
    nx, ny, nz = 2 * ex + 1, 2 * ey + 1, 2 * ez + 1
    xx, yy, zz = np.linspace(0, a, nx), np.linspace(0, b, ny), np.linspace(0, c, nz)
    Z, Y, X = np.meshgrid(zz, yy, xx, indexing="ij")            # node k = i + nx (j + ny k)
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    ci, cj, ck = np.meshgrid(np.arange(ex), np.arange(ey), np.arange(ez), indexing="ij")
    n0 = (2 * ci + nx * (2 * cj + ny * 2 * ck)).ravel()          # origin node of each cell

    def offset(i, j, k):                                        # grid steps within a cell
        return i + nx * (j + ny * k)

    if el.name == "hexahedron27":
        from aa540fem.core.shape_functions import HEX27_NODES

        offsets = np.array([offset(*(nd + 1).astype(int)) for nd in HEX27_NODES])
        conn = n0[:, None] + offsets[None, :]
        cells = {"hexahedron27": conn}
    else:
        from aa540fem.core.shape_functions import TET10_EDGES

        corner = {(i, j, k): offset(2 * i, 2 * j, 2 * k)
                  for i in (0, 1) for j in (0, 1) for k in (0, 1)}
        # Kuhn's six tetrahedra along the diagonal (0,0,0) -> (1,1,1)
        tets = []
        for order in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
            v = [np.zeros(3, dtype=int)]
            for axis in order:
                w = v[-1].copy()
                w[axis] = 1
                v.append(w)
            tets.append([corner[tuple(p)] for p in v])
        rows = []
        for t in tets:
            mids = [(t[p] + t[q]) // 2 for p, q in TET10_EDGES]      # grid midpoints
            rows.append(np.array(t + mids))
        conn = np.vstack([n0[:, None] + r[None, :] for r in rows])
        cells = {"tetra10": conn}
    mesh = Mesh(points, cells, {})
    mesh.check_orientation()
    faces = mesh.boundary_faces()
    nc = FACE_CORNERS[el.face_type]
    planes = {"left": (0, 0.0), "right": (0, a), "bottom": (1, 0.0), "top": (1, b),
              "front": (2, 0.0), "back": (2, c)}
    boundary = {}
    for tag, (axis, value) in planes.items():
        on = np.all(np.isclose(points[faces[:, :nc], axis], value), axis=1)
        boundary[tag] = faces[on]
    mesh.boundary = boundary
    return mesh


__all__ = ["Mesh", "geometry", "box", "SIDES", "BOX_SIDES", "ELEMENTS", "FACE_ELEMENTS",
           "FACE_TYPES"]
