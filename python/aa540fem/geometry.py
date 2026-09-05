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

from .elements import ELEMENTS, FACE_ELEMENTS, get_element

SIDES = ("top", "right", "left", "bottom")


@dataclass
class Mesh:
    """Nodes, connectivity and boundary information of a 2-D mesh.

    Attributes
    ----------
    points   : ``(n_nodes, 2)`` nodal coordinates.
    cells    : dict cell type -> ``(n_elems, nodes_per_element)`` 0-based
               connectivity in Gmsh/VTK local ordering (see
               :mod:`aa540fem.elements`).
    boundary : dict tag -> ``(n_edges, 2 or 3)`` boundary edges ordered
               ``(start, end[, mid])``.
    X, Y, nx, ny : grid form of the coordinates, only for structured meshes.
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
        if self.points.ndim != 2 or self.points.shape[1] != 2:
            raise ValueError("points must be an (n, 2) array")
        self.cells = {name: np.asarray(c, dtype=int) for name, c in self.cells.items()}
        for name, conn in self.cells.items():
            el = get_element(name)
            if conn.ndim != 2 or conn.shape[1] != el.n_nodes:
                raise ValueError(f"{name} connectivity must be (n, {el.n_nodes})")
        self.boundary = {tag: np.asarray(e, dtype=int).reshape(-1, np.asarray(e).shape[-1])
                         if np.asarray(e).size else np.zeros((0, 2), dtype=int)
                         for tag, e in self.boundary.items()}

    # ------------------------------------------------------------ basic info
    @property
    def x(self) -> np.ndarray:
        return self.points[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.points[:, 1]

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
        """dict tag -> sorted unique node indices on that boundary."""
        return {tag: np.unique(edges) for tag, edges in self.boundary.items()}

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
            _, dxi, deta = el.shape([el.centroid[0]], [el.centroid[1]])
            xe = self.x[conn]
            ye = self.y[conn]
            out[name] = ((xe @ dxi.T) * (ye @ deta.T) - (xe @ deta.T) * (ye @ dxi.T)).ravel()
        return out

    def check_orientation(self) -> int:
        """Reverse elements with a negative Jacobian so all are counter-clockwise.

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
        """Edges that belong to exactly one element, i.e. the outer boundary.

        Returned as an ``(n_edges, 2 or 3)`` array in the same ordering as
        ``boundary`` entries (traversed with the domain on the left).
        """
        face_types = {get_element(n).face_type for n in self.cells}
        if len(face_types) > 1:
            raise ValueError("cannot mix linear and quadratic cell blocks")
        faces = []
        for name, conn in self.cells.items():
            el = get_element(name)
            for f in el.faces:
                faces.append(conn[:, list(f)])
        faces = np.vstack(faces)
        key = np.sort(faces[:, :2], axis=1)
        _, idx, counts = np.unique(key, axis=0, return_index=True, return_counts=True)
        return faces[idx[counts == 1]]

    def triangulation(self) -> np.ndarray:
        """Corner-node triangles covering every element, for plotting."""
        tris = []
        for name, conn in self.cells.items():
            el = get_element(name)
            if el.family == "triangle":
                tris.append(conn[:, :3])
            else:
                tris.append(conn[:, [0, 1, 2]])
                tris.append(conn[:, [0, 2, 3]])
        return np.vstack(tris)

    def centroids(self) -> np.ndarray:
        """Physical centroid of every element, concatenated in block order."""
        out = []
        for name, conn in self.cells.items():
            el = get_element(name)
            phi, _, _ = el.shape([el.centroid[0]], [el.centroid[1]])
            out.append(np.column_stack([self.x[conn] @ phi[0], self.y[conn] @ phi[0]]))
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


__all__ = ["Mesh", "geometry", "SIDES", "ELEMENTS", "FACE_ELEMENTS"]
