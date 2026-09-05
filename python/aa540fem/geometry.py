"""Structured mesh generation on a rectangle.  Port of ``geometry.m``.

Node numbering is row-major with ``x`` varying fastest: node ``k`` sits at
``(xx[k % nx], yy[k // nx])``.  This differs from the MATLAB code, which
numbered nodes column-major (``y`` fastest) and hard-coded the 9-node
connectivity for a 100-element mesh; here the connectivity is derived from
the mesh size for every element type.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .shape_functions import NODES_PER_ELEMENT

SIDES = ("top", "right", "left", "bottom")


@dataclass
class Mesh:
    """Nodes, connectivity and boundary information of a structured mesh.

    Attributes
    ----------
    elem_type : 1 (3-node triangle), 2 (4-node quad) or 3 (9-node quad).
    x, y      : flat nodal coordinate arrays of length ``n_nodes``.
    X, Y      : the same coordinates as ``(ny, nx)`` grids (for plotting).
    conn      : ``(n_elems, nodes_per_element)`` integer connectivity,
                0-based, local ordering matching the shape functions.
    bc_nodes  : dict side -> 1-D array of node indices on that side.
    bc_edges  : dict side -> ``(n_edges, nodes_per_edge)`` array of the
                boundary edges on that side, ordered along the edge.
    nx, ny    : number of nodes along x and y.
    """

    elem_type: int
    x: np.ndarray
    y: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    conn: np.ndarray
    bc_nodes: dict
    bc_edges: dict
    nx: int
    ny: int

    @property
    def n_nodes(self) -> int:
        return self.x.size

    @property
    def n_elems(self) -> int:
        return self.conn.shape[0]

    @property
    def nodes_per_element(self) -> int:
        return NODES_PER_ELEMENT[self.elem_type]

    def to_grid(self, values: np.ndarray) -> np.ndarray:
        """Reshape a nodal vector to the ``(ny, nx)`` grid of ``X``/``Y``."""
        return np.asarray(values).reshape(self.ny, self.nx)


def geometry(a: float, b: float, elems: int, elem_type: int) -> Mesh:
    """Build a structured mesh of ``elems x elems`` cells on ``[0, a] x [0, b]``.

    Parameters
    ----------
    a, b      : domain size in the x and y directions.
    elems     : number of cells along each direction.  For ``elem_type == 1``
                every cell is split into two triangles (so there are
                ``2 * elems**2`` elements), otherwise there is one element per
                cell.
    elem_type : 1 --> 3-node linear triangles
                2 --> 4-node bilinear quadrilaterals
                3 --> 9-node biquadratic quadrilaterals
    """
    elems = int(elems)
    if elems < 1:
        raise ValueError("elems must be >= 1")
    if elem_type not in (1, 2, 3):
        raise ValueError(f"Unknown element type {elem_type!r}; expected 1, 2 or 3")

    # Nodes per cell edge: linear elements have 2, quadratic have 3.
    p = 2 if elem_type == 3 else 1
    nx = ny = p * elems + 1

    xx = np.linspace(0.0, a, nx)
    yy = np.linspace(0.0, b, ny)
    X, Y = np.meshgrid(xx, yy)          # shape (ny, nx), X[i, j] = xx[j]
    x = X.ravel()
    y = Y.ravel()

    # Global index of the lower-left node of each cell, cell (ex, ey)
    ex, ey = np.meshgrid(np.arange(elems), np.arange(elems))
    n0 = (p * ey * nx + p * ex).ravel()     # shape (elems**2,)

    if elem_type == 1:
        # Two triangles per cell, same split as the MATLAB code:
        #   lower-left triangle : origin, +eta (up), +xi (right)
        #   upper-right triangle: far corner, +eta (down), +xi (left)
        lower = np.column_stack([n0, n0 + nx, n0 + 1])
        upper = np.column_stack([n0 + nx + 1, n0 + 1, n0 + nx])
        conn = np.empty((2 * n0.size, 3), dtype=int)
        conn[0::2] = lower
        conn[1::2] = upper
    elif elem_type == 2:
        conn = np.column_stack([n0, n0 + 1, n0 + nx + 1, n0 + nx])
    else:
        offsets = np.array([0, 1, 2, nx, nx + 1, nx + 2, 2 * nx, 2 * nx + 1, 2 * nx + 2])
        conn = n0[:, None] + offsets[None, :]

    ix = np.arange(nx)
    iy = np.arange(ny)
    bc_nodes = {
        "top": (ny - 1) * nx + ix,
        "right": iy * nx + (nx - 1),
        "left": iy * nx,
        "bottom": ix,
    }

    # Boundary edges as consecutive groups of p+1 nodes along each side.
    def edges(nodes):
        starts = np.arange(0, nodes.size - 1, p)
        return np.stack([nodes[s:s + p + 1] for s in starts])

    bc_edges = {side: edges(nodes) for side, nodes in bc_nodes.items()}

    return Mesh(elem_type, x, y, X, Y, conn, bc_nodes, bc_edges, nx, ny)
