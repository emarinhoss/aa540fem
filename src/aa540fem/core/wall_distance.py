"""Distance from every node to the nearest wall, from the boundary edges."""

from __future__ import annotations

import numpy as np


def wall_distance(mesh, tags, chunk: int = 2000) -> np.ndarray:
    """Nodal distance to the union of the boundary edges of ``tags``.

    Quadratic edges are split into their two straight halves.  The distance
    of a node is the minimum over the segments of the point-segment distance
    (exact for straight walls; for curved walls the error is of the order of
    the sagitta of a half edge).
    """
    segments = []
    for tag in ([tags] if isinstance(tags, str) else tags):
        edges = mesh.boundary[tag]
        if edges.shape[1] == 3:
            segments.append(edges[:, [0, 2]])
            segments.append(edges[:, [2, 1]])
        else:
            segments.append(edges[:, :2])
    seg = np.vstack(segments)
    a = mesh.points[seg[:, 0]]                       # (m, 2)
    b = mesh.points[seg[:, 1]]
    ab = b - a
    ab2 = np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-300)
    d = np.empty(mesh.n_nodes)
    for s in range(0, mesh.n_nodes, chunk):
        P = mesh.points[s:s + chunk]                 # (c, 2)
        ap = P[:, None, :] - a[None, :, :]           # (c, m, 2)
        t = np.clip(np.einsum("cmj,mj->cm", ap, ab) / ab2[None, :], 0.0, 1.0)
        closest = a[None, :, :] + t[:, :, None] * ab[None, :, :]
        dist = np.linalg.norm(P[:, None, :] - closest, axis=2)
        d[s:s + chunk] = dist.min(axis=1)
    return d
