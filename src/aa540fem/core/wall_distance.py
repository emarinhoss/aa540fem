"""Distance from every node to the nearest wall, from the boundary faces (2-D and 3-D)."""

from __future__ import annotations

import numpy as np


def wall_distance(mesh, tags, chunk: int = 2000) -> np.ndarray:
    """Nodal distance to the union of the boundary faces of ``tags``.

    2-D: the faces are edges, quadratic edges are split into their two
    straight halves and the distance is the minimum point-segment distance.
    3-D: triangles are split into four flat sub-triangles through their
    mid-edge nodes (quadratic faces) and quadrilaterals into eight through
    the mid-edge and centre nodes, and the distance is the minimum
    point-triangle distance.  Exact for flat walls; for curved walls the
    error is of the order of the sagitta of a sub-face.
    """
    if mesh.dim == 3:
        return _wall_distance_3d(mesh, tags, chunk)
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


def _face_triangles(faces):
    """Flat sub-triangles (node index triples) of a block of boundary faces."""
    k = faces.shape[1]
    if k == 3:
        return faces
    if k == 6:                                       # triangle6: corners 0 1 2, mids 3 4 5
        return np.vstack([faces[:, [0, 3, 5]], faces[:, [3, 1, 4]], faces[:, [5, 4, 2]],
                          faces[:, [3, 4, 5]]])
    if k == 4:
        return np.vstack([faces[:, [0, 1, 2]], faces[:, [0, 2, 3]]])
    if k == 9:                                       # quad9: corners, mids 4-7, centre 8
        tris = []
        for i in range(4):
            c0, c1, m = i, (i + 1) % 4, 4 + i
            tris.append(faces[:, [c0, m, 8]])
            tris.append(faces[:, [m, c1, 8]])
        return np.vstack(tris)
    raise ValueError(f"unsupported boundary face with {k} nodes")


def point_triangle_distance(P, A, B, C):
    """Distance from points ``P (c, 3)`` to every triangle ``(A, B, C)`` ``(m, 3)``,
    returned as ``(c, m)`` (Ericson, Real-Time Collision Detection, 5.1.5)."""
    ab = B - A
    ac = C - A
    ap = P[:, None, :] - A[None, :, :]                      # (c, m, 3)
    d1 = np.einsum("cmj,mj->cm", ap, ab)
    d2 = np.einsum("cmj,mj->cm", ap, ac)
    bp = P[:, None, :] - B[None, :, :]
    d3 = np.einsum("cmj,mj->cm", bp, ab)
    d4 = np.einsum("cmj,mj->cm", bp, ac)
    cp = P[:, None, :] - C[None, :, :]
    d5 = np.einsum("cmj,mj->cm", cp, ab)
    d6 = np.einsum("cmj,mj->cm", cp, ac)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    denom = np.maximum(va + vb + vc, 1e-300)
    # barycentric coordinates of the projection, then clamp to the regions
    v = vb / denom
    w = vc / denom
    closest = A[None] + v[..., None] * ab[None] + w[..., None] * ac[None]
    # vertex regions
    with np.errstate(divide="ignore", invalid="ignore"):
        t_ab = np.clip(d1 / np.maximum(d1 - d3, 1e-300), 0.0, 1.0)
        t_ac = np.clip(d2 / np.maximum(d2 - d6, 1e-300), 0.0, 1.0)
        t_bc = np.clip((d4 - d3) / np.maximum((d4 - d3) + (d5 - d6), 1e-300), 0.0, 1.0)
    on_a = (d1 <= 0) & (d2 <= 0)
    on_b = (d3 >= 0) & (d4 <= d3)
    on_c = (d6 >= 0) & (d5 <= d6)
    on_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    on_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    on_bc = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
    closest = np.where(on_bc[..., None], B[None] + t_bc[..., None] * (C - B)[None], closest)
    closest = np.where(on_ac[..., None], A[None] + t_ac[..., None] * ac[None], closest)
    closest = np.where(on_ab[..., None], A[None] + t_ab[..., None] * ab[None], closest)
    closest = np.where(on_c[..., None], C[None], closest)
    closest = np.where(on_b[..., None], B[None], closest)
    closest = np.where(on_a[..., None], A[None], closest)
    return np.linalg.norm(P[:, None, :] - closest, axis=2)


def _wall_distance_3d(mesh, tags, chunk):
    tris = np.vstack([_face_triangles(block)
                      for tag in ([tags] if isinstance(tags, str) else tags)
                      for block in mesh.face_blocks(tag)])
    A, B, C = (mesh.points[tris[:, k]] for k in range(3))
    d = np.empty(mesh.n_nodes)
    # bound the (chunk x m) work arrays to a few hundred MB
    chunk = max(1, min(chunk, int(2e7 // max(tris.shape[0], 1))))
    for s in range(0, mesh.n_nodes, chunk):
        d[s:s + chunk] = point_triangle_distance(mesh.points[s:s + chunk], A, B, C).min(axis=1)
    return d
