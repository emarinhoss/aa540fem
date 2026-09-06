"""Mesh partitioning for domain decomposition.

:func:`partition_mesh` splits the elements of a :class:`Mesh` into ``n``
parts with METIS (``pymetis``; a coordinate bisection fallback keeps the
module usable without it) and builds, for every part, the owned nodes, the
one-layer halo of ghost nodes and the local numbering that a distributed
assembly needs: a rank assembles its own elements with the element kernels
of the serial code, the entries whose row is owned stay local and the
others are sent to the owner (PETSc does that in ``setValuesCOO``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from aa540fem.core.elements import get_element


@dataclass
class Partition:
    """One part of a partitioned mesh."""

    part: int
    elements: dict            # cell type -> indices of the elements of this part
    owned: np.ndarray         # global node numbers owned by this part (sorted)
    ghost: np.ndarray         # global node numbers used but owned elsewhere (sorted)
    owner: np.ndarray         # owner part of every global node
    local: np.ndarray         # global node -> local index (owned first, then ghosts; -1 elsewhere)

    @property
    def n_local(self) -> int:
        return self.owned.size + self.ghost.size

    def local_cells(self, mesh):
        """Connectivity of the part's elements in local numbering."""
        return {name: self.local[mesh.cells[name][idx]] for name, idx in self.elements.items()
                if idx.size}


def _element_graph(mesh):
    """Adjacency of the element graph: elements sharing a corner node."""
    conns = []
    for name, conn in mesh.cells.items():
        nc = get_element(name).n_corners
        conns.append(conn[:, :nc])
    n_elems = sum(c.shape[0] for c in conns)
    node_elems = [[] for _ in range(mesh.n_nodes)]
    offset = 0
    for c in conns:
        for e, nodes in enumerate(c):
            for node in nodes:
                node_elems[node].append(offset + e)
        offset += c.shape[0]
    adjacency = [set() for _ in range(n_elems)]
    for elems in node_elems:
        for e in elems:
            adjacency[e].update(elems)
    for e in range(n_elems):
        adjacency[e].discard(e)
    return [sorted(a) for a in adjacency]


def element_parts(mesh, n_parts: int) -> np.ndarray:
    """Part number of every element (blocks concatenated in ``mesh.cells`` order)."""
    n_elems = mesh.n_elems
    if n_parts <= 1:
        return np.zeros(n_elems, dtype=int)
    try:
        import pymetis

        _, parts = pymetis.part_graph(n_parts, adjacency=_element_graph(mesh))
        return np.asarray(parts, dtype=int)
    except ImportError:
        # fallback: sort the element centroids along the longer axis and cut into slabs
        centroids = mesh.centroids()
        axis = int(np.argmax(np.ptp(centroids, axis=0)))
        order = np.argsort(centroids[:, axis], kind="stable")
        parts = np.empty(n_elems, dtype=int)
        parts[order] = np.arange(n_elems) * n_parts // n_elems
        return parts


def partition_mesh(mesh, n_parts: int, parts=None) -> list[Partition]:
    """Partition ``mesh`` into ``n_parts`` :class:`Partition` objects.

    Node ownership goes to the lowest part number among the parts whose
    elements touch the node (a deterministic rule every rank can apply).
    """
    parts = element_parts(mesh, n_parts) if parts is None else np.asarray(parts, dtype=int)
    owner = np.full(mesh.n_nodes, n_parts, dtype=int)
    offset = 0
    blocks = {}
    for name, conn in mesh.cells.items():
        ne = conn.shape[0]
        block_parts = parts[offset:offset + ne]
        blocks[name] = block_parts
        for p in range(n_parts):
            nodes = conn[block_parts == p].ravel()
            np.minimum.at(owner, nodes, p)
        offset += ne
    result = []
    for p in range(n_parts):
        elements = {name: np.nonzero(bp == p)[0] for name, bp in blocks.items()}
        used = np.unique(np.concatenate([mesh.cells[name][idx].ravel()
                                         for name, idx in elements.items()]
                                        or [np.zeros(0, dtype=int)]))
        owned = used[owner[used] == p]
        ghost = used[owner[used] != p]
        local = -np.ones(mesh.n_nodes, dtype=int)
        local[owned] = np.arange(owned.size)
        local[ghost] = owned.size + np.arange(ghost.size)
        result.append(Partition(p, elements, owned, ghost, owner, local))
    return result
