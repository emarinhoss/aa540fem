"""Reading meshes from files and writing results for ParaView (via meshio).

``meshio`` is an optional dependency (``pip install meshio``).  Any format
it reads works, but Gmsh ``.msh`` files are the main target: physical
groups of dimension 1 become the boundary tags of the :class:`Mesh`.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from aa540fem.core.elements import ELEMENTS, FACE_ELEMENTS
from aa540fem.core.mesh import Mesh


def _meshio():
    try:
        import meshio
    except ImportError:  # pragma: no cover
        raise ImportError(
            "reading/writing mesh files requires meshio: pip install meshio") from None
    return meshio


def _compact(points, cells, boundary):
    """Drop nodes that no element uses and renumber the rest contiguously."""
    used = np.unique(np.concatenate([c.ravel() for c in cells.values()]))
    if used.size == points.shape[0] and np.array_equal(used, np.arange(points.shape[0])):
        return points, cells, boundary
    new = -np.ones(points.shape[0], dtype=int)
    new[used] = np.arange(used.size)
    cells = {n: new[c] for n, c in cells.items()}
    boundary = {t: new[e] for t, e in boundary.items()}
    for t, e in boundary.items():
        if (e < 0).any():
            raise ValueError(f"boundary {t!r} references nodes not used by any element")
    return points[used], cells, boundary


def from_meshio(m, orient: bool = True) -> Mesh:
    """Convert a ``meshio.Mesh`` to :class:`Mesh`."""
    points = np.asarray(m.points, dtype=float)[:, :2]

    blocks = defaultdict(list)
    for cb in m.cells:
        if cb.type in ELEMENTS:
            blocks[cb.type].append(np.asarray(cb.data, dtype=int))
    if not blocks:
        raise ValueError(f"no supported 2-D cells found; supported: {sorted(ELEMENTS)}")
    cells = {name: np.vstack(parts) for name, parts in blocks.items()}

    # Boundary edges from the 1-D cells, grouped by Gmsh physical group.
    names_by_tag = {int(v[0]): name for name, v in m.field_data.items() if int(v[1]) == 1}
    phys = m.cell_data.get("gmsh:physical")
    boundary = defaultdict(list)
    for i, cb in enumerate(m.cells):
        if cb.type not in FACE_ELEMENTS:
            continue
        data = np.asarray(cb.data, dtype=int)
        if phys is not None:
            tags = np.asarray(phys[i]).astype(int)
            for t in np.unique(tags):
                boundary[names_by_tag.get(t, str(t))].append(data[tags == t])
        else:
            boundary["boundary"].append(data)
    boundary = {t: np.vstack(parts) for t, parts in boundary.items()}

    points, cells, boundary = _compact(points, cells, boundary)
    mesh = Mesh(points, cells, boundary)
    if orient:
        mesh.check_orientation()
    if not boundary:
        mesh.boundary = {"boundary": mesh.boundary_faces()}
    return mesh


def read_mesh(path, orient: bool = True) -> Mesh:
    """Read a mesh file with meshio (Gmsh ``.msh``, VTK, XDMF, ...).

    Elements with a negative Jacobian are flipped unless ``orient=False``.
    Nodes not used by any element are removed.  If the file has no boundary
    line cells the outer boundary is detected and tagged ``"boundary"``.
    """
    return from_meshio(_meshio().read(path), orient=orient)


def to_meshio(mesh: Mesh, point_data=None, cell_data=None):
    """Convert :class:`Mesh` (plus optional fields) to a ``meshio.Mesh``.

    ``cell_data`` values are arrays over all elements in block order (as
    returned by the post-processing helpers) and are split per block.
    """
    meshio = _meshio()
    points = np.column_stack([mesh.points, np.zeros(mesh.n_nodes)])
    cells = [(name, conn) for name, conn in mesh.cells.items()]
    sizes = [conn.shape[0] for _, conn in cells]

    cd = None
    if cell_data:
        cd = {}
        for key, arr in cell_data.items():
            arr = np.asarray(arr)
            if arr.shape[0] != sum(sizes):
                raise ValueError(f"cell_data[{key!r}] must have one row per element")
            cd[key] = list(np.split(arr, np.cumsum(sizes)[:-1]))
    pd = {k: np.asarray(v) for k, v in (point_data or {}).items()}
    return meshio.Mesh(points, cells, point_data=pd, cell_data=cd)


def write_vtk(path, mesh: Mesh, point_data=None, cell_data=None):
    """Write the mesh and fields to a ParaView-readable file (``.vtu`` recommended)."""
    m = to_meshio(mesh, point_data, cell_data)
    m.write(str(path))
    return path
