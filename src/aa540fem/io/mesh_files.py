"""Reading meshes from files and writing results for ParaView (via meshio).

``meshio`` is an optional dependency (``pip install meshio``).  Any format
it reads works, but Gmsh ``.msh`` files are the main target: physical
groups of dimension ``d - 1`` (lines of a 2-D mesh, surfaces of a 3-D one)
become the boundary tags of the :class:`Mesh`; a file with 3-D cells is a
3-D mesh, otherwise a 2-D one.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np

from aa540fem.core.elements import ELEMENTS, FACE_TYPE_BY_NODES, FACE_TYPES
from aa540fem.core.mesh import Mesh

# Gmsh numbers the nodes of its 18-node prism differently from VTK (meshio passes
# them through unchanged): VTK node i is Gmsh node GMSH_WEDGE18[i]
GMSH_WEDGE18 = [0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11, 15, 17, 16]


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
    boundary = {t: ({k: new[v] for k, v in e.items()} if isinstance(e, dict) else new[e])
                for t, e in boundary.items()}
    for t, e in boundary.items():
        if any((b < 0).any() for b in (e.values() if isinstance(e, dict) else [e])):
            raise ValueError(f"boundary {t!r} references nodes not used by any element")
    return points[used], cells, boundary


def _face_group(parts):
    """One array per tag, or a dict face type -> array when the tag mixes
    triangles and quadrilaterals (3-D meshes with prisms)."""
    widths = sorted({p.shape[1] for p in parts})
    if len(widths) == 1:
        return np.vstack(parts)
    return {FACE_TYPE_BY_NODES[k]: np.vstack([p for p in parts if p.shape[1] == k])
            for k in widths}


def from_meshio(m, orient: bool = True, gmsh_order: bool = False) -> Mesh:
    """Convert a ``meshio.Mesh`` to :class:`Mesh` (``gmsh_order``: the file came
    from Gmsh, whose 18-node prisms meshio does not renumber to VTK order)."""
    dim = 3 if any(cb.type in ELEMENTS and ELEMENTS[cb.type].dim == 3 for cb in m.cells) else 2
    points = np.asarray(m.points, dtype=float)[:, :dim]

    blocks = defaultdict(list)
    for cb in m.cells:
        if cb.type in ELEMENTS and ELEMENTS[cb.type].dim == dim:
            data = np.asarray(cb.data, dtype=int)
            if gmsh_order and cb.type == "wedge18":
                data = data[:, GMSH_WEDGE18]
            blocks[cb.type].append(data)
    if not blocks:
        raise ValueError(f"no supported cells found; supported: {sorted(ELEMENTS)}")
    cells = {name: np.vstack(parts) for name, parts in blocks.items()}

    # Boundary faces from the (d-1)-dimensional cells, grouped by Gmsh physical group.
    names_by_tag = {int(v[0]): name for name, v in m.field_data.items() if int(v[1]) == dim - 1}
    phys = m.cell_data.get("gmsh:physical")
    boundary = defaultdict(list)
    for i, cb in enumerate(m.cells):
        if cb.type not in FACE_TYPES[dim]:
            continue
        data = np.asarray(cb.data, dtype=int)
        if phys is not None:
            tags = np.asarray(phys[i]).astype(int)
            for t in np.unique(tags):
                boundary[names_by_tag.get(t, str(t))].append(data[tags == t])
        else:
            boundary["boundary"].append(data)
    boundary = {t: _face_group(parts) for t, parts in boundary.items()}

    points, cells, boundary = _compact(points, cells, boundary)
    mesh = Mesh(points, cells, boundary)
    if orient:
        mesh.check_orientation()
    if not boundary:
        faces = mesh.boundary_faces()
        mesh.boundary = ({"boundary": faces} if not isinstance(faces, dict)
                         else {f"boundary_{k}": v for k, v in faces.items()})
    return mesh


def read_mesh(path, orient: bool = True) -> Mesh:
    """Read a mesh file with meshio (Gmsh ``.msh``, VTK, XDMF, ...).

    Elements with a negative Jacobian are flipped unless ``orient=False``.
    Nodes not used by any element are removed.  If the file has no boundary
    line cells the outer boundary is detected and tagged ``"boundary"``.
    """
    return from_meshio(_meshio().read(path), orient=orient,
                       gmsh_order=str(path).lower().endswith(".msh"))


def to_meshio(mesh: Mesh, point_data=None, cell_data=None):
    """Convert :class:`Mesh` (plus optional fields) to a ``meshio.Mesh``.

    ``cell_data`` values are arrays over all elements in block order (as
    returned by the post-processing helpers) and are split per block.
    """
    meshio = _meshio()
    points = (np.column_stack([mesh.points, np.zeros(mesh.n_nodes)]) if mesh.dim == 2
              else mesh.points)
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
