"""Generate the meshes used by the tests and examples with Gmsh.

Requires the ``gmsh`` Python package (``pip install gmsh``).  Writes the
annulus meshes ``annulus_tri.msh``, ``annulus_tri6.msh``, ``annulus_quad.msh``,
``annulus_quad9.msh``, the cylinder-in-channel mesh ``cylinder_tri6.msh`` (and
its boundary-layer variants ``cylinder_bl.msh`` and the twice-finer
``cylinder_bl_fine.msh``), the 3D-1Z cross-section ``cylinder3d1z_bl.msh``
(``cylinder3d1z.geo``: cylinder at x = 0.5 in the 2.5 channel, coarser, that
``examples/cylinder3d.py --extrude`` extrudes into prisms and hexahedra), the flat-plate meshes
``flat_plate_bl.msh`` and the NACA 0012 far-field mesh
``airfoil_naca0012_a5_tri6.msh`` and the small 3-D tetrahedral box
``box_tet10.msh`` next to this file.  The ``_bl`` meshes use
Gmsh's boundary-layer field, which extrudes quadrilaterals from the wall
(``quad9`` after the second-order pass) into an otherwise triangular mesh.
"""

from __future__ import annotations

import pathlib

import gmsh
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent

MESHES = {
    # name: (geometry file, element order, recombine into quads)
    "annulus_tri": ("annulus.geo", 1, False),
    "annulus_tri6": ("annulus.geo", 2, False),
    "annulus_quad": ("annulus.geo", 1, True),
    "annulus_quad9": ("annulus.geo", 2, True),
    "cylinder_tri6": ("cylinder.geo", 2, False),
}


def _boundary_layer(curves, size_wall, ratio, thickness, quads=True):
    """Gmsh boundary-layer field on ``curves`` (call after synchronize)."""
    f = gmsh.model.mesh.field
    bl = f.add("BoundaryLayer")
    f.setNumbers(bl, "CurvesList", curves)
    f.setNumber(bl, "Size", size_wall)
    f.setNumber(bl, "Ratio", ratio)
    f.setNumber(bl, "Thickness", thickness)
    f.setNumber(bl, "Quads", 1 if quads else 0)
    f.setAsBoundaryLayer(bl)
    return bl


def make(name: str, order: int, quads: bool, geo: pathlib.Path = HERE / "annulus.geo",
         boundary_layer=None, size_factor: float = 1.0):
    """Mesh a ``.geo`` file.  ``boundary_layer=(curves, size, ratio, thickness)``
    adds a quadrilateral boundary layer on those curve tags; ``size_factor``
    scales every mesh size of the file (``Mesh.MeshSizeFactor``), e.g. 0.5
    for a uniformly twice-finer mesh."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(geo))
        gmsh.option.setNumber("Mesh.MeshSizeFactor", size_factor)
        if boundary_layer is not None:
            _boundary_layer(*boundary_layer)
        if quads:
            gmsh.option.setNumber("Mesh.Algorithm", 8)          # Frontal-Delaunay for quads
            gmsh.option.setNumber("Mesh.RecombineAll", 1)
            gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
        gmsh.option.setNumber("Mesh.ElementOrder", order)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)  # quad9 / triangle6
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(2)
        out = HERE / f"{name}.msh"
        gmsh.write(str(out))
        return out
    finally:
        gmsh.finalize()


def make_flat_plate(name: str = "flat_plate_bl", order: int = 2, x0: float = -0.5,
                    x1: float = 2.0, plate: float = 1.5, height: float = 0.5,
                    lc: float = 0.05, lc_plate: float = 0.02, size_wall: float = 2e-3,
                    ratio: float = 1.15, thickness: float = 0.04,
                    lc_edge: float | None = None, edge_radius: float = 0.1):
    """Laminar flat-plate domain ``[x0, x1] x [0, height]``.

    The plate is ``y = 0, 0 <= x <= plate`` (tag ``plate``); upstream of it
    the bottom is a symmetry line (tag ``symmetry``); ``inlet``, ``outlet``
    and ``top`` are the other sides.  A quadrilateral boundary layer grows
    from the plate (first cell ``size_wall``, growth ``ratio``, total
    ``thickness``); the rest is triangles.  ``lc_edge`` refines the
    streamwise spacing towards the two ends of the plate (from ``lc_edge``
    at the ends to ``lc_plate`` at ``edge_radius``): the leading and
    trailing edges are singular points of the flow and cells a thousand
    times longer than high at those points defeat the Newton solver.
    """
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        geo = gmsh.model.geo
        p1 = geo.addPoint(x0, 0.0, 0.0, lc)
        p2 = geo.addPoint(0.0, 0.0, 0.0, lc_plate)
        p3 = geo.addPoint(plate, 0.0, 0.0, lc_plate)
        p4 = geo.addPoint(x1, 0.0, 0.0, lc)
        p5 = geo.addPoint(x1, height, 0.0, lc)
        p6 = geo.addPoint(x0, height, 0.0, lc)
        symmetry = geo.addLine(p1, p2)
        plate_line = geo.addLine(p2, p3)
        wake = geo.addLine(p3, p4)
        outlet = geo.addLine(p4, p5)
        top = geo.addLine(p5, p6)
        inlet = geo.addLine(p6, p1)
        loop = geo.addCurveLoop([symmetry, plate_line, wake, outlet, top, inlet])
        surface = geo.addPlaneSurface([loop])
        geo.synchronize()
        gmsh.model.addPhysicalGroup(1, [inlet], name="inlet")
        gmsh.model.addPhysicalGroup(1, [outlet], name="outlet")
        gmsh.model.addPhysicalGroup(1, [top], name="top")
        gmsh.model.addPhysicalGroup(1, [symmetry, wake], name="symmetry")
        gmsh.model.addPhysicalGroup(1, [plate_line], name="plate")
        gmsh.model.addPhysicalGroup(2, [surface], name="fluid")
        _boundary_layer([plate_line], size_wall, ratio, thickness)
        if lc_edge is not None:
            # size grows linearly with the distance to the plate ends, reaching
            # lc_plate at edge_radius and lc (the far-field size) further out;
            # Gmsh takes the minimum of this field and the point sizes, so the
            # plate keeps lc_plate and the far field lc
            f = gmsh.model.mesh.field
            dist = f.add("Distance")
            f.setNumbers(dist, "PointsList", [p2, p3])
            thr = f.add("Threshold")
            f.setNumber(thr, "InField", dist)
            f.setNumber(thr, "SizeMin", lc_edge)
            f.setNumber(thr, "SizeMax", lc)
            f.setNumber(thr, "DistMin", 0.0)
            f.setNumber(thr, "DistMax", edge_radius * (lc - lc_edge) / (lc_plate - lc_edge))
            f.setAsBackgroundMesh(thr)
        gmsh.option.setNumber("Mesh.ElementOrder", order)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(2)
        out = HERE / f"{name}.msh"
        gmsh.write(str(out))
        return out
    finally:
        gmsh.finalize()


def make_turbulent_flat_plate(name: str = "flat_plate_turb"):
    """Wall-resolved plate of the Spalart-Allmaras validation (``turbulent_flat_plate.py``).

    Domain ``[-0.5, 2.5] x [0, 1]``, plate ``0 <= x <= 2``; 30 quadrilateral
    layers from the wall, the first one 2e-5 thick (``y+`` about 1 at
    ``Re = 1e6`` per unit length), growth ratio 1.25, 6 % total thickness;
    the streamwise spacing shrinks to 2e-4 at the leading and trailing
    edges.
    """
    return make_flat_plate(name=name, x0=-0.5, x1=2.5, plate=2.0, height=1.0, lc=0.1,
                           lc_plate=0.02, size_wall=2e-5, ratio=1.25, thickness=0.06,
                           lc_edge=2e-4)


def make_turbulent_flat_plate_medium(name: str = "flat_plate_turb_medium"):
    """The turbulent plate with the same wall-normal layers (2e-5 first cell,
    ratio 1.25) and edge refinement but a coarser streamwise spacing (12.8k
    nodes against 23k): ``examples/turbulent_flat_plate.py --extrude`` extrudes
    it into a 3-D wall-resolved mesh that a direct solver still fits."""
    return make_flat_plate(name=name, x0=-0.5, x1=2.5, plate=2.0, height=1.0, lc=0.15,
                           lc_plate=0.035, size_wall=2e-5, ratio=1.25, thickness=0.06,
                           lc_edge=2e-4)


def naca4(code: str = "0012", n: int = 100, chord: float = 1.0):
    """Closed-trailing-edge NACA 4-digit profile, cosine spaced.

    Returns ``(upper, lower)`` arrays of ``(x, y)`` points ordered from the
    leading edge to the trailing edge.
    """
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    beta = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3
                  - 0.1036 * x ** 4)
    if m == 0:
        yc = np.zeros_like(x)
        dyc = np.zeros_like(x)
    else:
        front = x < p
        yc = np.where(front, m / p ** 2 * (2 * p * x - x ** 2),
                      m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x - x ** 2))
        dyc = np.where(front, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    th = np.arctan(dyc)
    upper = np.column_stack([x - yt * np.sin(th), yc + yt * np.cos(th)]) * chord
    lower = np.column_stack([x + yt * np.sin(th), yc - yt * np.cos(th)]) * chord
    return upper, lower


def make_airfoil(name: str = "airfoil_naca0012_a5_tri6", code: str = "0012",
                 alpha_deg: float = 5.0, n_points: int = 120, order: int = 2,
                 far=(-6.0, 12.0, -6.0, 6.0), lc_surface: float = 0.008, lc_far: float = 0.6,
                 lc_wake: float = 0.08):
    """Far-field mesh around a NACA 4-digit airfoil rotated by ``-alpha_deg``.

    The freestream is then along +x, so the x/y forces on the airfoil are
    drag and lift.  Physical curves: ``inlet`` (left), ``outlet`` (right),
    ``farfield`` (top and bottom), ``airfoil``.
    """
    upper, lower = naca4(code, n_points)
    a = np.radians(-alpha_deg)
    rot = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    pivot = np.array([0.25, 0.0])                       # quarter chord
    upper = (upper - pivot) @ rot.T + pivot
    lower = (lower - pivot) @ rot.T + pivot

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        geo = gmsh.model.geo
        # airfoil: trailing edge -> upper surface -> leading edge -> lower surface
        le = geo.addPoint(*upper[0], 0.0, lc_surface)
        te = geo.addPoint(*upper[-1], 0.0, lc_surface)
        up = [geo.addPoint(x, y, 0.0, lc_surface) for x, y in upper[1:-1]]
        lo = [geo.addPoint(x, y, 0.0, lc_surface) for x, y in lower[1:-1]]
        c_upper = geo.addSpline([te] + up[::-1] + [le])
        c_lower = geo.addSpline([le] + lo + [te])
        airfoil_loop = geo.addCurveLoop([c_upper, c_lower])

        x0, x1, y0, y1 = far
        p1 = geo.addPoint(x0, y0, 0.0, lc_far)
        p2 = geo.addPoint(x1, y0, 0.0, lc_far)
        p3 = geo.addPoint(x1, y1, 0.0, lc_far)
        p4 = geo.addPoint(x0, y1, 0.0, lc_far)
        bottom = geo.addLine(p1, p2)
        right = geo.addLine(p2, p3)
        top = geo.addLine(p3, p4)
        left = geo.addLine(p4, p1)
        outer_loop = geo.addCurveLoop([bottom, right, top, left])
        surface = geo.addPlaneSurface([outer_loop, airfoil_loop])
        geo.synchronize()

        gmsh.model.addPhysicalGroup(1, [left], name="inlet")
        gmsh.model.addPhysicalGroup(1, [right], name="outlet")
        gmsh.model.addPhysicalGroup(1, [top, bottom], name="farfield")
        gmsh.model.addPhysicalGroup(1, [c_upper, c_lower], name="airfoil")
        gmsh.model.addPhysicalGroup(2, [surface], name="fluid")

        # mesh size: fine at the surface, wake box, coarse far field
        f = gmsh.model.mesh.field
        dist = f.add("Distance")
        f.setNumbers(dist, "CurvesList", [c_upper, c_lower])
        f.setNumber(dist, "Sampling", 400)
        thr = f.add("Threshold")
        f.setNumber(thr, "InField", dist)
        f.setNumber(thr, "SizeMin", lc_surface)
        f.setNumber(thr, "SizeMax", lc_far)
        f.setNumber(thr, "DistMin", 0.05)
        f.setNumber(thr, "DistMax", 3.0)
        wake = f.add("Box")
        f.setNumber(wake, "VIn", lc_wake)
        f.setNumber(wake, "VOut", lc_far)
        f.setNumber(wake, "XMin", -0.3)
        f.setNumber(wake, "XMax", 4.0)
        f.setNumber(wake, "YMin", -0.5)
        f.setNumber(wake, "YMax", 0.5)
        f.setNumber(wake, "Thickness", 1.5)
        combined = f.add("Min")
        f.setNumbers(combined, "FieldsList", [thr, wake])
        f.setAsBackgroundMesh(combined)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)              # Frontal-Delaunay
        gmsh.option.setNumber("Mesh.ElementOrder", order)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(2)
        out = HERE / f"{name}.msh"
        gmsh.write(str(out))
        return out
    finally:
        gmsh.finalize()


def make_box_tet10(name: str = "box_tet10", size=(2.0, 1.0, 1.0), lc: float = 0.35):
    """Unstructured 10-node tetrahedral mesh of a box with the physical surfaces
    ``left``/``right`` (x), ``bottom``/``top`` (y), ``front``/``back`` (z); a
    small 3-D mesh for the tests and the 3-D channel example."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        a, b, c = size
        gmsh.model.occ.addBox(0, 0, 0, a, b, c)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
        planes = {"left": (0, 0.0), "right": (0, a), "bottom": (1, 0.0), "top": (1, b),
                  "front": (2, 0.0), "back": (2, c)}
        for tag, (axis, value) in planes.items():
            faces = [s for (_, s) in gmsh.model.getEntities(2)
                     if abs(gmsh.model.occ.getCenterOfMass(2, s)[axis] - value) < 1e-9]
            gmsh.model.addPhysicalGroup(2, faces, name=tag)
        gmsh.model.addPhysicalGroup(3, [1], name="fluid")
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(3)
        out = HERE / f"{name}.msh"
        gmsh.write(str(out))
        return out
    finally:
        gmsh.finalize()


def make_cylinder3d(name: str = "cylinder3d_tet10", lc: float = 0.05, lc_cyl: float = 0.01,
                    thickness: float = 0.06):
    """Schaefer-Turek 3D-1Z geometry: channel ``[0, 2.5] x [0, 0.41] x [0, 0.41]``
    with a cylinder of diameter 0.1 along ``z`` centred at ``(0.5, 0.2)``, in
    10-node tetrahedra graded from ``lc_cyl`` at the cylinder to ``lc`` at
    distance ``thickness``.  Physical surfaces: ``inlet`` (x = 0), ``outlet``
    (x = 2.5), ``walls`` (the four channel walls) and ``cylinder``."""
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add(name)
        L, H = 2.5, 0.41
        channel = gmsh.model.occ.addBox(0, 0, 0, L, H, H)
        cyl = gmsh.model.occ.addCylinder(0.5, 0.2, 0, 0, 0, H, 0.05)
        gmsh.model.occ.cut([(3, channel)], [(3, cyl)])
        gmsh.model.occ.synchronize()
        inlet, outlet, walls, cylinder = [], [], [], []
        for _, s in gmsh.model.getEntities(2):
            x, y, z = gmsh.model.occ.getCenterOfMass(2, s)
            xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, s)
            flat_x = (xmax - xmin) < 1e-6
            if flat_x and abs(x) < 1e-6:
                inlet.append(s)
            elif flat_x and abs(x - L) < 1e-6:
                outlet.append(s)
            elif abs(x - 0.5) < 0.06 and abs(y - 0.2) < 0.06 and (zmax - zmin) > 0.4:
                cylinder.append(s)                                  # the cylinder surface
            else:
                walls.append(s)
        for tag, faces in (("inlet", inlet), ("outlet", outlet), ("walls", walls),
                           ("cylinder", cylinder)):
            gmsh.model.addPhysicalGroup(2, faces, name=tag)
        gmsh.model.addPhysicalGroup(3, [v for _, v in gmsh.model.getEntities(3)], name="fluid")
        f = gmsh.model.mesh.field
        dist = f.add("Distance")
        f.setNumbers(dist, "SurfacesList", cylinder)
        f.setNumber(dist, "Sampling", 100)
        thr = f.add("Threshold")
        f.setNumber(thr, "InField", dist)
        f.setNumber(thr, "SizeMin", lc_cyl)
        f.setNumber(thr, "SizeMax", lc)
        f.setNumber(thr, "DistMin", 0.0)
        f.setNumber(thr, "DistMax", thickness)
        f.setAsBackgroundMesh(thr)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.ElementOrder", 2)
        gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(3)
        out = HERE / f"{name}.msh"
        gmsh.write(str(out))
        return out
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    for name, (geo, order, quads) in MESHES.items():
        print("wrote", make(name, order, quads, HERE / geo))
    print("wrote", make("cylinder_bl", 2, False, HERE / "cylinder.geo",
                        boundary_layer=([5, 6, 7, 8], 0.0015, 1.2, 0.012)))
    print("wrote", make("cylinder_bl_fine", 2, False, HERE / "cylinder.geo",
                        boundary_layer=([5, 6, 7, 8], 0.0008, 1.15, 0.012), size_factor=0.5))
    print("wrote", make("cylinder3d1z_bl", 2, False, HERE / "cylinder3d1z.geo",
                        boundary_layer=([5, 6, 7, 8], 0.002, 1.25, 0.012), size_factor=1.6))
    print("wrote", make_flat_plate())
    print("wrote", make_turbulent_flat_plate())
    print("wrote", make_turbulent_flat_plate_medium())
    print("wrote", make_airfoil())
    print("wrote", make_box_tet10())
    print("wrote", make_cylinder3d())
