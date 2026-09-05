"""Generate the meshes used by the tests and examples with Gmsh.

Requires the ``gmsh`` Python package (``pip install gmsh``).  Writes the
annulus meshes ``annulus_tri.msh``, ``annulus_tri6.msh``, ``annulus_quad.msh``,
``annulus_quad9.msh``, the cylinder-in-channel mesh ``cylinder_tri6.msh`` and
the NACA 0012 far-field mesh ``airfoil_naca0012_a5_tri6.msh`` next to this
file.
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


def make(name: str, order: int, quads: bool, geo: pathlib.Path = HERE / "annulus.geo"):
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(geo))
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


if __name__ == "__main__":
    for name, (geo, order, quads) in MESHES.items():
        print("wrote", make(name, order, quads, HERE / geo))
    print("wrote", make_airfoil())
