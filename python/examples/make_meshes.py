"""Generate the meshes used by the tests and examples with Gmsh.

Requires the ``gmsh`` Python package (``pip install gmsh``).  Writes the
annulus meshes ``annulus_tri.msh``, ``annulus_tri6.msh``, ``annulus_quad.msh``,
``annulus_quad9.msh`` and the cylinder-in-channel mesh ``cylinder_tri6.msh``
next to this file.
"""

from __future__ import annotations

import pathlib

import gmsh

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


if __name__ == "__main__":
    for name, (geo, order, quads) in MESHES.items():
        print("wrote", make(name, order, quads, HERE / geo))
