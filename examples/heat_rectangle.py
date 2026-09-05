"""Driver script.  Port of ``main.m``: edit the user inputs below and run

    python examples/heat_rectangle.py [--mesh FILE.msh] [--vtk OUT.vtu] [--method direct|cg]
                                      [--output temperature.png] [--show]

The only files that should need editing are this one and
``src/aa540fem/transport/material.py``.  With ``--mesh`` the boundary
dictionaries below must use the tags (Gmsh physical group names) of that
file instead of the four sides of the rectangle.  Giving the material
function a parameter named ``T`` makes the conductivity temperature
dependent; the problem is then solved with Newton's method.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import Problem, solve  # noqa: E402

# ---------------------------------------------------------------- User inputs
a = 4.0            # horizontal length of the domain
b = 6.0            # vertical length of the domain

elems = 50         # cells along x and y (must be equal)

elem_type = 3      # 1 --> 3-node linear triangles
                   # 2 --> 4-node linear quadrilaterals
                   # 3 --> 9-node quadratic quadrilaterals
                   # "triangle6" --> 6-node quadratic triangles

# Boundary condition types per boundary tag: 0 --> Dirichlet, 1 --> Neumann.
# Tags not listed are natural (zero flux).
bc_type = {
    "top": 0,
    "left": 1,
    "right": 1,
    "bottom": 0,
}

# Values at the boundaries: temperature for Dirichlet, normal flux for
# Neumann.  A value may also be a function of (x, y).
bc_val = {
    "top": 100.0,
    "left": 0.0,
    "right": 0.0,
    "bottom": 0.0,
}

# Order of the quadrature.  Quadrilaterals: Gauss-Legendre points per
# direction (1-4 in the original).  Triangles: 1, 3, 4, 6, 7, 9, 12 or 13
# Gauss points on the triangle.  None picks the lowest order that fully
# integrates the chosen element (1, 2, 3 for element types 1, 2, 3); the
# MATLAB default of 1 leaves 4- and 9-node elements rank deficient.
order = None
# -----------------------------------------------------------------------------


def build_problem(mesh_file=None) -> Problem:
    mesh = None
    if mesh_file:
        from aa540fem import read_mesh

        mesh = read_mesh(mesh_file)
    return Problem(a=a, b=b, elems=elems, elem_type=elem_type, order=order,
                   mesh=mesh, bc_type=dict(bc_type), bc_val=dict(bc_val))


def plot(solution, output=None, show=False, levels=20):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    mesh = solution.mesh
    tri = mtri.Triangulation(mesh.x, mesh.y, mesh.triangulation())
    fig, ax = plt.subplots()
    cs = ax.tricontour(tri, solution.T, levels)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Temperature")
    fig.colorbar(cs, ax=ax)
    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved {output}")
    if show:
        plt.show()
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", "-m", default=None,
                        help="mesh file to solve on instead of the rectangle (needs meshio)")
    parser.add_argument("--vtk", default=None, help="write the solution to this ParaView file")
    parser.add_argument("--method", default="direct", choices=["direct", "cg"],
                        help="linear solver (cg uses pyamg preconditioning if installed)")
    parser.add_argument("--output", "-o", default="temperature.png",
                        help="file for the contour plot ('' to skip)")
    parser.add_argument("--show", action="store_true", help="open an interactive window")
    args = parser.parse_args(argv)

    solution = solve(build_problem(args.mesh), verbose=True, method=args.method)
    print(f"T in [{solution.T.min():.6g}, {solution.T.max():.6g}]")
    if args.vtk:
        solution.save(args.vtk)
        print(f"Saved {args.vtk}")
    if args.output or args.show:
        plot(solution, args.output or None, args.show)
    return solution


if __name__ == "__main__":
    main()
