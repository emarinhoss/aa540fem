"""Heat conduction across an annulus read from a Gmsh mesh.

    python examples/annulus.py [--mesh examples/annulus_tri6.msh] [--method cg]
                               [--vtk annulus.vtu] [--no-plot]

T = 0 on the inner circle (r = 1), T = 1 on the outer circle (r = 2),
isotropic conductivity: the exact solution is T = ln(r) / ln(2).  Prints
the L2 error and writes a ParaView file.  The meshes are produced by
``make_meshes.py`` (linear and quadratic triangles and quadrilaterals).
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from aa540fem import DIRICHLET, Problem, error_norms, read_mesh, solve  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent


def exact(x, y):
    return np.log(np.hypot(x, y)) / np.log(2.0)


def exact_grad(x, y):
    r2 = x * x + y * y
    return x / r2 / np.log(2.0), y / r2 / np.log(2.0)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(HERE / "annulus_tri6.msh"))
    parser.add_argument("--method", default="direct", choices=["direct", "cg"])
    parser.add_argument("--vtk", default="annulus.vtu")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    mesh = read_mesh(args.mesh)
    problem = Problem(mesh=mesh,
                      bc_type={"inner": DIRICHLET, "outer": DIRICHLET},
                      bc_val={"inner": 0.0, "outer": 1.0})
    sol = solve(problem, verbose=True, method=args.method)

    errs = error_norms(mesh, sol.T, exact, exact_grad)
    print(f"{mesh.n_nodes} nodes, {mesh.n_elems} elements ({', '.join(mesh.cells)})")
    print(f"L2 error {errs['L2']:.3e}, H1 error {errs['H1']:.3e}")
    if args.vtk:
        sol.save(args.vtk)
        print(f"Saved {args.vtk}")

    if not args.no_plot:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri

        tri = mtri.Triangulation(mesh.x, mesh.y, mesh.triangulation())
        fig, ax = plt.subplots()
        cs = ax.tricontourf(tri, sol.T, 20)
        ax.triplot(tri, lw=0.2, color="k")
        ax.set_aspect("equal")
        fig.colorbar(cs, ax=ax, label="T")
        plt.show()
    return sol


if __name__ == "__main__":
    main()
