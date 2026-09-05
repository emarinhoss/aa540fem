"""Flow past a cylinder in a channel (Schaefer & Turek 1996 benchmark, case 2D-1).

    python examples/cylinder.py [--mesh examples/cylinder_tri6.msh] [--umax 0.3]
                                [--vtk cylinder.vtu] [--no-plot]

Channel 2.2 x 0.41 with a cylinder of diameter D = 0.1 at (0.2, 0.2),
parabolic inflow of maximum speed U_max (mean U = 2/3 U_max), nu = 1e-3,
rho = 1: Reynolds number Re = U D / nu = 20 for U_max = 0.3.  Reference
values: C_D = 5.5795, C_L = 0.0106, pressure difference between the front
and rear stagnation points dp = 0.1175.  Coefficients use
C = 2 F / (rho U^2 D).  P2/P1 triangles from ``cylinder.geo`` and Newton.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MESHES = HERE / "meshes"
REFERENCE = {"C_D": 5.57953523384, "C_L": 0.010618948146, "dp": 0.11752016697}


def run(mesh_file, umax=0.3, nu=1e-3, verbose=True):
    mesh = read_mesh(mesh_file)
    H = 0.41
    inflow = lambda x, y: 4.0 * umax * y * (H - y) / H ** 2
    prob = FlowProblem(mesh, mu=nu, rho=1.0,
                       bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0),
                           "cylinder": (0.0, 0.0), "outlet": "open"})
    sol = solve_flow(prob, verbose=verbose)

    umean = 2.0 / 3.0 * umax
    D = 0.1
    fx, fy = sol.forces("cylinder")
    scale = 2.0 / (1.0 * umean ** 2 * D)
    front = np.argmin(np.hypot(mesh.x - 0.15, mesh.y - 0.2))
    rear = np.argmin(np.hypot(mesh.x - 0.25, mesh.y - 0.2))
    p = sol.p_nodal
    results = {"C_D": scale * fx, "C_L": scale * fy, "dp": p[front] - p[rear],
               "Re": umean * D / nu}
    return sol, results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "cylinder_tri6.msh"))
    parser.add_argument("--umax", type=float, default=0.3)
    parser.add_argument("--vtk", default="cylinder.vtu")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    sol, res = run(args.mesh, args.umax)
    print(f"Re = {res['Re']:.1f}: {sol.mesh.n_nodes} nodes, {sol.mesh.n_elems} elements, "
          f"{sol.info['iterations']} Newton iterations")
    for key in ("C_D", "C_L", "dp"):
        line = f"{key:>4} = {res[key]:9.5f}"
        if abs(res["Re"] - 20) < 1e-6:
            line += f"   (reference {REFERENCE[key]:.5f}, " \
                    f"{100 * (res[key] - REFERENCE[key]) / REFERENCE[key]:+.2f} %)"
        print(line)
    if args.vtk:
        sol.save(args.vtk)
        print(f"Saved {args.vtk}")

    if not args.no_plot:
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri

        tri = mtri.Triangulation(sol.mesh.x, sol.mesh.y, sol.mesh.triangulation())
        fig, ax = plt.subplots(figsize=(11, 3))
        cs = ax.tricontourf(tri, sol.speed, 30)
        ax.set_aspect("equal")
        fig.colorbar(cs, ax=ax, label="|u|")
        plt.show()
    return sol, res


if __name__ == "__main__":
    main()
