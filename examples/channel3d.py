"""Three-dimensional channel flow: plane Poiseuille flow extruded in z.

    python examples/channel3d.py [--elem-type hexahedron27|tetra10] [--mesh box_tet10.msh]
                                 [--elems 4 3 2] [--re 20] [--no-prompt] [--no-plot]

Solves the steady Navier-Stokes equations in the box [0, 2] x [0, 1] x [0, 1]
with a parabolic inflow ``u = 4 U y (1 - y)``, no-slip walls at y = 0 and
y = 1, symmetry planes at z = 0 and z = 1 (only ``u_z`` fixed) and an open
outlet, on a structured mesh (``--elem-type``) or on the committed
unstructured tetrahedral mesh (``--mesh``).  The exact solution is the
inflow profile everywhere with the pressure ``8 mu U (2 - x)``; the script
reports the errors, the wall shear force and writes a ``.vtu`` for ParaView.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402
from aa540fem.core.mesh import box  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow  # noqa: E402

MESHES = pathlib.Path(__file__).resolve().parent / "meshes"
LENGTH, HEIGHT, DEPTH = 2.0, 1.0, 1.0


def problem(mesh, re: float, umax: float = 1.0):
    mu = umax * HEIGHT / re
    inflow = lambda x, y, z: 4.0 * umax * y * (HEIGHT - y) / HEIGHT ** 2
    return FlowProblem(mesh, mu=mu, rho=1.0, stabilisation=True,
                       bc={"left": (inflow, 0.0, 0.0), "bottom": (0.0, 0.0, 0.0),
                           "top": (0.0, 0.0, 0.0), "front": (None, None, 0.0),
                           "back": (None, None, 0.0), "right": "open"})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--elem-type", default="hexahedron27", choices=["hexahedron27", "tetra10"])
    parser.add_argument("--elems", type=int, nargs=3, default=[4, 3, 2])
    parser.add_argument("--mesh", default=None, help="an unstructured .msh (e.g. box_tet10.msh)")
    parser.add_argument("--re", type=float, default=20.0)
    parser.add_argument("--outdir", default="channel3d_out")
    parser.add_argument("--no-plot", action="store_true")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")

    if args.mesh:
        mesh = read_mesh(args.mesh if "/" in args.mesh else MESHES / args.mesh)
    else:
        mesh = box(LENGTH, HEIGHT, DEPTH, args.elems, args.elem_type)
    prob = problem(mesh, args.re)
    print(f"{mesh.n_nodes} nodes, {mesh.n_elems} elements ({', '.join(mesh.cells)}), "
          f"Re = {args.re:g}")
    t0 = time.time()
    sol = solve_flow(prob, verbose=True)
    print(f"steady solve: {sol.info['iterations']} Newton iterations, {time.time() - t0:.1f} s")

    u_exact = 4.0 * mesh.y * (HEIGHT - mesh.y) / HEIGHT ** 2
    p_exact = 8.0 * prob.mu * (LENGTH - mesh.x) / HEIGHT ** 2
    print(f"max |u - exact| = {np.abs(sol.u - u_exact).max():.2e}, "
          f"max |p - exact| = {np.abs(sol.p_nodal - p_exact).max():.2e}, "
          f"|div u| = {sol.divergence_norm():.1e}")
    fx, fy, fz = sol.forces("bottom")
    print(f"bottom wall force: Fx = {fx:.6f} (exact {4 * prob.mu * LENGTH * DEPTH / HEIGHT:.6f}), "
          f"Fy = {fy:.4f}, Fz = {fz:.1e}")

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        print("Saved", sol.save(outdir / "channel3d.vtu"))
    except ImportError:
        print("meshio not installed; no .vtu written")
    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        mid = np.isclose(mesh.z, 0.5 * DEPTH) & np.isclose(mesh.x, 0.5 * LENGTH)
        order = np.argsort(mesh.y[mid])
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot(sol.u[mid][order], mesh.y[mid][order], "o-", label="computed")
        ax.plot(u_exact[mid][order], mesh.y[mid][order], "k--", label="exact")
        ax.set_xlabel("u"), ax.set_ylabel("y"), ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "profile.png", dpi=120)
        print("Saved", outdir / "profile.png")


if __name__ == "__main__":
    main()
