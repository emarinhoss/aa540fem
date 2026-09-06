"""Three-dimensional lid-driven cavity: centreline profiles and mesh convergence.

    python examples/cavity3d.py [--re 100] [--elems 8 12] [--elem-type hexahedron27] [--no-prompt]

Unit cube, lid at y = 1 moving with ``u = 1`` in x, all other walls no-slip,
steady Navier-Stokes at the given Reynolds number (``mu = 1 / Re``).  For
every mesh the script reports the extrema of ``u`` on the vertical centreline
``x = z = 0.5`` and of ``v`` on the horizontal centreline ``y = z = 0.5``
(the quantities tabulated in the 3-D cavity literature: Ku, Hirsh and
Taylor 1987; Albensoeder and Kuhlmann 2005; Wong and Baker 2002) and writes
the profiles to ``<outdir>/centreline_<n>.csv``.  Successive meshes show the
convergence of the extrema; compare them with the published tables.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402
from aa540fem.core.mesh import box  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow  # noqa: E402


def problem(mesh, re: float):
    walls = {t: (0.0, 0.0, 0.0) for t in ("left", "right", "bottom", "front", "back")}
    lid = lambda x, y, z: np.ones_like(x)
    return FlowProblem(mesh, mu=1.0 / re, rho=1.0, stabilisation=True,
                       bc={**walls, "top": (lid, 0.0, 0.0)})


def centrelines(sol):
    mesh = sol.mesh
    vert = np.isclose(mesh.x, 0.5) & np.isclose(mesh.z, 0.5)
    horz = np.isclose(mesh.y, 0.5) & np.isclose(mesh.z, 0.5)
    iy, ix = np.argsort(mesh.y[vert]), np.argsort(mesh.x[horz])
    return ((mesh.y[vert][iy], sol.u[vert][iy]), (mesh.x[horz][ix], sol.v[horz][ix]))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--re", type=float, default=100.0)
    parser.add_argument("--elems", type=int, nargs="+", default=[8, 12])
    parser.add_argument("--elem-type", default="hexahedron27", choices=["hexahedron27", "tetra10"])
    parser.add_argument("--outdir", default="cavity3d_out")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    U0 = None
    rows = []
    for n in args.elems:
        mesh = box(1.0, 1.0, 1.0, n, args.elem_type)
        prob = problem(mesh, args.re)
        t0 = time.time()
        sol = solve_flow(prob, U0=None if U0 is None else U0(mesh), verbose=False)
        wall = time.time() - t0
        (y, u), (x, v) = centrelines(sol)
        row = {"n": n, "unknowns": sol.space.ndof, "iterations": sol.info["iterations"],
               "time": wall, "u_min": u.min(), "y_umin": y[u.argmin()],
               "v_max": v.max(), "x_vmax": x[v.argmax()], "v_min": v.min(),
               "x_vmin": x[v.argmin()]}
        rows.append(row)
        print(f"n = {n:3d} ({sol.space.ndof} unknowns, {sol.info['iterations']} Newton it., "
              f"{wall:.0f} s): u_min = {row['u_min']:.4f} at y = {row['y_umin']:.3f}, "
              f"v_max = {row['v_max']:.4f} at x = {row['x_vmax']:.3f}, "
              f"v_min = {row['v_min']:.4f} at x = {row['x_vmin']:.3f}")
        np.savetxt(outdir / f"centreline_{n}.csv", np.column_stack([y, u, x, v]), delimiter=",",
                   header="y,u(0.5,y,0.5),x,v(x,0.5,0.5)", comments="")
        # interpolate the coarse solution to the next mesh as the initial guess
        from scipy.interpolate import LinearNDInterpolator

        from aa540fem.incompressible.space import TaylorHoodSpace

        interp = [LinearNDInterpolator(mesh.points, comp, fill_value=0.0)
                  for comp in sol.velocity]

        def U0(next_mesh, interp=interp):
            vel = [f(next_mesh.points) for f in interp]
            return np.concatenate(vel + [np.zeros(TaylorHoodSpace(next_mesh).Np)])

    if len(rows) > 1:
        print("convergence of the extrema:")
        for key in ("u_min", "v_max", "v_min"):
            vals = [r[key] for r in rows]
            print(f"  {key}: " + ", ".join(f"{v:.4f}" for v in vals)
                  + f"   (change {abs(vals[-1] - vals[-2]):.4f} on the last refinement)")
    return rows


if __name__ == "__main__":
    main()
