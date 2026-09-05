"""Lid-driven cavity: the classic incompressible Navier-Stokes benchmark.

    python examples/cavity.py [--re 100] [--elems 32] [--vtk cavity.vtu] [--no-plot]

Unit square, no-slip walls, lid velocity 1, Q2/Q1 Taylor-Hood elements and
Newton's method.  Prints the centreline velocity extrema next to the
reference values of Ghia, Ghia & Shin (1982) for Re = 100, 400 and 1000
(higher Re is reached by continuation from the previous solution).
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from aa540fem import geometry  # noqa: E402
from aa540fem.flow import FlowProblem, solve_flow  # noqa: E402

GHIA = {   # Re: (u_min on x = 0.5, y of u_min, v_max on y = 0.5, v_min on y = 0.5)
    100: (-0.2109, 0.4531, 0.1753, -0.2453),
    400: (-0.3273, 0.2813, 0.3020, -0.4499),
    1000: (-0.3829, 0.1719, 0.3709, -0.5155),
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--re", type=float, default=100.0)
    parser.add_argument("--elems", type=int, default=32)
    parser.add_argument("--vtk", default="cavity.vtu")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    mesh = geometry(1.0, 1.0, args.elems, "quad9")
    walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
    U = None
    # continuation: Re 100 -> 400 -> 1000 -> target
    stages = [re for re in (100.0, 400.0, 1000.0) if re < args.re] + [args.re]
    for re in stages:
        prob = FlowProblem(mesh, mu=1.0 / re, rho=1.0, bc={**walls, "top": (1.0, 0.0)})
        sol = solve_flow(prob, U0=U, verbose=True)
        U = sol.U
        B = sol.assembler.Bx + sol.assembler.By
        print(f"Re = {re:g}: {sol.info['iterations']} Newton iterations, "
              f"discrete continuity |B u| = {np.linalg.norm(B @ sol.U):.2e}")

    mid = np.isclose(mesh.x, 0.5)
    u, y = sol.u[mid], mesh.y[mid]
    midy = np.isclose(mesh.y, 0.5)
    v, x = sol.v[midy], mesh.x[midy]
    print(f"u_min = {u.min():.4f} at y = {y[u.argmin()]:.4f};  "
          f"v_max = {v.max():.4f} at x = {x[v.argmax()]:.4f};  "
          f"v_min = {v.min():.4f} at x = {x[v.argmin()]:.4f}")
    if int(args.re) in GHIA:
        g = GHIA[int(args.re)]
        print(f"Ghia:   u_min = {g[0]:.4f} at y = {g[1]:.4f};  v_max = {g[2]:.4f};  "
              f"v_min = {g[3]:.4f}")
    if args.vtk:
        try:
            sol.save(args.vtk)
            print(f"Saved {args.vtk}")
        except ImportError as exc:
            print(exc)

    if not args.no_plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
        ax1.streamplot(mesh.X, mesh.Y, mesh.to_grid(sol.u), mesh.to_grid(sol.v), density=1.5,
                       color=mesh.to_grid(sol.speed), cmap="viridis")
        ax1.set_aspect("equal")
        ax1.set_title(f"streamlines, Re = {args.re:g}")
        ax2.plot(u, y, label="u on x = 0.5")
        ax2.plot(x, v, label="v on y = 0.5")
        ax2.legend()
        plt.show()
    return sol


if __name__ == "__main__":
    main()
