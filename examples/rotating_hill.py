"""Gaussian hill advected one full revolution by a rotating velocity field.

    python examples/rotating_hill.py [--elems 40] [--steps 200] [--sigma 0.15]
                                     [--scheme rk45|theta] [--series hill] [--no-plot]

Domain [-1, 1]^2, u = (-y, x), negligible diffusion, adaptive RK45 (default)
or Crank-Nicolson in time with SUPG in space.  After one revolution the exact solution equals the
initial hill, so the peak retention and L2 error measure the numerical
dissipation and dispersion.  Writes a ParaView time series (``hill.pvd``).
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import DIRICHLET, Problem, error_norms, geometry, solve_transient  # noqa: E402
from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--elems", type=int, default=40)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--sigma", type=float, default=0.15)
    parser.add_argument("--elem-type", default="quad")
    parser.add_argument("--scheme", default="rk45", choices=["rk45", "theta"])
    parser.add_argument("--series", default="hill", help="prefix of the .vtu/.pvd output")
    parser.add_argument("--no-plot", action="store_true")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")

    sigma = args.sigma
    hill = lambda x, y: np.exp(-((x - 0.5) ** 2 + y ** 2) / (2 * sigma ** 2))

    mesh = geometry(2.0, 2.0, args.elems, args.elem_type)
    for arr in (mesh.points, mesh.X, mesh.Y):
        arr -= 1.0                           # shift to [-1, 1]^2
    sides = ("top", "right", "left", "bottom")
    p = Problem(mesh=mesh, bc_type={s: DIRICHLET for s in sides},
                bc_val={s: 0.0 for s in sides},
                material=lambda x, y: (1e-6, 0.0, 0.0, 1e-6, 0.0),
                velocity=lambda x, y: (-y, x))

    period = 2 * np.pi
    if args.scheme == "rk45":
        sol = solve_transient(p, dt=period / args.steps, t_end=period, T0=hill, scheme="rk45",
                              output_interval=period / 20, verbose=True)
    else:
        sol = solve_transient(p, dt=period / args.steps, t_end=period, theta=0.5, T0=hill,
                              scheme="theta", store_every=max(1, args.steps // 20), verbose=True)
    err = error_norms(mesh, sol.T, hill)["L2"]
    print(f"peak retained: {sol.T.max():.3f} (initial 1), min T = {sol.T.min():+.4f}, "
          f"L2 error = {err:.3e}")
    if args.series:
        try:
            pvd = sol.save_series(args.series)
            print(f"Saved {pvd}")
        except ImportError as exc:
            print(exc)

    if not args.no_plot:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(9, 4))
        for ax, (title, T) in zip(axes, [("initial", sol.snapshots[0]), ("after one turn", sol.T)]):
            cs = ax.contourf(mesh.X, mesh.Y, mesh.to_grid(T), np.linspace(-0.05, 1.0, 22))
            ax.set_aspect("equal")
            ax.set_title(title)
        fig.colorbar(cs, ax=axes)
        plt.show()
    return sol


if __name__ == "__main__":
    main()
