"""Laminar flat plate: skin friction and velocity profiles against Blasius.

    python examples/flat_plate.py [--re 1e5] [--continuation auto] [--vtk plate.vtu] [--no-plot]

Domain [-0.5, 2] x [0, 0.5] with the plate on y = 0, 0 <= x <= 1.5 (no
slip), a symmetry line upstream and downstream of it, uniform inflow on the
left and top, do-nothing outflow.  Mesh ``flat_plate_bl.msh`` (Gmsh
boundary-layer quadrilaterals on the plate, triangles elsewhere,
``make_meshes.make_flat_plate``).  Reynolds number ``U L / nu`` with L = 1.
The steady solution uses SUPG/grad-div stabilisation and pseudo-transient
continuation.  Compared with the Blasius similarity solution:
``Cf = 0.664 / sqrt(Re_x)`` and ``u/U = f'(eta)``, ``eta = y sqrt(U/(nu x))``.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MESHES = HERE / "meshes"


def blasius():
    """Blasius boundary layer ``f''' + f f''/2 = 0``: returns ``f''(0)`` and ``f'(eta)``."""
    rhs = lambda e, y: [y[1], y[2], -0.5 * y[0] * y[2]]

    def shoot(fpp0):
        return solve_ivp(rhs, [0, 10], [0, 0, fpp0], rtol=1e-10, atol=1e-12).y[1, -1] - 1.0

    fpp0 = brentq(shoot, 0.2, 0.5)
    sol = solve_ivp(rhs, [0, 10], [0, 0, fpp0], rtol=1e-10, atol=1e-12, dense_output=True)
    return fpp0, lambda eta: sol.sol(np.minimum(np.asarray(eta, dtype=float), 10.0))[1]


def run(mesh_file=MESHES / "flat_plate_bl.msh", nu=1e-5, continuation="auto", verbose=True):
    mesh = read_mesh(mesh_file)
    prob = FlowProblem(mesh, mu=nu, rho=1.0, stabilisation=True,
                       bc={"inlet": (1.0, 0.0), "top": (1.0, 0.0), "symmetry": (None, 0.0),
                           "plate": (0.0, 0.0), "outlet": "open"})
    sol = solve_flow(prob, continuation=continuation, verbose=verbose)

    tr = sol.wall_traction("plate")
    rex = tr["x"] / nu
    cf = 2.0 * tr["tx"]
    cf_blasius = 0.664 / np.sqrt(np.maximum(rex, 1.0))
    _, fprime = blasius()
    profile_errors = {}
    profiles = {}
    for xs in (0.3, 0.5, 1.0):
        near = np.abs(mesh.x - xs) < 3e-3
        eta = mesh.y[near] * np.sqrt(1.0 / (nu * xs))
        keep = eta < 8
        u = sol.u[near][keep]
        profiles[xs] = (eta[keep], u)
        profile_errors[xs] = float(np.abs(u - fprime(eta[keep])).max())
    return sol, {"Re_x": rex, "Cf": cf, "Cf_blasius": cf_blasius, "profiles": profiles,
                 "profile_errors": profile_errors, "fprime": fprime}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "flat_plate_bl.msh"))
    parser.add_argument("--re", type=float, default=1e5, help="U L / nu with L = 1")
    parser.add_argument("--continuation", default="auto", choices=["auto", "newton", "ptc"])
    parser.add_argument("--vtk", default="")
    parser.add_argument("--no-plot", action="store_true")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")

    sol, res = run(args.mesh, nu=1.0 / args.re, continuation=args.continuation)
    print(f"{sol.info['continuation']}: {sol.info['iterations']} iterations, "
          f"converged = {sol.info['converged']}")
    rex, cf, cfb = res["Re_x"], res["Cf"], res["Cf_blasius"]
    sel = (rex > 1e4) & (rex < 1e5)
    dev = 100 * np.abs(cf[sel] / cfb[sel] - 1)
    print(f"Cf vs Blasius for 1e4 < Re_x < 1e5: max {dev.max():.1f} %, mean {dev.mean():.1f} %")
    for x, err in res["profile_errors"].items():
        print(f"profile at x = {x}: max |u/U - f'(eta)| = {err:.4f}")
    if args.vtk:
        sol.save(args.vtk)
        print(f"Saved {args.vtk}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        ax1.loglog(rex, cf, ".", ms=3, label="FEM")
        rr = np.logspace(3, np.log10(rex.max()), 100)
        ax1.loglog(rr, 0.664 / np.sqrt(rr), "k-", label="Blasius 0.664 / sqrt(Re_x)")
        ax1.set_xlabel("Re_x")
        ax1.set_ylabel("C_f")
        ax1.legend()
        ax1.grid(alpha=0.3, which="both")
        eta = np.linspace(0, 8, 200)
        ax2.plot(res["fprime"](eta), eta, "k-", label="Blasius")
        for x, (e, u) in res["profiles"].items():
            order = np.argsort(e)
            ax2.plot(u[order], e[order], "o", ms=3, label=f"x = {x}")
        ax2.set_xlabel("u / U")
        ax2.set_ylabel("eta")
        ax2.legend()
        ax2.grid(alpha=0.3)
        out = pathlib.Path("flat_plate.png")
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Saved {out}")
    return sol, res


if __name__ == "__main__":
    main()
