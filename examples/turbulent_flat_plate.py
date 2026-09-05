"""Turbulent flat plate with the Spalart-Allmaras model: skin friction and law of the wall.

    python examples/turbulent_flat_plate.py [--re 1e6] [--outdir turb_plate_out] [--no-plot]

Zero-pressure-gradient plate on ``y = 0, 0 <= x <= 2`` in the domain
``[-0.5, 2.5] x [0, 1]`` (mesh ``flat_plate_turb.msh``: 30 quadrilateral
boundary-layer cells from the wall, first cell 2e-5, inside a triangular
mesh), uniform inflow U = 1 on the inlet and the top, symmetry upstream and
downstream of the plate, do-nothing outflow.  ``--re`` is ``U / nu`` per
unit length (default 1e6, so Re_x runs to 2e6 at the trailing edge and the
first cell sits at y+ of about 1).

Compared with:
- the turbulent skin-friction correlation of White (Viscous Fluid Flow,
  3rd ed., eq. 6-78) ``Cf = 0.455 / ln^2(0.06 Re_x)`` and the 1/7-power law
  ``Cf = 0.0576 Re_x^(-1/5)``;
- the law of the wall ``u+ = y+`` (viscous sublayer) and
  ``u+ = ln(y+)/kappa + B`` with kappa = 0.41, B = 5.0 (log layer), at
  the station x = 1.5.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.incompressible import FlowProblem  # noqa: E402
from aa540fem.turbulence import solve_rans  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MESHES = HERE / "meshes"
KAPPA, B = 0.41, 5.0


def cf_white(rex):
    return 0.455 / np.log(0.06 * np.maximum(rex, 10.0)) ** 2


def cf_power(rex):
    return 0.0576 * np.maximum(rex, 1.0) ** (-0.2)


def law_of_the_wall(yplus):
    yplus = np.asarray(yplus, dtype=float)
    return np.where(yplus < 11.0, yplus, np.log(np.maximum(yplus, 1e-12)) / KAPPA + B)


def run(mesh_file=MESHES / "flat_plate_turb.msh", re=1e6, station=1.5, verbose=True, **kw):
    mesh = read_mesh(mesh_file)
    nu = 1.0 / re
    prob = FlowProblem(mesh, mu=nu, rho=1.0, stabilisation=True,
                       bc={"inlet": (1.0, 0.0), "top": (1.0, 0.0), "symmetry": (None, 0.0),
                           "plate": (0.0, 0.0), "outlet": "open"})
    rans = solve_rans(prob, wall_tags=["plate"], verbose=verbose, **kw)
    flow = rans.flow

    tr = flow.wall_traction("plate")
    rex = tr["x"] / nu
    cf = 2.0 * tr["tx"]

    # law of the wall at the station: friction velocity from the local wall shear
    i = np.argmin(np.abs(tr["x"] - station))
    u_tau = np.sqrt(max(tr["tx"][i], 1e-300))
    near = np.abs(mesh.x - station) < 5e-3
    order = np.argsort(mesh.y[near])
    yplus = mesh.y[near][order] * u_tau / nu
    uplus = flow.u[near][order] / u_tau
    return rans, {"Re_x": rex, "Cf": cf, "Cf_white": cf_white(rex), "Cf_power": cf_power(rex),
                  "yplus": yplus, "uplus": uplus, "u_tau": u_tau, "station": station,
                  "nu": nu}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "flat_plate_turb.msh"))
    parser.add_argument("--re", type=float, default=1e6, help="U / nu per unit length")
    parser.add_argument("--max-outer", type=int, default=40)
    parser.add_argument("--tol", type=float, default=1e-3)
    parser.add_argument("--viscosity-ramp", type=float, nargs="+", default=[100.0, 10.0, 1.0],
                        help="laminar-viscosity factors of the start-up continuation")
    parser.add_argument("--outdir", default="turb_plate_out")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    rans, res = run(args.mesh, args.re, max_outer=args.max_outer, tol=args.tol,
                    viscosity_ramp=tuple(args.viscosity_ramp))
    print(f"RANS: {len(rans.history)} outer iterations, converged = {rans.converged}, "
          f"{time.time() - t0:.0f} s")
    rex, cf = res["Re_x"], res["Cf"]
    sel = (rex > 3e5) & (rex < 2e6)
    dev = 100 * np.abs(cf[sel] / res["Cf_white"][sel] - 1)
    print(f"Cf vs White for 3e5 < Re_x < 2e6: max {dev.max():.1f} %, mean {dev.mean():.1f} %")
    for r in (3e5, 5e5, 1e6, 1.5e6, 2e6):
        i = np.argmin(np.abs(rex - r))
        print(f"  Re_x = {rex[i]:.2e}: Cf = {cf[i]:.5f}   White {res['Cf_white'][i]:.5f}   "
              f"1/7 law {res['Cf_power'][i]:.5f}")
    yp, up = res["yplus"], res["uplus"]
    log = (yp > 30) & (yp < 300)
    print(f"law of the wall at x = {res['station']}: first node y+ = {yp[1]:.2f}, "
          f"log layer max |u+ - law| = {np.abs(up[log] - law_of_the_wall(yp[log])).max():.2f} "
          f"({log.sum()} nodes)")
    rans.save(outdir / "turbulent_plate.vtu")
    print(f"Saved {outdir / 'turbulent_plate.vtu'}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        ax1.loglog(rex, cf, ".", ms=3, label="FEM + Spalart-Allmaras")
        rr = np.logspace(4, np.log10(rex.max()), 200)
        ax1.loglog(rr, cf_white(rr), "k-", label="White: 0.455 / ln^2(0.06 Re_x)")
        ax1.loglog(rr, cf_power(rr), "k--", label="1/7 power law")
        ax1.loglog(rr, 0.664 / np.sqrt(rr), "k:", label="laminar (Blasius)")
        ax1.set_xlabel("Re_x")
        ax1.set_ylabel("C_f")
        ax1.set_ylim(1e-3, 2e-2)
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3, which="both")
        yy = np.logspace(-1, 4, 300)
        ax2.semilogx(yy, law_of_the_wall(yy), "k-", label="u+ = y+ ;  ln(y+)/0.41 + 5")
        ax2.semilogx(yp[yp > 0], up[yp > 0], "o", ms=3, label=f"FEM, x = {res['station']}")
        ax2.set_xlabel("y+")
        ax2.set_ylabel("u+")
        ax2.set_xlim(0.1, 1e4)
        ax2.legend()
        ax2.grid(alpha=0.3, which="both")
        fig.savefig(outdir / "turbulent_plate.png", dpi=150, bbox_inches="tight")
        print(f"Saved {outdir / 'turbulent_plate.png'}")
    return rans, res


if __name__ == "__main__":
    main()
