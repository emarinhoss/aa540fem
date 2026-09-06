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
  the station x = 1.5;
- the Coles-Fernholz correlation ``Cf = 2 [ln(Re_theta)/0.384 + 4.127]^-2``
  at the momentum-thickness Reynolds number of the computed profile at
  that station (Nagib, Chauhan & Monkewitz 2007), which does not depend on
  where the boundary layer became turbulent.
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


def run(mesh_file=MESHES / "flat_plate_turb.msh", re=1e6, station=None, verbose=True,
        extrude=0.0, **kw):
    """``extrude > 0``: the 2-D mesh is extruded one layer over that depth and the
    same case is solved in 3-D between two symmetry planes (``u_z = 0``)."""
    mesh = read_mesh(mesh_file)
    if extrude:
        mesh = mesh.extrude(extrude, layers=1)
    nu = 1.0 / re
    bc = {"inlet": (1.0, 0.0), "top": (1.0, 0.0), "symmetry": (None, 0.0),
          "plate": (0.0, 0.0), "outlet": "open"}
    if mesh.dim == 3:
        bc = {tag: (spec if spec == "open" else spec + (0.0 if tag == "plate" else None,))
              for tag, spec in bc.items()}
        bc.update({"front": (None, None, 0.0), "back": (None, None, 0.0)})
    prob = FlowProblem(mesh, mu=nu, rho=1.0, stabilisation=True, bc=bc)
    # smooth initial profile (a boundary layer of thickness ~ 0.03 on the plate)
    # instead of the impulsive start, which the wall-resolved mesh cannot absorb
    def initial(x, y, z=None):
        delta = 0.03 * np.sqrt(np.maximum(x, 0.0) / 2.0 + 0.05)
        u = np.where(x > 0, 1.0 - np.exp(-y / delta), 1.0)
        zeros = np.zeros_like(y)
        return (u, zeros) if z is None else (u, zeros, zeros)

    rans = solve_rans(prob, wall_tags=["plate"], verbose=verbose, U0=initial, **kw)
    flow = rans.flow

    tr = flow.wall_traction("plate")
    rex = tr["x"] / nu
    cf = 2.0 * tr["tx"]
    if station is None:
        station = 0.75 * tr["x"].max()            # x = 1.5 on the plate of length 2

    # law of the wall at the station: friction velocity from the local wall shear
    i = np.argmin(np.abs(tr["x"] - station))
    u_tau = np.sqrt(max(tr["tx"][i], 1e-300))
    near = np.abs(mesh.x - station) < 5e-3
    if mesh.dim == 3:
        near &= np.isclose(mesh.z, mesh.z.min())              # the profile on one plane
    order = np.argsort(mesh.y[near])
    y, u = mesh.y[near][order], flow.u[near][order]
    yplus = y * u_tau / nu
    uplus = u / u_tau
    # momentum thickness of the profile and the Coles-Fernholz skin friction at
    # that Re_theta (Nagib, Chauhan & Monkewitz 2007: kappa 0.384, C 4.127),
    # the comparison that does not depend on where the boundary layer started
    inner = y < 0.1                               # the boundary layer, not the far field
    u_e = u[inner].max()                          # local edge velocity (displacement effect)
    theta = np.trapezoid(u[inner] / u_e * (1.0 - u[inner] / u_e), y[inner])
    re_theta = u_e * theta / nu
    cf_cf = 2.0 / (np.log(re_theta) / 0.384 + 4.127) ** 2
    return rans, {"Re_x": rex, "Cf": cf, "Cf_white": cf_white(rex), "Cf_power": cf_power(rex),
                  "yplus": yplus, "uplus": uplus, "u_tau": u_tau, "station": station,
                  "nu": nu, "Re_theta": re_theta, "Cf_station": 2.0 * u_tau ** 2,
                  "Cf_coles_fernholz": cf_cf}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "flat_plate_turb.msh"))
    parser.add_argument("--re", type=float, default=1e6, help="U / nu per unit length")
    parser.add_argument("--max-outer", type=int, default=40)
    parser.add_argument("--tol", type=float, default=1e-3)
    parser.add_argument("--viscosity-ramp", type=float, nargs="+", default=[100.0, 10.0, 1.0],
                        help="laminar-viscosity factors of the start-up continuation")
    parser.add_argument("--extrude", type=float, default=0.0, metavar="DEPTH",
                        help="solve the extruded 3-D case between symmetry planes")
    parser.add_argument("--outdir", default="turb_plate_out")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("-v", "--verbose", action="count", default=1,
                        help="-v: outer iterations (default), -vv: also the sub-solves")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    rans, res = run(args.mesh, args.re, verbose=args.verbose, max_outer=args.max_outer,
                    tol=args.tol, viscosity_ramp=tuple(args.viscosity_ramp),
                    extrude=args.extrude)
    print(f"RANS: {len(rans.history)} outer iterations, converged = {rans.converged}, "
          f"{time.time() - t0:.0f} s")
    rex, cf = res["Re_x"], res["Cf"]
    lo, hi = 0.15 * rex.max(), 0.98 * rex.max()      # 3e5 < Re_x < 2e6 at Re = 1e6
    sel = (rex > lo) & (rex < hi)
    for name, ref in (("White", res["Cf_white"]), ("the 1/7 law", res["Cf_power"])):
        dev = 100 * (cf[sel] / ref[sel] - 1)
        print(f"Cf vs {name} for {lo:.1e} < Re_x < {hi:.1e}: "
              f"{dev.min():+.1f} % to {dev.max():+.1f} %, mean {dev.mean():+.1f} %")
    for frac in (0.15, 0.25, 0.5, 0.75, 0.95):
        i = np.argmin(np.abs(rex - frac * rex.max()))
        print(f"  Re_x = {rex[i]:.2e}: Cf = {cf[i]:.5f}   White {res['Cf_white'][i]:.5f}   "
              f"1/7 law {res['Cf_power'][i]:.5f}")
    print(f"at x = {res['station']}: Re_theta = {res['Re_theta']:.0f}, "
          f"Cf = {res['Cf_station']:.5f}, "
          f"Coles-Fernholz Cf(Re_theta) = {res['Cf_coles_fernholz']:.5f} "
          f"({100 * (res['Cf_station'] / res['Cf_coles_fernholz'] - 1):+.1f} %)")
    yp, up = res["yplus"], res["uplus"]
    print(f"law of the wall at x = {res['station']}: first node y+ = {yp[1]:.2f}, "
          f"u_tau = {res['u_tau']:.4f}")
    for lo_p, hi_p in ((5, 30), (30, 300), (50, 300)):
        band = (yp > lo_p) & (yp < hi_p)
        if band.any():
            err = up[band] - law_of_the_wall(yp[band])
            print(f"  {lo_p} < y+ < {hi_p}: u+ - law in [{err.min():+.2f}, {err.max():+.2f}] "
                  f"({band.sum()} nodes)")
        else:
            print(f"  {lo_p} < y+ < {hi_p}: no nodes")
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
