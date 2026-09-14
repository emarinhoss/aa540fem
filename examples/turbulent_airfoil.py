"""Turbulent NACA 0012 with the Spalart-Allmaras model: lift, drag and Cp vs NASA TMR.

    python examples/turbulent_airfoil.py [--alpha 0] [--re 6e6] [--outdir turb_airfoil_out]
                                         [--no-plot]

Wall-resolved unit-chord NACA 0012 (mesh ``airfoil_naca0012_turb.msh``: 34
quadrilateral layers from the surface, first cell 4e-6 so y+ is about 1 at
Re = 6e6, fanned around the sharp trailing edge, inside triangles, far
field at 100 chords).  The mesh is unrotated: the angle of attack sets the
far-field velocity ``(cos a, sin a)``, and drag and lift are the force
components along and normal to the free stream.  Fully turbulent, freestream
``nu_tilde = 3 nu`` (the solve_rans default, NASA TMR practice); the pressure
is zero far away (do-nothing outflow), so ``Cp = 2 p``.

Compared with the NASA Turbulence Modeling Resource 2-D NACA 0012 validation
case (https://tmbwg.github.io/turbmodels/naca0012_val_sa.html), SA model,
Re = 6e6, 897x257 grid, at M = 0.15 -- weakly compressible, so about a
percent of difference in the pressure coefficients is expected:

    alpha      CL CFL3D / FUN3D        CD CFL3D / FUN3D
      0        0                       0.00819 / 0.00812
     10        1.0909 / 1.0983         0.01231 / 0.01242
     15        1.5461 / 1.5547         0.02124 / 0.02159

Only ``alpha = 0`` (the default) reproduces the reference here; the lifting
cases do not converge on this mesh, and no choice of ramp or tolerance
fixes that -- the residual norm is dominated by the far field while the
Newton direction stops descending inside the boundary layer, where the
frozen SUPG ``tau`` is least consistent.  Section 4b of
``docs/turbulence.md`` records the diagnosis and what would unblock it.

The sub-solves are bounded (``sa_max_steps``, ``max_ptc``): a capped flow
solve acts as pseudo-time marching and the outer iterations carry the state
to the attractor, so single ``NOT converged`` lines are part of the process,
not a failure.
"""

from __future__ import annotations

import argparse
import csv
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

# NASA TMR reference (SA, Re = 6e6, M = 0.15, 897x257): alpha -> (CL, CD)
TMR_CFL3D = {0.0: (0.0, 0.00819), 10.0: (1.0909, 0.01231), 15.0: (1.5461, 0.02124)}
TMR_FUN3D = {0.0: (0.0, 0.00812), 10.0: (1.0983, 0.01242), 15.0: (1.5547, 0.02159)}


def run(mesh_file=MESHES / "airfoil_naca0012_turb.msh", alpha=0.0, re=6e6,
        verbose=True, restart=None, **kw):
    mesh = read_mesh(mesh_file)
    nu = 1.0 / re
    a = np.radians(alpha)
    ux, uy = np.cos(a), np.sin(a)
    bc = {"inlet": (ux, uy), "farfield": (ux, uy), "airfoil": (0.0, 0.0), "outlet": "open"}
    prob = FlowProblem(mesh, mu=nu, rho=1.0, stabilisation=True, bc=bc)
    if restart is not None:
        # state of a previous run at another angle: rotate its velocity field
        # towards the new free stream, blended by wall distance (the wall keeps
        # the old, wall-aligned flow; the far field turns fully), so neither
        # the flow nor the SA solve faces the raw far-field jump.  ``alpha0``
        # in the file records the angle the state was converged at
        state = np.load(restart)
        U = np.array(state["U"], dtype=float)
        da = a - np.radians(float(state["alpha0"])) if "alpha0" in state else a
        if abs(da) > 1e-12:
            d = mesh.wall_distance("airfoil")
            theta = da * np.minimum(d, 1.0)           # 0 at the wall -> da at 1 chord
            n = mesh.n_nodes
            u0, v0 = U[:n].copy(), U[n:2 * n].copy()
            U[:n] = np.cos(theta) * u0 - np.sin(theta) * v0
            U[n:2 * n] = np.sin(theta) * u0 + np.cos(theta) * v0
        kw.update(U0=U, nu_tilde0=state["nu_tilde"])
    else:
        # smooth start: the free stream damped over ~ the trailing-edge boundary-
        # layer thickness, which the wall-resolved mesh can absorb (impulsive
        # starts cannot)
        d = mesh.wall_distance("airfoil")
        damp = 1.0 - np.exp(-d / 0.01)
        kw.setdefault("U0", lambda x, y: (ux * damp, uy * damp))  # evaluated at the nodes
    # before the eddy viscosity develops, the laminar trailing edge is weakly
    # unsteady and the steady residual stalls near 1e-4 of the start: accept
    # the flow there (the coupled outer iterations tighten it warm-started).
    # atol floors the warm-started sub-solves (whose rtol target, relative to
    # an already small start residual, is unreachable in a pseudo-time limit
    # cycle) and max_ptc bounds them like sa_max_steps bounds the SA solves
    # after a restart at a new angle the loose atol would accept the flow
    # before the circulation has adapted to the rotated far field, and the
    # SA solve then blows up on the inconsistent state: tighten it there
    # dtau0 well below 1: at incidence the start is shedding-unsteady and PTC
    # rejections (a full inner Newton each) are far dearer than a few extra
    # small accepted steps while SER grows the pseudo-time step back
    fo = {"rtol": 5e-4, "atol": 1e-5 if restart is None else 1e-8, "max_ptc": 60,
          "dtau0": 0.1}
    fo.update(kw.pop("flow_options", None) or {})
    kw["flow_options"] = fo
    if kw.get("viscosity_ramp") is None:
        kw["viscosity_ramp"] = (1.0,) if restart is not None else (100.0, 10.0, 1.0)
    # the transitional states of the viscosity ramp can stall the SA sub-solve
    # in a pseudo-time limit cycle: bound it and let the coupling continue
    kw.setdefault("sa_max_steps", 60)
    rans = solve_rans(prob, wall_tags=["airfoil"], verbose=verbose, **kw)
    flow = rans.flow

    fx, fy = flow.forces("airfoil")
    cd = 2.0 * (fx * ux + fy * uy)               # chord 1, U 1, rho 1: coefficients are 2 F
    cl = 2.0 * (fy * ux - fx * uy)
    tr = flow.wall_traction("airfoil")
    w, nx, ny, p = tr["weight"], tr["nx"], tr["ny"], tr["p"]
    # pressure part of the drag from the traction quadrature; the rest is friction
    cd_p = 2.0 * (np.sum(w * p * nx) * ux + np.sum(w * p * ny) * uy)
    # skin friction: the viscous part of the traction on the body (t - p n),
    # projected on the free stream (the pressure part dwarfs it near the nose)
    cf = 2.0 * ((tr["tx"] - p * nx) * ux + (tr["ty"] - p * ny) * uy)
    return rans, {"alpha": alpha, "CL": cl, "CD": cd, "CD_p": cd_p, "CD_f": cd - cd_p,
                  "x": tr["x"], "y": tr["y"], "Cp": 2.0 * tr["p"], "Cf": cf, "nu": nu}


def _read_tmr_zone(path, alpha):
    """``x, cp`` of one alpha zone of the TMR Tecplot file ``n0012cp_cfl3d_sa.dat``."""
    xs, cps, active = [], [], False
    for line in open(path):
        if line.lstrip().startswith(("#", "variables")):
            continue
        if line.lstrip().startswith("zone"):
            active = f"alpha={alpha:g}" in line.replace(" ", "").replace('"', "")
            continue
        if active:
            x, cp = line.split()
            xs.append(float(x))
            cps.append(float(cp))
    return np.array(xs), np.array(cps)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "airfoil_naca0012_turb.msh"))
    parser.add_argument("--alpha", type=float, default=0.0,
                        help="angle of attack (degrees); only 0 converges on this mesh, "
                             "see section 4b of docs/turbulence.md")
    parser.add_argument("--re", type=float, default=6e6, help="chord Reynolds number")
    parser.add_argument("--max-outer", type=int, default=12)
    parser.add_argument("--tol", type=float, default=5e-3,
                        help="outer tolerance on nu_t; below ~3e-3 the flow solves fight the "
                             "weakly unsteady trailing edge instead of converging further")
    parser.add_argument("--viscosity-ramp", type=float, nargs="+", default=None,
                        help="laminar-viscosity factors of the start-up continuation "
                             "(default 100 10 1, or just 1 after --restart)")
    parser.add_argument("--continuation", default="ptc", choices=["auto", "newton", "ptc"],
                        help="nonlinear strategy of the flow sub-solves (Newton only limps "
                             "on this case; PTC converges the start-up in a few steps)")
    parser.add_argument("--outdir", default="turb_airfoil_out")
    parser.add_argument("--restart", default=None, metavar="FILE",
                        help="state_a<alpha>.npz of a previous run: skip the viscosity ramp")
    parser.add_argument("--tmr-cp", default=None, metavar="FILE",
                        help="overlay the TMR Cp file n0012cp_cfl3d_sa.dat on the plot")
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
    rans, res = run(args.mesh, args.alpha, args.re, verbose=args.verbose,
                    max_outer=args.max_outer, tol=args.tol, restart=args.restart,
                    viscosity_ramp=(tuple(args.viscosity_ramp)
                                    if args.viscosity_ramp is not None else None),
                    flow_options={"continuation": args.continuation})
    print(f"RANS: {len(rans.history)} outer iterations, converged = {rans.converged}, "
          f"{time.time() - t0:.0f} s")
    suffix = f"a{args.alpha:g}"
    np.savez_compressed(outdir / f"state_{suffix}.npz",
                        U=rans.flow.U, nu_tilde=rans.nu_tilde, alpha0=args.alpha)
    print(f"Saved {outdir / f'state_{suffix}.npz'} (restart with --restart)")
    print(f"alpha = {args.alpha}: CL = {res['CL']:.4f}, CD = {res['CD']:.5f} "
          f"(pressure {res['CD_p']:.5f} + friction {res['CD_f']:.5f})")
    for name, table in (("CFL3D", TMR_CFL3D), ("FUN3D", TMR_FUN3D)):
        if args.alpha in table and args.re == 6e6:
            cl_ref, cd_ref = table[args.alpha]
            dcl = (f"dCL = {res['CL'] - cl_ref:+.4f}" if cl_ref == 0 else
                   f"dCL = {100 * (res['CL'] / cl_ref - 1):+.1f} %")
            print(f"  TMR {name}: CL = {cl_ref:.4f}, CD = {cd_ref:.5f}   "
                  f"{dcl}, dCD = {100 * (res['CD'] / cd_ref - 1):+.1f} %")

    with open(outdir / f"surface_{suffix}.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["x", "y", "Cp", "Cf"])
        writer.writerows(zip(res["x"], res["y"], res["Cp"], res["Cf"]))
    rans.save(outdir / f"turbulent_airfoil_{suffix}.vtu")
    print(f"Saved {outdir / f'surface_{suffix}.csv'} and "
          f"{outdir / f'turbulent_airfoil_{suffix}.vtu'}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        upper = res["y"] > 0.0                       # symmetric airfoil: camber line y = 0
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))
        for sel, label in ((upper, "upper"), (~upper, "lower")):
            ax1.plot(res["x"][sel], res["Cp"][sel], ".", ms=2, label=label)
            ax2.plot(res["x"][sel], res["Cf"][sel], ".", ms=2, label=label)
        if args.tmr_cp:
            xr, cpr = _read_tmr_zone(args.tmr_cp, args.alpha)
            ax1.plot(xr, cpr, "k-", lw=0.8, label="CFL3D (TMR)")
        ax1.invert_yaxis()
        ax1.set_xlabel("x/c")
        ax1.set_ylabel("C_p")
        ax1.legend()
        ax1.grid(alpha=0.3)
        ax2.set_xlabel("x/c")
        ax2.set_ylabel("C_f")
        ax2.legend()
        ax2.grid(alpha=0.3)
        fig.suptitle(f"NACA 0012, SA, Re = {args.re:g}, alpha = {args.alpha:g}")
        fig.savefig(outdir / f"turbulent_airfoil_{suffix}.png", dpi=150, bbox_inches="tight")
        print(f"Saved {outdir / f'turbulent_airfoil_{suffix}.png'}")
    return rans, res


if __name__ == "__main__":
    main()
