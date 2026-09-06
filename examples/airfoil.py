"""NACA 0012 airfoil at 5 degrees angle of attack: impulsive start, lift and drag.

    python examples/airfoil.py [--re 1000] [--scheme rk45|theta] [--dt 0.01] [--t-end 8]
                               [--output-interval 0.5] [--outdir airfoil_out]
                               [--steady/--no-steady] [--no-plot]

Laminar incompressible flow (chord c = 1, freestream U = 1, rho = 1,
mu = 1/Re) on the far-field mesh ``airfoil_naca0012_a5_tri6.msh`` from
``make_meshes.py``, in which the profile is rotated by -5 degrees so the
freestream is along +x and the x/y forces are drag and lift.  Boundary
conditions: uniform velocity on the inlet and the top/bottom far field,
no slip on the airfoil, do-nothing outflow.

The flow is started impulsively from the uniform freestream (projected
onto the divergence-free space) and integrated with the adaptive
Runge-Kutta 45 scheme (Dormand-Prince, the default) or, with
``--scheme theta``, with Crank-Nicolson at a fixed step after a few
backward-Euler start-up steps.  Lift and drag coefficients
C = 2 F / (rho U^2 c) are logged at every accepted step (``forces.csv``),
the fields are written at the fixed output times (multiples of
``--output-interval``) as a ParaView series (``airfoil.pvd``), and a steady
Newton solve started from the final transient state gives the settled
coefficients for comparison.
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
from aa540fem.incompressible import FlowProblem, solve_flow, solve_flow_transient  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MESHES = HERE / "meshes"
CHORD = 1.0
U_INF = 1.0
RHO = 1.0


def setup(mesh_file, re, stabilise=False):
    mesh = read_mesh(mesh_file)
    freestream = (U_INF, 0.0)
    prob = FlowProblem(mesh, mu=RHO * U_INF * CHORD / re, rho=RHO, stabilisation=stabilise,
                       bc={"inlet": freestream, "farfield": freestream,
                           "airfoil": (0.0, 0.0), "outlet": "open"})
    return prob


def coefficients(sol):
    """(C_L, C_D) from the traction integral over the airfoil."""
    fx, fy = sol.forces("airfoil")
    scale = 2.0 / (RHO * U_INF ** 2 * CHORD)
    return scale * fy, scale * fx


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "airfoil_naca0012_a5_tri6.msh"))
    parser.add_argument("--re", type=float, default=1000.0, help="chord Reynolds number")
    parser.add_argument("--scheme", default="rk45", choices=["rk45", "theta"])
    parser.add_argument("--dt", type=float, default=0.01,
                        help="initial (rk45) or fixed (theta) time step in c / U units")
    parser.add_argument("--t-end", type=float, default=8.0)
    parser.add_argument("--output-interval", type=float, default=0.5,
                        help="write the fields at multiples of this time")
    parser.add_argument("--rtol", type=float, default=1e-4, help="rk45 relative tolerance")
    parser.add_argument("--atol", type=float, default=1e-6, help="rk45 absolute tolerance")
    parser.add_argument("--stabilise", action="store_true",
                        help="SUPG/grad-div stabilisation (needed above Re ~ 1000 on this mesh)")
    parser.add_argument("--continuation", default="auto", choices=["auto", "newton", "ptc"],
                        help="steady solve: Newton, pseudo-transient continuation, or auto")
    parser.add_argument("--startup-steps", type=int, default=4,
                        help="theta: backward-Euler steps before Crank-Nicolson")
    parser.add_argument("--newton-rtol", type=float, default=1e-6, help="theta: Newton tolerance")
    parser.add_argument("--outdir", default="airfoil_out")
    parser.add_argument("--steady", action=argparse.BooleanOptionalAction, default=True,
                        help="finish with a steady Newton solve from the transient state")
    parser.add_argument("--no-plot", action="store_true")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prob = setup(args.mesh, args.re, args.stabilise)
    mesh = prob.mesh
    print(f"NACA 0012, alpha = 5 deg, Re = {args.re:g}: {mesh.n_nodes} nodes, "
          f"{mesh.n_elems} P2/P1 triangles")

    history = []
    next_report = [0.0]

    def log(step, t, sol):
        cl, cd = coefficients(sol)
        history.append((t, cl, cd))
        if t >= next_report[0] - 1e-9 or step == 1:
            print(f"  t = {t:7.3f}   C_L = {cl:8.4f}   C_D = {cd:8.4f}   "
                  f"dt = {sol.info.get('dt', args.dt):.2e}")
            next_report[0] = (np.floor(t / args.output_interval + 1e-9) + 1) * args.output_interval

    t0 = time.time()
    if args.scheme == "rk45":
        run = solve_flow_transient(prob, dt=args.dt, t_end=args.t_end, scheme="rk45",
                                   U0=lambda x, y: (U_INF, 0.0),
                                   output_interval=args.output_interval, callback=log,
                                   rtol=args.rtol, atol=args.atol)
        print(f"transient (RK45): {run.info['steps']} accepted steps, "
              f"{run.info['rejected']} rejected, dt in [{run.info['dt_min_used']:.2e}, "
              f"{run.info['dt_max_used']:.2e}], {time.time() - t0:.0f} s")
    else:
        store_every = max(1, int(round(args.output_interval / args.dt)))
        run = solve_flow_transient(prob, dt=args.dt, t_end=args.t_end, theta=0.5, scheme="theta",
                                   U0=lambda x, y: (U_INF, 0.0), store_every=store_every,
                                   callback=log, startup_steps=args.startup_steps,
                                   rtol=args.newton_rtol)
        print(f"transient (theta): {run.info['steps']} steps in {time.time() - t0:.0f} s, "
              f"Newton iterations per step: max {max(run.info['newton_iterations'])}")

    with open(outdir / "forces.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["t", "C_L", "C_D"])
        writer.writerows(history)
    pvd = run.save_series(outdir / "airfoil")
    print(f"Saved {outdir / 'forces.csv'} and {pvd} ({len(run.times)} time levels)")

    t_hist = np.array([h[0] for h in history])
    cl_hist = np.array([h[1] for h in history])
    cd_hist = np.array([h[2] for h in history])
    print(f"transient at t = {t_hist[-1]:g}:  C_L = {cl_hist[-1]:.4f}   C_D = {cd_hist[-1]:.4f}")

    steady = None
    if args.steady:
        t0 = time.time()
        steady = solve_flow(prob, U0=run.final.U, continuation=args.continuation)
        cl, cd = coefficients(steady)
        print(f"steady Newton from the final state: {steady.info['iterations']} iterations "
              f"in {time.time() - t0:.0f} s, converged = {steady.info['converged']}")
        print(f"steady:  C_L = {cl:.4f}   C_D = {cd:.4f}")
        steady.save(outdir / "airfoil_steady.vtu")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(t_hist, cl_hist, label="$C_L$")
        ax.plot(t_hist, cd_hist, label="$C_D$")
        if steady is not None:
            ax.axhline(cl, color="C0", ls=":", lw=1)
            ax.axhline(cd, color="C1", ls=":", lw=1)
        ax.set_xlabel("t U / c")
        ax.set_ylabel("coefficient")
        ax.set_ylim(0, 0.6)          # first-step spike of the impulsive start is off-scale
        ax.set_title(f"NACA 0012, alpha = 5 deg, Re = {args.re:g}, impulsive start ({args.scheme})")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.savefig(outdir / "forces.png", dpi=150, bbox_inches="tight")

        final = steady if steady is not None else run.final
        tri = mtri.Triangulation(mesh.x, mesh.y, mesh.triangulation())
        fig, ax = plt.subplots(figsize=(9, 4.5))
        cs = ax.tricontourf(tri, final.speed, np.linspace(0, 1.5, 31), extend="max")
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-0.75, 0.75)
        ax.set_aspect("equal")
        fig.colorbar(cs, ax=ax, label="|u| / U")
        ax.set_title("speed" + (" (steady)" if steady is not None else f" (t = {t_hist[-1]:g})"))
        fig.savefig(outdir / "speed.png", dpi=150, bbox_inches="tight")
        print(f"Saved {outdir / 'forces.png'} and {outdir / 'speed.png'}")
    return run, steady, history


if __name__ == "__main__":
    main()
