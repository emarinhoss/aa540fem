"""NACA 0012 airfoil at 5 degrees angle of attack: impulsive start, lift and drag.

    python examples/airfoil.py [--re 1000] [--dt 0.05] [--t-end 8] [--store-every 10]
                               [--outdir airfoil_out] [--steady/--no-steady] [--no-plot]

Laminar incompressible flow (chord c = 1, freestream U = 1, rho = 1,
mu = 1/Re) on the far-field mesh ``airfoil_naca0012_a5_tri6.msh`` from
``make_meshes.py``, in which the profile is rotated by -5 degrees so the
freestream is along +x and the x/y forces are drag and lift.  Boundary
conditions: uniform velocity on the inlet and the top/bottom far field,
no slip on the airfoil, do-nothing outflow.

The flow is started impulsively from the uniform freestream and integrated
at a fixed time step with Crank-Nicolson after a few backward-Euler start-up
steps (which damp the ringing of an impulsive start); lift and drag coefficients
C = 2 F / (rho U^2 c) are logged at every step (``forces.csv``), the fields
are written every ``store_every`` steps as a ParaView series
(``airfoil.pvd``), and a steady Newton solve started from the final
transient state gives the settled coefficients for comparison.
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.flow import FlowProblem, solve_flow, solve_flow_transient  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
CHORD = 1.0
U_INF = 1.0
RHO = 1.0


def setup(mesh_file, re):
    mesh = read_mesh(mesh_file)
    freestream = (U_INF, 0.0)
    prob = FlowProblem(mesh, mu=RHO * U_INF * CHORD / re, rho=RHO,
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
    parser.add_argument("--mesh", default=str(HERE / "airfoil_naca0012_a5_tri6.msh"))
    parser.add_argument("--re", type=float, default=1000.0, help="chord Reynolds number")
    parser.add_argument("--dt", type=float, default=0.05, help="time step (c / U units)")
    parser.add_argument("--t-end", type=float, default=8.0)
    parser.add_argument("--store-every", type=int, default=10,
                        help="write the fields every n steps")
    parser.add_argument("--startup-steps", type=int, default=4,
                        help="backward-Euler steps before Crank-Nicolson (impulsive start)")
    parser.add_argument("--newton-rtol", type=float, default=1e-6)
    parser.add_argument("--outdir", default="airfoil_out")
    parser.add_argument("--steady", action=argparse.BooleanOptionalAction, default=True,
                        help="finish with a steady Newton solve from the transient state")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prob = setup(args.mesh, args.re)
    mesh = prob.mesh
    print(f"NACA 0012, alpha = 5 deg, Re = {args.re:g}: {mesh.n_nodes} nodes, "
          f"{mesh.n_elems} P2/P1 triangles")

    history = []

    def log(step, t, sol):
        cl, cd = coefficients(sol)
        history.append((t, cl, cd))
        if step % args.store_every == 0 or step == 1:
            print(f"  t = {t:6.2f}   C_L = {cl:8.4f}   C_D = {cd:8.4f}")

    t0 = time.time()
    run = solve_flow_transient(prob, dt=args.dt, t_end=args.t_end, theta=0.5,
                               U0=lambda x, y: (U_INF, 0.0), store_every=args.store_every,
                               callback=log, startup_steps=args.startup_steps,
                               rtol=args.newton_rtol)
    print(f"transient: {run.info['steps']} steps in {time.time() - t0:.0f} s, "
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
        steady = solve_flow(prob, U0=run.final.U)
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
        ax.set_ylim(0, 0.6)                  # the first-step spike of the impulsive start is off-scale
        ax.set_title(f"NACA 0012, alpha = 5 deg, Re = {args.re:g}, impulsive start")
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
