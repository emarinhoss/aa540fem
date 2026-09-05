"""Vortex shedding behind a cylinder: the Schaefer-Turek 2D-2 benchmark (Re = 100).

    python examples/cylinder_shedding.py [--t-end 6] [--dt 0.01] [--scheme theta|rk45]
                                         [--outdir shedding_out] [--no-plot]

Same channel and cylinder as ``cylinder.py`` on the boundary-layer mesh
``cylinder_bl.msh``, parabolic inflow with U_max = 1.5 (mean 1), nu = 1e-3:
Re = U_mean D / nu = 100.  The steady Re = 20 flow (stabilised) is the
initial state; the time integration (Crank-Nicolson by default) runs to
``t_end`` while the lift and drag history is logged; the Strouhal number is
taken from the lift signal over the last periods.  Reference values:
St = 0.30, C_D,max = 3.23, C_L,max = 1.0 (Schaefer & Turek 1996).
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
from aa540fem.incompressible import FlowProblem, solve_flow, solve_flow_transient  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
MESHES = HERE / "meshes"
H, D, NU = 0.41, 0.1, 1e-3
REFERENCE = {"St": 0.2995, "C_D_max": 3.2298, "C_L_max": 1.0002}


def problem(mesh, umax):
    inflow = lambda x, y: 4.0 * umax * y * (H - y) / H ** 2
    # SUPG only: the grad-div term damps the shedding on this mesh
    return FlowProblem(mesh, mu=NU, rho=1.0, stabilisation=True, grad_div=False,
                       bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0),
                           "cylinder": (0.0, 0.0), "outlet": "open"})


def coefficients(sol, umean):
    fx, fy = sol.forces("cylinder")
    scale = 2.0 / (umean ** 2 * D)
    return scale * fx, scale * fy


def strouhal(t, cl, umean=1.0, last_fraction=0.5):
    """Shedding frequency from the zero crossings of the lift over the last part of the signal."""
    t, cl = np.asarray(t), np.asarray(cl)
    keep = t > t[-1] - last_fraction * (t[-1] - t[0])
    tt, cc = t[keep], cl[keep] - cl[keep].mean()
    up = np.where((cc[:-1] < 0) & (cc[1:] >= 0))[0]
    if up.size < 2:
        return float("nan")
    crossings = tt[up] - cc[up] * (tt[up + 1] - tt[up]) / (cc[up + 1] - cc[up])
    period = np.mean(np.diff(crossings))
    return D / (period * umean)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "cylinder_bl.msh"))
    parser.add_argument("--t-end", type=float, default=6.0)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--scheme", default="theta", choices=["theta", "rk45"])
    parser.add_argument("--output-interval", type=float, default=0.5)
    parser.add_argument("--outdir", default="shedding_out")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mesh = read_mesh(args.mesh)
    print(f"{mesh.n_nodes} nodes, {mesh.n_elems} elements ({', '.join(mesh.cells)})")

    t0 = time.time()
    start = solve_flow(problem(mesh, 0.3))          # Re = 20 steady state as the initial field
    print(f"steady Re = 20 start: {start.info['iterations']} iterations, {time.time() - t0:.0f} s")

    umean = 1.0
    prob = problem(mesh, 1.5)
    history = []
    next_report = [0.0]

    def log(step, t, sol):
        cd, cl = coefficients(sol, umean)
        history.append((t, cd, cl))
        if t >= next_report[0] - 1e-9:
            print(f"  t = {t:7.3f}   C_D = {cd:7.4f}   C_L = {cl:8.4f}")
            next_report[0] = (np.floor(t / args.output_interval + 1e-9) + 1) * args.output_interval

    t0 = time.time()
    if args.scheme == "theta":
        run = solve_flow_transient(prob, dt=args.dt, t_end=args.t_end, theta=0.5, scheme="theta",
                                   U0=start.U,
                                   store_every=max(1, int(round(args.output_interval / args.dt))),
                                   callback=log, startup_steps=2, rtol=1e-6)
    else:
        run = solve_flow_transient(prob, dt=args.dt, t_end=args.t_end, scheme="rk45", U0=start.U,
                                   output_interval=args.output_interval, callback=log)
    print(f"transient: {run.info['steps']} steps in {time.time() - t0:.0f} s")

    t = np.array([h[0] for h in history])
    cd = np.array([h[1] for h in history])
    cl = np.array([h[2] for h in history])
    last = t > t[-1] - 0.5 * (t[-1] - t[0])
    st = strouhal(t, cl, umean)
    results = {"St": st, "C_D_max": cd[last].max(), "C_L_max": cl[last].max()}
    for key, val in results.items():
        print(f"{key:>8} = {val:7.4f}   (reference {REFERENCE[key]:.4f})")
    with open(outdir / "forces.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t", "C_D", "C_L"])
        w.writerows(history)
    pvd = run.save_series(outdir / "shedding")
    print(f"Saved {outdir / 'forces.csv'} and {pvd}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(t, cd, label="$C_D$")
        ax.plot(t, cl, label="$C_L$")
        ax.set_xlabel("t U / D * D")
        ax.set_xlabel("t")
        ax.set_ylabel("coefficient")
        ax.set_title(f"cylinder in a channel, Re = 100: St = {st:.3f} (ref. 0.30)")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.savefig(outdir / "forces.png", dpi=150, bbox_inches="tight")
        print(f"Saved {outdir / 'forces.png'}")
    return run, results


if __name__ == "__main__":
    main()
