"""Timing of the distributed (domain-decomposed) flow solver.

    mpirun -n 4 python benchmarks/bench_mpi.py [--threads 1] [--cases steady theta]

Compares, on the cylinder mesh at Re 20 (steady) and Re 100 (20 theta
steps), the replicated path (every rank assembles the whole system, MUMPS
factorises it across the ranks: what ``solve_flow`` does under ``mpirun``)
with the distributed Newton of ``aa540fem.parallel.flow`` (every rank
assembles its own elements; distributed MUMPS or the parallel
fieldsplit/LSC solver), and checks that the states agree.  Rank 0 prints
the table; ``--record LABEL`` appends it to ``docs/performance.md``.  The
single-process numbers are in ``bench_flow.py``.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow, solve_flow_transient  # noqa: E402
from aa540fem.parallel.comm import world  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
MESHES = ROOT / "examples" / "meshes"


def cylinder(umax):
    mesh = read_mesh(MESHES / "cylinder_bl.msh")
    inflow = lambda x, y: 4.0 * umax * y * (0.41 - y) / 0.41 ** 2
    return FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=True,
                       bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0), "cylinder": (0.0, 0.0),
                           "outlet": "open"})


def timed(comm, fn):
    comm.Barrier()
    t0 = time.perf_counter()
    out = fn()
    comm.Barrier()
    return time.perf_counter() - t0, out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--cases", nargs="+", default=["steady", "theta"],
                        choices=["steady", "theta"])
    parser.add_argument("--record", default=None)
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args, interactive=False)
    warnings.simplefilter("ignore")
    comm = world()
    rank, size = comm.Get_rank(), comm.Get_size()
    from aa540fem.parallel.flow import DistributedFlowSystem

    lines = [f"{'case':30s} {'ranks':>5s} {'wall [s]':>9s}  note"]
    if "steady" in args.cases:
        prob = cylinder(0.3)
        t_rep, ref = timed(comm, lambda: solve_flow(prob))
        ref_U = ref.U
        lines.append(f"{'steady replicated':30s} {size:5d} {t_rep:9.2f}  "
                     f"{ref.info['iterations']} Newton it.")
        for method in ("direct", "fieldsplit"):
            t_dist, sol = timed(comm, lambda: DistributedFlowSystem(prob)
                                .solve_steady(method=method))
            err = np.abs(sol.U - ref_U).max() / np.abs(ref_U).max()
            note = f"{sol.info['iterations']} Newton it., |dU|/|U| = {err:.0e}"
            if method == "fieldsplit":
                note += f", {max(sol.info['linear_iterations'])} FGMRES it. max"
            lines.append(f"{'steady distributed ' + method:30s} {size:5d} {t_dist:9.2f}  {note}")
    if "theta" in args.cases:
        prob = cylinder(1.5)
        U0 = solve_flow(cylinder(0.3)).U
        t_rep, ref = timed(comm, lambda: solve_flow_transient(
            prob, dt=0.005, t_end=0.1, scheme="theta", U0=U0, startup_steps=2, rtol=1e-6))
        ref_U = ref.snapshots[-1]
        lines.append(f"{'theta 20 steps replicated':30s} {size:5d} {t_rep:9.2f}  "
                     f"{ref.info['factorisations']} factorisations")
        for method in ("direct", "fieldsplit"):
            t_dist, run = timed(comm, lambda: DistributedFlowSystem(prob).solve_transient(
                0.005, 0.1, U0=U0, startup_steps=2, rtol=1e-6, method=method))
            err = np.abs(run.snapshots[-1] - ref_U).max() / np.abs(ref_U).max()
            lines.append(f"{'theta distributed ' + method:30s} {size:5d} {t_dist:9.2f}  "
                         f"{run.info['factorisations']} factorisations, |dU|/|U| = {err:.0e}")
    if rank == 0:
        print("\n".join(lines))
        if args.record:
            doc = ROOT / "docs" / "performance.md"
            with open(doc, "a") as f:
                f.write(f"\n### {args.record}\n\n{config.describe()}\n\n```\n"
                        + "\n".join(lines) + "\n```\n")
            print(f"recorded in {doc}")


if __name__ == "__main__":
    main()
