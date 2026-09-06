"""Timing harness for the flow solver: where does the time go?

    python benchmarks/bench_flow.py [--cases steady theta rk45 rans] [--record LABEL]

Runs a few representative cases on the committed meshes, attributes the wall
time to element assembly, matrix factorisation, triangular solves and the
rest (Python overhead, sparse conversions, output), and prints a table.
``--record LABEL`` appends the table to ``docs/performance.md`` so that the
effect of every optimisation phase is on record.  Use the run-configuration
flags (``--threads``, ``--assembly``, ``--direct``) to compare backends.
"""

from __future__ import annotations

import argparse
import cProfile
import pathlib
import platform
import pstats
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow, solve_flow_transient  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
MESHES = ROOT / "examples" / "meshes"

# profile categories: (label, regular expression on "file:function")
CATEGORIES = (
    ("assembly", r"momentum_terms|_linear_matrices|residual_jacobian\b.*spalart|numba_kernels|"
                 r"numpy_kernels|ScatterPlan|scatter"),
    ("factorise", r"gstrf|splu|factorise|Factorisation|_petsc.*factor|setUp"),
    ("solve", r"SuperLU.*solve|_superlu_solve|\.solve\b.*(SuperLU|PETSc|Factor)|KSP.*solve"),
)


def cylinder_problem(umax=0.3, stabilisation=True):
    mesh = read_mesh(MESHES / "cylinder_bl.msh")
    H = 0.41
    inflow = lambda x, y: 4.0 * umax * y * (H - y) / H ** 2
    return FlowProblem(mesh, mu=1e-3, rho=1.0, stabilisation=stabilisation,
                       bc={"inlet": (inflow, 0.0), "walls": (0.0, 0.0),
                           "cylinder": (0.0, 0.0), "outlet": "open"})


def case_steady():
    prob = cylinder_problem(0.3)
    sol = solve_flow(prob)
    return f"{sol.info['iterations']} Newton it."


def case_theta():
    prob = cylinder_problem(1.5)
    run = solve_flow_transient(prob, dt=0.005, t_end=0.1, scheme="theta", startup_steps=2,
                               output_interval=0.1)
    return f"{run.info['steps']} steps"


def case_rk45():
    prob = cylinder_problem(1.5)
    run = solve_flow_transient(prob, dt=0.001, t_end=0.02, scheme="rk45", output_interval=0.02,
                               adaptive=False)
    return f"{run.info['steps']} steps"


def case_rans():
    from aa540fem.turbulence import solve_rans

    mesh = read_mesh(MESHES / "flat_plate_bl.msh")
    prob = FlowProblem(mesh, mu=1e-5, rho=1.0, stabilisation=True,
                       bc={"inlet": (1.0, 0.0), "top": (1.0, 0.0), "symmetry": (None, 0.0),
                           "plate": (0.0, 0.0), "outlet": "open"})
    rans = solve_rans(prob, wall_tags=["plate"], max_outer=2)
    return f"{len(rans.history)} outer it."


CASES = {"steady": case_steady, "theta": case_theta, "rk45": case_rk45, "rans": case_rans}


def attribute(stats: pstats.Stats):
    """Seconds per category from the profile (tottime, so categories are disjoint)."""
    import re

    totals = {label: 0.0 for label, _ in CATEGORIES}
    for (file, _, func), (_, _, tottime, cumtime, _) in stats.stats.items():
        key = f"{file}:{func}"
        for label, pattern in CATEGORIES:
            if re.search(pattern, key):
                # assembly is measured inclusively (its einsums are separate entries)
                totals[label] += cumtime if label == "assembly" and "momentum_terms" in key \
                    else tottime
                break
    return totals


def run_case(name):
    warnings.simplefilter("ignore")
    pr = cProfile.Profile()
    t0 = time.perf_counter()
    pr.enable()
    note = CASES[name]()
    pr.disable()
    wall = time.perf_counter() - t0
    parts = attribute(pstats.Stats(pr))
    parts["other"] = max(0.0, wall - sum(parts.values()))
    return wall, parts, note


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--cases", nargs="+", default=["steady", "theta", "rk45"],
                        choices=sorted(CASES))
    parser.add_argument("--record", default=None, help="label of the row in docs/performance.md")
    try:
        from aa540fem.cli import add_run_arguments, configure_from_args
    except ImportError:                                  # before the run configuration exists
        add_run_arguments = configure_from_args = None
    if add_run_arguments:
        add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args, interactive=False) if configure_from_args else None

    header = f"{'case':8s} {'wall [s]':>9s} {'assembly':>9s} {'factorise':>10s} {'solve':>7s} {'other':>7s}  note"
    lines = [header]
    for name in args.cases:
        wall, parts, note = run_case(name)
        lines.append(f"{name:8s} {wall:9.2f} {parts['assembly']:9.2f} {parts['factorise']:10.2f} "
                     f"{parts['solve']:7.2f} {parts['other']:7.2f}  {note}")
    print("\n".join(lines))

    if args.record:
        doc = ROOT / "docs" / "performance.md"
        machine = f"{platform.processor() or platform.machine()}, {np.__version__} numpy"
        cfg = f", {config.describe()}" if config is not None else ""
        with open(doc, "a") as f:
            f.write(f"\n### {args.record}\n\n{machine}{cfg}\n\n```\n" + "\n".join(lines) + "\n```\n")
        print(f"recorded in {doc}")


if __name__ == "__main__":
    main()
