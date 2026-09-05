"""Convergence studies.

    python scripts/convergence.py [--elems 4 8 16 32]      # spatial
    python scripts/convergence.py --transient              # temporal

Spatial: solves div(grad T) + f = 0 on [0, 2] x [0, 3] with homogeneous
Dirichlet conditions and f chosen so that T = sin(pi x / a) sin(pi y / b),
and prints the L2 and H1 errors and the observed convergence orders for
every element type.  Expected orders: L2 = p + 1 and H1 = p for degree-p
elements.

Temporal: integrates the decaying mode T = exp(-pi^2 t) sin(pi x) on
[0, 1] x [0, 0.5] with backward Euler and Crank-Nicolson and reports the
error against a fine-step reference for a sequence of time steps.
Expected orders: 1 (backward Euler) and 2 (Crank-Nicolson).
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from aa540fem import DIRICHLET, ELEMENTS, Problem, error_norms, solve, solve_transient  # noqa: E402

A, B = 2.0, 3.0


def exact(x, y):
    return np.sin(np.pi * x / A) * np.sin(np.pi * y / B)


def exact_grad(x, y):
    return (np.pi / A * np.cos(np.pi * x / A) * np.sin(np.pi * y / B),
            np.pi / B * np.sin(np.pi * x / A) * np.cos(np.pi * y / B))


def material(x, y):
    return 1.0, 0.0, 0.0, 1.0, np.pi ** 2 * (1 / A ** 2 + 1 / B ** 2) * exact(x, y)


def run(name, elems_list, method="direct"):
    rows = []
    for elems in elems_list:
        p = Problem(a=A, b=B, elems=elems, elem_type=name, material=material,
                    bc_type={s: DIRICHLET for s in ("top", "right", "left", "bottom")},
                    bc_val={s: 0.0 for s in ("top", "right", "left", "bottom")})
        sol = solve(p, method=method)
        e = error_norms(sol.mesh, sol.T, exact, exact_grad)
        rows.append((elems, sol.mesh.n_nodes, e["L2"], e["H1"]))
    return rows


def run_transient(dts, elems=8, elem_type="quad9", t_end=0.1):
    problem = Problem(a=1.0, b=0.5, elems=elems, elem_type=elem_type,
                      bc_type={"left": DIRICHLET, "right": DIRICHLET},
                      bc_val={"left": 0.0, "right": 0.0},
                      material=lambda x, y: (1.0, 0.0, 0.0, 1.0, 0.0))
    T0 = lambda x, y: np.sin(np.pi * x)
    ref = solve_transient(problem, dt=min(dts) / 50, t_end=t_end, theta=0.5, T0=T0,
                          store_every=10 ** 9).T
    print(f"\ndecaying mode, {elem_type} elements, {elems} cells, t_end = {t_end}")
    print(f"{'dt':>8} {'BE error':>11} {'rate':>6} {'CN error':>11} {'rate':>6}")
    prev = None
    for dt in dts:
        errs = [np.abs(solve_transient(problem, dt=dt, t_end=t_end, theta=theta, T0=T0,
                                       store_every=10 ** 9).T - ref).max()
                for theta in (1.0, 0.5)]
        if prev is None:
            print(f"{dt:8.4g} {errs[0]:11.3e} {'':>6} {errs[1]:11.3e}")
        else:
            r = np.log2(np.array(prev) / np.array(errs))
            print(f"{dt:8.4g} {errs[0]:11.3e} {r[0]:6.2f} {errs[1]:11.3e} {r[1]:6.2f}")
        prev = errs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--elems", type=int, nargs="+", default=[4, 8, 16, 32])
    parser.add_argument("--method", default="direct", choices=["direct", "cg"])
    parser.add_argument("--elements", nargs="+", default=sorted(ELEMENTS))
    parser.add_argument("--transient", action="store_true",
                        help="temporal convergence of backward Euler and Crank-Nicolson")
    parser.add_argument("--dts", type=float, nargs="+", default=[0.02, 0.01, 0.005, 0.0025])
    args = parser.parse_args(argv)

    if args.transient:
        run_transient(args.dts)
        return

    for name in args.elements:
        rows = run(name, args.elems, args.method)
        print(f"\n{name}")
        print(f"{'elems':>6} {'nodes':>7} {'L2':>11} {'rate':>6} {'H1':>11} {'rate':>6}")
        for i, (elems, nodes, l2, h1) in enumerate(rows):
            if i == 0:
                print(f"{elems:6d} {nodes:7d} {l2:11.3e} {'':>6} {h1:11.3e}")
            else:
                r2 = np.log2(rows[i - 1][2] / l2)
                r1 = np.log2(rows[i - 1][3] / h1)
                print(f"{elems:6d} {nodes:7d} {l2:11.3e} {r2:6.2f} {h1:11.3e} {r1:6.2f}")


if __name__ == "__main__":
    main()
