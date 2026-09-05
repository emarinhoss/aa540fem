"""Mesh-convergence study on a manufactured solution.

    python scripts/convergence.py [--elems 4 8 16 32]

Solves div(grad T) + f = 0 on [0, 2] x [0, 3] with homogeneous Dirichlet
conditions and f chosen so that T = sin(pi x / a) sin(pi y / b), and prints
the L2 and H1 errors and the observed convergence orders for every element
type.  Expected orders: L2 = p + 1 and H1 = p for degree-p elements.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from aa540fem import DIRICHLET, ELEMENTS, Problem, error_norms, solve  # noqa: E402

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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--elems", type=int, nargs="+", default=[4, 8, 16, 32])
    parser.add_argument("--method", default="direct", choices=["direct", "cg"])
    parser.add_argument("--elements", nargs="+", default=sorted(ELEMENTS))
    args = parser.parse_args(argv)

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
