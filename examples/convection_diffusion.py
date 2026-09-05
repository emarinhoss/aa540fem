"""Steady 1-D convection-diffusion boundary layer: Galerkin versus SUPG.

    python examples/convection_diffusion.py [--eps 0.01] [--elems 10]
                                            [--elem-type quad] [--no-plot]

u = (1, 0), T = 0 on the left, T = 1 on the right, insulated top and
bottom.  Exact solution T = (e^{x/eps} - 1) / (e^{1/eps} - 1).  With cell
Peclet number Pe = h / (2 eps) above one the Galerkin solution oscillates;
SUPG is nodally exact on uniform bilinear quadrilaterals.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import DIRICHLET, Problem, solve  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--eps", type=float, default=0.01)
    parser.add_argument("--elems", type=int, default=10)
    parser.add_argument("--elem-type", default="quad")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)

    eps = args.eps
    exact = lambda x, y: np.expm1(x / eps) / np.expm1(1.0 / eps)
    results = {}
    for supg in (False, True):
        p = Problem(a=1.0, b=1.0, elems=args.elems, elem_type=args.elem_type,
                    bc_type={"left": DIRICHLET, "right": DIRICHLET},
                    bc_val={"left": 0.0, "right": 1.0},
                    material=lambda x, y: (eps, 0.0, 0.0, eps, 0.0),
                    velocity=lambda x, y: (1.0, 0.0), supg=supg)
        sol = solve(p)
        results[supg] = sol
        err = np.abs(sol.T - exact(sol.mesh.x, sol.mesh.y)).max()
        label = "SUPG    " if supg else "Galerkin"
        print(f"{label}: min T = {sol.T.min():+.4f}, max nodal error = {err:.3e}")
    print(f"cell Peclet number = {1.0 / args.elems / (2 * eps):.2f}")

    if not args.no_plot:
        import matplotlib.pyplot as plt

        mesh = results[True].mesh
        row = mesh.ny // 2
        x = mesh.X[row]
        fig, ax = plt.subplots()
        xx = np.linspace(0, 1, 400)
        ax.plot(xx, exact(xx, 0), "k-", label="exact")
        ax.plot(x, results[False].grid[row], "o--", label="Galerkin")
        ax.plot(x, results[True].grid[row], "s-", label="SUPG")
        ax.set_xlabel("x")
        ax.set_ylabel("T")
        ax.legend()
        plt.show()
    return results


if __name__ == "__main__":
    main()
