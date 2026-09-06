"""Steady conduction with temperature-dependent conductivity, solved by Newton.

    python examples/nonlinear_conduction.py [--beta 1.0] [--elems 8]
                                            [--elem-type quad9] [--picard] [--no-plot]

kappa(T) = 1 + beta T on the unit square with T = 0 on the left and T = 1
on the right (insulated top and bottom) has the exact 1-D solution
T = (sqrt(1 + beta (2 + beta) x) - 1) / beta (Kirchhoff transform).
The script prints the Newton (or Picard) residual history and the error.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import DIRICHLET, Problem, error_norms, solve  # noqa: E402
from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--elems", type=int, default=8)
    parser.add_argument("--elem-type", default="quad9")
    parser.add_argument("--picard", action="store_true", help="fixed-point instead of Newton")
    parser.add_argument("--no-plot", action="store_true")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    print(f"run configuration: {config.describe()}")

    beta = args.beta
    exact = lambda x, y: (np.sqrt(1.0 + beta * (2.0 + beta) * x) - 1.0) / beta

    def material(x, y, T):
        k = 1.0 + beta * T
        return k, 0.0, 0.0, k, 0.0

    p = Problem(a=1.0, b=1.0, elems=args.elems, elem_type=args.elem_type,
                bc_type={"left": DIRICHLET, "right": DIRICHLET},
                bc_val={"left": 0.0, "right": 1.0}, material=material)
    sol = solve(p, verbose=True, newton=not args.picard, max_newton=200)
    errs = error_norms(sol.mesh, sol.T, exact)
    print(f"{sol.info['nonlinear']}: {sol.info['iterations']} iterations, "
          f"converged = {sol.info['converged']}, L2 error = {errs['L2']:.3e}")

    if not args.no_plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
        row = sol.mesh.ny // 2
        xx = np.linspace(0, 1, 200)
        ax1.plot(xx, exact(xx, 0), "k-", label="exact")
        ax1.plot(sol.mesh.X[row], sol.grid[row], "o", label="FEM")
        ax1.set_xlabel("x")
        ax1.set_ylabel("T")
        ax1.legend()
        ax2.semilogy(sol.info["residuals"], "s-")
        ax2.set_xlabel("iteration")
        ax2.set_ylabel("|R|")
        ax2.set_title(sol.info["nonlinear"])
        plt.show()
    return sol


if __name__ == "__main__":
    main()
