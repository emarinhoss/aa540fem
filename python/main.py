"""Driver script.  Port of ``main.m``: edit the user inputs below and run

    python main.py [--output temperature.png] [--show]

The only files that should need editing are this one and
``aa540fem/conductivity_and_forcing.py``.
"""

from __future__ import annotations

import argparse

from aa540fem import Problem, solve

# ---------------------------------------------------------------- User inputs
a = 4.0            # horizontal length of the domain
b = 6.0            # vertical length of the domain

elems = 50         # cells along x and y (must be equal)

elem_type = 3      # 1 --> 3-node linear
                   # 2 --> 4-node linear
                   # 3 --> 9-node quadratic

# Boundary condition types: 0 --> Dirichlet, 1 --> Neumann
top_bc_type = 0
left_bc_type = 1
right_bc_type = 1
bottom_bc_type = 0

# Values at the boundaries (temperature for Dirichlet, normal flux for Neumann)
top_bc_val = 100.0
left_bc_val = 0.0
right_bc_val = 0.0
bottom_bc_val = 0.0

# Order of the quadrature.  Quadrilaterals: Gauss-Legendre points per
# direction (1-4 in the original).  Triangles: 1, 3, 4, 6, 7, 9, 12 or 13
# Gauss points on the triangle.  None picks the lowest order that fully
# integrates the chosen element (1, 2, 3 for element types 1, 2, 3); the
# MATLAB default of 1 leaves 4- and 9-node elements rank deficient.
order = None
# -----------------------------------------------------------------------------


def build_problem() -> Problem:
    return Problem(
        a=a, b=b, elems=elems, elem_type=elem_type, order=order,
        bc_type={"top": top_bc_type, "right": right_bc_type,
                 "left": left_bc_type, "bottom": bottom_bc_type},
        bc_val={"top": top_bc_val, "right": right_bc_val,
                "left": left_bc_val, "bottom": bottom_bc_val},
    )


def plot(solution, output=None, show=False, levels=20):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mesh = solution.mesh
    fig, ax = plt.subplots()
    cs = ax.contour(mesh.X, mesh.Y, solution.grid, levels)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Temperature")
    fig.colorbar(cs, ax=ax)
    if output:
        fig.savefig(output, dpi=150, bbox_inches="tight")
        print(f"Saved {output}")
    if show:
        plt.show()
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", "-o", default="temperature.png",
                        help="file for the contour plot ('' to skip)")
    parser.add_argument("--show", action="store_true", help="open an interactive window")
    args = parser.parse_args(argv)

    solution = solve(build_problem(), verbose=True)
    print(f"T in [{solution.T.min():.6g}, {solution.T.max():.6g}]")
    if args.output or args.show:
        plot(solution, args.output or None, args.show)
    return solution


if __name__ == "__main__":
    main()
