"""Schaefer-Turek 3D-1Z benchmark: steady flow around a cylinder in a square channel.

    python examples/cylinder3d.py [--mesh cylinder3d_tet10.msh] [--linear petsc] [--no-prompt]

Channel ``[0, 2.5] x [0, 0.41] x [0, 0.41]``, cylinder of diameter D = 0.1
along z centred at (0.5, 0.2), inflow ``U(0, y, z) = 16 U_m y z (H - y)(H - z) / H^4``
with U_m = 0.45, nu = 1e-3, rho = 1, so the Reynolds number built on the mean
inflow velocity ``U_bar = 4 U_m / 9 = 0.2`` is 20.  Reference values (Schaefer
and Turek 1996, refined by Bayraktar, Mierka and Turek 2012):
``C_D = 6.185``, ``C_L = 0.0094``, ``dp = p(0.45, 0.2, 0.205) - p(0.55, 0.2, 0.205)
= 0.1710``; the benchmark's admissible intervals are ``[6.05, 6.25]``,
``[0.008, 0.01]`` and ``[0.165, 0.175]``.  The coefficients are
``C = 2 F / (rho U_bar^2 D H)``.  The meshes come from
``examples/meshes/make_meshes.py::make_cylinder3d``.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from aa540fem import read_mesh  # noqa: E402
from aa540fem.cli import add_run_arguments, configure_from_args  # noqa: E402
from aa540fem.incompressible import FlowProblem, solve_flow  # noqa: E402

MESHES = pathlib.Path(__file__).resolve().parent / "meshes"
H, D, NU, UM = 0.41, 0.1, 1e-3, 0.45
UBAR = 4.0 * UM / 9.0
REFERENCE = {"C_D": 6.185, "C_L": 0.0094, "dp": 0.1710}
INTERVALS = {"C_D": (6.05, 6.25), "C_L": (0.008, 0.01), "dp": (0.165, 0.175)}
PROBES = np.array([[0.45, 0.2, 0.205], [0.55, 0.2, 0.205]])


def problem(mesh, stabilisation=True):
    inflow = lambda x, y, z: 16.0 * UM * y * z * (H - y) * (H - z) / H ** 4
    return FlowProblem(mesh, mu=NU, rho=1.0, stabilisation=stabilisation,
                       bc={"inlet": (inflow, 0.0, 0.0), "walls": (0.0, 0.0, 0.0),
                           "cylinder": (0.0, 0.0, 0.0), "outlet": "open"})


def extruded_mesh(path, layers=5, grading=1.5):
    """The 3D-1Z channel from a 2-D cylinder mesh (boundary-layer quads around the
    cylinder, triangles elsewhere) extruded over the height ``H`` with layers graded
    towards the two end walls; the end planes join the ``walls`` tag."""
    m3 = read_mesh(path).extrude(H, layers, tags=("front", "back"), grading=grading)
    walls = {"quad9": [m3.boundary["walls"]], "triangle6": []}
    for tag in ("front", "back"):
        for block in m3.face_blocks(tag):
            walls["triangle6" if block.shape[1] == 6 else "quad9"].append(block)
        del m3.boundary[tag]
    m3.boundary["walls"] = {k: np.vstack(v) for k, v in walls.items() if v}
    return m3


def coefficients(sol):
    fx, fy, _ = sol.forces("cylinder")
    scale = 2.0 / (UBAR ** 2 * D * H)
    p = sol.pressure_at(PROBES)
    return {"C_D": fx * scale, "C_L": fy * scale, "dp": float(p[0] - p[1])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--mesh", default=str(MESHES / "cylinder3d_tet10.msh"))
    parser.add_argument("--extrude", metavar="MESH2D", default=None,
                        help="build the channel by extruding a 2-D cylinder mesh instead")
    parser.add_argument("--layers", type=int, default=5, help="cells across the height")
    parser.add_argument("--grading", type=float, default=1.5,
                        help="growth of the layer thickness away from the end walls")
    parser.add_argument("--outdir", default="cylinder3d_out")
    parser.add_argument("--no-save", action="store_true")
    parser.add_argument("--distributed", action="store_true",
                        help="use the MPI domain decomposition (run under mpirun)")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    config = configure_from_args(args)
    from aa540fem.parallel.comm import rank

    if rank() == 0:
        print(f"run configuration: {config.describe()}")
    mesh = (extruded_mesh(args.extrude, args.layers, args.grading) if args.extrude
            else read_mesh(args.mesh))
    prob = problem(mesh)
    if rank() == 0:
        print(f"{mesh.n_nodes} nodes, {mesh.n_elems} elements ({', '.join(mesh.cells)})")
    t0 = time.time()
    if args.distributed:
        from aa540fem.parallel.flow import DistributedFlowSystem

        method = "fieldsplit" if config.linear_backend.startswith("petsc") else "direct"
        sol = DistributedFlowSystem(prob).solve_steady(method=method, verbose=rank() == 0)
    else:
        sol = solve_flow(prob, verbose=True)
    wall = time.time() - t0
    if rank() != 0:
        return sol
    print(f"steady solve: {sol.info['iterations']} Newton iterations, {wall:.0f} s"
          + (f", {max(sol.info['linear_iterations'])} FGMRES it. max"
             if sol.info.get("linear_iterations") else ""))
    coef = coefficients(sol)
    for key, val in coef.items():
        lo, hi = INTERVALS[key]
        ok = "ok" if lo <= val <= hi else "OUTSIDE"
        print(f"{key:>4} = {val:8.4f}   (reference {REFERENCE[key]:.4f}, "
              f"interval [{lo}, {hi}] {ok})")
    if not args.no_save:
        outdir = pathlib.Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        print("Saved", sol.save(outdir / "cylinder3d.vtu"))
    return sol


if __name__ == "__main__":
    main()
