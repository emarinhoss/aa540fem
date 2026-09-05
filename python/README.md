# Python port of the AA 540 anisotropic heat conduction FEM code

This directory is a NumPy/SciPy port of the MATLAB code in `../project`,
extended into a small but general 2-D finite element framework for the
scalar transport equation

    rho_c dT/dt + u . grad T - div( kappa(x, y) . grad T ) = f(x, y)

Without a velocity and time derivative this is the original
`div(kappa grad T) + f = 0`.  `kappa` is a full 2x2 conductivity tensor,
`u` an optional velocity field (with SUPG stabilisation), and the
boundaries carry Dirichlet (`T = T0`) or Neumann (`n . kappa grad T = q_n`)
conditions by name.  Elements are 3- and 6-node triangles or 4- and 9-node
quadrilaterals, on the built-in structured rectangle or on any mesh read
through meshio (Gmsh `.msh` files in particular).  Time integration uses
the theta-method (backward Euler or Crank-Nicolson).

## Installation

```bash
pip install -e ".[mesh,amg,dev]"     # from this directory
# or just the runtime pieces:
pip install -r requirements.txt
```

Extras: `mesh` (meshio, for reading mesh files and writing ParaView output),
`amg` (pyamg, algebraic multigrid preconditioner for the CG solver),
`meshgen` (gmsh, only needed to regenerate `examples/*.msh`), `dev`
(pytest, ruff).

## Running

```bash
python main.py                          # rectangle from main.m, writes temperature.png
python main.py --show                   # interactive window
python main.py --method cg --vtk out.vtu
python examples/annulus.py              # curved annulus mesh, exact solution ln(r)/ln(2)
python examples/convection_diffusion.py # Galerkin vs SUPG on a boundary layer
python examples/rotating_hill.py        # transient advection, writes a ParaView series
python scripts/convergence.py           # mesh-convergence study, prints observed orders
python scripts/convergence.py --transient   # temporal orders of backward Euler / Crank-Nicolson
python -m pytest                        # verification suite
ruff check .                            # lint
```

As in the MATLAB version, the files meant to be edited are `main.py`
(domain, mesh, element type, boundary conditions, quadrature order) and
`aa540fem/conductivity_and_forcing.py` (conductivity tensor and heat
source).

### From Python

```python
from aa540fem import Problem, solve, read_mesh

# Structured rectangle, same inputs as main.m
problem = Problem(a=4, b=6, elems=50, elem_type=3,
                  bc_type={"top": 0, "right": 1, "left": 1, "bottom": 0},
                  bc_val={"top": 100, "right": 0, "left": 0, "bottom": 0})
solution = solve(problem)
solution.T            # nodal temperatures
solution.grid         # temperatures on the (ny, nx) grid of solution.mesh.X / .Y

# Unstructured mesh with Gmsh physical groups "inner" and "outer"
mesh = read_mesh("examples/annulus_tri6.msh")
problem = Problem(mesh=mesh, bc_type={"inner": 0, "outer": 0},
                  bc_val={"inner": 0.0, "outer": 1.0})
solution = solve(problem, method="cg")
solution.flux()       # element-centroid gradients and heat fluxes
solution.save("annulus.vtu")   # open in ParaView
```

Boundary values and the optional `material(x, y)` callable may be functions
of position.  Boundary tags that are not listed in `bc_type` get the natural
zero-flux condition.

### Convection and time stepping

```python
import numpy as np
from aa540fem import Problem, solve, solve_transient

# Steady boundary layer with u = (1, 0): SUPG keeps the solution monotone
problem = Problem(a=1, b=1, elems=10, elem_type="quad",
                  bc_type={"left": 0, "right": 0}, bc_val={"left": 0.0, "right": 1.0},
                  material=lambda x, y: (0.01, 0, 0, 0.01, 0.0),   # kxx, kxy, kyx, kyy, f
                  velocity=lambda x, y: (1.0, 0.0), supg=True)
steady = solve(problem, method="gmres")          # non-symmetric: gmres, not cg

# Transient: rho_c dT/dt + ..., Crank-Nicolson from T0 to t = 1
result = solve_transient(problem, dt=0.01, t_end=1.0, theta=0.5,
                         T0=lambda x, y: np.zeros_like(x), store_every=10)
result.T                    # final field
result.snapshot(0.5)        # stored field closest to t = 0.5
result.save_series("out/run")   # out/run_0000.vtu ... + out/run.pvd for ParaView
```

`material`, `velocity`, `rho_c` and the boundary values may take a third
argument `t` (`lambda x, y, t: ...`); the solver then re-evaluates them each
step.  With time-independent coefficients the system matrix is factorised
once.  The SUPG parameter is the classic `tau = h/(2|u|) (coth Pe - 1/Pe)`
for steady problems and the transient form
`tau = [(2/dt)^2 + (2|u|/h)^2 + (4 kappa/h^2)^2]^(-1/2)` in time stepping;
the diffusion term of the residual is omitted, which is exact for linear
elements.  On uniform bilinear quadrilaterals SUPG is nodally exact for the
1-D boundary layer (`examples/convection_diffusion.py`).

### Meshing with Gmsh

Tag the boundary curves with `Physical Curve("name")` (see
`examples/annulus.geo`); those names become the keys of `bc_type` and
`bc_val`.  Generate with `gmsh -2 file.geo` or through the Python API as in
`examples/make_meshes.py`, which also shows how to get quadrilaterals
(`Mesh.RecombineAll`) and curved second-order elements
(`Mesh.ElementOrder = 2`).  Elements with negative Jacobians are flipped on
read and unused nodes are dropped.

## Layout

| Module                              | Contents                                              | MATLAB origin |
|-------------------------------------|-------------------------------------------------------|---------------|
| `aa540fem/elements.py`              | `ReferenceElement` registry: nodes, faces, quadrature | (new)         |
| `aa540fem/shape_functions.py`       | Lagrange shape functions, Gmsh node ordering          | `interpfunc_*.m` |
| `aa540fem/quadrature.py`            | Gauss-Legendre and triangle rules                     | `gauss_legendre_quad.m`, `gauss_trgl.m` |
| `aa540fem/geometry.py`              | `Mesh` data structure, structured rectangle           | `geometry.m`  |
| `aa540fem/mesh_io.py`               | meshio import/export, VTK output                      | (new)         |
| `aa540fem/element.py`               | Jacobians, element diffusion/convection/mass, SUPG    | `elem_eqn.m`  |
| `aa540fem/boundary.py`              | Dirichlet elimination, Neumann edge integrals         | `dirichlet.m` |
| `aa540fem/solver.py`                | `Problem`, operator assembly, direct / CG / GMRES     | `main.m`      |
| `aa540fem/transient.py`             | theta-method time stepping, ParaView series output    | (new)         |
| `aa540fem/util.py`                  | time-argument detection for user callables            | (new)         |
| `aa540fem/postprocess.py`           | Centroid gradients and fluxes, L2/H1 error norms      | (new)         |
| `aa540fem/conductivity_and_forcing.py` | User material and source                           | `conductivity_and_forcing.m` |
| `main.py`                           | Driver with the user inputs                           | `main.m`      |
| `examples/`, `scripts/`, `tests/`   | Annulus case, convergence study, verification suite   |               |

## Verification

`tests/` checks the shape functions (nodal property, partition of unity,
finite-difference derivatives), quadrature tables, mesh orientation and
boundary detection, exact reproduction of linear and quadratic fields with
anisotropic conductivity and mixed boundary conditions, the Gmsh annulus
with curved elements, CG/GMRES against the direct solver, VTK round trips,
mass and convection operators, SUPG (nodally exact boundary layer, bounded
skew advection), a decaying mode and a transient manufactured solution with
convection and time-dependent Dirichlet data.
`scripts/convergence.py` on the manufactured solution
`T = sin(pi x/a) sin(pi y/b)` gives the expected orders:

| element   | L2 order | H1 order |
|-----------|----------|----------|
| triangle  | 2.0      | 1.0      |
| quad      | 2.0      | 1.0      |
| triangle6 | 3.0      | 2.0      |
| quad9     | 3.0      | 2.0      |

and `--transient` on the decaying mode `exp(-pi^2 t) sin(pi x)` gives
temporal orders 1.0 (backward Euler) and 2.0 (Crank-Nicolson).

## Differences from the MATLAB code

The port keeps the user interface of the original but fixes several
problems found while translating it:

* **Shape-function derivatives.** `interpfunc_4.m` had the wrong sign on
  half of the derivatives and `interpfunc_9.m` had a wrong `dPhi/dxi` for
  node 1.  The Python versions are built as tensor products of 1-D Lagrange
  polynomials and are checked against finite differences.
* **Node ordering.** Elements use the Gmsh/VTK local ordering
  (counter-clockwise corners, then edge mid-nodes, then centre) so meshes
  from files work directly.  The MATLAB triangle had nodes 2 and 3 swapped
  and the 9-node element was row-major.
* **Quadrature.** The MATLAB code only sampled the diagonal `xi = eta` of
  the quadrilateral and used the 1-D rule for triangles too.  The port uses
  the full tensor-product Gauss-Legendre rule on quadrilaterals and the
  FSELIB `gauss_trgl` rule on triangles.  The MATLAB default `order = 1`
  under-integrates 4- and 9-node elements and makes the system singular, so
  `order = None` now selects each element's `full_order` (exact stiffness,
  mass and convection matrices for constant coefficients) and a warning is
  issued below the minimum order.
* **9-node mesh.** The connectivity was hard-coded for a 100-element mesh
  and elements overlapped; it is now derived from the mesh size.
* **`elems` for triangles** is the number of cells per side (each split in
  two triangles), the same meaning as for quadrilaterals.
* **Neumann conditions** are implemented (the MATLAB branch was empty) as
  the boundary integral of the prescribed flux.  Where two Dirichlet tags
  share a node the last listed tag wins consistently (the MATLAB code mixed
  both values at corners).
* The global system is a SciPy sparse matrix and can be solved iteratively.
* **Convection and time** are new: the MATLAB code was steady diffusion only.

## Roadmap toward flow simulation

Done: the infrastructure (elements, unstructured meshes, output, solvers,
CI), transient conduction (mass matrix, theta-method) and
convection-diffusion with SUPG stabilisation, which is the scalar prototype
of a flow solver.  Next, in order: nonlinear conductivity `kappa(T)`
(Newton), incompressible Navier-Stokes (Taylor-Hood or PSPG), then
compressible flow, where a finite-volume or discontinuous Galerkin
discretisation replaces continuous Galerkin.  Aircraft-scale RANS cases are
better run in an established solver such as SU2; this code is the place to
understand what such a solver does.
