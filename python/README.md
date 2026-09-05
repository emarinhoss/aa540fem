# Python port of the AA 540 anisotropic heat conduction FEM code

This directory is a NumPy/SciPy port of the MATLAB code in `../project`,
extended into a small but general 2-D finite element framework.  It solves
the scalar transport equation

    rho_c dT/dt + u . grad T - div( kappa(x, y) . grad T ) = f(x, y)

Without a velocity and time derivative this is the original
`div(kappa grad T) + f = 0`.  `kappa` is a full 2x2 conductivity tensor
that may depend on `T` (solved with Newton's method), `u` an optional
velocity field (with SUPG stabilisation), and the
boundaries carry Dirichlet (`T = T0`) or Neumann (`n . kappa grad T = q_n`)
conditions by name.  Elements are 3- and 6-node triangles or 4- and 9-node
quadrilaterals, on the built-in structured rectangle or on any mesh read
through meshio (Gmsh `.msh` files in particular).  Time integration uses
the theta-method (backward Euler or Crank-Nicolson).

It also solves the incompressible Navier-Stokes equations

    rho ( du/dt + (u . grad) u ) - mu lap(u) + grad p = rho f,   div u = 0

with Taylor-Hood elements (quadratic velocity, linear pressure) and
Newton's method; see "Incompressible flow" below.

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
python examples/nonlinear_conduction.py # kappa(T) = 1 + beta T, Newton vs Picard
python examples/cavity.py               # lid-driven cavity vs Ghia et al.
python examples/cylinder.py             # flow past a cylinder, Schaefer-Turek drag/lift
python examples/airfoil.py              # NACA 0012 at 5 deg: impulsive start, C_L / C_D history
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

### Temperature-dependent coefficients

```python
def material(x, y, T):                  # a parameter named T makes the problem nonlinear
    k = 1.0 + 0.5 * T
    return k, 0.0, 0.0, k, 0.0

problem = Problem(a=1, b=1, elems=8, elem_type="quad9",
                  bc_type={"left": 0, "right": 0}, bc_val={"left": 0.0, "right": 1.0},
                  material=material)
sol = solve(problem)                    # Newton; sol.info["residuals"] holds |R| per iteration
sol = solve(problem, newton=False)      # Picard (fixed point) iteration
result = solve_transient(problem, dt=0.1, t_end=2.0, theta=0.5)   # Newton at every step
```

The Jacobian carries the derivatives of `kappa` and `f` with respect to
`T` (central finite differences at the quadrature points), so convergence
is quadratic; the source may depend on `T` too (`f = g - T**3`).  A
temperature-dependent `velocity` or `rho_c` is evaluated at the current
`T` but not differentiated.  Newton systems are non-symmetric: use
`method="direct"` or `"gmres"`.

`material`, `velocity`, `rho_c` and the boundary values may take a third
argument `t` (`lambda x, y, t: ...`); parameters are matched by name, so
`def material(x, y, T, t)` gets both.  The solver re-evaluates
time-dependent data each step.  With time-independent coefficients the system matrix is factorised
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

### Incompressible flow

```python
from aa540fem import geometry, read_mesh
from aa540fem.flow import FlowProblem, solve_flow, solve_flow_transient

# Lid-driven cavity, Re = 100 (mu = 1/Re, rho = 1, lid speed 1)
mesh = geometry(1, 1, 32, "quad9")           # Q2 velocity, Q1 pressure on the corners
walls = {s: (0.0, 0.0) for s in ("left", "right", "bottom")}
cavity = FlowProblem(mesh, mu=0.01, rho=1.0, bc={**walls, "top": (1.0, 0.0)})
sol = solve_flow(cavity)                      # Newton from the Stokes solution
sol.u, sol.v, sol.p_nodal, sol.speed          # nodal fields
sol.save("cavity.vtu")                        # velocity vector + pressure for ParaView

# Channel with a body: Gmsh mesh with tags inlet / outlet / walls / cylinder
mesh = read_mesh("examples/cylinder_tri6.msh")   # P2/P1 triangles
flow = FlowProblem(mesh, mu=1e-3, rho=1.0,
                   bc={"inlet": (lambda x, y: 4 * 0.3 * y * (0.41 - y) / 0.41**2, 0.0),
                       "walls": (0.0, 0.0), "cylinder": (0.0, 0.0), "outlet": "open"})
sol = solve_flow(flow)
fx, fy = sol.forces("cylinder")               # traction integral of -p I + mu (grad u + grad u^T)

# Time dependent (theta = 0.5 Crank-Nicolson), Dirichlet values may take t
run = solve_flow_transient(flow, dt=0.01, t_end=1.0, theta=0.5, store_every=10,
                           startup_steps=4,                 # backward Euler first (impulsive start)
                           callback=lambda n, t, sol: print(t, sol.forces("cylinder")))
run.final.speed
run.save_series("out/flow")                   # .vtu per step + .pvd
```

The transient solver takes one Newton factorisation per step and reuses it
while the iteration contracts well (modified Newton), refreshing it
otherwise.

`examples/airfoil.py` is the aerodynamic case: a NACA 0012 at 5 degrees
angle of attack (the profile is rotated in the mesh so the freestream is
along x), chord Reynolds number 1000, impulsively started and integrated at
a fixed time step with fields written every `--store-every` steps and lift
and drag coefficients logged every step to `forces.csv`, followed by a
steady Newton solve from the final state.  `make_meshes.py` builds the
far-field mesh with Gmsh (`naca4()` generates any 4-digit profile, and
`make_airfoil(alpha_deg=...)` any angle of attack).

Boundary values are `(ux, uy)` pairs (constants, callables of `(x, y[, t])`,
or `None` for a free component, e.g. `(None, 0.0)` on a symmetry line) or
`"open"` for the do-nothing outflow `mu du/dn - p n = 0`.  When no boundary
is open the pressure is pinned at one node (`pin_value`).  The Jacobian is
an indefinite saddle-point matrix, so the direct solver is used; a
continuation in Reynolds number is done by passing a previous `sol.U` as
`U0` (see `examples/cavity.py`).

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
| `aa540fem/nonlinear.py`             | damped Newton / Picard iteration                      | (new)         |
| `aa540fem/flow.py`                  | Navier-Stokes: Taylor-Hood space, Newton, theta-method, forces | (new) |
| `aa540fem/util.py`                  | `t` / `T` argument matching for user callables        | (new)         |
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
skew advection), a decaying mode, a transient manufactured solution with
convection and time-dependent Dirichlet data, and nonlinear cases: the
Kirchhoff problem `kappa = 1 + T` against its exact solution (Newton
converging quadratically in a handful of iterations, Picard needing more),
manufactured solutions with `kappa = 1 + T^2` and a `T^3` source, with
convection, and in time.  For the flow solver: Poiseuille flow reproduced to
machine precision (velocity, pressure, wall forces, on quads and
triangles), the Kovasznay solution at Re = 40 converging at order 3 in
velocity and 2 in pressure, the lid-driven cavity at Re = 100 against Ghia,
Ghia & Shin, the decaying Taylor-Green vortex in time, and the
Schaefer-Turek cylinder benchmark at Re = 20:

| quantity | computed | reference |
|----------|----------|-----------|
| C_D      | 5.5791   | 5.5795    |
| dp       | 0.11748  | 0.11752   |
| C_L      | 0.0067   | 0.0106    |

(lift is two orders of magnitude smaller than drag and needs a finer mesh
around the cylinder to converge).  The airfoil case at Re = 1000 has no
exact reference; `examples/airfoil.py` on its 15k-node mesh (dt = 0.05,
160 steps, 4 minutes) gives

| NACA 0012, alpha = 5 deg, Re = 1000 | C_L    | C_D    |
|-------------------------------------|--------|--------|
| transient at t U / c = 8            | 0.2565 | 0.1289 |
| steady Newton from that state       | 0.2481 | 0.1285 |

Published laminar computations for this case (e.g. Kurtulus 2015) report
C_L of roughly 0.25-0.3 and C_D of roughly 0.13; the slow drift of C_L
after t = 2 is the wake and the separated region on the suction side
settling, which take many chord times at this Reynolds number.
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
* **Convection, time, nonlinear conductivity and Navier-Stokes** are new:
  the MATLAB code was steady linear diffusion only.

## Roadmap toward flow simulation

Done: the infrastructure (elements, unstructured meshes, output, solvers,
CI), transient conduction, convection-diffusion with SUPG, nonlinear
coefficients with Newton's method, and incompressible Navier-Stokes with
Taylor-Hood elements including body forces on a boundary (drag and lift).
What the flow solver still lacks for aerodynamic work, roughly in order of
usefulness: SUPG/PSPG stabilisation for higher Reynolds numbers on coarser
meshes, an iterative saddle-point solver
(block preconditioning) to go beyond ~10^5 unknowns, a turbulence model,
and finally compressibility, where a finite-volume or discontinuous
Galerkin discretisation replaces continuous Galerkin.  Aircraft-scale RANS
cases are better run in an established solver such as SU2; this code is
the place to understand what such a solver does.
