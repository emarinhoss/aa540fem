# aa540fem: finite element solvers for transport and incompressible flow

A NumPy/SciPy finite element code that grew out of the AA 540 (University
of Washington) anisotropic heat conduction project.  The original MATLAB
code is kept in `matlab/`, its report in `report/`; everything else is the
Python framework.  It solves the scalar transport equation

    rho_c dT/dt + u . grad T - div( kappa(x, y) . grad T ) = f(x, y)

Without a velocity and time derivative this is the original
`div(kappa grad T) + f = 0`.  `kappa` is a full 2x2 conductivity tensor
that may depend on `T` (solved with Newton's method), `u` an optional
velocity field (with SUPG stabilisation), and the
boundaries carry Dirichlet (`T = T0`) or Neumann (`n . kappa grad T = q_n`)
conditions by name.  Elements are 3- and 6-node triangles or 4- and 9-node
quadrilaterals, on the built-in structured rectangle or on any mesh read
through meshio (Gmsh `.msh` files in particular).  Time integration uses
the adaptive explicit Runge-Kutta 45 (Dormand-Prince) by default or the
implicit theta-method (backward Euler or Crank-Nicolson).

It also solves the incompressible Navier-Stokes equations

    rho ( du/dt + (u . grad) u ) - mu lap(u) + grad p = rho f,   div u = 0

with Taylor-Hood elements (quadratic velocity, linear pressure) and
Newton's method; see "Incompressible flow" below.

## Repository layout

```
src/aa540fem/          the package (src layout, `pip install -e .`)
  core/                meshes, reference elements, shape functions, quadrature, helpers
  io/                  mesh files through meshio, VTK and .pvd time-series output
  linalg/              Dirichlet elimination, direct/Krylov solvers, Newton iteration
  timestepping/        explicit Runge-Kutta 45 (the theta-method lives with each physics)
  transport/           scalar transport: material, element operators, boundary conditions,
                       steady/nonlinear/transient solvers, post-processing
  incompressible/      Navier-Stokes: problem, Taylor-Hood space, assembler, solution
                       and forces, steady and transient solvers
examples/              runnable cases (meshes and Gmsh geometries in examples/meshes/)
scripts/               convergence studies
tests/                 verification suite (pytest)
docs/architecture.md   how the pieces fit and where 3-D, turbulence, compressible go
matlab/, report/       the original course project
```

## Installation

```bash
pip install -e ".[mesh,amg,dev]"     # from the repository root
# or just the runtime pieces:
pip install -r requirements.txt
```

Extras: `mesh` (meshio, for reading mesh files and writing ParaView output),
`amg` (pyamg, algebraic multigrid preconditioner for the CG solver),
`meshgen` (gmsh, only needed to regenerate `examples/*.msh`), `dev`
(pytest, ruff).

## Running

```bash
python examples/heat_rectangle.py       # rectangle from main.m, writes temperature.png
python examples/heat_rectangle.py --show --method cg --vtk out.vtu
python examples/annulus.py              # curved annulus mesh, exact solution ln(r)/ln(2)
python examples/convection_diffusion.py # Galerkin vs SUPG on a boundary layer
python examples/rotating_hill.py        # transient advection, writes a ParaView series
python examples/nonlinear_conduction.py # kappa(T) = 1 + beta T, Newton vs Picard
python examples/cavity.py               # lid-driven cavity vs Ghia et al.
python examples/cylinder.py             # flow past a cylinder, Schaefer-Turek drag/lift
python examples/airfoil.py              # NACA 0012 at 5 deg: impulsive start, C_L / C_D history
python examples/flat_plate.py           # laminar plate at Re 1e5 vs Blasius (stabilised, PTC)
python examples/cylinder_shedding.py    # Re 100 vortex shedding, Strouhal number (Schaefer-Turek 2D-2)
python examples/turbulent_flat_plate.py # Spalart-Allmaras RANS plate at Re 1e6 vs law of the wall
python scripts/convergence.py           # mesh-convergence study, prints observed orders
python scripts/convergence.py --transient   # temporal orders of backward Euler / Crank-Nicolson
python -m pytest                        # verification suite
ruff check .                            # lint
```

As in the MATLAB version, the files meant to be edited for a quick run are
`examples/heat_rectangle.py` (domain, mesh, element type, boundary
conditions, quadrature order) and `src/aa540fem/transport/material.py`
(conductivity tensor and heat source); anything beyond that is done through
the `Problem` / `FlowProblem` objects.

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
mesh = read_mesh("examples/meshes/annulus_tri6.msh")
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
`examples/meshes/annulus.geo`); those names become the keys of `bc_type` and
`bc_val`.  Generate with `gmsh -2 file.geo` or through the Python API as in
`examples/meshes/make_meshes.py`, which also shows how to get quadrilaterals
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
mesh = read_mesh("examples/meshes/cylinder_tri6.msh")   # P2/P1 triangles
flow = FlowProblem(mesh, mu=1e-3, rho=1.0,
                   bc={"inlet": (lambda x, y: 4 * 0.3 * y * (0.41 - y) / 0.41**2, 0.0),
                       "walls": (0.0, 0.0), "cylinder": (0.0, 0.0), "outlet": "open"})
sol = solve_flow(flow)
fx, fy = sol.forces("cylinder")               # traction integral of -p I + mu (grad u + grad u^T)

# Time dependent: adaptive RK45 (default) with the fields stored every 0.1 time units
run = solve_flow_transient(flow, dt=0.005, t_end=1.0, output_interval=0.1,
                           callback=lambda n, t, sol: print(t, sol.forces("cylinder")))
# ... or Crank-Nicolson at a fixed step with backward-Euler start-up (impulsive start)
run = solve_flow_transient(flow, dt=0.01, t_end=1.0, scheme="theta", theta=0.5,
                           store_every=10, startup_steps=4)
run.final.speed
run.save_series("out/flow")                   # .vtu per step + .pvd
```

### Three dimensions

The same solver runs in 3-D on 27-node hexahedra (Q2/Q1), 10-node
tetrahedra (P2/P1) and 18-node prisms (the triangle6 x line3 tensor
product, with a 6-node linear prism for the pressure); every routine takes
the dimension from the mesh, so a
`FlowProblem` on a 3-D mesh has three velocity components in its boundary
conditions, callables receive `(x, y, z[, t])`, forces come back as
`(Fx, Fy, Fz)` and `wall_traction` integrates over the boundary faces.

```python
from aa540fem.core.mesh import box                      # structured [0,a]x[0,b]x[0,c]
mesh = box(2.0, 1.0, 1.0, (4, 3, 2), "hexahedron27")   # or "tetra10" (6 per cell)
mesh = read_mesh("examples/meshes/box_tet10.msh")       # Gmsh: physical surfaces are the tags
inflow = lambda x, y, z: 4 * y * (1 - y)
prob = FlowProblem(mesh, mu=0.05, stabilisation=True,
                   bc={"left": (inflow, 0.0, 0.0), "bottom": (0.0, 0.0, 0.0),
                       "top": (0.0, 0.0, 0.0), "front": (None, None, 0.0),
                       "back": (None, None, 0.0), "right": "open"})
sol = solve_flow(prob)
sol.forces("bottom")                                    # (Fx, Fy, Fz)
```

`examples/channel3d.py` runs this extruded Poiseuille flow (Re 20, either
mesh); the computed velocity, pressure and wall force match the exact
solution to 1e-14, and the tests (`tests/test_flow3d.py`,
`tests/test_elements3d.py`) check the same on both element types with and
without stabilisation, the stabilised Jacobian against finite differences,
the theta and RK45 integrators, the fieldsplit solver and the MPI
domain decomposition in 3-D.  Tetrahedra use a conical-product Gauss-Jacobi
rule (27 points for the P2/P1 pair).  The Spalart-Allmaras model and the
RANS coupling also run in 3-D (`|curl u|`, wall distance to faces,
`Mesh.extrude` for 2-D/3-D consistency tests, see `docs/turbulence.md`
section 3.8).  `Mesh.extrude(depth, layers, grading=...)` turns any 2-D
quadratic mesh into hexahedra and prisms, with the layers graded towards
the two end planes, so the committed boundary-layer meshes (quadrilaterals
along the walls, triangles elsewhere) serve in 3-D as well:
`examples/cylinder3d.py --extrude` and `examples/turbulent_flat_plate.py
--extrude` below.

**Schaefer-Turek 3D-1Z** (`examples/cylinder3d.py`: cylinder of diameter
0.1 along z in the 2.5 x 0.41 x 0.41 channel, Re 20 on the mean inflow
velocity 0.2, stabilised P2/P1 tetrahedra from
`make_meshes.py::make_cylinder3d`).  Reference values 6.185 / 0.0094 /
0.1710 (John 2002, Bayraktar, Mierka and Turek 2012) and the benchmark
intervals of Schaefer and Turek (1996):

| mesh (unknowns) | solver | wall | C_D | C_L | dp |
|---|---|---|---|---|---|
| `cyl3d_coarse.msh` (81k) | MUMPS / fieldsplit-LSC / 4-rank domain decomposition | 268 s / 235 s / 97 s | 6.240 | -0.016 | 0.1758 |
| `cylinder3d_tet10.msh` (132k) | MUMPS, 6 GB | 592 s | 6.227 | 0.014 | 0.1742 |
| `cyl3d_fine2.msh` (227k) | fieldsplit-LSC, 8 GB | 992 s | 6.207 | 0.014 | 0.1718 |
| `cylinder3d1z_bl.msh` extruded, 5 graded layers (164k) | 4-rank domain decomposition, 2-D start | 364 s | 6.170 | 0.0112 | 0.1656 |
| `cylinder3d1z_bl.msh` extruded, 7 graded layers (222k) | fieldsplit-LSC, 2-D start | 1055 s | 6.181 | 0.0105 | 0.1685 |
| reference (interval) | | | 6.185 (6.05-6.25) | 0.0094 (0.008-0.010) | 0.1710 (0.165-0.175) |

On the tetrahedral meshes drag and pressure difference converge
monotonically towards the reference (0.4 % and 0.5 % away on the finest
mesh) and lie inside the intervals from 132k unknowns on; the lift, 0.15 %
of the drag, has the right sign from 132k unknowns but stays at 0.014,
which needs a boundary-layer resolution of the cylinder that isotropic
tetrahedra do not have (refining the cylinder surface alone while
coarsening the channel, `cyl3d_bl.msh`, makes every number worse: at this
blockage the channel walls and the wake set the drag).  The extruded
boundary-layer mesh (`--extrude`: the 2-D cross-section with quadrilateral
layers on the cylinder, extruded into prisms and hexahedra with layers
graded towards the end walls) brings the lift to 0.0112 at 164k unknowns
with drag and pressure difference inside the intervals, and to 0.0105
with seven layers (222k unknowns, drag 0.06 % and pressure difference
1.5 % from the reference), converging with the spanwise resolution; the
next refinement no longer fits the 15 GB of this machine.
The direct solver and the Krylov solver give identical numbers; the Krylov
path is the one that fits the 227k case in memory.

**Lid-driven cube** (`examples/cavity3d.py`, Re 100, hexahedron27, lid
``u = 1`` at y = 1).  Extrema of the centreline velocities on successive
meshes (nodal values, so the positions are quantised to the node spacing):

| cells per side (unknowns) | u_min on x = z = 0.5 (at y) | v_max on y = z = 0.5 (at x) | v_min (at x) | wall |
|---|---|---|---|---|
| 8 (15k) | -0.1663 (0.500) | 0.1217 (0.250) | -0.2053 (0.812) | 37 s |
| 12 (49k) | -0.1834 (0.458) | 0.1337 (0.208) | -0.2162 (0.833) | 248 s |
| 16 (113k) | -0.1917 (0.469) | 0.1384 (0.219) | -0.2263 (0.812) | 1298 s |

The changes halve with each refinement (second order in the cell size), so
the converged extrema are about -0.20, 0.14 and -0.24.  The published
tables of Ku, Hirsh and Taylor (1987), Wong and Baker (2002) and
Albensoeder and Kuhlmann (2005) could not be retrieved from this
environment, so the quantitative comparison for this case is still open;
the profiles are written to `centreline_<n>.csv` for it.

### High Reynolds numbers: stabilisation and continuation

```python
prob = FlowProblem(mesh, mu=1e-5, rho=1.0, bc=..., stabilisation=True)   # SUPG + grad-div
sol = solve_flow(prob, continuation="auto")   # Newton, then pseudo-transient continuation if needed
tr = sol.wall_traction("plate")               # x, y, tx, ty, nx, ny, p, weight along the wall
cf = 2 * tr["tx"] / (rho * U**2)              # skin friction distribution
```

`stabilisation=True` adds residual-based SUPG to the momentum equations
(the streamline weight applied to the full residual, including the viscous
term through the shape-function Hessians, so it is consistent: Poiseuille
flow stays exact to round-off) and grad-div stabilisation; `pspg=True` adds
pressure stabilisation, which Taylor-Hood does not need.  It is off by
default.  The parameters use the element metric tensor
(`tau = [(2/dt)^2 + u.G u + nu^2 G:G/2]^(-1/2)`, Shakib 1991 / Bazilevs et
al. 2007), which is smooth in the velocity; Tezduyar's flow-direction
element length is available as `element_length="streamline"`, but on
strongly stretched boundary-layer cells it jumps between the cell length
and its height with the slightest rotation of the velocity and stalls
Newton.  The Jacobian is consistent, including the derivatives of the
stabilisation parameters (a finite-difference check to 1e-9 is in the
tests).
`continuation="ptc"` (or the `"auto"` fallback, which restarts from the
best Newton iterate) solves backward-Euler pseudo-time steps with a few
Newton iterations each and grows the step by switched evolution
relaxation; the step is local (scaled by the convective time scale of the
surrounding cells, so `dtau0` is a CFL number) and starts from the divergence-free projection of the
initial velocity.  It is what makes the laminar flat plate at Re 1e5
converge from rest (10 pseudo-time steps), where plain Newton diverges.  Boundary-layer
meshes (quadrilaterals extruded from the wall inside a triangular mesh)
come from `make_meshes.make_flat_plate()` and the `boundary_layer=` option
of `make_meshes.make()`.

### Turbulence (RANS)

```python
from aa540fem.turbulence import solve_rans

prob = FlowProblem(mesh, mu=1e-6, rho=1.0, stabilisation=True, bc=...)   # laminar viscosity
rans = solve_rans(prob, wall_tags=["plate"])     # Spalart-Allmaras, segregated coupling
rans.flow.wall_traction("plate")                 # turbulent skin friction
rans.nu_t, rans.nu_tilde, rans.distance          # eddy viscosity, working variable, wall distance
rans.save("rans.vtu")
```

The Spalart-Allmaras one-equation model (negative variant, no trip term)
is documented with its equations, references and implementation in
[`docs/turbulence.md`](docs/turbulence.md); it needs wall-resolved meshes
(first cell at y+ of about 1, see `make_meshes.make_turbulent_flat_plate`).
`examples/turbulent_flat_plate.py` validates it on a zero-pressure-gradient
plate at Re 1e6 per unit length (Re_x up to 2e6; 4 minutes with numba and
MUMPS, 9 with NumPy and SuperLU):

| turbulent flat plate, x = 1.5 (Re_theta = 2884) | computed | reference |
|-------------------------------------------------|----------|-----------|
| Cf                                              | 0.00327  | 0.00323 (Coles-Fernholz) |
| u+ - (ln(y+)/0.41 + 5) for 30 < y+ < 300        | 0.15 to 0.87 | 0 |
| Cf vs White's Re_x correlation, 3e5 < Re_x < 2e6 | 4.6 to 12.8 % low | (turbulent from x = 0) |

For the incompressible system the pressure is a constraint multiplier, not
an ODE unknown, so the RK45 scheme is applied to the velocity with a
pressure projection at every stage: each stage solves the constant
saddle-point system `[[M, -B^T], [B, 0]]` (factorised once for the whole
run), which keeps the discrete velocity exactly divergence-free and needs
no Newton iteration.  Its step is limited by the viscous and convective
stability limits of the finest cells.  The theta scheme takes one Newton
factorisation per step and reuses it while the iteration contracts well
(modified Newton), refreshing it otherwise; it allows much larger steps.

`examples/airfoil.py` is the aerodynamic case: a NACA 0012 at 5 degrees
angle of attack (the profile is rotated in the mesh so the freestream is
along x), chord Reynolds number 1000, impulsively started and integrated
with RK45 (or `--scheme theta`) with fields written at fixed output times
(`--output-interval`) and lift and drag coefficients logged at every step
to `forces.csv`, followed by a steady Newton solve from the final state.
`examples/meshes/make_meshes.py` builds the far-field mesh with Gmsh
(`naca4()` generates any 4-digit profile, and `make_airfoil(alpha_deg=...)`
any angle of attack).

Boundary values are `(ux, uy)` pairs (constants, callables of `(x, y[, t])`,
or `None` for a free component, e.g. `(None, 0.0)` on a symmetry line) or
`"open"` for the do-nothing outflow `mu du/dn - p n = 0`.  When no boundary
is open the pressure is pinned at one node (`pin_value`).  The Jacobian is
an indefinite saddle-point matrix, so the direct solver is used; a
continuation in Reynolds number is done by passing a previous `sol.U` as
`U0` (see `examples/cavity.py`).

### Running faster

Every example takes the run-configuration flags of `aa540fem.cli`
(`--machine {50,100}`, `--threads N`, `--assembly {numpy,numba}`,
`--direct {superlu,petsc}`, `--linear {scipy,petsc,petsc-cuda}`,
`--memory-budget GB`, `--no-prompt`); on a terminal the script prints what
the machine offers and asks once whether to use 50 % or 100 % of it.
`python -m aa540fem.cli` prints the probe, and `AA540FEM_*` environment
variables set the same choices for batch jobs.

```
pip install -e ".[numba]"     # threaded element assembly (numba)
pip install -e ".[petsc]"     # petsc4py + mpi4py: MUMPS factorisations, MPI, Krylov solver
pip install -e ".[fast]"      # both, plus psutil/threadpoolctl for the probe
```

The solver assembles every matrix on one fixed sparsity pattern, assembles
residuals only where no Jacobian is needed (Runge-Kutta stages, Newton
trial points), runs the element kernels in numba threads and factorises
with MUMPS when PETSc is installed.  Wall times on a 4-core machine
(details in [`docs/performance.md`](docs/performance.md)):

| case | before | numba + SuperLU | numba + PETSc/MUMPS |
|---|---|---|---|
| steady Newton, cylinder mesh, 6 iterations | 9.1 s | 9.6 s | 4.8 s |
| theta scheme, 20 steps, cylinder mesh | 73 s | 29 s | 18 s |
| RK45, 20 steps, cylinder mesh | 44 s | 7.2 s | 3.2 s |
| RANS coupling, 2 outer iterations, plate mesh | 33 s | 26 s | 17 s |

[`docs/parallel.md`](docs/parallel.md) explains which option fits which
machine, and what the MPI (`aa540fem.parallel`: replicated assembly with a
distributed MUMPS factorisation, or the METIS domain decomposition of
`parallel/flow.py` in which every rank assembles its own elements; the
steady cylinder takes 1.2 s on four ranks against 4.8 s in one process,
`mpirun -n 4 python benchmarks/bench_mpi.py`) and GPU (PETSc Krylov path: FGMRES with a fieldsplit
preconditioner, least-squares-commutator Schur complement, `--linear petsc`
or `petsc-cuda`) routes are.

## Modules

| Module | Contents | MATLAB origin |
|--------|----------|---------------|
| `core/elements.py` | `ReferenceElement` registry (2-D and 3-D): nodes, faces, edges, quadrature, Taylor-Hood pairs | (new) |
| `core/shape_functions.py` | Lagrange shape functions in meshio/VTK node ordering (triangles, quads, tetrahedra, hexahedra) | `interpfunc_*.m` |
| `core/quadrature.py` | Gauss-Legendre, triangle, hexahedral and tetrahedral (conical product) rules | `gauss_legendre_quad.m`, `gauss_trgl.m` |
| `core/mesh.py` | `Mesh` data structure (2-D / 3-D), structured `geometry()` and `box()` | `geometry.m` |
| `core/util.py` | argument matching (`t`, `T`) for user callables | (new) |
| `io/mesh_files.py`, `io/series.py` | meshio import/export, VTK and `.pvd` output | (new) |
| `linalg/dirichlet.py` | symmetric elimination of prescribed values | `dirichlet.m` |
| `linalg/solvers.py`, `linalg/newton.py` | direct / CG / GMRES, damped Newton | (new) |
| `timestepping/rk.py` | Dormand-Prince RK45 with step control | (new) |
| `transport/material.py` | user material and source | `conductivity_and_forcing.m` |
| `transport/element.py` | Jacobians, diffusion/convection/mass operators, SUPG, Newton terms | `elem_eqn.m` |
| `transport/boundary.py` | Dirichlet wrapper, Neumann edge integrals | `dirichlet.m` |
| `transport/problem.py` | `Problem`, operator assembly, boundary conditions, steady solve | `main.m` |
| `transport/nonlinear.py`, `transport/transient.py` | Newton steady solve, RK45 / theta time stepping | (new) |
| `transport/postprocess.py` | centroid gradients and fluxes, L2/H1 error norms | (new) |
| `incompressible/problem.py`, `space.py`, `assembler.py` | `FlowProblem`, Taylor-Hood dofs, Navier-Stokes residual/Jacobian with SUPG / grad-div / PSPG | (new) |
| `incompressible/solution.py`, `steady.py`, `transient.py` | fields, wall traction and forces, Newton / PTC, RK45 / theta | (new) |
| `turbulence/spalart_allmaras.py`, `turbulence/rans.py` | Spalart-Allmaras model and its coupling to the flow | (new) |
| `core/wall_distance.py`, `linalg/continuation.py` | wall distance, pseudo-transient continuation | (new) |
| `examples/heat_rectangle.py` | driver with the user inputs of the original | `main.m` |

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

(lift is two orders of magnitude smaller than drag; on the boundary-layer
mesh `cylinder_bl.msh` with stabilisation the same case gives
C_D = 5.5795, C_L = 0.0106, both on the reference).  The laminar flat plate
at Re_L = 1e5 (`examples/flat_plate.py`, stabilised, pseudo-transient
continuation from rest in 10 steps, 30 s) gives a skin friction within 5 %
of Blasius for 1e4 < Re_x < 1e5 (mean 2.7 %) and velocity profiles within
0.014 of the similarity solution at three stations; near the leading edge
the Navier-Stokes skin friction exceeds Blasius, as it should.  The cylinder at
Re = 100 (`examples/cylinder_shedding.py`, stabilised Crank-Nicolson)
sheds vortices with, on the standard mesh (`cylinder_bl.msh`, 10.8k nodes,
dt = 0.005: 19 minutes with numba and MUMPS, over an hour with NumPy and
SuperLU), with the time step halved (`--dt 0.0025`, 29 minutes) and on a
twice-finer mesh (`--mesh examples/meshes/cylinder_bl_fine.msh`, 41k nodes,
39 minutes),

| Schaefer-Turek 2D-2 | 10.8k nodes, dt 0.005 | 10.8k nodes, dt 0.0025 | 41k nodes, dt 0.005 | reference |
|---------------------|--------|--------|--------|--------|
| Strouhal number     | 0.3013 | 0.3016 | 0.3014 | 0.2995 |
| C_D,max             | 3.2303 | 3.2271 | 3.2320 | 3.2298 |
| C_L,max             | 0.9958 | 0.9885 | 0.9972 | 1.0002 |

All three runs are inside the intervals of Schaefer and Turek (St 0.295 to
0.305, C_D,max 3.22 to 3.24, C_L,max 0.99 to 1.01) and periodic from
t = 6 on.  Halving the time step changes the force maxima by 0.3 to 0.7 %
and the mesh refinement by 0.05 to 0.15 %, while the Strouhal number moves
by less than 0.1 % under either, so the 0.6 % between the computed St and
the reference is a property of the discretisation as run (the
stabilisation terms and the open outlet), not of the resolution.

The airfoil case at Re = 1000 has no
exact reference; `examples/airfoil.py` on its 15k-node mesh (dt = 0.05,
160 steps, 4 minutes) gives

| NACA 0012, alpha = 5 deg, Re = 1000 | C_L    | C_D    | steps | wall time |
|-------------------------------------|--------|--------|-------|-----------|
| RK45 (rtol 1e-4), t U / c = 8       | 0.2562 | 0.1289 | 3956 accepted, 203 rejected, dt 1e-4 to 3e-3 | 11.5 min |
| theta (Crank-Nicolson, dt = 0.05)   | 0.2565 | 0.1289 | 160 | 4 min |
| steady Newton from the final state  | 0.2481 | 0.1285 | 4 Newton iterations | 5 s |

The two time integrators agree to three decimals along the whole history;
the explicit scheme needs about 25 times more steps (its step is set by
the viscous stability limit of the smallest cells at the leading edge) but
each step is cheap (no Newton, one factorisation for the whole run).

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
temporal orders 1.0 (backward Euler), 2.0 (Crank-Nicolson) and 5.0
(fixed-step RK45).

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
Also done: SUPG/grad-div stabilisation, pseudo-transient continuation,
boundary-layer meshes and wall-shear output, validated on the Blasius plate
and the shedding cylinder; the Spalart-Allmaras RANS model
(`docs/turbulence.md`); threaded assembly, MUMPS, the fieldsplit/LSC
Krylov solver and MPI domain decomposition (`docs/parallel.md`); and the
3-D discretisation (hexahedra and tetrahedra, verified on extruded exact
flows).  What remains for aerodynamic work: 3-D validation cases with
reference data (Schaefer-Turek 3D-1Z, the 3-D cavity), the turbulence
model in 3-D, and finally compressibility, where a finite-volume or
discontinuous Galerkin discretisation replaces continuous Galerkin.  Aircraft-scale RANS
cases are better run in an established solver such as SU2; this code is
the place to understand what such a solver does.
