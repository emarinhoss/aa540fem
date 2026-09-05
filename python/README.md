# Python port of the AA 540 anisotropic heat conduction FEM code

This directory is a NumPy/SciPy port of the MATLAB code in `../project`,
extended into a small but general 2-D finite element framework.  It solves

    div( kappa(x, y) . grad T ) + f(x, y) = 0

where `kappa` is a full 2x2 conductivity tensor, with Dirichlet
(`T = T0`) or Neumann (`n . kappa grad T = q_n`) conditions on named
boundaries, using 3- and 6-node triangles or 4- and 9-node quadrilaterals,
on the built-in structured rectangle or on any mesh read through meshio
(Gmsh `.msh` files in particular).

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
python scripts/convergence.py           # mesh-convergence study, prints observed orders
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
| `aa540fem/element.py`               | Jacobians, element stiffness and load (vectorised)    | `elem_eqn.m`  |
| `aa540fem/boundary.py`              | Dirichlet elimination, Neumann edge integrals         | `dirichlet.m` |
| `aa540fem/solver.py`                | `Problem`, assembly, direct / CG+AMG solve            | `main.m`      |
| `aa540fem/postprocess.py`           | Centroid gradients and fluxes, L2/H1 error norms      | (new)         |
| `aa540fem/conductivity_and_forcing.py` | User material and source                           | `conductivity_and_forcing.m` |
| `main.py`                           | Driver with the user inputs                           | `main.m`      |
| `examples/`, `scripts/`, `tests/`   | Annulus case, convergence study, verification suite   |               |

## Verification

`tests/` checks the shape functions (nodal property, partition of unity,
finite-difference derivatives), quadrature tables, mesh orientation and
boundary detection, exact reproduction of linear and quadratic fields with
anisotropic conductivity and mixed boundary conditions, the Gmsh annulus
with curved elements, CG against the direct solver, and VTK round trips.
`scripts/convergence.py` on the manufactured solution
`T = sin(pi x/a) sin(pi y/b)` gives the expected orders:

| element   | L2 order | H1 order |
|-----------|----------|----------|
| triangle  | 2.0      | 1.0      |
| quad      | 2.0      | 1.0      |
| triangle6 | 3.0      | 2.0      |
| quad9     | 3.0      | 2.0      |

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
  `order = None` now selects the lowest order that fully integrates each
  element and a warning is issued for anything lower.
* **9-node mesh.** The connectivity was hard-coded for a 100-element mesh
  and elements overlapped; it is now derived from the mesh size.
* **`elems` for triangles** is the number of cells per side (each split in
  two triangles), the same meaning as for quadrilaterals.
* **Neumann conditions** are implemented (the MATLAB branch was empty) as
  the boundary integral of the prescribed flux.
* The global system is a SciPy sparse matrix and can be solved iteratively.

## Roadmap toward flow simulation

This milestone covers the infrastructure (elements, unstructured meshes,
output, solvers, CI).  The physics steps that build on it, in order:
transient conduction (mass matrix, implicit time stepping), nonlinear
conductivity (Newton), convection-diffusion with SUPG stabilisation,
incompressible Navier-Stokes (Taylor-Hood or PSPG), then compressible flow,
where a finite-volume or discontinuous Galerkin discretisation replaces
continuous Galerkin.  Aircraft-scale RANS cases are better run in an
established solver such as SU2; this code is the place to understand what
such a solver does.
