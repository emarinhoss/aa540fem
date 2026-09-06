# Architecture

## Data flow

1. **Mesh** (`core/mesh.py`): nodes, one connectivity block per element
   type, named boundary edge sets.  Built by `geometry()` (structured
   rectangle) or read from a file (`io/mesh_files.py`, Gmsh physical groups
   become the tag names).
2. **Reference elements** (`core/elements.py`): each element type owns its
   shape functions (`core/shape_functions.py`), quadrature
   (`core/quadrature.py`), edge faces, orientation and, for Taylor-Hood,
   its pressure partner.  Everything downstream loops over
   `mesh.cells.items()` and asks `get_element(name)` for what it needs.
3. **Physics** assemble per cell block with vectorised element routines:
   - `transport/element.py` returns `ElementMatrices(K, C, M, f, dA)` for
     the scalar equation (diffusion, convection with SUPG, mass, source,
     Newton terms);
   - `incompressible/assembler.py` builds the viscous, mass and divergence
     matrices once and, per evaluation, the convective residual with its
     Jacobian plus (optionally) the SUPG, grad-div and PSPG stabilisation
     terms, which need the strong momentum residual at the quadrature
     points: the Laplacian of the discrete velocity comes from the
     shape-function Hessians (`core/shape_functions.py`, mapped in
     `transport/element.physical_laplacian`), and the stabilisation
     parameters measure the cell with the element metric tensor
     (`transport/element.element_metric`, stored per quadrature block),
     which keeps them smooth in the velocity on stretched cells;
   - `turbulence/spalart_allmaras.py` assembles the one-equation model on
     the same quadrature blocks and `turbulence/rans.py` couples it to the
     flow solver.
   Element matrices are scattered into global COO/CSR matrices with the
   connectivity; the `TaylorHoodSpace` maps nodes to `[u_x, u_y, p]` dofs.
4. **Boundary conditions**: Neumann fluxes are edge integrals added to the
   load vector; Dirichlet values are eliminated symmetrically by
   `linalg/dirichlet.py` (matrix rows/cols zeroed, unit diagonal, values
   moved to the right-hand side) so the reduced matrix can be reused for
   many right-hand sides.
5. **Solvers**: `linalg/solvers.py` (direct LU factorised once, CG, GMRES
   with pyamg preconditioning), `linalg/newton.py` (damped Newton on any
   `residual_jacobian(T)` callable, optional frozen Jacobian and divergence
   abort), `incompressible/steady.pseudo_transient` (backward-Euler
   pseudo-time steps with inner Newton and switched evolution relaxation,
   the robust path at high Reynolds number), `timestepping/rk.py`
   (Dormand-Prince RK45 on any `rhs(t, y)`).  The
   physics modules only provide residual/Jacobian or right-hand-side
   callables; the theta-method loops live with the physics because their
   Newton residual mixes the mass matrix with the physics operators.
6. **Output**: `io/mesh_files.write_vtk` (meshio) and `io/series.write_series`
   (`.pvd` collections); post-processing in `transport/postprocess.py` and
   `incompressible/solution.py` (traction forces).

## Conventions

- Local node ordering follows Gmsh / VTK for every element (corners
  counter-clockwise, then mid-edge nodes, then centre).
- Boundary edges are stored `(start, end[, mid])`; element faces walk the
  boundary with the element on the left, which fixes outward normals.
- Coefficient callables are matched by parameter name: `(x, y)`, plus `t`
  for time and `T` for temperature (`core/util.call_coeff`).
- Dirichlet tags listed later win at shared nodes.

## Where future capabilities go

**3-D.**  Add tetrahedra/hexahedra to `core/elements.py` (shape functions,
quadrature, faces as triangles/quads instead of edges) and let
`core/mesh.py` carry `(n, 3)` points and face sets for the boundary.  The
Jacobian and gradient helpers in `transport/element.py` need a 3x3 inverse,
the Neumann integral and the traction integral need surface quadrature, and
`TaylorHoodSpace` a third velocity block; assembly, solvers, Newton, RK45
and the output layer do not change.  Meshio already reads 3-D Gmsh files.

**Turbulence.**  A `turbulence/` subpackage with one module per closure
(Spalart-Allmaras, k-omega SST): each adds transported scalar(s) whose
operators reuse `transport/element.py` with the eddy viscosity feeding back
into the momentum viscosity of `incompressible/assembler.py`.  Wall
distance is a mesh utility (`core/`), wall functions a boundary-condition
type.  The coupling loop (segregated or monolithic Newton) belongs in
`incompressible/steady.py` / `transient.py`.

**Turbulence (done: Spalart-Allmaras, `turbulence/`).**  The closure is
one transported scalar assembled on the flow solver's quadrature blocks
(`turbulence/spalart_allmaras.py`), the wall distance is a mesh utility
(`core/wall_distance.py`), the eddy viscosity enters the momentum
equations as a nodal field (`FlowProblem.eddy_viscosity`), and
`turbulence/rans.py` runs the segregated coupling.  A second model (k-omega
SST) would follow the same pattern with two transported scalars; see
`docs/turbulence.md`.

**Compressible flow.**  Continuous Galerkin handles shocks poorly; a
`compressible/` subpackage would more naturally be a finite-volume or
discontinuous Galerkin discretisation on the same `Mesh`, `elements` and
`timestepping` layers (explicit RK is the standard choice there).

**Iterative solvers for large cases.**  Block preconditioning for the
saddle-point system (pressure Schur complement approximations) is what
lifts the direct-solver limit of roughly 10^5 unknowns; `linalg/krylov.py`
holds the PETSc fieldsplit version (least-squares-commutator Schur
complement, ILU or multigrid on the velocity block).

**Performance layers** (see `docs/parallel.md`): `backends/pattern.py`
(fixed sparsity pattern, deterministic scatter), `backends/numpy_kernels.py`
and `backends/numba_kernels.py` / `numba_sa.py` (element kernels, the
NumPy ones being the reference), `linalg/direct.py` (SuperLU / PETSc-MUMPS
factorisations), `hardware.py` and `cli.py` (probe, run configuration, the
50 % / 100 % prompt), `parallel/` (MPI).  Physics modules never import a
backend; they take element-local arrays from a kernel and scatter them, so
a GPU or distributed backend plugs in below them.
