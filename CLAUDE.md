# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`aa540fem` is a Python finite-element package that started as a port of a MATLAB
AA 540 heat-conduction code (the original lives in `matlab/`, the report in
`report/`) and grew into a 2D incompressible Navier-Stokes solver (Taylor-Hood
Q2/Q1 and P2/P1, SUPG/PSPG/grad-div stabilisation, Newton and pseudo-transient
continuation, theta and RK45 time stepping, Spalart-Allmaras RANS) with
threaded assembly, MUMPS/PETSc direct solvers, a fieldsplit Krylov solver and
an MPI domain-decomposition path. Development happens on the branch
`python-implementation`. Everything is validated against reference cases
(Poiseuille, Kovasznay, Ghia cavity, Schaefer-Turek cylinder, Blasius and
Coles-Fernholz plates); the numbers in `README.md` and `docs/performance.md`
are the regression gate for any change to the solver.

## Commands

```
pip install -e ".[mesh,amg,dev]"        # base: NumPy kernels, SuperLU; meshio for .msh files
pip install -e ".[numba]"               # threaded element kernels (the default when importable)
pip install -e ".[petsc]"               # petsc4py + mpi4py: MUMPS, MPI, Krylov (docs/parallel.md §6 has the build recipe)
ruff check .                            # lint (line length 100, isort); CI fails on any finding
python -m pytest -q                     # full suite (~1.5-2.5 min); optional-dependency tests skip when missing
python -m pytest -q tests/test_flow.py::test_name
python -m pytest -q tests/test_backends.py tests/test_krylov.py      # numba vs numpy, PETSc solver
mpirun --oversubscribe -n 4 python -m pytest -q tests/test_mpi.py    # add --allow-run-as-root in containers
python examples/cylinder_shedding.py --no-prompt --no-plot --outdir /tmp/out   # any example; --help lists flags
python benchmarks/bench_flow.py [--cases steady theta rk45 rans krylov] [--record LABEL]
mpirun -n 4 python benchmarks/bench_mpi.py --threads 1 [--record LABEL]
python -m aa540fem.cli                  # print the hardware probe and resolved run configuration
python examples/meshes/make_meshes.py   # regenerate every committed mesh (needs gmsh); call make(...) for one
```

`--record LABEL` appends a benchmark table to `docs/performance.md`; keep the
prose there in sync when the numbers change. Examples take the run-configuration
flags from `aa540fem.cli.add_run_arguments` (`--machine {50,100}`, `--threads`,
`--assembly`, `--direct`, `--linear`, `--memory-budget`, `--no-prompt`). Without
`--no-prompt` or `--machine` an example asks on a TTY; `tests/conftest.py` sets
`AA540FEM_NO_PROMPT=1`. Long validation runs (shedding benchmark, turbulent
plate) take 15-75 minutes: run them detached with a log file.

## Architecture

**Layers** (`docs/architecture.md`, `docs/parallel.md`): `core/` (mesh,
reference elements, quadrature, wall distance) -> `transport/` (scalar
convection-diffusion, the original heat code) and `incompressible/`
(Navier-Stokes) -> `turbulence/` (Spalart-Allmaras coupled through
`FlowProblem.eddy_viscosity`) -> `linalg/` and `timestepping/` -> `parallel/`.
Physics modules never import a backend or a solver library directly; they go
through `aa540fem.backends` (assembly) and `aa540fem.linalg.direct.factorise`
(factorisation), whose defaults come from `aa540fem.hardware.get_config()`.

**Dof layout.** `TaylorHoodSpace` numbers the unknowns `[u_x (N), u_y (N),
p (Np corner nodes)]`; element-local arrays use the same order
`[u_x nodes, u_y nodes, p corners]` (L = 22 for quad9, 15 for triangle6), and
`_Block` holds the precomputed quadrature data of one cell type.

**Fixed sparsity pattern.** `FlowAssembler` builds one `SparsityPattern`
(`backends/pattern.py`) per space; every matrix (K, M, B, B^T, the Jacobian
parts JN and JS, pseudo-time and theta modifications) is a data vector on it
(`asm.K_data` etc.), combined with vector arithmetic and wrapped as a
`PatternMatrix`. Element kernels return element-local `(ne, L, L)` arrays that a
`ScatterPlan` sums in a fixed order, so assembly is bitwise identical for any
thread count. Dirichlet elimination on a pattern is `PatternDirichlet`
(`linalg/dirichlet.py`, cached masks). This element-local-then-scatter split is
what the numba, MPI and GPU paths all consume, so new physics terms belong in
the kernels, not in a bespoke sparse assembly.

**Kernels.** `backends/numpy_kernels.py` (einsum) is the reference,
`backends/numba_kernels.py` / `numba_sa.py` the threaded twins; `tests/test_backends.py`
checks they agree to 1e-12 over the option matrix. `momentum_terms(...,
jacobian=False)` is the residual-only path; the Newton/PTC loops take a
`residual=` callable so trial points never assemble a Jacobian.

**Nonlinear solves.** `linalg/newton.py::newton_iterate` (damped, optional
frozen Jacobian, accepts a `solver=` from a previous call so time steps share a
factorisation) and `linalg/continuation.py::pseudo_transient`. `method` may be
`"direct"`, `"fieldsplit"` or a callable `A -> solver`;
`incompressible/steady.py::linear_method` resolves it, including `--linear petsc`
from the run configuration. The theta scheme carries its factorisation across
steps (`reuse_jacobian=True`) and refactorises only when modified Newton stops
contracting by 3x per iteration; `info["factorisations"]` reports the count.

**Linear solvers.** `linalg/direct.py::factorise` picks SuperLU or PETSc LU with
MUMPS (`AA540FEM_DIRECT`), falls back on a residual check, and under `mpirun`
uses `parallel/distributed.py::DistributedFactorisation` (replicated matrix,
distributed MUMPS). `linalg/krylov.py::FieldSplitSolver` is FGMRES + PCFIELDSPLIT
(Schur, lower) with the least-squares-commutator preconditioner
`LSCPreconditioner` (algebraic, handles SUPG/grad-div/convection and pinned
pressure dofs), ILU(1) or GAMG on the velocity block, restart 200. PETSc's options
database is global: every solver instance sets its own options prefix, and any
new PETSc solver code must do the same or tests contaminate each other.

**Hardware and run configuration.** `hardware.py` probes cores, memory, GPUs
and MPI; `RunConfig` holds threads, memory budget and backend choices with the
precedence explicit arguments > `AA540FEM_*` environment variables > defaults
(100 % of usable cores, 80 % of free RAM), and applies thread counts to numba,
BLAS and OpenMP so pools never oversubscribe.

**MPI.** `parallel/partition.py` (METIS parts, one-layer halo, local numbering),
`parallel/flow.py::DistributedFlowSystem` (each rank assembles its own sub-mesh
with the ordinary `FlowAssembler`, inserts through a local-to-global map into a
distributed PETSc matrix, halo exchange by VecScatter, Dirichlet by
`MatZeroRowsColumns`, distributed Newton and theta scheme, state gathered back in
the original numbering so post-processing stays serial). Every rank still reads
the whole mesh. Under `mpirun`, code that only rank 0 runs must not call
`factorise` or anything else collective, or the ranks deadlock.

**Meshes.** `examples/meshes/*.msh` are committed and generated by
`make_meshes.py` from the `.geo` files (`cylinder.geo` sizes come only from its
background field; `make(..., size_factor=)` scales them). Tests and benchmarks
read `cylinder_bl.msh` and `flat_plate_bl.msh`; the `_bl` meshes mix quad9
boundary-layer cells with triangle6.

## Conventions that matter here

- Optional dependencies (numba, petsc4py, mpi4py, pymetis, meshio, gmsh, psutil)
  are always import-guarded; tests for them skip when missing and CI runs three
  jobs (base, numba, MPI on 4 ranks).
- Solver changes must reproduce the validated numbers: run the relevant
  example or benchmark and update `README.md` / `docs/performance.md` with the
  measured values, never estimated ones.
- Kill stray background runs by PID or with anchored `pkill -f "^python ..."`
  patterns; an unanchored pattern matches the shell that launched it.
